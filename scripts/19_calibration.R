#!/usr/bin/env Rscript
# Calibration of the across-stratum combination, not of CAMERA itself.
#
# CAMERA is a published competitive test and is not re-validated here. What this project
# adds to it is the step afterwards: one-sided p-values from k strata combined by
# Stouffer's method. That step assumes each stratum's p is uniform under the competitive
# null and that the strata are independent. Both assumptions are checkable by simulation,
# and both can fail in ways that matter for this corpus.
#
# The null simulated is COMPETITIVE, which is the only null the claim is about: the RFO
# genes and the background genes both respond to drought, the RFO genes are co-expressed
# with one another, and there is NO extra RFO-specific displacement. A calibrated test
# rejects such data at the nominal rate. A test that rejects more often is detecting
# co-expression or a shared response, not pathway specificity.
#
# Scenarios, each a separate column of the output:
#   independent          - k independent strata, RFO co-expressed, shared drought response
#   dependent_pair       - as above, but two strata are re-sequencings of the same samples
#   direction_missing    - as above, but genes that fall are likelier to drop out
suppressMessages({library(limma); library(edgeR)})

args <- commandArgs(trailingOnly = TRUE)
getf <- function(flag, default) {
  i <- match(flag, args); if (is.na(i)) default else args[i + 1]
}
n_rep    <- as.integer(getf("--reps", "400"))
k_strata <- as.integer(getf("--strata", "8"))
n_genes  <- as.integer(getf("--genes", "6000"))
n_set    <- as.integer(getf("--set-size", "7"))
n_per    <- as.integer(getf("--per-arm", "6"))
rho_set  <- as.numeric(getf("--set-correlation", "0.55"))
seed     <- as.integer(getf("--seed", "1"))
out_dir  <- getf("--out", "calibration")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(seed)

# One stratum of counts under the competitive null. The set's members share a latent factor
# (so they co-vary at roughly rho_set) and the drought response is drawn from the SAME
# distribution for set and background, which is what makes the null competitive.
sim_stratum <- function(drop_by_direction = FALSE) {
  mu <- 2^runif(n_genes, 2, 11)
  lfc <- rnorm(n_genes, 0, 0.8)                 # shared response, set and background alike
  idx <- seq_len(n_set)                          # the set occupies the first rows
  grp <- rep(c(0, 1), each = n_per)
  latent <- rnorm(2 * n_per)                     # library-level factor the set loads on
  lmu <- outer(log2(mu), rep(1, 2 * n_per)) + outer(lfc, grp)
  noise <- matrix(rnorm(n_genes * 2 * n_per, 0, 0.6), n_genes)
  lmu[idx, ] <- lmu[idx, ] +
    sqrt(rho_set) * matrix(rep(latent, each = n_set), n_set) * 0.6 +
    sqrt(1 - rho_set) * noise[idx, ]
  lmu[-idx, ] <- lmu[-idx, ] + noise[-idx, ]
  cnt <- matrix(rpois(length(lmu), pmax(2^lmu, 0.1)), n_genes)
  rownames(cnt) <- sprintf("g%05d", seq_len(n_genes))
  if (drop_by_direction) {
    # genes that fall are likelier to be unmeasured — the missingness depends on the
    # response, which is exactly the case the competitive null has no defence against
    p_drop <- plogis(-2 - 1.5 * lfc)
    cnt[runif(n_genes) < p_drop, ] <- 0L
  }
  list(counts = cnt, grp = factor(grp), set = rownames(cnt)[idx])
}

one_p <- function(s) {
  d <- DGEList(s$counts)
  keep <- filterByExpr(d, group = s$grp)
  d <- d[keep, , keep.lib.sizes = FALSE]
  members <- intersect(s$set, rownames(d))
  if (length(members) < 3) return(NA_real_)
  d <- normLibSizes(d)
  design <- model.matrix(~ s$grp)
  v <- voom(d, design)
  cam <- camera(v, list(set = which(rownames(v) %in% members)), design,
                contrast = 2, inter.gene.cor = NA)
  if (cam$Direction[1] == "Up") cam$PValue[1] / 2 else 1 - cam$PValue[1] / 2
}

stouffer <- function(p) {
  p <- p[!is.na(p)]
  if (!length(p)) return(NA_real_)
  p <- pmin(pmax(p, 1e-12), 1 - 1e-12)
  pnorm(sum(qnorm(1 - p)) / sqrt(length(p)), lower.tail = FALSE)
}

scenarios <- c("independent", "dependent_pair", "direction_missing")
res <- matrix(NA_real_, n_rep, length(scenarios), dimnames = list(NULL, scenarios))
for (r in seq_len(n_rep)) {
  strata <- lapply(seq_len(k_strata), function(i) sim_stratum())
  res[r, "independent"] <- stouffer(vapply(strata, one_p, numeric(1)))
  # a dependent pair: stratum 2 is the same material re-sequenced, so its p is near-identical
  dup <- strata; dup[[2]] <- strata[[1]]
  res[r, "dependent_pair"] <- stouffer(vapply(dup, one_p, numeric(1)))
  miss <- lapply(seq_len(k_strata), function(i) sim_stratum(drop_by_direction = TRUE))
  res[r, "direction_missing"] <- stouffer(vapply(miss, one_p, numeric(1)))
  if (r %% 25 == 0) cat(sprintf("  %d/%d reps\n", r, n_rep))
}

summ <- do.call(rbind, lapply(scenarios, function(s) {
  p <- res[, s]; p <- p[!is.na(p)]
  hit <- sum(p < 0.05)
  ci <- binom.test(hit, length(p))$conf.int
  data.frame(scenario = s, n_usable = length(p), false_positive_rate = hit / length(p),
             ci_low = ci[1], ci_high = ci[2],
             ks_p_uniform = suppressWarnings(ks.test(p, "punif")$p.value),
             stringsAsFactors = FALSE)
}))
summ$nominal <- 0.05
summ$reps <- n_rep; summ$strata <- k_strata; summ$genes <- n_genes
summ$set_size <- n_set; summ$per_arm <- n_per; summ$set_correlation <- rho_set; summ$seed <- seed
write.csv(summ, file.path(out_dir, "PUBLIC_calibration_summary.csv"), row.names = FALSE)
write.csv(as.data.frame(res), file.path(out_dir, "calibration_pvalues.csv"), row.names = FALSE)
print(summ, row.names = FALSE)
