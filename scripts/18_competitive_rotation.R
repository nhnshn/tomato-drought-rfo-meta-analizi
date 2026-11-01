#!/usr/bin/env Rscript
# Sample-level competitive and self-contained gene-set tests, one stratum at a time.
#
# Why this exists
# ---------------
# 11_geneset_test.py builds its null by resampling genes, which treats them as independent,
# and then widens that null by a variance-inflation factor estimated from 13 pooled effect
# estimates. Both steps are approximations of the same thing: the correlation between the
# set's members. CAMERA does it properly, estimating the inter-gene correlation from the
# libraries themselves inside each study, where there are tens of residual degrees of
# freedom rather than thirteen points. fry adds the self-contained question (did the set
# move at all?) using sample rotation, which preserves the correlation by construction.
#
# The two questions are different and both are reported. A set can move (fry) without
# moving more than the rest of the transcriptome (camera); only camera speaks to the
# claim this project makes.
#
# Strata are read from the design table the differential-expression step wrote, so the
# libraries and blocking terms are exactly those that were fitted there rather than a
# reconstruction from the samplesheet.
suppressMessages({library(limma); library(edgeR)})

args <- commandArgs(trailingOnly = TRUE)
get <- function(flag, default = NA) {
  i <- match(flag, args)
  if (is.na(i)) {
    if (is.na(default)) stop(sprintf("%s is required", flag))
    return(default)
  }
  args[i + 1]
}
design_file <- get("--design")
matrix_dir  <- get("--matrix-dir")
sets_file   <- get("--sets")
out_dir     <- get("--out")
estimand    <- get("--estimand", "concurrent_control")
tissue_tag  <- get("--tissue", "leafshoot")
dependent   <- get("--dependent-group", "")
drop_nontarget <- "--drop-nontarget-panel" %in% args
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

dsn  <- read.csv(design_file)
sets <- read.csv(sets_file)
# Every gene declared in any RFO set; the GT8 control set is not a panel member
# and stays in the background under both background definitions.
panel_genes <- unique(sets$gene[sets$set != "gt8_control"])
dsn  <- dsn[dsn$estimand == estimand & grepl(paste0("__", tissue_tag, "$"), dsn$stratum), ]
strata <- sort(unique(dsn$stratum))
cat(sprintf("%d strata with estimand '%s' and tissue '%s'\n", length(strata), estimand, tissue_tag))

set_names <- sort(unique(sets$set))
rows <- list()

for (st in strata) {
  meta <- dsn[dsn$stratum == st, ]
  study <- meta$study[1]
  f <- file.path(matrix_dir, paste0(study, "_counts.csv.gz"))
  if (!file.exists(f)) { cat(sprintf("  %s: no count matrix, skipped\n", st)); next }
  cnt <- read.csv(f, row.names = 1, check.names = FALSE)
  meta <- meta[meta$sample_id %in% colnames(cnt), ]
  if (length(unique(meta$treatment_class)) < 2) {
    cat(sprintf("  %s: one arm only, skipped\n", st)); next
  }
  cnt <- cnt[, meta$sample_id, drop = FALSE]

  grp <- factor(meta$treatment_class, levels = c("control", "drought"))
  # Blocking terms are whatever the fitted design carried, kept only where they still vary
  # and leave the model estimable in this subset.
  terms <- setdiff(intersect(c("subtissue", "timepoint", "genotype"), colnames(meta)), NA)
  terms <- terms[sapply(terms, function(t) length(unique(meta[[t]])) > 1)]
  mm <- NULL
  repeat {
    form <- if (length(terms)) paste("~", paste(c(terms, "grp"), collapse = " + ")) else "~ grp"
    md <- cbind(meta, grp = grp)
    cand <- model.matrix(as.formula(form), data = md)
    if (qr(cand)$rank == ncol(cand) && ncol(cand) < nrow(md)) { mm <- cand; break }
    if (!length(terms)) { mm <- model.matrix(~ grp); break }
    terms <- terms[-length(terms)]
  }
  coef_i <- which(colnames(mm) == "grpdrought")

  d <- DGEList(round(as.matrix(cnt)))
  keep <- filterByExpr(d, group = grp)
  d <- d[keep, , keep.lib.sizes = FALSE]
  d <- normLibSizes(d)
  v <- voom(d, mm)

  for (sn in set_names) {
    members <- intersect(sets$gene[sets$set == sn], rownames(v))
    if (length(members) < 3) next
    # camera's background is every other measured gene. Non-target panel members
    # therefore sit in the background, where they drag the comparison toward the
    # set and make the test conservative. --drop-nontarget-panel removes them so
    # the background matches the one the stratum-level contrast uses (panel genes
    # excluded), which is what the replication pre-registration specifies.
    vv <- v
    if (drop_nontarget) {
      others <- setdiff(panel_genes, members)
      k <- !(rownames(vv) %in% others)
      vv <- vv[k, ]
    }
    idx <- list(which(rownames(vv) %in% members)); names(idx) <- sn
    cam <- camera(vv, idx, mm, contrast = coef_i, inter.gene.cor = NA)
    fr  <- fry(vv, idx, mm, contrast = coef_i)
    # camera's PValue is two-sided; convert to a one-sided p for induction so that
    # strata can be combined on a common, directional scale.
    p_up <- if (cam$Direction[1] == "Up") cam$PValue[1] / 2 else 1 - cam$PValue[1] / 2
    rows[[length(rows) + 1]] <- data.frame(
      stratum = st, gene_set = sn, n_genes = length(members),
      n_genes_measured = nrow(v), inter_gene_cor = cam$Correlation[1],
      direction = cam$Direction[1], camera_p_two_sided = cam$PValue[1],
      camera_p_up = p_up, fry_p = fr$PValue[1], fry_p_mixed = fr$PValue.Mixed[1],
      stringsAsFactors = FALSE)
  }
  cat(sprintf("  %s: %d libraries, %d genes, design %s\n",
              st, ncol(v), nrow(v), paste(colnames(mm)[-1], collapse = "+")))
}

per_stratum <- do.call(rbind, rows)
write.csv(per_stratum, file.path(out_dir, "PUBLIC_rotation_per_stratum.csv"), row.names = FALSE)

# Combine across strata with Stouffer's method on the one-sided camera p-values. This
# assumes the strata are independent; strata named in --dependent-group are collapsed to
# the first of them in a declared sensitivity rather than counted twice.
dep <- if (nchar(dependent)) strsplit(dependent, ",")[[1]] else character(0)
combine <- function(df, label) {
  do.call(rbind, lapply(split(df, df$gene_set), function(g) {
    p <- pmin(pmax(g$camera_p_up, 1e-12), 1 - 1e-12)
    z <- sum(qnorm(1 - p)) / sqrt(nrow(g))
    data.frame(cohort = label, gene_set = g$gene_set[1], k = nrow(g),
               n_up = sum(g$direction == "Up"),
               median_inter_gene_cor = median(g$inter_gene_cor),
               stouffer_z = z, combined_p_up = pnorm(z, lower.tail = FALSE),
               stringsAsFactors = FALSE)
  }))
}
res <- combine(per_stratum, "all strata")
# vapply, not sapply: on an empty `dep` sapply returns list(), which cannot index a
# character vector, and the whole combined table is lost after the per-stratum fits.
present <- dep[vapply(dep, function(x) any(grepl(x, per_stratum$stratum)), logical(1))]
if (length(present) > 1) {
  drop <- present[-1]
  kept <- per_stratum[!Reduce(`|`, lapply(drop, function(x) grepl(x, per_stratum$stratum))), ]
  res <- rbind(res, combine(kept, "independent samples"))
}
res$bh_within_cohort <- ave(res$combined_p_up, res$cohort, FUN = function(p) p.adjust(p, "BH"))
write.csv(res, file.path(out_dir, "PUBLIC_rotation_combined.csv"), row.names = FALSE)
cat("\n"); print(res[order(res$cohort, res$combined_p_up), ], row.names = FALSE)
