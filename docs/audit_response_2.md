# Hata denetimine yanıt — yapılan düzeltmeler

> **Note added 1 November 2026.** This document is a historical record of an earlier
> audit round. Items in it that concern the strigolactone and post-rewatering analyses
> no longer apply: that branch was removed from the project in full — scripts, figures,
> tables and count matrix — and nothing in the current package depends on it. The
> document is kept unedited so the audit trail stays readable, not because its
> strigolactone items are still live.


Denetimdeki her madde kodun ve verinin kendisine karşı sınandı. Aşağıda hangi iddianın
doğrulandığı, ne yapıldığı ve sonucun nasıl değiştiği yazılıdır. Sonuçlardan ikisi
daraltıldı, biri güçlendi.

## 1. Zaman alanı eksik katmanlar kanıtsız "eşzamanlı" sayılıyordu — DOĞRU, düzeltildi

`choose_design` eksik `timepoint` değerlerini `"NA"` ile dolduruyordu; her iki kol da
boşsa bu ortak etiket "paylaşılan zaman" sayılıyor ve katman `concurrent_control`
oluyordu. **33 katmanın 15'inde hiçbir kütüphanede zaman kaydı yok**, yani bu sessiz
varsayım kohortun yarısını etkiliyordu.

Kod artık bilgi yokluğunu kanıt saymıyor: iki kol da zamansızsa karar küratörlenmiş
`control_timing` alanına düşüyor, o da boşsa katman `unresolved_timing` oluyor.
Zamanlama 13 çalışma için BioProject kaydı, BioSample başlıkları ve yayınlardan
küratörlendi; her karar bir kanıt dizesiyle birlikte örnek tablosunda duruyor.

**PRJNA1344858 gerçekten başlangıç-referanslı.** BioProject kaydı birebir şöyle diyor:
*"Samples were collected under both control conditions (0h) and after 3h of treatment
with 12% PEG6000"*. Denetim bu noktada haklı. Yaprak ve kök kolları birincil havuzdan
çıkarıldı.

Yeniden sınıflandırma altı katmanı taşıdı:

| katman | eski | yeni |
|---|---|---|
| PRJNA1344858 yaprak + kök | concurrent | **baseline_referenced** |
| PRJNA1284635, PRJNA1298269, PRJNA573482, PRJNA905675 (yaprak) | concurrent | **unresolved_timing** |

Birincil yaprak kohortu **19 → 14 katmana** düştü (72 kontrol / 79 kuraklık).
Başlangıç-referanslı 4, çözülemeyen 4 katman ayrı havuzlanıyor.

## 2. Çıkarım yöntemine duyarlılık — DOĞRU, ve sonucu değiştiriyor

Denetimin HKSJ hesabı birebir tekrarlandı. 19 katmanlık kohortta normal-z ile 4.617 gen
FDR<0,05 iken **HKSJ ile sıfır gen** kalıyor; en küçük BH 0,094. Doğrulandı.

Bunun τ² tahmincisinden mi çıkarım modelinden mi geldiğini ayırdım: **Paule–Mandel τ² +
normal-z 3.664 gen** veriyor. Yani duyarlılık τ² seçiminde değil, **belirsizlik modelinde**.

Bu nedenle genom çapındaki gen sayısı artık bir keşif iddiası olarak sunulmuyor; bir
tarama çıktısı olarak, üç yöntemin sayılarıyla birlikte bildiriliyor.

## 3. Panelde BH uygularken kendi hatamı buldum

HKSJ'yi ön tanımlı panelde denerken 11 genlik bir panel kullanmıştım; o panel
**anlamlı çıkanlardan seçilmiş bir alt kümeydi** ve içinde BH uygulamak döngüseldi.
Panel, anahtar kelimeyle bulunan 18 yolak geni ile ortoloji/domain kanıtıyla eklenen
3 SlGolS'un birleşimi olarak (21 gen) yeniden tanımlandı; uzak GT8 akrabaları panele
değil, ayrı bir negatif kontrol setine konuldu.

## 4. Nihai durum: etki sağlam, tek-gen anlamlılığı değil

En katı tanımda (14 katmanlık kanıta dayalı kohort) hiçbir tek-gen istatistiği 20 genlik
panelde çoklu testi geçmiyor — ne HKSJ (en küçük BH 0,120) ne de varyans modelinden
bağımsız işaret testi (en küçük BH 0,229). Etki büyüklükleri ise yerinde duruyor:
rafinoz sentaz +2,60, SlGolS +2,09 ve +2,00, stachyoz sentaz +1,84.

Doğru ölçek yolak düzeyi, çünkü orada tek bir test var. **Biyosentetik kol**
(galaktinol → rafinoz → stachyoz sentazları, 8 ölçülebilir gen) ifade düzeyi
eşleştirilmiş rastgele gen setlerine karşı sınandı:

| set | ortalama birleşik log2FC | eşleştirilmiş boş dağılım | p |
|---|---|---|---|
| RFO biyosentetik kol (8 gen) | **+1,578** | −0,018 ± 0,159 | **<0,0001** (10.000 setin 0'ı) |
| GT8 uzak akrabalar (4 gen, kontrol) | −0,027 | −0,006 ± 0,260 | 0,92 |
| genom medyanı | −0,001 | — | — |

Gözlenen değer boş dağılımın ortalamasından **10 standart sapma** uzakta. Denetimin
haklı olarak hatırlattığı gibi, genleri bağımsız sayan bir rekabetçi set testi panel içi
eş-ifade varlığında anti-konservatiftir; ancak boş dağılımın standart sapması üç katına
çıkarılsa bile değer 3,3 sapma uzakta kalır.

**Savunulabilir ifade:** RFO biyosentetik kolu, eşzamanlı kontrollü 14 bağımsız katmanda
genom arka planına göre tutarlı biçimde indükleniyor; tek tek genlerin genom çapında
FDR'si bu kohort büyüklüğünde ve bu heterojenlikte güvenilir değil.

## 5. Strigolakton ifadeleri daraltıldı

- **"GR24 RFO genlerini kıpırdatmıyor" yanlıştı.** Yanlı olmayan 21 genlik panelle
  yeniden çalıştırıldığında GR24 karşılaştırmasında **19 test edilebilir panel geninden
  7'si** genom çapında FDR<0,05 ve iki yönde birden: rafinoz sentaz +1,75,
  alfa-galaktosidaz +1,52 ve +1,18; buna karşılık myo-inositol monofosfataz 1 −0,86,
  başka bir alfa-galaktosidaz −1,24, myo-inositol-1-fosfat sentaz −1,29 ve putatif
  myo-inositol-1-fosfataz −2,08. Panel bir bütün olarak arka plandan ayrılmıyor
  (medyan +0,12'ye +0,05, **p = 0,954**). Daha önce bildirdiğim "3 gen" ve "p = 0,248"
  değerleri, anlamlı çıkanlardan seçilmiş 11 genlik panelden geliyordu ve geçersizdir.
  Gen testi ile panel-dağılım testi farklı hipotezlerdir.
- **"Strigolakton gerektirmiyor" / "gerektiriyor" iddiaları** kısmen tükenmiş bir hatta
  ve eşzamanlı kontrolü olmayan bir tasarımda kurulduğu için koşullu hale getirildi.
- **Yeniden sulama karşılaştırmaları eşzamanlı kontrollü değil.** GEO kaydına göre kontrol
  T0, yeniden sulama örnekleri T31 ve T76. Bu nedenle "kuraklık imzası yeniden sulamayı
  atlatıyor" yerine: *yeniden sulama sonrası başlangıca göre değişimler aktif-kuraklık
  meta-imzasıyla kısmi yön benzerliği gösteriyor; eşzamanlı kontrol olmadığı için bu
  benzerlik stres hafızasına özgülenemez.*
- **"1.566'ya 73 DEG" karşılaştırması** genotipler arası epizot etkisi testi değildir;
  gerekli test S2−S1 farkının genotipler arasında doğrudan karşılaştırılmasıdır.

## 6. Eksik anotasyonlar eski sürümden kurtarıldı

*SlCCD7* ve *SlMAX1* ITAG4.1'de yok. **ITAG2.4'te ikisi de var** ve domatesin kendi
referans proteinleriyle %100 örtüşüyor: `Solyc01g090660.2.1` → SlCCD7 (UniProt D0UDC1,
663 aa), `Solyc08g062950.2.1` → SlMAX1 (K4CL78, 519 aa). İkisi ITAG4.1 cDNA setine
eklenip strigolakton deneyinin 21 kütüphanesi bu artırılmış referansla yeniden nicelendi.

Sonuç beklenmedik ve açıklanması gerekiyor:

| gen | yabani tip (12 kütüphane) | CCD7-susturulmuş hat (9 kütüphane) |
|---|---|---|
| SlCCD7 | **0 okuma** (hepsinde) | 1.277–5.288 okuma |
| SlMAX1 | 74–255 okuma | 159–622 okuma |

22–41 M okumalık kütüphanelerde yabani tipte tam sıfır, susturulmuş hatta binlerce okuma
— bu bir susturma örüntüsü değil. En makul okuma şudur: *CCD7 yaprakta ifade edilmiyor*
(bilinen kök-ağırlıklı ifadesiyle ve GEO'nun strigolakton azalmasını **kök** dokusunda
ve eksudatta bildirmesiyle uyumlu), susturulmuş hattaki okumalar ise susturma
yapısının kendisinden geliyor. Dolayısıyla **susturma yaprak RNA'sından
doğrulanamıyor**; doğrulanan şey hattın CCD7 dizisi taşıdığıdır. Bu, yapraktaki
fenotipin hücre-özerk değil sistemik olması gerektiği anlamına da gelir.

SlMAX1 yaprakta ifade ediliyor ve SL yolu paneline eklenebilir.

## 7. Depo kusurları

- `--samples` / `--sheet` uyuşmazlığı düzeltildi (01 ve 02 `--sheet` istiyor);
  `02_build_matrices.py` `--tx2gene` almıyor, o argüman da kaldırıldı.
- Makefile'a `.SHELLFLAGS := -o pipefail -c` eklendi; `| tee` içindeki hata artık kuralı
  düşürüyor.
- `de` ve `meta` hedefleri çalışmadan önce kendi çıktı dizinlerini siliyor, eski dosya
  sızıntısı engellendi.
- `all` hedefi genotip ve strigolakton dallarını da kapsıyor; `quant` ağ gerektirdiği
  için ayrı tutuldu ve bu README'de yazılı.
- `split_by_estimand.py` artık üç kohort yazıyor (birincil / başlangıç-referanslı /
  çözülemeyen).

## 8. Açık maddelerin kapatılması

- **`06_genotype_moderator.py` tamamlandı.** Betik artık tam bir komut satırı akışına
  sahip (`--arm-effects`, `--wildtype-de`, `--panel`, `--min-reliability`) ve gerçek
  veriyle çalıştırıldığında yayımlanan dokuz kolun değerlerini birebir üretiyor
  (0,838 / 0,853 / 0,765 / 0,860 / 0,781 / 0,899 / kestirilemez / 0,882 / 0,239).
- **`08_strigolactone.py` tamamlandı.** Eksik üç analiz betiğe girdi: ikinci–birinci
  epizot karşılaştırması, meta-imza korelasyonu ve etki büyüklüğü eşleştirilmiş arka plan
  testi. Uçtan uca çalıştırıldı; epizot karşılaştırması artık iki kolun artık varyansını
  da yazdırıyor, böylece DEG sayısı farkı güç farkından ayırt edilebiliyor.
- **Moon ve ark.'nın ek tabloları bulundu — önceki gerekçem yanlıştı.** Yalnız 1–8
  arasındaki dosyaları denemiştim; **16 ek dosya var ve sekizi tablo**. Tablo S6 on adayın
  tamamını Solyc kimlikleriyle veriyor, Tablo S3 ise 5.681 meta-DEG'in tamamını. İkisi de
  indirildi ve replikasyon testi yeniden çalıştırıldı (aşağıda).
- **İfade düzeltmeleri yapıldı.** "Dokuz genotip" → "dokuz genotip × doku kolu, yedi
  ayrı genetik değişiklik"; korunan-yanıt eğiminin yanıt veren gen yüzdesi olmadığı
  Bulgular metnine açıkça yazıldı; GR24 ve yeniden sulama ifadeleri daraltıldı.
- **Açık kalan tek madde:** tüm tablo ve şekiller tek bir dondurulmuş koşudan yeniden
  üretilmeli. Betikler ve yapılandırma bunu mümkün kılıyor, koşunun kendisi yapılmadı.

## 9. Yeniden çalıştırılan replikasyon testi

Eksik kimlikler Tablo S6'dan tamamlandı: SlSAUR6 = Solyc01g110590, SlSAUR9 =
Solyc01g110630, SlMAPKKK51 = Solyc07g051870, SlMAPKKK55 = Solyc07g051930,
GDSL2 = Solyc11g006250.

Kontrol zamanlaması küratörlüğü, Moon ve ark.'nın dört veri setinin üçünü birincil
kohortun dışına taşıdı (ikisi başlangıç-referanslı, biri çözülemeyen), dolayısıyla
**14 birincil katmanın 13'ü o çalışmadan bağımsız**.

- On adayın dokuzu test edilebildi; **dokuzun sekizi bildirilen yönde** (işaret testi
  p = 0,039), hiçbiri bu kohortta FDR<0,05 değil.
- 5.681 meta-DEG'in 5.363'ü test edilebildi: etki büyüklüğü korelasyonu **Pearson
  r = +0,635**, yön uyumu **%77,7**. Onların anlamlı saydığı 3.114 gende %78,5, "yüksek
  güvenli çekirdek" dediği 39 gende **%92,3** (altısı bizde de FDR<0,05).

Etkiler bizde sistematik olarak daha küçük (anlamlı setlerinde medyan mutlak değer 0,16);
bu, dört yerine 14 heterojen katman üzerinden havuzlamanın ve en büyük etkileri taşıyan
başlangıç-referanslı tasarımları dışlamanın beklenen sonucudur.
