# =============================================================================
#  Genetic diversity and population structure of Cucumis metuliferus
#  (African horned cucumber) from DArT-seq SNP genotyping
#
#  Reproducible, single-file analysis accompanying:
#    Koopa et al. (2026) "DArT-based characterisation of genetic diversity in
#    genotypes of Cucumis metuliferus collected from SACU." PLOS ONE.
#
#  What this script does, in order (section numbers match the code below):
#    1. Read the DArT SNP report and the sample metadata; populations = the
#       9 geographic sampling locations; exclude MP1 (contaminated library)
#    2. Quality-filter the SNPs (-> 119 individuals) and report
#       descriptive statistics of the retained dataset
#    3. Diversity indices (Ho, He, FIS, PIC), overall and per location
#    4. PCA, plus a check that skipping the MAF filter does not bias it
#    5. Genomic relationship matrix (identity by descent among genotypes)
#    6. Dendrogram of genetic distances (Czekanowski)
#    7. Population structure with fastSTRUCTURE (K = 1..10, 10 replicates);
#       K chosen by the mean marginal likelihood
#    8. AMOVA (partitioning of variation among vs within locations)
#    9. Selection and validation of a "core" set of the most informative
#       genotypes
#   Figures 1-10 of the manuscript are written to figures/ along the way;
#   tables to outputs/.
#
#  Everything uses the dartRverse suite: https://github.com/green-striped-gecko
#
#  To run, from THIS folder (so the relative data/ and outputs/ paths resolve):
#      setwd("path/to/cucumis-metuliferus-dartseq")
#      source("analysis.R")
# =============================================================================


# -----------------------------------------------------------------------------
# 0. Packages and settings
# -----------------------------------------------------------------------------
# dartRverse loads the sub-packages used below (dartR.base, dartR.popgen).
# One-time install:  install.packages("dartRverse"); dartRverse::dartRverse_install()
library(dartRverse)   # read/filter genlight, PCA, distances, heterozygosity
library(poppr)        # AMOVA
library(rrBLUP)       # A.mat(): additive genomic relationship matrix
library(ggplot2)      # figures 1, 4, 5

# Optional packages: a figure that needs one of them is skipped if it is absent.
has <- function(p) requireNamespace(p, quietly = TRUE)
for (p in c("pheatmap", "dendextend", "reshape2", "patchwork", "hierfstat",
            "sf", "rnaturalearth", "ggspatial"))
  if (!has(p)) message("optional package '", p, "' not installed; the figure(s) using it will be skipped")

# fastSTRUCTURE (section 7) is an external program, not an R package. Point
# these at the fastStructure and plink executables (or set the env vars).
FASTSTRUCTURE <- Sys.getenv("FASTSTRUCTURE", path.expand("~/programs/fastStructure"))
PLINK         <- Sys.getenv("PLINK",         path.expand("~/programs/plink"))

set.seed(1)                                   # reproducible random steps (section 9)
dir.create("outputs", showWarnings = FALSE)   # tables and intermediate results
dir.create("figures", showWarnings = FALSE)   # manuscript figures (300 dpi PNG)

# colour-blind-safe palette for the 9 sampling locations
PAL9 <- c("#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE",
          "#AA3377", "#BBBBBB", "#EE8866", "#000000")
# Figures follow the PLOS ONE specification: at most 7.5 in wide and 8.75 in
# high, 300 dpi, Arial 8-12 pt, multi-panel labels (A), (B).
FONT <- "Arial"
png300 <- function(name, w = 7.5, h = 6) {
  png(file.path("figures", name), width = w, height = h, units = "in", res = 300)
  par(family = FONT, cex = 0.9)
}
theme_set(theme_bw(base_family = FONT, base_size = 10))


# -----------------------------------------------------------------------------
# Supporting functions  (used only in section 9 - safe to skip on a first read)
# -----------------------------------------------------------------------------
# The main analysis flow below just calls top_ind() and pairwise_identity();
# you do not need to read these functions to follow the analysis.

# Hill-number diversity of a set of counts (alpha diversity, order q):
#   q = 0 -> richness,  q = 1 -> Shannon,  q = 2 -> Simpson.  (Chao et al. 2014.)
d.chao <- function(counts, q) {
  p <- counts / sum(counts)          # relative frequencies
  p <- p[p > 0]                      # drop zeros (log is undefined at 0)
  if (q == 1) exp(-sum(p * log(p)))  # Shannon index (the q -> 1 limit)
  else        (sum(p^q))^(1 / (1 - q))
}

# Number of polymorphic loci in a genlight (loci with both alleles present).
n_poly <- function(g) { a <- colMeans(as.matrix(g), na.rm = TRUE) / 2; sum(a > 0 & a < 1, na.rm = TRUE) }

# Percentage of identical genotype calls for every pair of individuals,
# counting only loci called in both. Returns an individuals x individuals matrix.
pairwise_identity <- function(X) {
  n <- nrow(X); out <- matrix(NA_real_, n, n, dimnames = list(rownames(X), rownames(X)))
  for (i in seq_len(n)) {
    both <- !is.na(X[i, ]) & !is.na(t(X))            # loci called in i and in each other row
    same <- (X[i, ] == t(X)) & both
    out[i, ] <- 100 * colSums(same) / colSums(both)
  }
  out
}

# Core-subset selection, exactly as described in the Methods:
#   score every genotype on four complementary "informativeness" metrics,
#     AR   - number of SNPs at which it carries a rare allele (MAC < MAC_min)
#     Het  - individual observed heterozygosity
#     Sha  - alpha diversity of its own genotype profile (Hill number of
#            order q; q = 0, the default, is the number of SNPs at which the
#            genotype carries the alternate allele, as used in Table 4)
#     Diff - mean Czekanowski distance to all other genotypes
#   shortlist the top `top.ind` genotypes of EACH metric (union), then drop,
#   from every pair sharing > rep_threshold identical calls, the member with
#   more missing data (gl.report.replicates). What remains is the core.
# Returns a list: table (one row per genotype, values, ranks, flags) and
# replicates (the > rep_threshold pairs).
top_ind <- function(x, MAC_min = 5, top.ind = 20, rep_threshold = 0.98, q = 0) {
  pop(x) <- rep("pop", nInd(x))            # treat all individuals as one group
  x <- gl.filter.allna(x, verbose = 0)
  x <- gl.recalc.metrics(x, verbose = 0)
  X <- as.matrix(x)                        # individuals x SNPs, coded 0/1/2/NA

  # (AR) rare-allele carriers -------------------------------------------------
  thr <- MAC_min / (nInd(x) * 2)           # minor-allele-count expressed as a frequency
  af  <- gl.alf(x); af$loc <- rownames(af) # per-locus allele frequencies
  rare_ref <- af$loc[af$alf1 <= thr]       # loci whose reference allele is rare
  rare_alt <- af$loc[af$alf2 <= thr]       # loci whose alternate allele is rare
  carriers <- c(
    unlist(lapply(rare_ref, function(l) rownames(X)[which(X[, l] %in% c(0, 1))])),
    unlist(lapply(rare_alt, function(l) rownames(X)[which(X[, l] %in% c(1, 2))])))
  AR <- as.data.frame(table(factor(carriers, levels = indNames(x))))  # 0 if none
  names(AR) <- c("ind.name", "AR_value"); AR$ind.name <- as.character(AR$ind.name)

  # (Het) individual observed heterozygosity ----------------------------------
  het <- gl.report.heterozygosity(x, method = "ind", plot.display = FALSE, verbose = 0)
  het <- data.frame(ind.name = het[[1]], Het_value = het$Ho)

  # (Sha) alpha diversity of each genotype profile ----------------------------
  sha <- data.frame(
    ind.name  = rownames(X),
    Sha_value = apply(X, 1, function(g) d.chao(g[!is.na(g) & g > 0], q = q)))

  # (Diff) mean genetic distance to every other genotype ----------------------
  d   <- as.matrix(gl.dist.ind(x, method = "Czekanowski", plot.display = FALSE, verbose = 0))
  dif <- data.frame(ind.name = rownames(d), Diff_value = rowMeans(d, na.rm = TRUE))

  # combine and convert each metric to a rank (1 = most informative) ----------
  tab <- Reduce(function(a, b) merge(a, b, by = "ind.name"), list(AR, het, sha, dif))
  tab$AR   <- rank(-tab$AR_value,   ties.method = "min")
  tab$Het  <- rank(-tab$Het_value,  ties.method = "min")
  tab$Sha  <- rank(-tab$Sha_value,  ties.method = "min")
  tab$Diff <- rank(-tab$Diff_value, ties.method = "min")
  tab$mean_rank <- rowMeans(tab[, c("AR", "Het", "Sha", "Diff")])

  # shortlist = top `top.ind` of any metric; drop near-identical replicates ---
  tab$shortlisted <- tab$AR <= top.ind | tab$Het <= top.ind | tab$Sha <= top.ind | tab$Diff <= top.ind
  reps <- gl.report.replicates(gl.keep.ind(x, ind.list = tab$ind.name[tab$shortlisted], verbose = 0),
                               perc_geno = rep_threshold, plot.out = FALSE, verbose = 0)
  tab$replicate_dropped <- tab$ind.name %in% reps$ind.list.drop
  tab$select <- tab$shortlisted & !tab$replicate_dropped
  tab <- tab[order(!tab$select, tab$mean_rank), ]     # selected first, then by mean rank
  rownames(tab) <- NULL
  list(table = tab, replicates = reps$table.rep)
}


# -----------------------------------------------------------------------------
# 1. Read the data and assign populations
# -----------------------------------------------------------------------------
# gl.read.dart() reads the DArT two-row SNP report into a "genlight" object and
# attaches the per-sample metadata (id, sampling location, coordinates).
gl <- gl.read.dart(
  filename     = "data/Report_DCu24-9392_SNP_mapping_2.csv",
  ind.metafile = "data/metadata_with_decimal_coordinates.csv",
  verbose = 1)

gl <- gl.sort(gl, sort.by = "ind", order.by = indNames(gl))  # order-independent

# MP1 is excluded. Its reported genotype came from a technical-replicate
# library that was contaminated during library preparation with DNA from
# outside this project (read counts of the two MP1 libraries in DArT order
# DCu24-9392: the foreign alleles occur in one library only). The DNA sample
# itself is sound, but the reported calls are not; 119 genotypes remain, and
# loci that were polymorphic only because of MP1 are removed.
gl <- gl.drop.ind(gl, ind.list = "MP1", verbose = 0)
n_before <- nLoc(gl)
gl <- gl.filter.monomorphs(gl, verbose = 0)
cat(sprintf("Excluded MP1: %d genotypes; %d loci monomorphic without MP1 removed, %d SNPs remain\n",
            nInd(gl), n_before - nLoc(gl), nLoc(gl)))

# Population = geographic SAMPLING LOCATION (metadata column "pop2"; 9 sites).
# This is the grouping used for the AMOVA and for colouring the plots.
pop(gl) <- gl$other$ind.metrics$pop2
popf <- factor(pop(gl)); locs <- levels(popf); names(PAL9) <- locs
write.csv(as.data.frame(table(location = popf)), "outputs/genotypes_per_location.csv", row.names = FALSE)

# Figure 1: map of the sampling sites, one point per accession (public-domain
# Natural Earth basemap).
if (all(sapply(c("sf", "rnaturalearth", "ggspatial"), has))) {
  im <- gl$other$ind.metrics
  xy <- data.frame(loc = as.character(popf),
                   lon = suppressWarnings(as.numeric(im$lon)),
                   lat = suppressWarnings(as.numeric(im$lat)))
  xy <- xy[is.finite(xy$lon) & is.finite(xy$lat), ]
  if (nrow(xy) > 0) {
    world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
    p1 <- ggplot(world) + geom_sf(fill = "grey96", colour = "grey70", linewidth = 0.2) +   # one point per accession
      coord_sf(xlim = range(xy$lon) + c(-3, 3), ylim = range(xy$lat) + c(-3, 3), expand = FALSE) +
      ggspatial::annotation_scale(location = "bl") +
      ggspatial::annotation_north_arrow(location = "tr", style = ggspatial::north_arrow_fancy_orienteering()) +
      geom_point(data = xy, aes(lon, lat, fill = loc), shape = 21, size = 2.6, colour = "black", alpha = 0.8) +
      scale_fill_manual(values = PAL9, name = "Sampling location") +
      labs(x = NULL, y = NULL) + theme_bw()
    ggsave("figures/Figure1.png", p1, width = 7.5, height = 5.8, dpi = 300)
  } else message("Figure 1 skipped: no usable coordinates in the metadata")
} else message("Figure 1 skipped: needs sf, rnaturalearth and ggspatial")


# -----------------------------------------------------------------------------
# 2. Quality filtering and descriptive statistics
# -----------------------------------------------------------------------------
# Three standard DArT-seq filters, applied in sequence (Table 1):
gl <- gl.filter.callrate(gl, threshold = 0.80, verbose = 1)         # >= 80% called
gl <- gl.filter.rdepth(gl, lower = 5, upper = 1000, verbose = 1)    # read depth 5-1000
gl <- gl.filter.reproducibility(gl, threshold = 0.99, verbose = 1)  # >= 99% reproducible

# NOTE: we deliberately do NOT apply a minor-allele-frequency (MAF) filter or LD
# pruning. Rare alleles are informative for the diversity and core-subset aims;
# section 4 checks that retaining them does not change the PCA structure.

cat(sprintf("Filtered data: %d individuals x %d SNPs\n", nInd(gl), nLoc(gl)))
saveRDS(gl, "outputs/genlight_filtered.rds")

# Descriptive statistics of the retained SNPs (call rate, missingness, read
# depth, MAF, PIC), as reported in the Methods.
M    <- as.matrix(gl)                                  # individuals x SNPs, 0/1/2/NA
lm   <- gl$other$loc.metrics
af   <- colMeans(M, na.rm = TRUE) / 2                  # alternate-allele frequency per locus
maf  <- pmin(af, 1 - af)
miss <- colMeans(is.na(M))
desc <- data.frame(
  metric = c("mean_call_rate", "mean_missing_per_locus", "max_missing_per_locus",
             "mean_missing_per_individual", "max_missing_per_individual",
             "median_read_depth", "min_read_depth", "max_read_depth",
             "mean_MAF", "median_MAF", "prop_loci_MAF_below_0.05", "mean_PIC"),
  value  = round(c(mean(1 - miss), mean(miss), max(miss),
                   mean(rowMeans(is.na(M))), max(rowMeans(is.na(M))),
                   median(lm$rdepth), min(lm$rdepth), max(lm$rdepth),
                   mean(maf), median(maf), mean(maf < 0.05), mean(lm$AvgPIC, na.rm = TRUE)), 4))
write.csv(desc, "outputs/descriptive_stats.csv", row.names = FALSE)
print(desc)

# Alignment of the marker sequences to the two reference genomes carried in
# the DArT report (cucumber Chinese Long v3, melon DHL92 v4). Annotation
# only: no analysis below uses genome position.
aligned <- function(g) { ch <- lm[[paste0("Chrom_", g)]]; !(is.na(ch) | ch == "" | ch == "0") }
cu <- aligned("Cucumber_ChineseLong_v3"); me <- aligned("Cucumis_melo_DHL92_v4")
map <- data.frame(reference = c("cucumber_ChineseLong_v3", "melo_DHL92_v4", "either", "both"),
                  prop_filtered_SNPs_aligned = round(c(mean(cu), mean(me), mean(cu | me), mean(cu & me)), 4))
write.csv(map, "outputs/mapping_summary.csv", row.names = FALSE)
print(map)


# -----------------------------------------------------------------------------
# 3. Diversity indices
# -----------------------------------------------------------------------------
He <- mean(2 * af * (1 - af), na.rm = TRUE)               # expected heterozygosity
Ho <- mean(colMeans(M == 1, na.rm = TRUE), na.rm = TRUE)  # observed heterozygosity
cat(sprintf("Overall  Ho = %.4f   He = %.4f   FIS = %.4f   PIC = %.4f\n",
            Ho, He, 1 - Ho / He, mean(lm$AvgPIC, na.rm = TRUE)))

# Per-location Ho, He and FIS, written to a table for reference.
het_pop <- gl.report.heterozygosity(gl, method = "pop", plot.display = FALSE, verbose = 0)
write.csv(het_pop, "outputs/diversity_by_location.csv", row.names = FALSE)


# -----------------------------------------------------------------------------
# 4. Principal component analysis, and its robustness to MAF filtering
# -----------------------------------------------------------------------------
pca <- gl.pcoa(gl, nfactors = 5, plot.out = FALSE, verbose = 0)
ve  <- round(100 * pca$eig[1:3] / sum(pca$eig[pca$eig > 0]), 1)  # % variance explained
cat(sprintf("PCA variance: PC1 %.1f%%  PC2 %.1f%%  PC3 %.1f%%\n", ve[1], ve[2], ve[3]))

# Figure 2: (a) PC1 vs PC2 and (b) PC1 vs PC3, coloured by sampling location.
pc_plot <- function(i, j, col, pch = 19, main = "") {
  plot(pca$scores[, i], pca$scores[, j], col = col, pch = pch, cex = 1.3, main = main,
       xlab = sprintf("PC%d (%.1f%%)", i, ve[i]), ylab = sprintf("PC%d (%.1f%%)", j, ve[j]))
}
png300("Figure2.png", w = 7.5, h = 3.9)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 0.5), oma = c(0, 0, 0, 5.5), xpd = NA, family = FONT, cex = 0.75)
pc_plot(1, 2, PAL9[popf], main = "(A)"); pc_plot(1, 3, PAL9[popf], main = "(B)")
legend(par("usr")[2] * 1.04, par("usr")[4], legend = locs, col = PAL9, pch = 19, bty = "n", title = "Location", cex = 0.9)
dev.off()

# Robustness: repeat the PCA after filtering to MAF >= 0.05 and correlate the
# individual scores with those of the full dataset (r ~ 1 => structure unchanged).
gl_maf  <- gl.filter.maf(gl, threshold = 0.05, plot.display = FALSE, verbose = 0)
pca_maf <- gl.pcoa(gl_maf, nfactors = 5, plot.out = FALSE, verbose = 0)
rob <- data.frame(PC = 1:3, n_SNPs_MAF_filtered = nLoc(gl_maf),
  r_full_vs_MAF_filtered = round(abs(sapply(1:3, function(k) cor(pca$scores[, k], pca_maf$scores[, k]))), 3))
write.csv(rob, "outputs/pca_maf_robustness.csv", row.names = FALSE)
cat(sprintf("MAF >= 0.05 keeps %d SNPs; |r| of PC scores vs full data: PC1 %.2f  PC2 %.2f  PC3 %.2f\n",
            nLoc(gl_maf), rob$r[1], rob$r[2], rob$r[3]))


# -----------------------------------------------------------------------------
# 5. Genomic relationship matrix (identity by descent among genotypes)
# -----------------------------------------------------------------------------
# Additive relationship matrix of Endelman (2011) from rrBLUP; genotypes are
# recoded to -1/0/1 and missing calls mean-imputed, as A.mat() expects.
G <- rrBLUP::A.mat(M - 1, impute.method = "mean", return.imputed = FALSE)
dimnames(G) <- list(indNames(gl), indNames(gl))
saveRDS(G, "outputs/grm.rds")

# Figure 3: heatmap of the relationship matrix, annotated by sampling location.
# Colour scale centred on 0 (blue = less related than the sample average,
# red = more related) and clipped at +/- 4, as in the original figure.
HEAT_COL <- c(colorRampPalette(c("#0000FF", "#00FFFF"))(50), colorRampPalette(c("#FFFF00", "#FF0000"))(50))  # boundary at 0
HEAT_BRK <- seq(-4, 4, length.out = 101)
if (has("pheatmap")) {
  ann <- data.frame(Location = popf); rownames(ann) <- indNames(gl)
  pheatmap::pheatmap(pmin(pmax(G, -4), 4), annotation_row = ann, annotation_col = ann,
                     annotation_colors = list(Location = PAL9), color = HEAT_COL, breaks = HEAT_BRK,
                     border_color = NA, show_rownames = TRUE, show_colnames = TRUE, fontsize_row = 3, fontsize_col = 3,
                     fontfamily = FONT, fontsize = 8, main = "Probability of identity by descent",
                     filename = "figures/Figure3.png", width = 7.5, height = 6.7)
} else message("Figure 3 skipped: needs pheatmap")


# -----------------------------------------------------------------------------
# 6. Dendrogram of genetic distances between samples (Czekanowski)
# -----------------------------------------------------------------------------
D  <- gl.dist.ind(gl, method = "Czekanowski", plot.display = FALSE, verbose = 0)
hc <- hclust(D, method = "average")                   # UPGMA
png300("Figure6.png", w = 7.5, h = 4.5)
if (has("dendextend")) {
  dend <- dendextend::set(as.dendrogram(hc), "labels_cex", 0.35)
  dendextend::labels_colors(dend) <- PAL9[popf][hc$order]
  plot(dend, ylab = "Czekanowski distance")
} else plot(hc, cex = 0.35, main = "", xlab = "", sub = "", ylab = "Czekanowski distance")
dev.off()


# -----------------------------------------------------------------------------
# 7. Population structure with fastSTRUCTURE
# -----------------------------------------------------------------------------
# K = 1..10, ten replicate runs per K with different seeds; the optimal K is
# the one maximising the mean marginal likelihood across replicates.
# The dartRverse wrapper gl.run.faststructure() uses a single seed for all
# replicates (which would make them identical), so the binary is driven
# directly. Skipped if the fastStructure / plink executables are not found.
fs_ok <- file.exists(FASTSTRUCTURE) && file.exists(PLINK)
if (!fs_ok) {
  message("fastSTRUCTURE skipped: set FASTSTRUCTURE and PLINK (section 0) to the executables and re-run")
} else {
  fsdir <- "outputs/faststructure"; dir.create(fsdir, showWarnings = FALSE)
  gl2plink(gl, outfile = "gl_plink", outpath = fsdir, verbose = 0)   # .ped/.map
  system2(PLINK, c("--file", file.path(fsdir, "gl_plink"), "--make-bed",
                   "--allow-extra-chr", "--out", file.path(fsdir, "gl_plink")), stdout = FALSE)
  Ks <- 1:10; reps <- 1:10
  for (K in Ks) for (r in reps) if (!file.exists(file.path(fsdir, sprintf("rep_%d.%d.log", r, K))))   # skip runs already done
    system2(FASTSTRUCTURE, c(sprintf("-K %d", K),
                             sprintf("--input=%s",  file.path(fsdir, "gl_plink")),
                             sprintf("--output=%s", file.path(fsdir, sprintf("rep_%d", r))),
                             "--tol=1e-5", "--prior=simple", sprintf("--seed=%d", r)),
            stdout = FALSE, stderr = FALSE)

  # marginal likelihood of every run -> mean and SD per K -> choose K
  ml_runs <- do.call(rbind, lapply(Ks, function(K) data.frame(K = K, rep = reps,
    marginal_likelihood = sapply(reps, function(r) {
      L <- readLines(file.path(fsdir, sprintf("rep_%d.%d.log", r, K)))
      as.numeric(sub(".*= ", "", tail(grep("^Marginal Likelihood =", L, value = TRUE), 1))) }))))
  ml <- do.call(rbind, lapply(split(ml_runs, ml_runs$K), function(d)
    data.frame(K = d$K[1], mean_marginal_likelihood = mean(d$marginal_likelihood),
               sd = sd(d$marginal_likelihood), best_rep = d$rep[which.max(d$marginal_likelihood)])))
  write.csv(ml_runs, "outputs/faststructure_marginal_likelihood_runs.csv", row.names = FALSE)
  write.csv(ml,      "outputs/faststructure_marginal_likelihood.csv",      row.names = FALSE)
  print(ml)
  K_best <- ml$K[which.max(ml$mean_marginal_likelihood)]
  # fastSTRUCTURE's second criterion: number of components that together carry
  # 99.99% of the ancestry (chooseK.py "model components"), at each K
  used <- sapply(Ks, function(K) {
    Q <- as.matrix(read.table(file.path(fsdir, sprintf("rep_%d.%d.meanQ", ml$best_rep[ml$K == K], K))))
    m <- sort(colMeans(Q), decreasing = TRUE); which(cumsum(m) >= 0.9999)[1] })
  ml$components_used <- used
  write.csv(ml, "outputs/faststructure_marginal_likelihood.csv", row.names = FALSE)
  cat(sprintf("fastSTRUCTURE: K maximising the mean marginal likelihood = %d; components used at K >= %d: %d\n",
              K_best, K_best, used[K_best]))

  # Figure 4: mean marginal likelihood (+/- SD over replicates) against K.
  p4 <- ggplot(ml, aes(K, mean_marginal_likelihood)) + geom_line() + geom_point(size = 2) +
    geom_errorbar(aes(ymin = mean_marginal_likelihood - sd, ymax = mean_marginal_likelihood + sd), width = 0.15) +
    scale_x_continuous(breaks = Ks) +
    labs(x = "Number of clusters (K)", y = "Mean marginal likelihood") + theme_bw()
  ggsave("figures/Figure4.png", p4, width = 6.5, height = 4.5, dpi = 300)

  # Figure 5: ancestry proportions for K = 2 to 6. Replicate runs of each K
  # are aligned with CLUMPP and grouped into modes with the CLUMPAK method
  # (Kopelman et al. 2015) by dartR.popgen::gl.plot.faststructure(); each row
  # is one K.mode (e.g. 4.1, 4.2). Individuals are ordered as in the
  # Czekanowski dendrogram of Figure 6, so related genotypes sit together.
  Kfig <- 2:6
  q_list <- lapply(Kfig, function(K) lapply(reps, function(r) {
    Q <- read.table(file.path(fsdir, sprintf("rep_%d.%d.meanQ", r, K))); colnames(Q) <- paste0("cluster", seq_len(K))
    data.frame(id = indNames(gl), orig.pop = as.character(pop(gl)), Q) }))
  names(q_list) <- Kfig; q_list <- lapply(q_list, function(y) { names(y) <- reps; y })
  pdf(NULL)                                            # the function draws; we only need its output
  modes <- gl.plot.faststructure(list(q_list = q_list), k.range = Kfig, den = FALSE, ind_name = FALSE)
  dev.off()
  if (has("reshape2")) {
    d <- do.call(rbind, lapply(modes, function(m) {
      long <- reshape2::melt(m[, c("Label", "K", "ord", grep("^cluster", names(m), value = TRUE))],
                             id.vars = c("Label", "K", "ord"), variable.name = "cluster", value.name = "Q")
      long }))
    d$Label <- factor(d$Label, levels = indNames(gl)[hc$order])   # same order as the Fig 6 dendrogram
    d$K <- factor(d$K, levels = unique(sapply(modes, function(m) m$K[1])))
    p5 <- ggplot(d, aes(Label, Q, fill = cluster)) + geom_col(width = 1, colour = "black", linewidth = 0.05) +
      facet_grid(K ~ ., switch = "y") + scale_fill_manual(values = unname(PAL9)[c(4, 2, 3, 5, 6, 1)]) +
      scale_y_continuous(expand = c(0, 0)) + labs(x = NULL, y = NULL) +
      theme_minimal(base_family = FONT, base_size = 8) +
      theme(legend.position = "none", panel.grid = element_blank(), panel.spacing = unit(0.5, "mm"),
            axis.text.y = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 3.5),
            strip.text.y.left = element_text(angle = 0, size = 8), strip.placement = "outside")
    ggsave("figures/Figure5.png", p5, width = 7.5, height = 6, dpi = 300)
    write.csv(d, "outputs/faststructure_modes_K2-6.csv", row.names = FALSE)
  } else message("Figure 5 skipped: needs reshape2")
}


# -----------------------------------------------------------------------------
# 8. AMOVA  (how genetic variation partitions among vs within locations)
# -----------------------------------------------------------------------------
gi <- gl2gi(gl, verbose = 0)                 # convert genlight -> genind for poppr
strata(gi) <- data.frame(pop = pop(gi))
amova <- poppr.amova(gi, ~pop, nperm = 9999, method = "ade4")
print(amova$results)
print(amova$componentsofcovariance)
capture.output(amova$results, amova$componentsofcovariance, amova$statphi, file = "outputs/amova.txt")

# Sampling location is a sampling stratum, not a demonstrated biological
# population, and several locations hold genotypes from more than one genetic
# cluster. So, when fastSTRUCTURE ran (section 7), tabulate location against
# the best-K cluster (S3 Table) and repeat the AMOVA with cluster as the
# grouping, to see how much of the within-location variance is the
# co-occurrence of distinct clusters at one site.
if (file.exists("outputs/faststructure_Q_best.csv")) {
  Q  <- read.csv("outputs/faststructure_Q_best.csv")
  comp <- apply(Q[, grep("^V", names(Q))], 1, which.max)            # component of max membership
  size_rank <- rank(-table(comp), ties.method = "first")             # renumber clusters by size (1 = largest)
  Q$cluster <- paste0("C", size_rank[as.character(comp)])
  write.csv(Q[, c("ind.name", "location", "cluster")], "outputs/cluster_assignment.csv", row.names = FALSE)
  # Diversity within each cluster: with the between-cluster (Wahlund) component
  # removed, FIS shows how much of the homozygote excess is inbreeding.
  cl_ind <- Q$cluster[match(indNames(gl), Q$ind.name)]
  div_cl <- do.call(rbind, lapply(c(sort(unique(cl_ind)), "all"), function(k) {
    Mk <- if (k == "all") M else M[cl_ind == k, , drop = FALSE]
    afk <- colMeans(Mk, na.rm = TRUE) / 2; Hek <- mean(2 * afk * (1 - afk), na.rm = TRUE); Hok <- mean(colMeans(Mk == 1, na.rm = TRUE), na.rm = TRUE)
    data.frame(cluster = k, n = nrow(Mk), Ho = round(Hok, 4), He = round(Hek, 4), FIS = round(1 - Hok / Hek, 3)) }))
  print(div_cl); write.csv(div_cl, "outputs/diversity_by_cluster.csv", row.names = FALSE)
  He_cl <- setNames(div_cl$He[div_cl$cluster != "all"], div_cl$cluster[div_cl$cluster != "all"])
  F_ind <- 1 - rowMeans(M == 1, na.rm = TRUE) / He_cl[cl_ind]           # individual F relative to own cluster
  write.csv(data.frame(ind.name = indNames(gl), cluster = cl_ind, F_ind = round(F_ind, 3)), "outputs/inbreeding_by_individual.csv", row.names = FALSE)
  cat(sprintf("Individual F vs own cluster: median %.2f; below 0.5: %s\n", median(F_ind), paste(indNames(gl)[F_ind < 0.5], collapse = ", ")))
  xt <- table(location = Q$location, cluster = Q$cluster)
  print(xt); write.csv(as.data.frame.matrix(xt), "outputs/location_by_cluster.csv")
  strata(gi) <- data.frame(cluster = Q$cluster[match(indNames(gl), Q$ind.name)])
  amova_cl <- poppr.amova(gi, ~cluster, nperm = 9999, method = "ade4")
  print(amova_cl$componentsofcovariance); print(amova_cl$statphi)
  capture.output(amova_cl$results, amova_cl$componentsofcovariance, amova_cl$statphi,
                 file = "outputs/amova_by_cluster.txt")
}


# -----------------------------------------------------------------------------
# 9. Core subset of the most informative genotypes
# -----------------------------------------------------------------------------
# Score every genotype on four diversity metrics, shortlist the top 20 of each
# (43 unique genotypes), drop the member with more missing data from every
# pair with > 98% identical calls (8 genotypes), and keep the rest (35).
sel  <- top_ind(gl, MAC_min = 5, top.ind = 20, rep_threshold = 0.98)
tab  <- sel$table
core <- tab$ind.name[tab$select]
write.csv(tab,            "outputs/core_selection_table.csv",    row.names = FALSE)  # Table 4
write.csv(sel$replicates, "outputs/core_replicate_pairs.csv",    row.names = FALSE)  # Table 5
writeLines(core,          "outputs/core_genotypes.txt")
cat(sprintf("Core selection: %d shortlisted, %d dropped as replicates, %d selected\n",
            sum(tab$shortlisted), sum(tab$replicate_dropped), length(core)))

# Validation: polymorphic loci retained by the core vs 1,000 random subsets of
# the same size.
full_poly <- n_poly(gl)
core_poly <- n_poly(gl.keep.ind(gl, ind.list = core, verbose = 0))
rand_poly <- replicate(1000,
  n_poly(gl.keep.ind(gl, ind.list = sample(indNames(gl), length(core)), verbose = 0)))
cat(sprintf("Core retains %.1f%% of polymorphic loci; beats %.1f%% of random subsets (random mean %.1f%%)\n",
            100 * core_poly / full_poly, 100 * mean(core_poly > rand_poly), 100 * mean(rand_poly) / full_poly))

# Figure 7: distribution of pairwise % identical genotypes, all pairs (top) and
# pairs among the shortlisted genotypes (bottom); the 98% threshold is marked.
PI <- pairwise_identity(M)
short <- tab$ind.name[tab$shortlisted]
png300("Figure7.png", w = 6.5, h = 7)
par(mfrow = c(2, 1), mar = c(4.5, 4.5, 2, 1))
hist(PI[upper.tri(PI)], breaks = 50, col = "grey80", border = "white", main = "(A) All pairs",
     xlab = "Identical genotype calls (%)"); abline(v = 98, lty = 2, col = "red")
PIs <- PI[short, short]
hist(PIs[upper.tri(PIs)], breaks = 50, col = "grey80", border = "white", main = "(B) Shortlisted genotypes",
     xlab = "Identical genotype calls (%)"); abline(v = 98, lty = 2, col = "red")
dev.off()

# Figure 8: accumulated allelic richness and number of polymorphic loci as the
# selected genotypes are added one at a time, from the top-ranked downwards.
core_ord <- core                              # already ordered by mean rank (see top_ind)
steps <- 2:length(core_ord)
acc <- t(sapply(steps, function(k) {
  g  <- gl.filter.allna(gl.keep.ind(gl, ind.list = core_ord[1:k], verbose = 0), verbose = 0)
  ar <- if (has("hierfstat"))
    mean(hierfstat::allelic.richness(hierfstat::genind2hierfstat(gl2gi(g, verbose = 0)))$Ar[, 1], na.rm = TRUE)
  else NA
  c(n = k, allelic_richness = ar, polymorphic_loci = n_poly(g)) }))
write.csv(acc, "outputs/core_accumulation.csv", row.names = FALSE)
png300("Figure8.png", w = 6, h = 7.5)
par(mfrow = c(2, 1), mar = c(4.5, 4.5, 1, 1))
plot(acc[, "n"], acc[, "allelic_richness"], type = "b", pch = 19, col = "deeppink",
     xlab = "Number of genotypes", ylab = "Allelic richness", main = "(A)")
plot(acc[, "n"], acc[, "polymorphic_loci"], type = "b", pch = 19, col = "deeppink",
     xlab = "Number of genotypes", ylab = "Polymorphic loci", main = "(B)")
dev.off()

# Figure 9: PCA highlighting the selected (blue) and non-selected (green) genotypes.
is_core <- indNames(gl) %in% core
png300("Figure9.png", w = 7.5, h = 5.6)
par(mar = c(5, 5, 2, 8), xpd = TRUE)
pc_plot(1, 2, col = ifelse(is_core, "#4477AA", "#228833"), pch = ifelse(is_core, 19, 1))
legend("topright", inset = c(-0.28, 0), legend = c("Selected", "Not selected"),
       col = c("#4477AA", "#228833"), pch = c(19, 1), bty = "n")
dev.off()

# Figure 10: dendrogram + heatmap of pairwise relatedness (section 5 matrix),
# annotated by selected / not selected. (The manuscript's kinship estimates
# come from the external program EMIBD9; the additive relationship matrix is
# used here so the figure is reproducible from R alone.)
if (has("pheatmap")) {
  ann <- data.frame(Subset = factor(ifelse(is_core, "Selected", "Not selected"))); rownames(ann) <- indNames(gl)
  pheatmap::pheatmap(pmin(pmax(G, -4), 4), annotation_row = ann, annotation_col = ann,
                     annotation_colors = list(Subset = c(Selected = "#EE8866", "Not selected" = "#4477AA")),
                     color = HEAT_COL, breaks = HEAT_BRK, border_color = NA,
                     show_rownames = TRUE, show_colnames = TRUE, fontsize_row = 3, fontsize_col = 3,
                     fontfamily = FONT, fontsize = 8, main = "Pairwise relatedness of selected and non-selected genotypes",
                     filename = "figures/Figure10.png", width = 7.5, height = 6.7)
} else message("Figure 10 skipped: needs pheatmap")

cat("\nDone. Tables written to outputs/, figures to figures/.\n")
