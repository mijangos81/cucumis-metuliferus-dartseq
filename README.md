# Genetic diversity of *Cucumis metuliferus* (DArT-seq)

Reproducible analysis of DArT-seq SNP genotyping of the African horned cucumber
(*Cucumis metuliferus*), accompanying:

> Koopa *et al.* (2026). *DArT-based characterisation of genetic diversity in
> genotypes of Cucumis metuliferus collected from SACU.* PLOS ONE.

All analyses are in a **single, commented R script** (`analysis.R`) and run on the
[dartRverse](https://github.com/green-striped-gecko/dartRverse) suite. The script
reproduces every number and figure in the manuscript from the raw DArT report.

## Contents

```
├── analysis.R    # the complete analysis, top to bottom (read this)
├── data/
│   ├── Report_DCu24-9392_SNP_mapping_2.csv     # raw DArT SNP report
│   └── metadata_with_decimal_coordinates.csv   # sample id, location, coordinates
├── outputs/      # created when you run the script (tables, fastSTRUCTURE runs)
├── figures/      # created when you run the script (Figures 1-10, 300 dpi PNG)
└── README.md
```

## What the script does

1. Read the DArT SNP report + metadata into a `genlight` object; populations =
   the 9 geographic **sampling locations** (Figure 1)
2. Quality-filter the SNPs (call rate ≥ 80%, read depth 5–1000, reproducibility ≥ 99%)
   → **120 individuals × 5,875 SNPs**, and descriptive statistics of the retained
   SNPs (call rate, missingness, read depth, MAF, PIC)
3. Diversity indices (Ho, He, FIS, PIC), overall and per location
4. PCA (Figure 2), plus a check that the PCA is unchanged after filtering to
   MAF ≥ 0.05 (correlation of PC scores)
5. Genomic relationship matrix, rrBLUP `A.mat()` (Figure 3)
6. Dendrogram of Czekanowski distances (Figure 6)
7. Population structure with fastSTRUCTURE, K = 1–10 with ten seeded replicates
   per K; optimal K by mean marginal likelihood (**K = 3**; Figures 4 and 5)
8. AMOVA among / within locations, 9,999 permutations
9. Core subset: top 20 genotypes on each of four diversity metrics (43 unique),
   minus 8 near-identical replicates (> 98% identical calls) → **35 genotypes**;
   validated against 1,000 random subsets of the same size (Figures 7–10)

## Requirements

- **R ≥ 4.2** with `dartRverse`, `poppr`, `rrBLUP` and `ggplot2`:

  ```r
  install.packages(c("dartRverse", "poppr", "rrBLUP", "ggplot2"))
  library(dartRverse); dartRverse_install()   # installs the dartR sub-packages
  ```

- Optional, for individual figures (skipped with a message if absent):
  `pheatmap`, `dendextend`, `reshape2`, `patchwork`, `hierfstat`, `sf`,
  `rnaturalearth`, `ggspatial`.

- **fastSTRUCTURE** and **plink** (step 7 only) — external programs, not R
  packages: <https://rajanil.github.io/fastStructure/>,
  <https://www.cog-genomics.org/plink/>. Set the paths at the top of
  `analysis.R` (or the `FASTSTRUCTURE` / `PLINK` environment variables). If
  they are not found, step 7 is skipped automatically and every other step
  still runs. The wrapper `gl.run.faststructure()` is not used because it
  gives every replicate the same seed.

## Running it

```r
setwd("path/to/cucumis-metuliferus-dartseq")   # this folder
source("analysis.R")
```

A full run takes about 40 minutes on a laptop, almost all of it the 100
fastSTRUCTURE runs. Tables (CSV) go to `outputs/`, figures to `figures/`, and the
key numbers (dataset size, diversity indices, PCA variance, optimal K, AMOVA,
core-subset retention) are printed to the console.

## Notes

- **No MAF/LD filtering** is applied on purpose — rare alleles are informative
  for the diversity and core-subset aims. Step 4 shows that removing them does
  not change the PCA structure (2,136 SNPs at MAF ≥ 0.05; PC1–3 score
  correlations ≥ 0.99).
- The CLUMPAK grouping of fastSTRUCTURE replicate runs into modes is done on the
  CLUMPAK web server with the `outputs/faststructure/rep_*.<K>.meanQ` files; it
  is not reproduced in R.
- The manuscript's kinship estimates for Figure 10 come from the external
  program EMIBD9; the script draws the same figure from the additive
  relationship matrix of step 5 so it is reproducible from R alone.
- Coordinates in the metadata are decimal degrees, converted from the original
  field GPS readings.
- `set.seed(1)` at the top makes the random steps (core-vs-random comparison)
  reproducible.

## Citation

Please cite the manuscript above and the dartRverse suite (Mijangos *et al.*,
*Methods in Ecology and Evolution*; <https://github.com/green-striped-gecko/dartRverse>).
