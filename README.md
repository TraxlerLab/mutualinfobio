![](assets/mibio.png)

## mutualinfobio

**mutualinfobio** is an R package for computing mutual information (MI) and coevolutionary scores from protein multiple sequence alignments (MSAs), with built-in support for FASTA and plain-text input formats. The package includes visualization tools such as MI heatmaps.

------------------------------------------------------------------------

## 🔧 Installation

To install the development version directly from GitHub:

``` r
install.packages("remotes")  # if not already installed
remotes::install_github("TraxlerLab/mutualinfobio")
```

------------------------------------------------------------------------

## 📦 Example usage

``` r
library(mutualinfobio)

# Load an example alignment from the Cell 2008 paper by Weigt et al.
path <- system.file("extdata", "cell3925mmc4.fasta", package = "mutualinfobio")

# Construct an Alignment object
aln <- Alignment(alignment_path = path)

# Compute mutual information-based covariance scores
mi_scores <- getMICovarianceScore(aln, cores = 8)

# Visualize as a heatmap
# Visualize as a heatmap
draw_MI_heatmap(
  data = mi_scores,
  alignment_labels = c("HK", "RR"),
  alignment_label_locations = c(50, 200),
  alignment_domain_residues = c(100, 308),
  heatmap_lines = FALSE
)
```

> ✨ The file `cell3925mmc4.fasta` is included as example data under a [Creative Commons Attribution (CC BY)](https://creativecommons.org/licenses/by/4.0/) license, from the supplementary material of Wkerker det al. (2008), *Cell*. DOI: [10.1016/j.cell.2008.04.040](https://doi.org/10.1016/j.cell.2008.04.040)

------------------------------------------------------------------------

## 🧬 Features

-   Compute MI-based covariation from aligned sequences
-   Heatmap visualizations of MI scores
-   Support for FASTA and plain text alignments
-   Parallel computation support in UNIX systems(with CRAN-safe defaults)
-   Cross-Platform (Parallel computation not supported in Windows)

------------------------------------------------------------------------

## 📖 Citation

This work was based on the methods described by Skerker et. al:

> Skerker, J.M., Perchuk, B.S., Siryaporn, A., Lubin, E.A., Ashenberg, O., Goulian, M., & Laub, M.T. (2008).\
> Rewiring the specificity of two-component signal transduction systems. *Cell*, 133(6), 1043–1054.\
> <https://doi.org/10.1016/j.cell.2008.04.040>

------------------------------------------------------------------------

## 📄 License

This package is released under the MIT License. See `LICENSE` file for details.

The included alignment file (`cell3925mmc4.fasta`) is licensed separately under CC BY and is used with attribution.

------------------------------------------------------------------------
