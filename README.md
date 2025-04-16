# mutualinfobio

Tools for calculating mutual information and covariation scores from multiple sequence alignments, with heatmap visualization and support for FASTA or plain text input.

## Installation

You can install the development version of `mutualinfobio` from GitHub:

```r
install.packages("remotes")
remotes::install_github("TraxlerLab/mutualinfobio")

library(mutualinfobio)

# Load example alignment (FASTA file)
path <- system.file("extdata", "cell3925mmc4.fasta", package = "mutualinfobio")
aln <- Alignment(alignment_path = path)


# Compute mutual information
mi_scores <- getMICovarianceScore(aln)

# Visualize
draw_MI_heatmap(mi_scores, alignment_labels = c("GeneA", "GeneB"), alignment_label_locations = c(1, 50), alignment_domain_residues = c(50, 100), heatmap_lines = TRUE)
```
