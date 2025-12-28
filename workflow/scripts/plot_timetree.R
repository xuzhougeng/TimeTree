#!/usr/bin/env Rscript
# Plot time-calibrated phylogenetic tree using ggtree
# Visualize divergence times with 95% HPD confidence intervals

suppressPackageStartupMessages({
  library(ggtree)
  library(treeio)
  library(ggplot2)
  library(dplyr)
  library(scales)
})

# Get snakemake variables
input_tree <- snakemake@input[["figtree"]]
output_pdf <- snakemake@output[["pdf"]]
output_png <- snakemake@output[["png"]]
output_svg <- snakemake@output[["svg"]]
time_unit <- snakemake@params[["time_unit"]]

if (is.null(time_unit) || time_unit == "") {
  time_unit <- "Ma"
}

# Read the MCMCtree FigTree output
cat("Reading tree from:", input_tree, "\n")
tree <- read.beast(input_tree)

# Extract tree data
tree_data <- as_tibble(tree)

# Get node ages and HPD intervals
# MCMCtree annotates nodes with 'age', '95%HPD' or 'height_0.95_HPD'
node_data <- tree_data %>%
  filter(!is.na(node)) %>%
  mutate(
    # Try different column names for node age
    node_age = case_when(
      !is.na(age) ~ as.numeric(age),
      !is.na(height) ~ as.numeric(height),
      TRUE ~ NA_real_
    )
  )

# Check for HPD columns (different MCMCtree versions use different names)
hpd_cols <- grep("HPD|CI", names(tree_data), value = TRUE, ignore.case = TRUE)
cat("Found HPD columns:", paste(hpd_cols, collapse = ", "), "\n")

# Base ggtree plot
p <- ggtree(tree, ladderize = TRUE, right = TRUE) +
  # Time axis (reversed so present is at 0)
  theme_tree2() +
  scale_x_continuous(
    labels = function(x) abs(x),
    name = paste0("Time (", time_unit, ")")
  )

# Add tip labels (species names in italics for Latin names)
p <- p + geom_tiplab(
  fontface = "italic",
  size = 3.5,
  align = TRUE,
  linesize = 0.3,
  offset = 0.5
)

# Add node age labels at internal nodes
# Get internal node data
internal_nodes <- tree_data %>%
  filter(is.na(label) | label == "" | !isTip)

# Check if we have age/height information
if ("age" %in% names(tree_data) || "height" %in% names(tree_data)) {
  age_col <- if ("age" %in% names(tree_data)) "age" else "height"

  # Add node age labels
  p <- p + geom_nodelab(
    aes(label = round(!!sym(age_col), 1)),
    size = 2.5,
    hjust = -0.1,
    vjust = -0.5,
    color = "darkblue"
  )
}

# Add HPD bars if available
# MCMCtree typically outputs height_0.95_HPD or 95%HPD
hpd_95_col <- grep("95.*HPD|HPD.*95|height_0.95_HPD", names(tree_data),
                   value = TRUE, ignore.case = TRUE)[1]

if (!is.na(hpd_95_col) && hpd_95_col %in% names(tree_data)) {
  cat("Using HPD column:", hpd_95_col, "\n")

  # Extract HPD range
  p <- p + geom_range(
    range = hpd_95_col,
    color = "steelblue",
    alpha = 0.6,
    size = 2,
    center = "age"
  )
}

# Add node points
p <- p + geom_nodepoint(
  size = 1.5,
  color = "darkred",
  alpha = 0.7
)

# Add scale bar
p <- p + geom_treescale(
  x = 0, y = -1,
  fontsize = 3,
  linesize = 0.5,
  offset = 0.5
)

# Adjust plot margins for tip labels
max_label_length <- max(nchar(as.character(tree_data$label[!is.na(tree_data$label)])), na.rm = TRUE)
right_margin <- max_label_length * 0.15 + 2

p <- p + theme(
  plot.margin = margin(10, right_margin, 10, 10, unit = "mm"),
  axis.text.x = element_text(size = 10),
  axis.title.x = element_text(size = 12)
)

# Add title
p <- p + ggtitle("Time-Calibrated Phylogenetic Tree") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

# Calculate plot dimensions based on number of tips
n_tips <- length(tree@phylo$tip.label)
plot_height <- max(6, n_tips * 0.4)
plot_width <- max(10, 8 + right_margin * 0.3)

cat("Number of tips:", n_tips, "\n")
cat("Plot dimensions:", plot_width, "x", plot_height, "\n")

# Save outputs
cat("Saving PDF:", output_pdf, "\n")
ggsave(output_pdf, p, width = plot_width, height = plot_height, units = "in", dpi = 300)

cat("Saving PNG:", output_png, "\n")
ggsave(output_png, p, width = plot_width, height = plot_height, units = "in", dpi = 300)

cat("Saving SVG:", output_svg, "\n")
ggsave(output_svg, p, width = plot_width, height = plot_height, units = "in")

cat("Done!\n")
