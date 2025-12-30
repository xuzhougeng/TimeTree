#!/usr/bin/env Rscript
# Plot phylogenetic tree (rooted or unrooted) using ggtree
# For visualizing intermediate trees (species tree, consensus tree, etc.)

suppressPackageStartupMessages({
  library(ggtree)
  library(treeio)
  library(ggplot2)
  library(ape)
})

# Get snakemake variables
input_tree <- snakemake@input[["tree"]]
output_pdf <- snakemake@output[["pdf"]]
output_png <- snakemake@output[["png"]]
output_svg <- snakemake@output[["svg"]]
title <- snakemake@params[["title"]]
show_support <- snakemake@params[["show_support"]]

if (is.null(title) || title == "") {
  title <- "Phylogenetic Tree"
}

if (is.null(show_support)) {
  show_support <- TRUE
}

cat("Reading tree from:", input_tree, "\n")

# Read tree
tree <- tryCatch({
  read.tree(input_tree)
}, error = function(e) {
  cat("read.tree failed, trying read.nexus...\n")
  read.nexus(input_tree)
})

cat("Tree parsed successfully\n")
cat("Number of tips:", length(tree$tip.label), "\n")
cat("Is rooted:", is.rooted(tree), "\n")

n_tips <- length(tree$tip.label)

# Check if tree has node labels (bootstrap/support values)
has_support <- !is.null(tree$node.label) && any(!is.na(tree$node.label) & tree$node.label != "")

# Base ggtree plot
if (is.rooted(tree)) {
  p <- ggtree(tree, ladderize = TRUE, right = TRUE)
} else {
  # For unrooted trees, use daylight layout
  p <- ggtree(tree, layout = "daylight")
}

# Add tip labels
p <- p + geom_tiplab(
  fontface = "italic",
  size = 3.5,
  offset = 0.01
)

# Add support values if available and requested
if (has_support && show_support) {
  cat("Adding support values\n")
  p <- p + geom_nodelab(
    aes(label = label),
    size = 2.5,
    hjust = -0.2,
    vjust = -0.5,
    color = "darkgreen"
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
  fontsize = 3,
  linesize = 0.5,
  offset = 0.3
)

# Calculate plot dimensions
max_label_length <- max(nchar(tree$tip.label))
right_margin_mm <- max_label_length * 2.5 + 10

p <- p + theme_tree2() +
  coord_cartesian(clip = "off") +
  theme(
    plot.margin = margin(10, right_margin_mm, 10, 10, unit = "mm"),
    axis.text.x = element_text(size = 10)
  )

# Add title
p <- p + ggtitle(title) +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

# Calculate plot dimensions based on number of tips
plot_height <- max(5, n_tips * 0.5)
plot_width <- max(8, 6 + max_label_length * 0.1)

cat("Plot dimensions:", plot_width, "x", plot_height, "\n")

# Create output directory if needed
dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)

# Save outputs
cat("Saving PDF:", output_pdf, "\n")
ggsave(output_pdf, p, width = plot_width, height = plot_height, units = "in", dpi = 300)

cat("Saving PNG:", output_png, "\n")
ggsave(output_png, p, width = plot_width, height = plot_height, units = "in", dpi = 300)

cat("Saving SVG:", output_svg, "\n")
ggsave(output_svg, p, width = plot_width, height = plot_height, units = "in")

cat("Done!\n")
