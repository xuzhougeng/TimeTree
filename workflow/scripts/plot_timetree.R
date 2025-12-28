#!/usr/bin/env Rscript
# Plot time-calibrated phylogenetic tree using ggtree
# Visualize divergence times with 95% HPD confidence intervals

suppressPackageStartupMessages({
  library(ggtree)
  library(treeio)
  library(ggplot2)
  library(ape)
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

cat("Reading tree from:", input_tree, "\n")

# Read the file content
file_content <- readLines(input_tree)
cat("File content:\n")
cat(file_content, sep = "\n")
cat("\n")

# Find the tree line and fix UTREE -> TREE for ape compatibility
tree_line_idx <- grep("UTREE|TREE", file_content, ignore.case = TRUE)
if (length(tree_line_idx) > 0) {
  # Replace UTREE with TREE
  file_content <- gsub("UTREE", "TREE", file_content, ignore.case = TRUE)
}

# Write to temp file
temp_file <- tempfile(fileext = ".tre")
writeLines(file_content, temp_file)

# Extract HPD annotations before removing them
tree_text <- paste(file_content, collapse = "\n")

# Parse HPD annotations: [&95%HPD={min, max}]
hpd_pattern <- "\\[&95%HPD=\\{([0-9.]+),\\s*([0-9.]+)\\}\\]"
hpd_matches <- gregexpr(hpd_pattern, tree_text, perl = TRUE)
hpd_strings <- regmatches(tree_text, hpd_matches)[[1]]

hpd_data <- data.frame(
  hpd_min = numeric(0),
  hpd_max = numeric(0)
)

if (length(hpd_strings) > 0) {
  for (h in hpd_strings) {
    nums <- as.numeric(regmatches(h, gregexpr("[0-9.]+", h))[[1]])
    if (length(nums) >= 2) {
      hpd_data <- rbind(hpd_data, data.frame(hpd_min = nums[1], hpd_max = nums[2]))
    }
  }
}

cat("Found", nrow(hpd_data), "HPD annotations\n")

# Remove HPD annotations for ape parsing
clean_content <- gsub("\\s*\\[&95%HPD=\\{[^}]+\\}\\]", "", file_content)
clean_file <- tempfile(fileext = ".tre")
writeLines(clean_content, clean_file)

# Read tree with ape
tree <- tryCatch({
  read.nexus(clean_file)
}, error = function(e) {
  cat("read.nexus failed, trying read.tree...\n")
  # Extract just the newick string
  newick_line <- grep("^\\s*\\(", file_content, value = TRUE)
  if (length(newick_line) == 0) {
    # Try to extract from TREE line
    tree_line <- grep("TREE", clean_content, value = TRUE, ignore.case = TRUE)
    if (length(tree_line) > 0) {
      newick_str <- sub(".*=\\s*", "", tree_line[1])
      newick_str <- trimws(newick_str)
      if (!grepl(";$", newick_str)) newick_str <- paste0(newick_str, ";")
      read.tree(text = newick_str)
    } else {
      stop("Could not parse tree file")
    }
  } else {
    read.tree(text = newick_line[1])
  }
})

cat("Tree parsed successfully\n")
cat("Number of tips:", length(tree$tip.label), "\n")
cat("Tip labels:", paste(tree$tip.label, collapse = ", "), "\n")

# Get node information
n_tips <- length(tree$tip.label)
n_nodes <- tree$Nnode

# Calculate node ages from branch lengths (tree is ultrametric, time from root)
# For ultrametric trees, we need to calculate node heights
node_depths <- node.depth.edgelength(tree)
max_depth <- max(node_depths)
node_ages <- max_depth - node_depths

cat("Root age (max depth):", max_depth, "\n")

# Create a data frame for node annotations
node_info <- data.frame(
  node = 1:(n_tips + n_nodes),
  age = node_ages
)

# Add HPD data to internal nodes if available
internal_nodes <- (n_tips + 1):(n_tips + n_nodes)
if (nrow(hpd_data) > 0 && nrow(hpd_data) == length(internal_nodes)) {
  node_info$hpd_min <- NA
  node_info$hpd_max <- NA
  node_info$hpd_min[internal_nodes] <- hpd_data$hpd_min
  node_info$hpd_max[internal_nodes] <- hpd_data$hpd_max
}

# Base ggtree plot - use branch.length for time scale
p <- ggtree(tree, ladderize = TRUE, right = TRUE)

# Reverse x-axis so time goes from past (left) to present (right=0)
p <- p + scale_x_continuous(
  labels = function(x) round(max_depth - x, 1),
  name = paste0("Time (", time_unit, ")")
) + theme_tree2()

# Add tip labels (species names in italics)
p <- p + geom_tiplab(
  fontface = "italic",
  size = 3.5,
  offset = max_depth * 0.02
)

# Add node age labels at internal nodes
p <- p %<+% node_info

# Add age labels for internal nodes
p <- p + geom_nodelab(
  aes(x = x, label = ifelse(!isTip, round(max_depth - x, 1), "")),
  size = 2.8,
  hjust = -0.1,
  vjust = -0.5,
  color = "darkblue"
)

# Add HPD bars if available
if ("hpd_min" %in% names(node_info) && any(!is.na(node_info$hpd_min))) {
  cat("Adding HPD bars\n")

  # Get internal node positions - convert to regular data.frame
  tree_data <- as.data.frame(p$data)

  # Merge HPD data
  hpd_subset <- node_info[, c("node", "hpd_min", "hpd_max")]
  tree_data <- merge(tree_data, hpd_subset, by = "node", all.x = TRUE)

  # Add HPD error bars for internal nodes
  # Convert HPD (ages) to x-coordinates
  hpd_plot_data <- tree_data[!is.na(tree_data$hpd_min), ]
  if (nrow(hpd_plot_data) > 0) {
    hpd_plot_data$x_min <- max_depth - hpd_plot_data$hpd_max  # older age = smaller x
    hpd_plot_data$x_max <- max_depth - hpd_plot_data$hpd_min  # younger age = larger x

    p <- p + geom_errorbarh(
      data = hpd_plot_data,
      aes(xmin = x_min, xmax = x_max, y = y),
      height = 0.3,
      color = "steelblue",
      linewidth = 1.5,
      alpha = 0.7
    )
  }
}

# Add node points for internal nodes
p <- p + geom_nodepoint(
  size = 2,
  color = "darkred",
  alpha = 0.8
)

# Add scale bar
p <- p + geom_treescale(
  x = 0,
  y = -0.5,
  fontsize = 3,
  linesize = 0.5,
  offset = 0.3
)

# Calculate plot dimensions
max_label_length <- max(nchar(tree$tip.label))

# Use coord_cartesian with clip="off" to allow labels to extend beyond plot area
# This avoids expanding the x-axis unnecessarily
p <- p + coord_cartesian(clip = "off")

# Calculate right margin based on label length (in mm)
right_margin_mm <- max_label_length * 3 + 10

p <- p + theme(
  plot.margin = margin(10, right_margin_mm, 10, 10, unit = "mm"),
  axis.text.x = element_text(size = 10),
  axis.title.x = element_text(size = 12)
)

# Add title
p <- p + ggtitle("Time-Calibrated Phylogenetic Tree") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

# Calculate plot dimensions based on number of tips and label length
plot_height <- max(5, n_tips * 0.7)
plot_width <- max(10, 8 + max_label_length * 0.12)

cat("Plot dimensions:", plot_width, "x", plot_height, "\n")

# Save outputs
cat("Saving PDF:", output_pdf, "\n")
ggsave(output_pdf, p, width = plot_width, height = plot_height, units = "in", dpi = 300)

cat("Saving PNG:", output_png, "\n")
ggsave(output_png, p, width = plot_width, height = plot_height, units = "in", dpi = 300)

cat("Saving SVG:", output_svg, "\n")
ggsave(output_svg, p, width = plot_width, height = plot_height, units = "in")

# Clean up temp files
unlink(temp_file)
unlink(clean_file)

cat("Done!\n")
