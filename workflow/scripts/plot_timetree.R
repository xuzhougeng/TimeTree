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

# ============================================================================
# Parse tree with HPD annotations using recursive descent parser
# ============================================================================

# Extract the Newick string with annotations
tree_text <- paste(file_content, collapse = "\n")

# Find the actual tree definition line (contains "= (" pattern)
# This avoids matching "BEGIN TREES;"
tree_line <- grep("=\\s*\\(", file_content, value = TRUE)
cat("Found", length(tree_line), "tree definition line(s)\n")

if (length(tree_line) > 0) {
  # Extract everything after the first "= " followed by "("
  # Use a specific pattern to match "TREE name = " prefix
  newick_with_hpd <- sub("^[^=]+=\\s*", "", tree_line[1])
  newick_with_hpd <- trimws(newick_with_hpd)
  cat("Extracted Newick string (first 200 chars):\n")
  cat(substr(newick_with_hpd, 1, 200), "...\n")
} else {
  # Fallback: try to find any line starting with parenthesis
  newick_line <- grep("^\\s*\\(", file_content, value = TRUE)
  if (length(newick_line) > 0) {
    newick_with_hpd <- trimws(newick_line[1])
  } else {
    stop("Could not find tree definition in file")
  }
}

cat("Parsing tree with HPD annotations...\n")

# Function to parse Newick with HPD and track node HPDs
parse_tree_with_hpd <- function(newick_str) {
  # First, get clean newick for parsing
  clean_newick <- gsub("\\s*\\[&95%HPD=\\{[^}]+\\}\\]", "", newick_str)

  # Parse tree to get structure
  tree <- read.tree(text = clean_newick)

  n_tips <- length(tree$tip.label)
  n_internal <- tree$Nnode

  # Find all HPD annotations with their positions
  hpd_pattern <- "\\[&95%HPD=\\{([0-9.]+),\\s*([0-9.]+)\\}\\]"
  matches <- gregexpr(hpd_pattern, newick_str, perl = TRUE)

  if (matches[[1]][1] == -1) {
    cat("No HPD annotations found in tree\n")
    return(list(tree = tree, hpd_data = NULL))
  }

  # Extract HPD values
  match_starts <- matches[[1]]
  match_lengths <- attr(matches[[1]], "match.length")

  hpd_values <- list()
  for (i in seq_along(match_starts)) {
    hpd_str <- substr(newick_str, match_starts[i], match_starts[i] + match_lengths[i] - 1)
    # Extract only the numbers inside the curly braces {min, max}
    inner_match <- regmatches(hpd_str, regexpr("\\{[^}]+\\}", hpd_str))
    if (length(inner_match) > 0) {
      nums <- as.numeric(regmatches(inner_match, gregexpr("[0-9.]+", inner_match))[[1]])
      if (length(nums) >= 2) {
        hpd_values[[i]] <- list(pos = match_starts[i], min = nums[1], max = nums[2])
        cat(sprintf("  HPD %d: [%.2f, %.2f]\n", i, nums[1], nums[2]))
      }
    }
  }

  cat("Found", length(hpd_values), "HPD annotations\n")

  # Calculate node ages
  node_depths <- node.depth.edgelength(tree)
  max_depth <- max(node_depths)
  node_ages <- max_depth - node_depths

  internal_nodes <- (n_tips + 1):(n_tips + n_internal)
  internal_ages <- node_ages[internal_nodes]

  # Create HPD data frame with min/max values
  hpd_df <- data.frame(
    order = 1:length(hpd_values),
    hpd_min = sapply(hpd_values, function(x) x$min),
    hpd_max = sapply(hpd_values, function(x) x$max)
  )

  # Initialize node to HPD mapping
  node_to_hpd <- data.frame(
    node = internal_nodes,
    age = internal_ages,
    hpd_min = NA_real_,
    hpd_max = NA_real_
  )

  # Match each internal node to its HPD by finding the HPD interval
  # that contains the node's age (point estimate should be within CI)
  for (i in seq_len(nrow(node_to_hpd))) {
    node_age <- node_to_hpd$age[i]

    # Find all HPD intervals that contain this age
    for (j in seq_len(nrow(hpd_df))) {
      hpd_min <- hpd_df$hpd_min[j]
      hpd_max <- hpd_df$hpd_max[j]

      # Check if node age falls within this HPD interval (with small tolerance)
      if (node_age >= hpd_min - 1 && node_age <= hpd_max + 1) {
        # Check if this HPD hasn't been assigned yet
        already_used <- any(node_to_hpd$hpd_min == hpd_min &
                           node_to_hpd$hpd_max == hpd_max, na.rm = TRUE)
        if (!already_used) {
          node_to_hpd$hpd_min[i] <- hpd_min
          node_to_hpd$hpd_max[i] <- hpd_max
          break
        }
      }
    }
  }

  # Report matching results
  matched <- sum(!is.na(node_to_hpd$hpd_min))
  cat("Matched", matched, "HPD annotations to", nrow(node_to_hpd), "internal nodes\n")

  # Debug: print matching details
  cat("\nNode-HPD matching details:\n")
  for (i in seq_len(nrow(node_to_hpd))) {
    if (!is.na(node_to_hpd$hpd_min[i])) {
      cat(sprintf("  Node %d: age=%.1f -> HPD=[%.1f, %.1f]\n",
                  node_to_hpd$node[i],
                  node_to_hpd$age[i],
                  node_to_hpd$hpd_min[i],
                  node_to_hpd$hpd_max[i]))
    } else {
      cat(sprintf("  Node %d: age=%.1f -> NO MATCH\n",
                  node_to_hpd$node[i],
                  node_to_hpd$age[i]))
    }
  }
  cat("\n")

  return(list(tree = tree, hpd_data = node_to_hpd, max_depth = max_depth))
}

# Parse tree
result <- parse_tree_with_hpd(newick_with_hpd)
tree <- result$tree
hpd_node_data <- result$hpd_data
max_depth <- result$max_depth

cat("Tree parsed successfully\n")
cat("Number of tips:", length(tree$tip.label), "\n")
cat("Tip labels:", paste(tree$tip.label, collapse = ", "), "\n")

# Get node information
n_tips <- length(tree$tip.label)
n_nodes <- tree$Nnode

# Calculate node ages
node_depths <- node.depth.edgelength(tree)
if (is.null(max_depth)) {
  max_depth <- max(node_depths)
}
node_ages <- max_depth - node_depths

cat("Root age (max depth):", max_depth, "\n")

# Create a data frame for node annotations
node_info <- data.frame(
  node = 1:(n_tips + n_nodes),
  age = node_ages,
  hpd_min = NA_real_,
  hpd_max = NA_real_
)

# Add HPD data from parsed results
if (!is.null(hpd_node_data) && nrow(hpd_node_data) > 0) {
  for (i in seq_len(nrow(hpd_node_data))) {
    node_idx <- hpd_node_data$node[i]
    node_info$hpd_min[node_idx] <- hpd_node_data$hpd_min[i]
    node_info$hpd_max[node_idx] <- hpd_node_data$hpd_max[i]
  }
  cat("HPD data added to node_info\n")
  cat("Nodes with HPD:", sum(!is.na(node_info$hpd_min)), "\n")
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

# Add node info with HPD data
p <- p %<+% node_info

# Check if HPD data is available
has_hpd <- any(!is.na(node_info$hpd_min))

# Add HPD bars FIRST (before node points) so they appear behind
if (has_hpd) {
  cat("Adding HPD confidence interval bars\n")

  # Get tree data with node positions
  tree_data <- as.data.frame(p$data)

  # Merge HPD data
  hpd_subset <- node_info[, c("node", "hpd_min", "hpd_max")]
  tree_data <- merge(tree_data, hpd_subset, by = "node", all.x = TRUE)

  # Filter to nodes with HPD data
  hpd_plot_data <- tree_data[!is.na(tree_data$hpd_min), ]

  if (nrow(hpd_plot_data) > 0) {
    # Convert HPD ages to x-coordinates
    # older age = smaller x (further left on reversed axis)
    hpd_plot_data$x_min <- max_depth - hpd_plot_data$hpd_max
    hpd_plot_data$x_max <- max_depth - hpd_plot_data$hpd_min

    cat("HPD bars to plot:", nrow(hpd_plot_data), "\n")
    for (i in seq_len(nrow(hpd_plot_data))) {
      cat(sprintf("  Node %d: age=%.1f, HPD=[%.1f, %.1f]\n",
                  hpd_plot_data$node[i],
                  max_depth - hpd_plot_data$x[i],
                  hpd_plot_data$hpd_min[i],
                  hpd_plot_data$hpd_max[i]))
    }

    # Add horizontal bars for 95% HPD intervals
    p <- p + geom_segment(
      data = hpd_plot_data,
      aes(x = x_min, xend = x_max, y = y, yend = y),
      color = "#3182bd",
      linewidth = 3,
      alpha = 0.5,
      lineend = "round"
    )

    # Add vertical caps at the ends of HPD bars
    cap_height <- 0.15
    p <- p + geom_segment(
      data = hpd_plot_data,
      aes(x = x_min, xend = x_min, y = y - cap_height, yend = y + cap_height),
      color = "#3182bd",
      linewidth = 1,
      alpha = 0.7
    ) + geom_segment(
      data = hpd_plot_data,
      aes(x = x_max, xend = x_max, y = y - cap_height, yend = y + cap_height),
      color = "#3182bd",
      linewidth = 1,
      alpha = 0.7
    )
  }
}

# Add node points for internal nodes (on top of HPD bars)
p <- p + geom_nodepoint(
  size = 2.5,
  color = "darkred",
  alpha = 0.9
)

# Add age labels with HPD for internal nodes
if (has_hpd) {
  # Create label with age and HPD range
  node_info$label_text <- ifelse(
    !is.na(node_info$hpd_min),
    sprintf("%.1f\n[%.1f-%.1f]",
            node_info$age,
            node_info$hpd_min,
            node_info$hpd_max),
    ifelse(node_info$node > n_tips,
           sprintf("%.1f", node_info$age),
           "")
  )

  # Re-add node info with labels
  p <- p %<+% node_info[, c("node", "label_text")]

  p <- p + geom_nodelab(
    aes(label = label_text),
    size = 2.2,
    hjust = -0.05,
    vjust = 0.3,
    color = "darkblue",
    lineheight = 0.8
  )
} else {
  # Simple age labels without HPD
  p <- p + geom_nodelab(
    aes(x = x, label = ifelse(!isTip, round(max_depth - x, 1), "")),
    size = 2.8,
    hjust = -0.1,
    vjust = -0.5,
    color = "darkblue"
  )
}

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

# Add title and subtitle
if (has_hpd) {
  p <- p + labs(
    title = "Time-Calibrated Phylogenetic Tree",
    subtitle = "Blue bars: 95% HPD confidence intervals"
  ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "#3182bd")
    )
} else {
  p <- p + ggtitle("Time-Calibrated Phylogenetic Tree") +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
}

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

cat("Done!\n")
