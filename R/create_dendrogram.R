create_primary_dendrogram <- function(cluster_result, output_dir) {
  tree <- cluster_result$specifications$main_27_hcs$tree
  labels <- cluster_result$specifications$main_27_hcs$data$country
  labels[labels == "Turkey"] <- "Türkiye"
  draw <- function() {
    old <- graphics::par(mar = c(4.2, 4.2, 2.0, 1.0), family = "sans")
    on.exit(graphics::par(old), add = TRUE)
    graphics::plot(tree, labels = labels, hang = -1, cex = .72,
      main = "", xlab = "Country", sub = "",
      ylab = "Gower dissimilarity")
    stats::rect.hclust(tree, k = 4,
      border = c("black", "grey35", "grey55", "grey70"))
  }
  pdf_file <- file.path(output_dir, "07a_supplementary_primary_dendrogram.pdf")
  png_file <- file.path(output_dir, "07a_supplementary_primary_dendrogram.png")
  grDevices::pdf(pdf_file, width = 8.5, height = 6.2, useDingbats = FALSE)
  draw(); grDevices::dev.off()
  grDevices::png(png_file, width = 8.5, height = 6.2, units = "in",
                 res = 600, type = "cairo-png", bg = "white")
  draw(); grDevices::dev.off()
  invisible(c(pdf = pdf_file, png = png_file))
}
