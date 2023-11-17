save_plot_and_data <- function(plot_object, width, height, file_name, file_type = "jpg", dpi = 500) {
  # Save the plot using ggsave
  ggsave(paste0(file_name, ".", file_type), plot_object, width = width, height = height, dpi = dpi)
  
  # Save the plot object using saveRDS
  saveRDS(plot_object, file = paste0(file_name, ".rds"))
}

