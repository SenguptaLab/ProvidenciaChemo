files <- fs::dir_ls(path = '/Volumes/4TB_NGS/Providencia_files/barcoding.tasks.lima-0/JUb39_assembly_5kmin/ideel_preds/lengths/', regexp = "chromosome")
lenghts <- purrr::map_df(files, read_tsv, col_names = FALSE, .id = "filename")

lenghts %>%
  mutate(dataset = basename(filename)) %>%
  ggplot() + geom_histogram(aes(x=X1/X2, fill = dataset),
                            position = "identity", alpha = 0.5) +
  coord_cartesian(ylim = c(0,500))
