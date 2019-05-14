filepath <- here::here("extdata/tdcKO_plateChemo.csv")

tdcKO <- read.csv(filepath) %>%
  format_AvoidData(day.correct = "OP50") %>%
  mutate(strain = factor(strain, levels = c("OP50", "JUb39", "JUb39; tdcDel::cmR delAADC")),
         plate = factor(seq(1:nrow(.))))

p1 <- tdcKO %>%
  #filter(date %in% c("4_18_17", "7_22_18")) %>%
  #filter(date %in% c("4_18_17", "5_31_17","5_12_17","6_12_17", "7")) %>%
  plot_plasticityIndex(xvar = strain) +
  facet_grid(.~genotype) +
  scale_color_plot(palette = "grey-blue", drop = TRUE) +
  guides(colour = FALSE) +
  theme(strip.text.x = element_text(size = 12)) #+
scale_y_continuous(limits = c(-2.5,2.5))
p1

