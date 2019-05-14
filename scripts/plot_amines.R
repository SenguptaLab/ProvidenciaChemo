filepath <- here::here("extdata/Figure_2A.csv")

enzymes <- read.csv(filepath) %>%
  dplyr::filter(date %in% c("4_18_17", "6_10_17","5_31_17","5_12_17","6_12_17", "5_10_16", "6_7_16", "6_14_16", "6_16_18", "2019_05_07"), # worms were contaminated after June 2018, bleached for subsequent expreiments
                genotype %in% c("N2", "tbh-1", "tdc-1", "tph-1", "cat-2")) %>%
  format_AvoidData(day.correct = FALSE) %>%
  mutate(strain = factor(strain, levels = c("OP50", "JUb39")),
         genotype = fct_relevel(genotype, c("N2", "tdc-1", "tbh-1", "tph-1", "cat-2")), #"cat-2", "tph-1")),
         plate = factor(seq(1:nrow(.))))

p1 <- enzymes %>%
  #filter(date %in% c("4_18_17", "7_22_18")) %>%
  #filter(date %in% c("4_18_17", "5_31_17","5_12_17","6_12_17", "7")) %>%
  plot_plasticityIndex(xvar = strain) +
  facet_grid(.~genotype) +
  scale_color_plot(palette = "grey-blue", drop = TRUE) +
  guides(colour = FALSE) +
  theme(strip.text.x = element_text(size = 12)) #+
  scale_y_continuous(limits = c(-2.5,2.5))
p1

