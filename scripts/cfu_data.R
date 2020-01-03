colonies <- read_csv('data/Colony_counts.csv') %>%
  mutate(supe_total = n_sup * 10, # plated 1/10 volume of last wash after bleaching to get n cells left in supe
         cells_per_worm = ( (n_colonies * dilut) - supe_total ) / n_worms) # normalize based on # cells left over after bleaching

colonies %>% ggplot(aes(x = food, y = cells_per_worm)) +
  geom_bardots(fillvar = food) +
  scale_fill_plot(palette = "grey-blue", drop = TRUE) +
  labs(y = "cfu / worm")



