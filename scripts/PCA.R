library(tidyverse)
wormpost <- read_csv(here::here('data/wormpost_strains.csv')) %>%
  separate(source, into = c("worm_shorthand", "replicate"), sep = c("_")) %>%
  mutate(Genus = interaction("g", Genus, sep = "__"),
         Species = "s__")

taxonomy <- read_tsv('/Volumes/4TB_NGS/Providencia_files/taxonomy/taxonomy_r86_August2018.tsv',
                     col_names = c("ID", "taxonomy")) %>%
                     separate(taxonomy, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";") %>%
  select(3:7) %>%
  distinct() %>%
  filter(Genus %in% wormpost$Genus) %>%
  mutate(OTU = interaction("OTU", 1:nrow(.), sep = "_"))

top_family <- wormpost %>%
  filter(!is.na(My_Family)) %>%
  group_by(My_Family) %>%
  tally() %>%
  filter(n > 4)

wormpost %>%
  filter(My_Family %in% top_family$My_Family) %>%
  #group_by(worm_genotype, replicate, genus) %>%
  filter(!is.na(My_Family), !My_Family %in% c("Rhodotorula.yeast", "Trichosporonaceae")) %>%
  #tally() %>%
  ggplot(aes(x = worm_strain, fill = My_Family)) +
  geom_bar(position = "fill") +
  #facet_grid(.~worm_strain) +
  scale_fill_brewer(palette = "Dark2")

MO_pseudomonas <- wormpost %>% filter(Genus == "Pseudomonas", set == "MO")

family_by_group <- wormpost %>%
  filter(!is.na(Family)) %>%
  group_by(worm_strain, Family) %>%
  tally()

sample_OTU_all <- left_join(wormpost, taxonomy) %>%
  select(worm_strain, replicate, OTU, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
  mutate(sample_ID = interaction(worm_strain, replicate, sep = "_")) %>%
  filter(!is.na(sample_ID))

OTU_ids <- sample_OTU_all %>%
  select(OTU, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
  distinct()

OTU_table <- sample_OTU_all %>%
  group_by(OTU, sample_ID) %>%
  tally() %>%
  pivot_wider(names_from = sample_ID,
              values_from = n,
              values_fill = list(n = 0)) %>%
  left_join(., OTU_ids) %>%
  filter(!is.na(OTU))

metadata = tibble(SampleID = unique(sample_OTU_all$sample_ID),
                  type = c(rep("JU322", 3),
                           rep("N2", 3),
                           rep("PX178", 3),
                           "compost",
                           rep("fruit", 4),
                           rep("mollusc", 3)))

myotutable <- ampvis2::amp_load(otutable = OTU_table,
                                metadata = metadata)

ampvis2::amp_ordinate(
  data = myotutable,
  type = "PCA",
  constrain = "type",
  distmeasure = "hellinger",
  sample_color_by = "type",
  sample_colorframe = TRUE,
  sample_colorframe_label = "SampleID",
  species_plotly = TRUE)


