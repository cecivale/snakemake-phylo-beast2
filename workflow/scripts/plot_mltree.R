# ------------------------------------------------------------------------------
#          ---        
#        / o o \    Project:  cov-armee Phylodynamics
#        V\ Y /V    Plot phylogenetic trees with metadata
#    (\   / - \     
#     )) /    |     
#     ((/__) ||     Code by Ceci VA 
# ------------------------------------------------------------------------------


library(tidyverse)
library(lubridate)
library(treedataverse)
library(ggsci)
library(scales)

tree_file <- "/Users/ceciliav/Projects/2105-cov-armee/8.final-analysis/results/analysis/old/iqtree/asymp_symp_phylogenetics/asymp_armydem.1.treefile"
#tree_file <- "/Users/ceciliav/Projects/2105-cov-armee/8.final-analysis/results/analysis/old/iqtree/asymp_symp_phylogenetics/sw_community_ns.treefile"
tree_file <- snakemake@input[["tree"]]
tree_ref <- treeio::read.tree(tree_file)
tree <- drop.tip(tree_ref, snakemake@params[["outgroup"]])

# TODO generalize metadata input
metadata_asymp <- read_tsv("/Users/ceciliav/Projects/2105-cov-armee/8.final-analysis/results/data/screening/metadata.tsv") %>%
  mutate(deme = "asymp")
metadata_symp <- read_tsv("/Users/ceciliav/Projects/2105-cov-armee/8.final-analysis/results/data/community/metadata.tsv") %>%
  mutate(deme = "symp")
  
metadata <- bind_rows(metadata_symp, metadata_asymp) 

d <- as_tibble(tree) %>% left_join(metadata, by = c("label" = "sample_id")) %>%
  mutate(cat = case_when(deme == "asymp" ~ "asymp",
                army_dem2 ~ "symp_armydem",
                !army_dem2 ~ "symp",
                TRUE ~ deme),
         age_cat = case_when(#altersjahr <= 18 ~ "0-18",
                             altersjahr <= 30 ~ "0-30",
                             altersjahr <= 60 ~ "31-60",
                             #altersjahr <= 70 ~ "51-70",
                             altersjahr > 60 ~ "+60"))

col_deme <-  c("#ffd900", "#ffa600", "#1464b5")
col_deme <-  c("black", "black", "#00A1FF")
names(col_deme) <- c("symp", "symp_armydem", "asymp")

p <- ggtree(as.treedata(d), size = 0.4, color="grey60")  +# aes(color = cat))  +
            #aes(color =  pangoLineage == "B.1.160.16"))  +
  geom_tippoint(aes(color = cat), size = 1, alpha = 0.8) +
  geom_nodepoint(aes(subset = as.numeric(label) >= 95),shape = 15, size=0.8, color = "grey30") +
  scale_size_continuous(range = c(0, 2), name = "bootstrap") +
  scale_alpha(name = "bootstrap") +
  scale_color_manual(values = col_deme, na.value = "grey30") +
  geom_treescale(x = 0, y = 150, label = "subs/site", fontsize = 3) 

ann <- d %>% filter(!is.na(deme)) %>% 
  mutate(screening_week = case_when(isoweek(date) == 2 ~ 1,
                                    isoweek(date) == 3 ~ 2,
                                    isoweek(date) == 6 ~ 3,
                                    TRUE ~ 3),
         sex = ifelse(is.na(sex) & deme == "asymp", "Männlich", sex)) %>%
  # select(label, cat, screening_week, division, kanton, altersjahr, 
  #        age_cat, sex, 
  #        ct, nextstrainClade, 
  #        ) %>%
  select(label, cat, division, altersjahr, 
         ct, nextstrainClade, 
  ) %>%
  column_to_rownames("label")


col_week <- pal_jama()(length(unique(ann$screening_week)))
names(col_week) <- unique(ann$screening_week)

col_div <- pal_igv()(length( unique(ann$division)))
names(col_div) <- unique(ann$division)

col_kanton <- pal_igv()(length(unique(ann$kanton)))
names(col_kanton) <- unique(ann$kanton)

col_age <- viridis_pal()(length(unique(ann$altersjahr)))
names(col_age) <- unique(ann$altersjahr)

col_age_cat <-  viridis_pal()(length(unique(ann$age_cat)))
names(col_age_cat) <- unique(ann$age_cat)

col_sex <- pal_npg()(length(unique(ann$sex)))
names(col_sex) <- unique(ann$sex)

col_ct <- viridis_pal()(length(unique(ann$ct)))
names(col_ct) <- unique(ann$ct)

col_clade <- pal_aaas()(length(unique(ann$nextstrainClade)))
names(col_clade) <- unique(ann$nextstrainClade)

my_colors <- c(col_deme,
               #col_week,
               col_div,
               #col_kanton,
               col_age,
               #col_age_cat,
               #col_sex,
               col_ct,
               col_clade
            )

ph <- gheatmap(p, ann, 
         colnames_angle = 90, colnames_offset_y = -15, font.size = 3) +
  theme(legend.position = "none") +
  scale_fill_manual(values = my_colors, na.value = "grey40") +
  vexpand(direction = -1, .1)


ggsave(ph, filename = snakemake@output[["fig"]], width = 297, heigh = 210, units = "mm")

# Trying things
# p <- ggtree(tree_ref)
# 
# p <- p %<+% d + geom_tippoint(aes(color=cat))
# 
# age_data <-  d %>% filter(!is.na(deme)) %>%
#   select(id = label, altersjahr, cat, nextstrainClade, pangoLineage)
# p  %<+% d + geom_facet(panel = "Age", data = age_data, geom = geom_col,
#              aes(x = altersjahr, color = cat,
#                  fill = cat))
# facet_plot(p, panel = "Age", data = age_data, geom = geom_col,
#                aes(x = altersjahr,#color = cat,
#                    fill = nextstrainClade), orientation = 'y')
# 
# ggplot(ann %>% filter(cat == "symp")) +
#   geom_boxplot(aes(nextstrainClade, altersjahr))
# 
# 
# 
# c = groupClade(as_tibble(tree), .node = 249) %>% filter(group == 1)
# d %>% filter(!is.na(deme)) %>%
#   mutate(clade = (label %in% c$label)) %>%
#   ggplot() +
#   geom_violin(aes(x = altersjahr, y = clade))


ggtree(as.treedata(d), size = 0.4, color = "grey60")  +
  #aes(color =  pangoLineage == "B.1.160.16"))  +
  geom_tippoint(aes(color = cat), size = 3, alpha = 0.9) +
  geom_nodepoint(aes(subset = as.numeric(label) >= 95, alpha = as.numeric(label)), shape = 15, color = "grey30") +
  scale_alpha(name = "bootstrap") +
  scale_color_manual(values = c(symp = "#4B5938", asymp="#BDBC86"), na.value = "grey30") +
  geom_treescale(label = "subs/site", fontsize = 3) +
  theme(legend.position = "none")

