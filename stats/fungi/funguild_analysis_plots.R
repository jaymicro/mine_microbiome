
# Import packages ---------------------------------------------------------------------------------------


library(tidyverse)
library(vegan)
library(FUNGuildR)
library(janitor)
library(pscl)
library(performance)
library(patchwork)



# Import data -------------------------------------------------------------------------------------------


df_asv_count <- read.csv("../../data/fungi/ASVs_counts.tsv", sep = "\t") %>% 
  select(-HVC11) %>% 
  clean_names(., "all_caps") %>% 
  column_to_rownames(., var = "X") %>% 
  select_if(colSums(.)>2000) %>% 
  rownames_to_column(., var = "ASVs") %>% 
  adorn_totals(where = "col") %>% 
  filter(Total > 152) %>% 
  # adorn_percentages(denominator = "col") %>% 
  select(-Total) %>% 
  column_to_rownames(., var = "ASVs") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(., var = "sample_id")

cleaned_asv <- df_asv_count %>% 
  select(-sample_id) %>% 
  colnames() 

tax <- read.csv("../../data/fungi/ASVs_taxonomy.tsv", sep = "\t") %>% 
  filter(X %in% cleaned_asv)

meta_df <- read.csv("../../data/soil/meta2.csv") %>% 
  janitor::clean_names() %>% 
  filter(site == "HVC" & site_code != "TBL2_CTRL") %>% 
 # select(-p_h) %>% 
  mutate(
    # al_ppm = al_ppm/100,
    #      b_ppm = b_ppm/100,
    #      cu_ppm = cu_ppm/100,
    #      mn_ppm = mn_ppm/100,
    #      mo_ppm = mo_ppm/100,
    #      na_ppm = na_ppm/100,
    #      zn_ppm = zn_ppm/100,
    #      fe_ppm = fe_ppm/100,
    
         sample_id = str_c(site_code, "_", replicate),
         sample_id = str_remove(sample_id, pattern = "_$"),
         age = recode(age, `#VALUE!` = "reference"))

names(meta_df) <- gsub(x = names(meta_df), pattern = "ppm", replacement = "percent")  

setdiff(df_asv_count$sample_id, meta_df$sample_id)

filt <- df_asv_count$sample_id %>% 
  unique()

df_meta <- meta_df %>% 
  filter(sample_id %in% filt)



setdiff(filt,meta_df$sample_id)
setdiff(meta_df$sample_id,filt)



count_meta_df <- left_join(df_asv_count, df_meta, by = "sample_id") %>% 
  column_to_rownames(., var = "sample_id")




df_count <- count_meta_df %>% 
  select_at(vars(contains("ASV"))) %>% 
  rownames_to_column(., var = "sample_id") %>% 
  janitor::adorn_percentages() %>% 
  column_to_rownames(., var = "sample_id")


meta_data <- count_meta_df %>% 
  select_at(vars(!contains("ASV")))



# Fungal guild analysis ---------------------------------------------------------------------------------

fung <- get_funguild_db()


tax_fun <- tax %>% 
  drop_na(Genus) %>% 
  replace(is.na(.), "") %>% 
  mutate(Taxonomy = paste(Kingdom, Phylum, Class, Order, Family, Genus, Species, sep = ','))


sample_guilds <- funguild_assign(tax_fun)


guilds <- sample_guilds %>% 
  drop_na(guild) %>% 
  filter(confidenceRanking %in% c("Probable", "Highly Probable"))

guilds_filtered <- guilds %>% 
  filter(trophicMode %in% c("Pathotroph", "Saprotroph", "Symbiotroph"))


pathotroph <- guilds_filtered %>% 
  filter(trophicMode == "Pathotroph") %>% 
  pull(X)

saprotroph <- guilds_filtered %>% 
  filter(trophicMode == "Saprotroph") %>% 
  pull(X)

symbiotroph <- guilds_filtered %>% 
  filter(trophicMode == "Symbiotroph") %>% 
  pull(X)


# Observed richness of fungal guilds --------------------------------------------------------------------


# Pathotrophs -------------------------------------------------------------------------------------------


path_rich <- df_count %>% 
  select(any_of(pathotroph)) %>% 
  mutate(rich = specnumber(.,)) %>% 
  select(rich) %>% 
  rownames_to_column(., var = "sample_id") %>% 
  mutate(age = factor(sample_id, levels = df_meta$sample_id, labels = df_meta$age),
         reclamation = factor(sample_id, levels = df_meta$sample_id, labels = df_meta$reclamation_method),
         site_code = str_remove(sample_id, pattern = "_..$")) %>% 
  filter(age != "reference") %>% 
  mutate(age = as.numeric(paste(age)))


m3 <-glm(rich ~ poly(age,2) + reclamation, family = "poisson",  path_rich)
summary(m3)
car::Anova(m3, test.statistic = "F")
r2(m3)
plot(resid(m3))


patho_plot <- ggplot(path_rich, aes(y = rich, x = age, fill = age)) +
  geom_point(pch = 21, size =3) +
  geom_smooth(method = "lm", color = "black") +
  scale_fill_viridis_c(direction = -1) +
  theme_bw()+
  theme(legend.position = "none",
        legend.title = element_text(size =18),
        legend.text = element_text(size = 16),
        legend.key.width = unit(1.25, "cm"),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16, color = "black")) +
  labs(x = NULL,
       y = "Richness",
       title = "Pathotrophs") +
  scale_x_continuous(breaks = c(seq(0, 30, by =5))) +
  annotate(geom = "text", x = 20, y=30,
           label =  expression(paste(R^2 == 0.326, "; pval < 0.005")), size =4,
           parse = T) 
patho_plot


# Saprotrophs -------------------------------------------------------------------------------------------


sap_rich <- df_count %>% 
  dplyr::select(any_of(saprotroph))%>% 
  mutate(rich = specnumber(.,)) %>% 
  select(rich) %>% 
  rownames_to_column(., var = "sample_id") %>% 
  mutate(age = factor(sample_id, levels = df_meta$sample_id, labels = df_meta$age),
         reclamation = factor(sample_id, levels = df_meta$sample_id, labels = df_meta$reclamation_method),
         site_code = str_remove(sample_id, pattern = "_..$")) %>% 
  filter(age != "reference") %>% 
  mutate(age = as.numeric(paste(age)))


m4 <- glm(rich ~ poly(age,2) + reclamation, family = "poisson",  sap_rich)
summary(m4)
car::Anova(m4, test.statistic = "F")
r2(m4)

sapro_plot <- ggplot(sap_rich, aes(y = rich, x = age, fill = age)) +
  geom_point(pch = 21, size =3) +
  geom_smooth(method = "lm",
              formula = y ~ x + I(x^2),
              color = "black") +
  scale_fill_viridis_c(name = "Reclamation age (yr)",direction = -1) +
  theme_bw()+
  theme(legend.position = "bottom",
        legend.title = element_text(size =18),
        legend.text = element_text(size = 16),
        legend.key.width = unit(1.25, "cm"),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16, color = "black"))+
  labs(x = NULL,
       y = NULL,
       title = "Saprotrophs") +
  scale_x_continuous(breaks = c(seq(0, 30, by =5))) +
  annotate(geom = "text", x = 20, y=120,
           label =  expression(paste(R^2 ==  0.819, "; pval < 0.01")), size =4,
           parse = T) 
sapro_plot


# Symbiotrophs ------------------------------------------------------------------------------------------


symbio_rich <- df_count %>% 
  dplyr::select(any_of(symbiotroph)) %>% 
  mutate(rich = specnumber(.,)) %>% 
  select(rich) %>% 
  rownames_to_column(., var = "sample_id") %>% 
  mutate(age = factor(sample_id, levels = df_meta$sample_id, labels = df_meta$age),
         reclamation = factor(sample_id, levels = df_meta$sample_id, labels = df_meta$reclamation_method),
         site_code = str_remove(sample_id, pattern = "_..$")) %>% 
  filter(age != "reference") %>% 
  mutate(age = as.numeric(paste(age)))


m5 <- hurdle(rich ~ poly(age,2) + reclamation,dist = "poisson", zero.dist = "binomial", link = "logit", symbio_rich)
summary(m5)
car::Anova(m5, test.statistic = "F")
r2(m5)


symbio_plot <- ggplot(symbio_rich, aes(y = rich, x = age, fill = age)) +
  geom_point(pch = 21, size =3) +
  geom_smooth(method = "lm",
              formula = y ~ x + I(x^2),
              color = "black") +
  scale_fill_viridis_c(direction = -1) +
  theme_bw()+
  theme(legend.position = "none",
        legend.title = element_text(size =18),
        legend.text = element_text(size = 16),
        legend.key.width = unit(1.25, "cm"),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16, color = "black")) +
  labs(x = NULL,
       y = NULL,
       title = "Symbiotrophs") +
  scale_x_continuous(breaks = c(seq(0, 30, by =5))) +
  annotate(geom = "text", x = 10, y=40,
           label =  expression(paste(R^2 ==   0.78, "; pval < 0.001")), size =4,
           parse = T) 
symbio_plot




# Combining all plots -----------------------------------------------------------------------------------

funguild_plot <- patho_plot + sapro_plot + symbio_plot+ plot_annotation(tag_levels = 'A') #+
  #plot_layout(guides = "collect")  & theme(legend.position = 'bottom')
funguild_plot 


#ggsave(funguild_plot, filename = "../../plot/funguild_plot.jpg", width = 14, height =6, dpi = 1000)


#save.image("funguild.RData")

