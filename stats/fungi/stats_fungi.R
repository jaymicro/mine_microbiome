
# Import packages ---------------------------------------------------------


library(tidyverse)
library(janitor)
library(vegan)
library(ape)
library(ggpubr)
library(glue)
library(lme4)
library(lmerTest)
library(patchwork)
library(scales)
library(iCAMP)









# Import and wrangle data -------------------------------------------------


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
#  select(-p_h) %>% 
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
  select_at(vars(contains("ASV")))


meta_data <- count_meta_df %>% 
  select_at(vars(!contains("ASV")))


# Bray Curtis dissimilarity on hellinger transformed data -----------------


bray_dist <- df_count %>% 
  decostand(., method = "hellinger") %>% 
  vegdist(., method = "bray")


# Principle  coordinate analysis ------------------------------------------

pc <- ape::pcoa(bray_dist)
pc$values
biplot(pc)

# Variance of the axes ----------------------------------------------------

variance <- round(pc$values[1:2, 2] *100, 1 )



# Permanova ---------------------------------------------------------------

set.seed(111)
mod_adonis <- adonis2(bray_dist ~  meta_data$reclamation_method + meta_data$age,
        by= "margin", permutations = 999)
mod_adonis
reclam <-  round(c(mod_adonis$F[1], mod_adonis$`Pr(>F)`[1]),3)
age <- round(c(mod_adonis$F[2], mod_adonis$`Pr(>F)`[2]),3)


# Plot PCoa ---------------------------------------------------------------

plot_df <- data.frame(meta_data,
                      PC1= pc$vectors[,1],
                      PC2= pc$vectors[,2]) %>% 
  mutate(age = factor(age, levels = c("3", "4", "5","6","10","14","17","19","20","24","26","reference")))

p1 <- plot_df %>% 
  ggplot(., aes(x = PC1, y = PC2, fill =age)) +
  geom_point(pch= 21, size = 5, alpha = 0.75) +
  scale_fill_viridis_d(name = "Reclamation age (yr)", direction = -1) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_text(size =16),
        legend.text = element_text(size =16),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16, color = "black")) +
  guides(fill = guide_legend(nrow = 2)) +
  labs(x = glue("Axis1 ({variance[1]} %)"),
       y = glue("Axis2 ({variance[2]} %)"),
       title = "Fungal community") +
  annotate(geom = "text", x = c(-0.25,-0.25), y = c(-0.45,-0.35),
           label = c(glue("Reclamation method:
                           pseudoF = {reclam[1]}, p-val < {reclam[2]}"),
                     glue("Reclamation age:
                           pseudoF = {age[1]}, p-val < {age[2]}")))

p1




# Convert bray-curtis dissimilarity matrix to long form and similarity --------

bc_mat <- bray_dist %>% 
  as.matrix(.) %>% 
  as.data.frame(.) %>% 
  rownames_to_column(., var = "name1") %>% 
  pivot_longer(-name1) %>% 
  filter(name1 != name) %>% 
  mutate(name1_id = as.numeric(paste(factor(name1, levels = row.names(meta_data), labels = meta_data$age))),
         name_id = factor(name, levels = row.names(meta_data), labels = meta_data$age),
         name1_reclamation = factor(name1, levels = row.names(meta_data), labels = meta_data$reclamation_method),
         name_reclamation = factor(name, levels = row.names(meta_data), labels = meta_data$reclamation_method),
         sim = (1- value)*100,
         sampling_yr = factor(name1, levels = row.names(meta_data), labels = meta_data$sampling_year),
         site_code = factor(name1, levels = row.names(meta_data), labels = meta_data$site_code)) %>% 
  filter( name_id == "reference",
          name1_reclamation != "control") %>% 
  mutate(age = as.numeric(paste(name1_id)))



# model for similarity against age and reclamation treatment --------------

m2 <- lmer(sim ~ age + name1_reclamation + (1|name1) ,
          data = bc_mat,
          na.action = na.omit)
summary(m2)
car::Anova(m2, test.statistic = "F")
round(MuMIn::r.squaredGLMM(m2)[1],3)
emmeans::emmeans(m2, pairwise ~ name1_reclamation)

# Plot similarity ---------------------------------------------------------

p2 <- ggplot(bc_mat, aes(x = age, y = sim)) +
  geom_point(pch = 21, size = 2,alpha = 0.75,  aes(fill = age),show.legend = F)+
  geom_smooth(method = "lm", se = T, color = "black", size = 1.25) +
  scale_x_continuous(breaks = seq(0,30,by = 5)) +
  scale_y_continuous(labels = function(x) paste0(x, "%"))+
  scale_fill_viridis_c( direction = -1) +
  theme_bw() +
  theme(legend.title = element_text(size =16),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16, color = "black")) +
  labs(x = "Reclamation age (yr)",
       y = "Similarity with reference site") +
  annotate(geom = "text", x = 7, y=30,
           label =  expression(atop(paste(R[m]^2 == "0.193", "; pval < 0.001"))), size =4,
  parse = T)
  
p2

fungi_plot <- p1 + p2 + plot_annotation(tag_levels = 'A')  + plot_layout(guides = "collect")  & theme(legend.position = 'bottom')

fungi_plot

#ggsave(fungi_plot, filename = "../../../plot/fungi.jpg", width = 13, height = 7, dpi = 1000)





# Alpha diversity ---------------------------------------------------------------------------------------



alpha_div <- data.frame(
  shann = diversity(df_count, index = "shannon"),
  simp = diversity(df_count, index = "simpson"),
  richness = specnumber(df_count)
) %>% 
  rownames_to_column(., var = "sample_id")

alpha_df <- meta_data %>% 
  rownames_to_column(., var = "sample_id") %>% 
  inner_join(alpha_div, by = "sample_id") %>% 
  filter(age !=  "reference") %>% 
  mutate(age = as.numeric(paste(age)))

str(alpha_df)

#not significant
m3 <- lmer(richness ~ age + reclamation_method + (1|site_code), 
           data = alpha_df)
summary(m3)
performance::r2(m3)

#not significant
m4 <- lmer(simp ~ age + reclamation_method + (1|site_code), 
           data = alpha_df)
summary(m4)
performance::r2(m4)

#not significant
m5 <- lmer(simp ~ age + reclamation_method + (1|site_code), 
           data = alpha_df)
summary(m5)
performance::r2(m5)





