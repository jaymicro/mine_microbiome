
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
library(broom.mixed)







# Import and wrangle data -------------------------------------------------


df_asv_count <- read.csv("../../data/bacteria/ASVs_counts.tsv", sep = "\t") %>% 
  dplyr::select(-HVC11) %>% 
  janitor::clean_names(., "all_caps") %>% 
  column_to_rownames(., var = "X") %>% 
  select_if(colSums(.)>2000) %>% 
  rownames_to_column(., var = "ASVs") %>% 
  adorn_totals(where = "col") %>% 
  filter(Total > 52) %>% 
 # adorn_percentages(denominator = "col") %>% 
  dplyr::select(-Total) %>% 
  column_to_rownames(., var = "ASVs") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(., var = "sample_id")

tax <- read.csv("../../data/bacteria/ASVs_taxonomy.tsv", sep = "\t")


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



# Variance of the axes ----------------------------------------------------

variance <- round(pc$values[1:2, 2] *100, 1 )



# Permanova ---------------------------------------------------------------

set.seed(111)
mod_adonis <- adonis2(bray_dist ~  meta_data$reclamation_method + meta_data$age,
        by= "margin",  permutations = 999)
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
       title = "Bacterial community") +
  annotate(geom = "text", x = c(-0.3,-0.3), y = c(-0.30,-0.40),
           label = c(glue("Reclamation method:
                           pseudoF = {reclam[1]}, p-val = {reclam[2]}"),
                     glue("Reclamation age:
                           pseudoF = {age[1]}, p-val = {age[2]}")))

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

m2 <- lmer(sim ~ age + name1_reclamation + (1|site_code)  ,
          data = bc_mat,
          na.action = na.omit, REML = FALSE)
summary(m2)
anova(m2)
MuMIn::r.squaredGLMM(m2)
emmeans::emmeans(m2, pairwise ~ name1_reclamation)



# Plot similarity ---------------------------------------------------------

p2 <- ggplot(bc_mat, aes(x = age, y = sim)) +
  geom_point(pch = 21, size = 2,alpha = 0.75,  aes(fill = age),show.legend = F)+
  geom_smooth(method = "lm", se = T, color = "black", size = 1.25) +
  scale_x_continuous(breaks = seq(0,30,by = 5)) +
  scale_y_continuous(labels = function(x) paste0(x, "%"))+
  scale_fill_viridis_c(name = NULL,
                       direction = -1) +
  theme_bw() +
  theme(
        # legend.title = element_text(size =16),
        # legend.position = "bottom",
        # legend.text = element_text(size =16),
        # legend.key.width = unit(2, "cm"),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16, color = "black")) +
  labs(x = "Reclamation age (yr)",
       y = "Similarity with reference site") +
  annotate(geom = "text", x = 7, y=35,
           label =  expression(atop(paste(R[m]^2 == 0.167, "; pval < 0.001"))), size =4,
  parse = T)
  
p2

bac_plot <- p1 + p2 + plot_annotation(tag_levels = 'A') + plot_layout(guides = "collect")  & theme(legend.position = 'bottom')

bac_plot


#ggsave(bac_plot, filename = "../../../plot/bac.jpg", width = 13, height = 7, dpi = 1000)


# BNTI bacteria -------------------------------------------------------------------------------

bnti_bac <- read.csv("../../data/bacteria/bnti.csv") 



# Logistic regression -------------------------------------------------------------------------
bnti_bac$name_reclamation <- as.factor(bnti_bac$name_reclamation)
bnti_bac$name_reclamation <- relevel(bnti_bac$name_reclamation, ref = "no_biosolids")

mod_log <-  glmer(assembly ~ age + name_reclamation +  (1|name1/name) ,
             family =binomial(link = "logit"), bnti_bac)

car::Anova(mod_log)




# Plot odd ratio ------------------------------------------------------------------------------

p3 <- tidy(mod_log,conf.int=TRUE,exponentiate=TRUE,effects="fixed") %>% 
  filter(term != "(Intercept)") %>% 
  ggplot(., aes(x = estimate, y  = term)) +
  geom_point(aes(y = term, x = estimate), pch= 23, fill = "black", size = 2, show.legend = F) +
  geom_linerange(aes(xmin = conf.low, xmax = conf.high), show.legend = F) +
  geom_vline(aes(xintercept = 1), linetype = "dashed")+
  expand_limits(x = 3)+
  scale_x_log10() +
  annotation_logticks(sides = "b") +
  scale_y_discrete(name = NULL,
                   labels = c("name_reclamationbiosolids" = "Biosolids\nTreatment",
                              "age" = "Reclamation\nAge")) +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12, color = "black"),
        axis.text.y = element_text(hjust = 0.5)) +
  annotate(geom = "text", x = c(0.916, 0.394), y = c(1.05,2.05), label = c("**", "*"))+ 
  labs(x = "Odds ratio")

p3
#ggsave(p3, filename = "../../../plot/bac_bnti.jpg", width = 5, height = 3, dpi = 1000)

p4 <- ggeffects::ggpredict(mod_log, terms = "age") %>%
  plot() +
  theme_bw() +
  theme(axis.text = element_text(size = 16, color = "black"),
        axis.title = element_text(size = 16, color = "black"))+
  scale_x_continuous(breaks = seq(0,30,by = 5)) +
  labs(title = "Predicted probablities of\nDeterministic Assembly",
       y = NULL,
       x = "Reclamation age")
  
p4

p5 <- ggeffects::ggpredict(mod_log, terms = "name_reclamation") %>%
  plot() +
  theme_bw() +
  theme(axis.text = element_text(size = 16, color = "black"),
        axis.title = element_text(size = 16, color = "black"))+
   labs(title = NULL,
       y = NULL,
       x = "Reclamation method")

p5


assembly_plot <- p4 + p5 + plot_annotation(tag_levels = "A")
assembly_plot

#ggsave(assembly_plot, filename = "../../../plot/assembly_plot.jpg", width = 12, height = 7, dpi = 1000)


# Alpha diversity ---------------------------------------------------------------------------------------
#no significant differences in alpha diversity

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

#Not significant
m3 <- lm(shann ~ age + reclamation_method , 
           data = alpha_df)
summary(m3)
anova(m3)


m4 <- glm(richness ~ age + reclamation_method , 
           family = "poisson",data = alpha_df)
summary(m4)
car::Anova(m4, test.statistic = "F")


m5 <- lm(simp ~ age + reclamation_method , 
            data = alpha_df)
summary(m5)
anova(m5)


#alpha diversity plot
bac_shann <- ggplot(alpha_df, aes(x = age, y = shann)) +
  geom_point(pch = 21, size = 2,alpha = 0.75,  aes(fill = age),show.legend = F)+
  geom_smooth(method = "lm", se = T, color = "black", size = 1.25) +
  scale_fill_viridis_c(name = "Reclamation age (yr)", direction = -1) +
  theme_bw() +
  theme(legend.title = element_text(size =16),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16, color = "black")) +
  scale_x_continuous(breaks = seq(0,35,by = 5)) +
  labs(x = NULL,
       y = "Shannon-Weiner index")
bac_shann

bac_richness <- ggplot(alpha_df, aes(x = age, y = richness)) +
  geom_point(pch = 21, size = 2,alpha = 0.75,  aes(fill = age),show.legend = F)+
  geom_smooth(method = "lm", se = T, color = "black", size = 1.25) +
  scale_fill_viridis_c(name = "Reclamation age (yr)", direction = -1) +
  theme_bw() +
  theme(legend.title = element_text(size =16),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16, color = "black")) +
  scale_x_continuous(breaks = seq(0,35,by = 5)) +
  labs(x = NULL,
       y = "Observed richness")
bac_richness 


bac_simp <- ggplot(alpha_df, aes(x = age, y = simp)) +
  geom_point(pch = 21, size = 2,alpha = 0.75,  aes(fill = age),show.legend = T)+
  geom_smooth(method = "lm", se = T, color = "black", size = 1.25) +
  scale_fill_viridis_c(name = "Reclamation age (yr)", direction = -1) +
  theme_bw() +
  theme(legend.title = element_text(size =16),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16, color = "black"),
        legend.position = 'bottom') +
  scale_x_continuous(breaks = seq(0,35,by = 5)) +
  labs(x = NULL,
       y = "Simpson index")
bac_simp


alpha_div <- bac_shann + bac_simp + bac_richness +
  plot_annotation(tag_levels = "A") 

alpha_div

#ggsave(alpha_div, filename = "../../../plot/bac_alpha_div.jpg", width = 16, height = 7, dpi = 1000)

div <- c("Observed\nrichness", "Shannon-Weiner\nindex", "Simpson\nindex")
names(div) <- c("richness", "shann", "simp")

comp <- list(c("biosolids", "no_biosolids"))

bac_alpha_div_reclam_method <- alpha_df %>% 
  select(reclamation_method, shann, simp, richness) %>% 
  pivot_longer(-reclamation_method) %>% 
  ggplot(.,aes(x = reclamation_method, y = value, fill = reclamation_method)) +
  geom_violin(alpha = 0.6, show.legend = F) +
  geom_boxplot(width = 0.3, color = "black")+
  facet_wrap(~ name, scales = "free", labeller = labeller(name = div))+
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_fill_viridis_d(name = "Reclamation method") +
  ggpubr::stat_compare_means(comparisons = comp, label = "p.signif") +
  labs(x = NULL,
       y = "Alpha diversity")

bac_alpha_div_reclam_method

#ggsave(bac_alpha_div_reclam_method, filename = "../../../plot/bac_alpha_div_reclam_method.jpg", width = 16, height = 7, dpi = 1000)




