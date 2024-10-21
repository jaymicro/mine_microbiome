
# Import packages ---------------------------------------------------------------------------------------


library(tidyverse)
library(vegan)
library(lme4)
library(lmerTest)
library(janitor)
library(glue)
library(patchwork)



# Import gene abundance ----------------------------------------------------------------------------------------


gene <- read_tsv("../../data/shotgun/relab_genefamilies.tsv", show_col_types = FALSE) %>% 
  janitor::clean_names() %>% 
  filter(!str_detect(gene_family, pattern = "\\|"))

metadata <- read.csv("../../data/shotgun/meta_data.csv")



# Bray curtis dissimilarity of gene abundance -----------------------------------------------------------


bray_dist <- gene %>% 
  column_to_rownames(., var = "gene_family") %>% 
  t() %>% 
  as.data.frame() %>% 
 # decostand(., method = "hellinger") %>% 
  vegdist(., method = "bray")


# permanova ---------------------------------------------------------------------------------------------


set.seed(121)
mod_adonis<- adonis2(bray_dist ~ metadata$age,  permutations = 999)
mod_adonis
age <- round(c(mod_adonis$F[1], mod_adonis$`Pr(>F)`[1]),3)

pc <- ape::pcoa(bray_dist)
variance <- round(pc$values$Relative_eig[c(1,2)]*100,2)



# long data of bray-curtis dissimilarity ----------------------------------------------------------------


bd <- bray_dist %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  rownames_to_column(., var = "sample") %>% 
  pivot_longer(-sample) %>% 
  mutate(sample = str_remove(sample, pattern = "_merged_abundance_rp_ks"),
         name = str_remove(name, pattern = "_merged_abundance_rp_ks"),
         sample_age = factor(sample, levels = metadata$sample, label = metadata$age),
         name_age = factor(name, levels = metadata$sample, label = metadata$age)) %>% 
  filter(sample_age == "control",
         name_age != "control") %>% 
  mutate(age = as.numeric(paste(name_age)),
         site = str_remove(name, pattern = "_([A-z][0-9]|[0-9][A-z])_s[0-9]"),
         sim = 100*(1-value))

str(bd)

m2 <- lmer(sim ~ age + (1|sample), bd)
summary(m2)
car::Anova(m2, test.statistic = "F")
round(MuMIn::r.squaredGLMM(m2)[[1]],4)



# plot --------------------------------------------------------------------------------------------------


p1 <- data.frame(pc$vectors[,1:2]) %>% 
  rownames_to_column(., var = "sample") %>% 
  mutate(sample = str_remove(sample, pattern = "_merged_abundance_rp_ks")) %>% 
  inner_join(., metadata, by = "sample") %>% 
  mutate(age = factor(age, levels = c( "4", "10", "24","control" )),
         age = recode(age, control="reference")) %>% 
  ggplot(., aes ( x = Axis.1, y = Axis.2, fill = age)) +
  geom_point(pch = 21, size = 5) +
  scale_fill_manual(name = "Reclamation Age (yr)",
                       values = c("#fde725", "#7fd250", "#346c8d", "#440053" ) ) +
  theme_bw()+
  theme(legend.position = "bottom",
        legend.title = element_text(size =16),
        legend.text = element_text(size =16),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16, color = "black")) +
  guides(fill = guide_legend(nrow = 1)) +
  labs(x = glue("Axis1 ({variance[1]} %)"),
       y = glue("Axis2 ({variance[2]} %)"),
       title = "Functional potential") +
  annotate(geom = "text", x =-0.2 , y = 0.18,
           label = c(glue("Reclamation age:
                           pseudoF = {age[1]}, pval < {age[2]}")))

p1 



p2 <- ggplot(bd, aes(x = age, y = sim)) +
  geom_point(pch = 21, size = 4,alpha = 1,  aes(fill = age),show.legend = F)+
  geom_smooth(method = "lm", se = T, color = "black", linewidth = 1.25) +
  scale_x_continuous(breaks = seq(0,30,by = 5)) +
  scale_y_continuous(labels = function(x) paste0(x, "%"))+
  scale_fill_viridis_c( direction = -1, limits = c(4, 31)) +
  theme_bw() +
  theme(legend.title = element_text(size =16),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16, color = "black")) +
  labs(x = "Reclamation age (yr)",
       y = "Similarity with reference site") +
  annotate(geom = "text", x = 9, y=20,
           label =  expression(paste(R[m]^2 == 0.643, "; pval < 0.001")), size =4,
           parse = T)

p2


func_plot <- p1 + p2 + plot_annotation(tag_levels = "A")+ plot_layout(guides = "collect")  & theme(legend.position = 'bottom')

func_plot

#ggsave(func_plot, filename = "../../plot/function.jpg", width = 12, height = 6.5, dpi = 1000)

