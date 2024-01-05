# Import library ----------------------------------------------------------------------------------------

library(tidyverse)
library(vegan)
library(lme4)
library(lmerTest)
library(patchwork)
library(easystats)



# Import data -------------------------------------------------------------------------------------------

meta <- read.csv("../../data/soil/meta1.csv") %>% 
  janitor::clean_names()

plant <- read.csv("../../data/plant/veg_combined.csv") %>% 
  janitor::clean_names() %>% 
  replace(is.na(.), 0) %>% 
  inner_join(., meta, by = c("site_code", "replicate")) %>% 
  mutate(row_name = paste0(site_code, replicate)) %>% 
  column_to_rownames(., var = "row_name") %>% 
  mutate(age = gsub(age, pattern = "#VALUE!", replacement = "reference"))


df_plant <- plant %>% 
  select(western_stickseed:timothy)
df_litter <- plant %>% 
  select(litter) %>% 
  rownames_to_column(., var = "id") 




# plant alpha div ---------------------------------------------------------

alpha_div <- data.frame(
  rich = specnumber(df_plant),
  shann = diversity(df_plant, index = "shannon"),
  simp = diversity(df_plant, index = "simpson")
) %>% 
  rownames_to_column(., var = "id")

div_df <- plant %>% 
  rownames_to_column(., var = "id") %>% 
  inner_join(., alpha_div, by = "id") %>% 
  filter(reclamation_method != "control") %>% 
  mutate(age = as.numeric(paste(age))) %>% 
droplevels()
str(div_df)

div_df$reclamation_method <- as.factor(div_df$reclamation_method) %>% 
  stats::relevel(div_df$reclamation_method, ref = "no_biosolids")


m1 <- glmer(rich ~ poly(age,2) + reclamation_method + (1|site_code), family = "poisson",
             data = div_df)
summary(m1)
car::Anova(m1)
emmeans::emmeans(m1, pairwise ~ reclamation_method)$emmeans

rich1 <- ggplot(div_df, aes(x = age, y = rich)) +
  geom_boxplot(data = div_df, aes(x = age, y = rich,fill = age, group = age), inherit.aes = F)+
  geom_point(aes(fill = age),pch = 21, size =2, width = 0.2,height = 0.01, show.legend = T) +
  geom_smooth(method= "lm", color = "black", se =T) +
  scale_fill_viridis_c( name = "Reclamation age (yr)",direction = -1) +
  theme_bw() +
  theme(legend.position = "bottom",
       legend.title = element_text(size =16),
       legend.text = element_text(size =16),
       legend.key.width = unit(1.25, 'cm'),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 16, color = "black")) +
  scale_x_continuous(breaks = c(seq(0, 30, by =5))) + 
  scale_y_continuous(breaks = c(seq(0, 10, by = 2))) +#,
 #                    labels = c("0","5","10","15","20","25","30")) +
  annotate(geom = "text", x = 7, y=9,
           label =  expression(paste(R[m]^2 == 0.131, "; pval < 0.05")), size =4,
           parse = T) +
  labs(x = NULL,
       y = "Plant Species richness") +
  expand_limits(y = c(0,10))
  #guides(fill = guide_legend(nrow = 2, override.aes = list(size=4))) 
rich1

rich2 <- ggplot(div_df, aes(x = reclamation_method, y = rich))+
  geom_violin(width = 0.3, outlier.alpha = 0,
              scale = "count",  alpha = 0.65) +
#  geom_jitter(width = 0.1, pch = 21, size =2.5, color = "#440154", height = 0L) +
  stat_summary(geom = "pointrange", fun.data = "mean_cl_boot", size = 0.85, linewidth = 1.5)+
  ggpubr::geom_bracket(xmin = "biosolids", xmax = "no_biosolids", y.position = 9,
                       label = "*") +
  scale_x_discrete(name = NULL,
                   labels = c("No\nbiosolids", "Biosolids")) +
  scale_y_continuous(name  = NULL, breaks = seq(0,10,by = 2)) +
  theme_bw() +
  theme(legend.title = element_text(size =16),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 16, color = "black")) +
  expand_limits(y = c(0,10))

rich2
plnt_rich <- rich1 + rich2+ plot_layout(guides = "collect")  & theme(legend.position = 'bottom')
plnt_rich 

#ggsave(plnt_rich, filename = "../../../plot/plant_richness.jpg", width = 10, height = 6, dpi = 1000)


# Litter ------------------------------------------------------------------------------------------------

m2 <- lmer(scale(litter) ~ poly(age,2) + reclamation_method + (1|site_code), 
           data = div_df)

summary(m2)
car::Anova(m2, test.statistic = "F")
MuMIn::r.squaredGLMM(m2)
emmeans::emmeans(m2, pairwise ~ reclamation_method)$emmeans



lit_plot <- ggplot(div_df, aes(x = age, y = litter)) +
  #geom_boxplot(data = div_df, aes(x = age, y = rich,fill = age, group = age), inherit.aes = F)+
  geom_point(aes(fill = age),pch = 21, size =3, width = 0.2,height = 0.01, show.legend = T) +
  #geom_smooth(method= "lm", color = "black", se =T) +
  geom_smooth(method= "lm",formula = y ~ x + I(x^2),
              color = "black", se =T) +
  scale_fill_viridis_c( name = "Reclamation age (yr)",direction = -1) +
  scale_x_continuous(breaks = c(seq(0, 30, by =5))) +
  scale_y_continuous(labels = scales::percent_format(scale = 1))+
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_text(size =16),
        legend.text = element_text(size =16),
        legend.key.width = unit(1.25, 'cm'),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 16, color = "black")) +
  labs(x = NULL,
       y = "Litter cover") +
  annotate(geom = "text", x = 6, y=100,
           label =  expression(paste(R[m]^2 == 0.317, "; pval < 0.001")), size =4)
lit_plot

#ggsave(lit_plot, filename = "../../plot/lit_plot.jpg", width = 7, height = 5, dpi = 1000)

