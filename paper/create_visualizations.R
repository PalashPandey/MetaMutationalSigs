setwd("C:\\Users\\pande\\OneDrive - Drexel University\\Documents\\Fall-2021\\Coop\\CGC\\SanjeeVCFFiles/PLOS_review_paper/metaSignatures/test_results")


library(ggplot2)
# library(tidyverse)
library(dplyr)
# remotes::install_github("tylermorganwall/rayshader")
# library(rayshader)

decon_legacy_df = as_tibble(read.csv("deconstructsigs_results/legacy_sample_errors.csv")) %>% 
  rename( sample = X,  error = x) %>%
  mutate(toolname = "deconstructsigs")

mut_legacy_df = as_tibble(read.csv("mutational_patterns_results/legacy_sample_errors.csv") )%>% 
  rename( sample = X,  error = x)%>%
  mutate(toolname = "mutationalPatterns")
  # %>%
  # mutate(sample__ = "")%>%
  # unite(sample , sample__, sample_, sep = "", remove = FALSE)%>% 
  # select(-c(sample_, sample__) )

sigfit_legacy_df = as_tibble(read.csv("sigfit_results/sample_errors_legacy.csv"))%>% 
  rename( sample = X,  error = errors_per_sample)%>%
  mutate(toolname = "sigfit")

sigflow_legacy_df = as_tibble(read.csv("sigflow/legacy_fitting_reconstruction_errors.csv"))%>%
  mutate(toolname = "sigflow")

final_legacy_df = rbind(rbind(rbind(decon_legacy_df, mut_legacy_df), sigfit_legacy_df), sigflow_legacy_df) %>% 
  mutate(mode = "Legacy SBS")



decon_sbs_df = as_tibble(read.csv("deconstructsigs_results/sbs_sample_errors.csv")) %>% 
  rename( sample = X,  error = x) %>%
  mutate(toolname = "deconstructsigs")

mut_sbs_df = as_tibble(read.csv("mutational_patterns_results/sample_errors.csv") )%>% 
  rename( sample = X,  error = x)%>%
  mutate(toolname = "mutationalPatterns")
  # %>%
  # mutate(sample__ = "")%>%
  # unite(sample , sample__, sample_, sep = "", remove = FALSE)%>% 
  # select(-c(sample_, sample__) )

sigfit_sbs_df = as_tibble(read.csv("sigfit_results/sample_errors_sbs.csv"))%>% 
  rename( sample = X,  error = errors_per_sample_sbs)%>%
  mutate(toolname = "sigfit")

sigflow_sbs_df = as_tibble(read.csv("sigflow/SBS_fitting_reconstruction_errors.csv"))%>%
  mutate(toolname = "sigflow")

final_sbs_df = rbind(rbind(rbind(decon_sbs_df, mut_sbs_df), sigfit_sbs_df), sigflow_sbs_df) %>%
  mutate(mode = "V3 SBS")


decon_dbs_df = as_tibble(read.csv("deconstructsigs_results/dbs_sample_errors.csv")) %>% 
  rename( sample = X,  error = x) %>%
  mutate(toolname = "deconstructsigs")

mut_dbs_df = as_tibble(read.csv("mutational_patterns_results/dbs_sample_errors.csv") )%>% 
  rename( sample = X,  error = x)%>%
  mutate(toolname = "mutationalPatterns")
  # %>%
  # mutate(sample__ = "")%>%
  # unite(sample , sample__, sample_, sep = "", remove = FALSE)%>% 
  # select(-c(sample_, sample__) )

sigfit_dbs_df = as_tibble(read.csv("sigfit_results/sample_errors_dbs.csv"))%>% 
  rename( sample = X,  error = errors_per_sample_dbs)%>%
  mutate(toolname = "sigfit")

sigflow_dbs_df = as_tibble(read.csv("sigflow/DBS_fitting_reconstruction_errors.csv"))%>%
  mutate(toolname = "sigflow")

final_dbs_df = rbind(rbind(rbind(decon_dbs_df, mut_dbs_df), sigfit_dbs_df), sigflow_dbs_df) %>%
  mutate(mode = "dbs")


decon_id_df = as_tibble(read.csv("deconstructsigs_results/indel_sample_errors.csv")) %>% 
  rename( sample = X,  error = x) %>%
  mutate(toolname = "deconstructsigs")

mut_id_df = as_tibble(read.csv("mutational_patterns_results/id_sample_errors.csv") )%>% 
  rename( sample = X,  error = x)%>%
  mutate(toolname = "mutationalPatterns")
# %>%
#   mutate(sample__ = "")%>%
#   unite(sample , sample__, sample_, sep = "", remove = FALSE)%>% 
#   select(-c(sample_, sample__) )

sigfit_id_df = as_tibble(read.csv("sigfit_results/sample_errors_indel.csv"))%>% 
  rename( sample = X,  error = errors_per_sample)%>%
  mutate(toolname = "sigfit")

sigflow_id_df = as_tibble(read.csv("sigflow/ID_fitting_reconstruction_errors.csv"))%>%
  mutate(toolname = "sigflow")

final_id_df = rbind(rbind(rbind(decon_id_df, mut_id_df), sigfit_id_df), sigflow_id_df) %>%
  mutate(mode = "ID")



final_df =  rbind(rbind(final_legacy_df, final_sbs_df) , final_id_df)
final_df[nrow(final_df) + 1,] = list(NA, 0, "sigflow", "DBS")

final_df_all_modes =  rbind(rbind(rbind(final_legacy_df, final_sbs_df) , final_id_df), final_dbs_df)
final_df_all_modes =  rbind(rbind(final_legacy_df, final_sbs_df) , final_dbs_df) 

setwd("C:\\Users\\pande\\OneDrive - Drexel University\\Documents\\Fall-2021\\Coop\\CGC\\SanjeeVCFFiles/PLOS_review_paper/metaSignatures/paper/")

# devtools::install_github("kassambara/ggpubr")
library(ggpubr)
library(ggbeeswarm)

ggplot(final_df, aes(x = toolname, y = error, color= mode )) +
  geom_boxplot(outlier.shape = NA ) +
  geom_beeswarm(groupOnX=FALSE)
  
unique(final_df$mode)

ggboxplot(final_df, x = "toolname",
                     y = "error",
                   combine = TRUE,
                     color = "mode", palette = "jco",
                     ylab = "error", 
                     add = "jitter",                              # Add jittered points
                     add.params = list(size = 0.1, jitter = 0.2)  # Point size and the amount of jittering
           )


ggplot(final_df_all_modes, aes(x = toolname, y = error, color = toolname, shape = mode )) +
  # geom_blank() + 
    # geom_violin(trim = FALSE) +
  geom_boxplot( ) +
  geom_jitter(position=position_jitter(0.2))


final_df_all_modes$mode

ggplot(final_df_all_modes, aes(x = toolname, y = error, color= mode )) +
  geom_boxplot( ) +
  geom_jitter(position=position_jitter(0.2)) 

ggsave("figure_d.pdf")

# + 
#   scale_colour_discrete(limits = c("Indel", "Legacy", "SBS", "DBS") , drop = FALSE )


# "circle", "triangle", "square", "star"



# ggplot(final_df, aes(x = toolname, y = error , color = mode)) +
#   geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size=4)
# 
# p + geom_jitter(shape=16, position=position_jitter(0.2))