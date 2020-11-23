#################
# Example script to reproduce 
# main figures from:
# Hertzog et al Model-based data integration
# problems can be reported by mailing:
# lionel.hertzog@thuenen.de
#############

# set working directory to the local file path to the repo
setwd("")

# load libraries
library(tidyverse)
library(purrr)
library(slider)
library(cvms)

# set path to write the figures
path_out <- "figures"

# load species names
mhb <- read.csv("data/example_mhb.csv")
sp <- unique(mhb$Artname)

# go through the species and load the files already there and plot it
keep_all <- NULL
cor_all <- NULL
trend_all <- NULL

for(species in sp){

  # collect files with rmse infos
  file_rmse <- list.files(path = paste0("model_output/", species), pattern = "rmse")
  file_rmse %>%
    map_df(~read.csv(paste0("model_output/", species, "/", .))) %>%
    mutate(model = rep(substr(file_rmse, 11, 11), each = 2)) -> rmse_dat
  
  # check if trim results is there and load
  rmse_trim <- tryCatch(read.csv(paste0("model_output/",species, "/r2_modeltrim.csv")),
                        error = NULL, warning = function(w) {})
  
  if(!is.null(rmse_trim)){
    rmse_trim$mad <- NA
    rmse_trim$low.ci <- NA
    rmse_trim$high.ci <- NA
    names(rmse_trim)[2] <- "median"
    rmse_trim$model <- "trim"
    rmse_dat <- rbind(rmse_dat, rmse_trim[,c(2:6, 1)])
  }
  
  rmse_dat$model <- factor(rmse_dat$model, 
                           levels = c("trim", as.character(1:5)),
                           labels = c("trim", "CBBS", "CBBS+ornitho.list",
                                      "CBBS+ornitho.all", "CBBS+all.list",
                                      "CBBS+all.list+all.unstr"))
  
  # load the totals
  file_total <- list.files(path = paste0("model_output/", species), pattern = "totals_model[1-6].csv")
  model_id <- substr(file_total, 13, 13)
  model_id <- c(rep(model_id[model_id == 1], 14),
                rep(model_id[model_id != 1], each = 15))
  
  file_total %>%
    map_df(~read.csv(paste0("model_output/", species, "/", .))) %>%
    mutate(model = model_id) -> total_dat
  
  # add trim
  total_trim <- tryCatch(read.csv(paste0("model_output/", species, "/totals_modeltrim.csv")),
                         error = NULL, warning = function(w) {})
  if(!is.null(total_trim)){
    total_trim$model <- "trim"
    names(total_trim) <- names(total_dat)
    total_dat <- rbind(total_dat, total_trim)
  }
  
  total_dat2 <- total_dat

  total_dat$model <- factor(total_dat$model, 
                            levels = c("trim", as.character(1:5)),
                            labels = c("trim", "CBBS", "CBBS+ornitho.list",
                                       "CBBS+ornitho.all", "CBBS+all.list","CBBS+all.list+all.unstr"))
  
  # first graph temporal trend based on the diff models
  gg_1 <- ggplot(total_dat, aes(x=Jahr, y=M, ymin=lo, ymax=hi)) +
    geom_vline(xintercept = 2012, linetype = "dashed") +
    geom_ribbon(alpha = 0.2, aes(fill = model)) +
    geom_ribbon(data = subset(total_dat, Jahr > 2011 & Jahr < 2019), alpha = 0.7, aes(fill = model)) +
    geom_linerange(color = "white") +
    geom_point(aes(color = model)) +
    facet_wrap(~model) +
    guides(color=FALSE,fill=FALSE)+
    theme_bw() +
    labs(x = "Year",
         y = "Predicted total farmland bird abundance with 95% C(r)I",
         title = species)
  
  ggsave(paste0(path_out, "/temporal_trend_", species, ".png"), gg_1)
  # 
  # for the second figure compute the precision from the totals
  total_dat %>%
    filter(Jahr > 2011 & Jahr < 2019) %>%
    mutate(precision = (hi - lo) / M) %>%
    group_by(model) %>%
    summarise(precision_m = mean(precision),
              precision_sd = sd(precision)) -> total_df
  
  # put rmse / accuracy and precision together
  rmse_dat %>%
    filter(metric == "rmse") %>%
    select(median, mad, model) %>%
    right_join(total_df, by = "model") -> df2
  
  df2$Species <- species
  # # save this for later use
  keep_all <- rbind(keep_all, df2)
  
  # third figure, get correlation between TRIM and the other models
  total_dat %>%
    filter(Jahr > 2011 & Jahr < 2019) %>%
    pivot_wider(id_cols = Jahr, names_from = model, values_from = M) %>%
    summarise(across(2:6, ~ cor(., trim, method = "spearman"))) %>%
    mutate(species = species) %>%
    pivot_longer(cols = 1:5, names_to = "model", values_to = "cor") -> cor_trim
  
 
  cor_all <- rbind(cor_all, cor_trim)
  
  # fourth, check the long-term trend for each model
  # first grab the 95% confidence band for the first year from TRIM
  total_dat %>%
    group_by(model) %>%
    summarise(trend = ifelse(M[14] < lo[1], "decreasing",
                             ifelse(M[14] > hi[1], "increasing", "stable"))) %>%
    mutate(species = species) -> tr
  
  trend_all <- rbind(trend_all, tr)
  
 
  print(species)
}

## figure 1: trend classification comparison

# using cvms
  trend_all %>%
  pivot_wider(id_cols = species, names_from = model, values_from = trend) %>%
  pivot_longer(3:7, names_to = "model") %>%
  mutate(trim = as.character(factor(trim, levels = c("decreasing", "stable", "increasing"),
                                    labels = c("c0","c1","c2"))),
         value = as.character(factor(value, levels = c("decreasing", "stable", "increasing"),
                                     labels = c("c0","c1","c2")))) %>%
  group_by(model) %>%
  evaluate(target_col = "trim", prediction_cols = "value", type = "multinomial") -> tt_dd

# the rows represent the different models
# and are somehow disorganized

# confusion matrix for model 2
m2<-plot_confusion_matrix(tt_dd[1,"Confusion Matrix"],
                          add_col_percentages = FALSE,
                          add_row_percentages = FALSE,
                          place_x_axis_above = FALSE)

m2.1 <- m2 +
  scale_x_discrete(labels=c("decreasing", "stable", "increasing"),
                   name = "Model 1") +
  scale_y_discrete(labels=c("decreasing", "stable", "increasing"),
                   name = "Model 2") +
  labs(title = "a)")

m5<-plot_confusion_matrix(tt_dd[2,"Confusion Matrix"],
                          add_col_percentages = FALSE,
                          add_row_percentages = FALSE,
                          place_x_axis_above = FALSE)

m5.1 <- m5 +
  scale_x_discrete(labels=c("decreasing", "stable", "increasing"),
                   name = "Model 1") +
  scale_y_discrete(labels=c("decreasing", "stable", "increasing"),
                   name = "Model 5") +
  labs(title = "d)")

m6<-plot_confusion_matrix(tt_dd[3,"Confusion Matrix"],
                          add_col_percentages = FALSE,
                          add_row_percentages = FALSE,
                          place_x_axis_above = FALSE)

m6.1 <- m5 +
  scale_x_discrete(labels=c("decreasing", "stable", "increasing"),
                   name = "Model 1") +
  scale_y_discrete(labels=c("decreasing", "stable", "increasing"),
                   name = "Model 6") +
  labs(title = "e)")

m4<-plot_confusion_matrix(tt_dd[4,"Confusion Matrix"],
                          add_col_percentages = FALSE,
                          add_row_percentages = FALSE,
                          place_x_axis_above = FALSE)

m4.1 <- m4 +
  scale_x_discrete(labels=c("decreasing", "stable", "increasing"),
                   name = "Model 1") +
  scale_y_discrete(labels=c("decreasing", "stable", "increasing"),
                   name = "Model 4") +
  labs(title = "c)")

m3<-plot_confusion_matrix(tt_dd[5,"Confusion Matrix"],
                          add_col_percentages = FALSE,
                          add_row_percentages = FALSE,
                          place_x_axis_above = FALSE)

m3.1 <- m3 +
  scale_x_discrete(labels=c("decreasing", "stable", "increasing"),
                   name = "Model 1") +
  scale_y_discrete(labels=c("decreasing", "stable", "increasing"),
                   name = "Model 3") +
  labs(title = "b)")

ma <- grid.arrange(m2.1, m3.1,m4.1,m5.1, m6.1,ncol=2)
ggsave(paste0(path_out, "/figures/confusion_matrix.png"), ma)


# figure 2: correlation between the m1-m6 and TRIM across species
cor_all$model <- factor(cor_all$model,
                        levels = c("CBBS", "CBBS+ornitho.list",
                                   "CBBS+ornitho.all", "CBBS+all.list", "CBBS+all.list+all.unstr"),
                        labels = paste("Model", 2:6, sep=" "))


gg_cor <- ggplot(cor_all, aes(x=model, y=cor, color = model)) +
  geom_jitter(width = 0.1, height = 0) +
  stat_summary(fun.data = "mean_cl_boot", color = "red") +
  labs(y = "Correlation coefficient",
       x = "") +
  guides(color = FALSE) +
  theme(text = element_text(size=25))

ggsave(paste0(path_out, "figures/correlation_with_trim.png"), gg_cor)


# figure 3: compute percent change in precision / accuracy from trim

# first for accuracy  
keep_all %>%
  pivot_wider(id_cols = "Species", names_from = "model", 
              values_from = "median") %>%
  mutate(across(2:6,  ~ ((trim - .) / abs(trim)) * 100)) %>%
  pivot_longer(cols = 2:6, names_to = "model",
               values_to = "pchange") -> acc_dd

acc_dd$model <- factor(acc_dd$model,
                       levels = c("CBBS", "CBBS+ornitho.list",
                                  "CBBS+ornitho.all", "CBBS+all.list",
                                  "CBBS+all.list+all.unstr"),
                       labels = paste("Model", 2:6, sep=" "))


# now for precision
keep_all %>%
  pivot_wider(id_cols = "Species", names_from = "model", 
              values_from = "precision_m") %>%
  mutate(across(2:6,  ~ ((trim - .) / abs(trim)) * 100)) %>%
  pivot_longer(cols = 2:6, names_to = "model",
               values_to = "pchange") -> pre_dd


pre_dd$model <- factor(pre_dd$model,
                       levels = c("CBBS", "CBBS+ornitho.list",
                                  "CBBS+ornitho.all", "CBBS+all.list",
                                  "CBBS+all.list+all.unstr"),
                       labels = paste("Model", 2:6, sep=" "))

# now put the two together for plotting 
names(acc_dd)[4] <- "accuracy"
names(pre_dd)[4] <- "precision"

all <- cbind(acc_dd[,c(1, 3, 4)], pre_dd[,4])
all_d <- pivot_longer(all, 3:4, names_to = "type", values_to = "pchange")

# remove some weird results by m5
# all_d <- all_d[all_d$pchange > -20 & !is.na(all_d$pchange),]

gg_pa <- ggplot(all_d, aes(x=model, y=pchange, color = model)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_jitter(width = 0.1) +
  facet_grid(~type) +
  guides(color=FALSE)+
  stat_summary(fun.data = "mean_cl_boot",color = "red") +
  labs(y = "Change compared to model 1 (%)",
       x="") +
  theme(text = element_text(size=20))

ggsave(paste0(path_out, "figures/accuracy_precision.png"), gg_pa)

