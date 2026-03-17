# ROI SIM CODE
library(dplyr)
library(tidyr)
library(ggplot2)
RNGkind("L'Ecuyer-CMRG")
library(data.table)
library(Cardinal)
library(broom)
source("Utilities.R")
source("SimulationUtilities.R")
##### Helper Functions #####


##### Generate base images ##### 

baseline <- 80
set.seed(2025)
pixel_sd_base <- 30
pixel_sd_base <- 50
pixel_sd <- rnorm(3000, mean = pixel_sd_base, sd = 7) + 5

# set.seed(2025)
# pixel_sd_base <- 40
# pixel_sd <- rnbinom(3000, mu = pixel_sd_base, size = 4)
# pixel_sd <- -pixel_sd + max(pixel_sd) + 5
# hist(pixel_sd)
# mean(pixel_sd)


if (any(pixel_sd < 5)) {
  warning("Some pixel_sd values are below 5", call. = F, immediate. = T)
}

pdx <- PositionDataFrame(expand_grid(x = 1:100, y = 1:100))
pdx$all <- TRUE
fdx <- MassDataFrame(1:3000)
fdx$all <- baseline


samples <- data.frame(subject = c(1:8,1:8),
                      subtissue = c(rep(0,8),rep(1,8)),
                      condition = c(rep(0,4),rep(1,4),rep(0,4),rep(1,4))) %>%
  mutate(id = paste0(subject,"_",subtissue,"_",condition))


load("Simulation-Data/sim_list_original.RData")

##### Add in effects and run models #####

delta_in_roi <- 30
gamma <- sqrt(35)
resid_gamma <- sqrt(35)
resid <- sqrt(70)


ROI_def_delta <- 40
set.seed(2024)
effects_normal_subtissue <- tibble(mz = mz(sim_list_original[[1]]),
                  subtissue = c(rep(0,2250), rep(delta_in_roi, 250), rep(0, 250), rep(delta_in_roi, 250)),
                  #condition = c(rep(0,2500), rep(delta_in_roi,500)),
                  condition = c(rep(0,2500), rep(-1*delta_in_roi,250), rep(delta_in_roi,250)),
                  ROI_defining_features = c(rep(ROI_def_delta,50),rep(-ROI_def_delta,50), rep(0,2900)))

set.seed(2024)
effects_shuffled_subtissue <- tibble(mz = mz(sim_list_original[[1]]),
                  subtissue = c(rep(0,100), sample(rep(c(0, delta_in_roi),2900/2))),
                  condition = c(rep(0,2500), rep(delta_in_roi,500)),
                  ROI_defining_features = c(rep(ROI_def_delta,50),rep(-ROI_def_delta,50), rep(0,2900)))

set.seed(2025)
gamma_mat_cad <- matrix(rnorm(4*nrow(effects_normal_subtissue), sd = gamma), nrow = 4)
gamma_mat_oa <- matrix(rnorm(4*nrow(effects_normal_subtissue), sd = gamma), nrow = 4)
gamma_mat_cad <- scale(gamma_mat_cad, center = TRUE, scale = FALSE)
gamma_mat_oa <- scale(gamma_mat_oa, center = TRUE, scale = FALSE)
gamma_mat <- rbind(gamma_mat_oa, gamma_mat_cad, gamma_mat_oa, gamma_mat_cad)

set.seed(2026)
residual_mat <- matrix(rnorm(16*nrow(effects_normal_subtissue), sd = resid), nrow = 16)
residual_mat <- scale(residual_mat, center = TRUE, scale = FALSE)

set.seed(2026)
residual_mat_gamma <- matrix(rnorm(16*nrow(effects_normal_subtissue), sd = resid_gamma), nrow = 16)
residual_mat_gamma <- scale(residual_mat_gamma, center = TRUE, scale = FALSE)



#### effects without gamma ####
effects_matrix_c2 <- sapply(1:16, function(i) {
  #samples[[i, "subtissue"]]*effects_normal_subtissue$subtissue + 
    samples[[i, "condition"]]*effects_normal_subtissue$condition +
    1*effects_normal_subtissue$ROI_defining_features +
    #gamma_mat[i,] +
    residual_mat[i,] +
    baseline
})

effects_matrix_c2_subtissue_normal <- sapply(1:16, function(i) {
  samples[[i, "subtissue"]]*effects_normal_subtissue$subtissue + 
  samples[[i, "condition"]]*effects_normal_subtissue$condition +
    1*effects_normal_subtissue$ROI_defining_features +
    #gamma_mat[i,] +
    residual_mat[i,] +
    baseline
})

effects_matrix_c2_subtissue_shuffled <- sapply(1:16, function(i) {
  samples[[i, "subtissue"]]*effects_shuffled_subtissue$subtissue + 
  samples[[i, "condition"]]*effects_shuffled_subtissue$condition +
    1*effects_shuffled_subtissue$ROI_defining_features +
    #gamma_mat[i,] +
    residual_mat[i,] +
    baseline
})

effects_matrix_bg <- sapply(1:16, function(i) {
  -1*effects_normal_subtissue$ROI_defining_features +
    #gamma_mat[i,] +
    residual_mat[i,] +
    baseline
})

i <- 1
for (eff in c(effects_matrix_bg, 
              effects_matrix_c2_subtissue_shuffled, 
              effects_matrix_c2_subtissue_normal, 
              effects_matrix_c2)) {
  
  if (any(eff < 0)){
    stop(paste0("Some NON GAMMA effects in ", i, " have negative means."))
  }
  i <- i + 1
}

#### No subtissue, no gamma ####

sim_list_adj <- list()
gc(verbose = T)
sim_list_adj <- bplapply(seq_along(sim_list_original), FUN = function(i){
  m <- sim_list_original[[i]]
  runNames(m) <- paste0("run_", i)
  
  intensity(m) <- scale_row_sd_mean_subset(as.matrix(intensity(m)),
                                           target_sd = pixel_sd,
                                           target_mean = effects_matrix_bg[,i])
  
  
  pData(m) <- addShape(addShape(pData(m), 
                                c(50,50), 20, name = "circle_1"), 
                       c(50,50), 40, name = "circle_2")
  
  pData(m)$circle_2 <- ifelse(pData(m)$circle_2 & pData(m)$circle_1,
                              FALSE,
                              pData(m)$circle_2)
  
  pData(m)$bg <- ifelse(pData(m)$circle_2,
                        FALSE,
                        TRUE)
  
  intensity(m) <- scale_row_sd_mean_subset(as.matrix(intensity(m)), 
                                           target_sd = pixel_sd,
                                           columns_to_scale = pData(m)$circle_2,
                                           target_mean = effects_matrix_c2[,i])
  
  summarizeFeatures(summarizeFeatures(m, stat = c("mean", "sd"), na.rm = T, verbose = F, groups = ifelse(pData(m)$circle_2, 
                                                                                                         "oracle_donut", 
                                                                                                         "other")), 
                    stat = c("mean", "sd"), na.rm = T, verbose = F)
}, BPPARAM = MulticoreParam(4))


gc(verbose = T)
sscl <- bplapply(sim_list_adj, FUN = function(m){
  set.seed(123456789, kind = "L'Ecuyer-CMRG")
  spatialShrunkenCentroids(m, k = 4, s = c(6,9,12,15), weights = "adaptive")
}, BPPARAM = MulticoreParam(4))

results_list <- analyze_spatial_clusters(sscl,
                                         sim_list_adj, 
                                         samples, 
                                         n_top_per_condition = 4,
                                         n_top_features = 20)

mse_frame_to_test <- extract_data_for_modeling_alt(results_list[["sim_to_test"]], results_list[["samples_to_test"]])

run_models_to_test_ssc_roi<- run_models_no_me(mse_frame_to_test, 
                                              formu = "ssc_roi.mean ~ condition")

run_models_to_test_oracle_roi<- run_models_no_me(mse_frame_to_test, 
                                                 formu = "oracle_donut.mean ~ condition")


gc(verbose = T)
skml_7c <- bplapply(sim_list_adj, FUN = function(m){
  set.seed(123456789, kind = "L'Ecuyer-CMRG")
  spatialKMeans(m, k = 7, weights = "adaptive")
}, BPPARAM = MulticoreParam(4))

results_list_skm_7c <- analyze_spatial_clusters_skm(skml_7c,
                                                    sim_list_adj, 
                                                    samples, 
                                                    n_top_per_condition = 4,
                                                    n_top_features = 20)

mse_frame_to_test_skm_7c <- extract_data_for_modeling_alt(results_list_skm_7c[["sim_to_test"]], results_list_skm_7c[["samples_to_test"]])

run_models_to_test_skm_roi_7c<- run_models_no_me(mse_frame_to_test_skm_7c, 
                                                 formu = "skm_roi.mean ~ condition")

results_list[[3]]

results_list_skm_7c[[3]]

intersect(results_list[[3]] %>% filter(condition == 1) %>%  pull(id), 
                                       results_list_skm_7c[[3]] %>% filter(condition == 1) %>%  pull(id))

#### for figures ####

trt_sample_id <- intersect(results_list[[3]] %>% filter(condition == 1) %>%  pull(id), 
          results_list_skm_7c[[3]] %>% filter(condition == 1) %>%  pull(id))[1]

control_sample_id <- intersect(results_list[[3]] %>% filter(condition == 0) %>%  pull(id), 
                            results_list_skm_7c[[3]] %>% filter(condition == 0) %>%  pull(id))[1]

m <- cbind(sim_list_adj[[which(samples$id == control_sample_id)]], 
           sim_list_adj[[which(samples$id == trt_sample_id)]])
png(filename = "ROI_sim_sii_ground_truth_roi_defining.png", width = 600, height = 1200)
print(image(m,
      i = c(1),
      layout = c(2,1),
      scale = T,
      superpose = F,
      key = F,
      grid = F,
      xaxt = "n",
      yaxt = "n",
      ylab = NA,
      xlab = NA))
dev.off()


png(filename = "ROI_sim_sii_ground_truth_no_eff.png", width = 600, height = 1200)
print(image(m,
      i = c(1000),
      layout = c(2,1),
      scale = T,
      superpose = F,
      key = F,
      grid = F,
      xaxt = "n",
      yaxt = "n",
      ylab = NA,
      xlab = NA))
dev.off()

png(filename = "ROI_sim_sii_ground_truth_diff_abun.png", width = 600, height = 1200)
print(image(m,
      i = c(2700),
      layout = c(2,1),
      scale = T,
      superpose = F,
      key = F,
      grid = F,
      xaxt = "n",
      yaxt = "n",
      ylab = NA,
      xlab = NA))
dev.off()


png(filename = "ROI_sim_sii_segmentation_oracle.png", width = 600, height = 1200)
print(image(m,
      "circle_2",
      layout = c(2,1),
      scale = T,
      superpose = F,
      key = F,
      grid = F,
      xaxt = "n",
      yaxt = "n",
      ylab = NA,
      xlab = NA))
dev.off()

m <- cbind(results_list[[2]][[which(results_list[[3]]$id == control_sample_id)]], 
           results_list[[2]][[which(results_list[[3]]$id == trt_sample_id)]])
png(filename = "ROI_sim_sii_segmentation_ssc.png", width = 600, height = 1200)
print(image(m,
      "roi_ssc",
      layout = c(2,1),
      scale = T,
      superpose = F,
      key = F,
      grid = F,
      xaxt = "n",
      yaxt = "n",
      ylab = NA,
      xlab = NA))
dev.off()

m <- cbind(results_list_skm_7c[[2]][[which(results_list_skm_7c[[3]]$id == control_sample_id)]], 
           results_list_skm_7c[[2]][[which(results_list_skm_7c[[3]]$id == trt_sample_id)]])
png(filename = "ROI_sim_sii_segmentation_skm.png", width = 600, height = 1200)
print(image(m,
      "roi_skm",
      layout = c(2,1),
      scale = T,
      superpose = F,
      key = F,
      grid = F,
      xaxt = "n",
      yaxt = "n",
      ylab = NA,
      xlab = NA))
dev.off()


##### charting #####
res_labels <- tibble(
  isPos = c(TRUE, FALSE, TRUE, FALSE),
  truePos = c(TRUE, FALSE, FALSE, TRUE),
  label = c("True Positive", "True Negative", "False Positive", "False Negative")
)

#ssc 
accuracy_frame_ssc <- run_models_to_test_ssc_roi[[2]] %>% 
  mutate(isPos = p.adjust( p.value_condition1, method = "fdr") < 0.05,
         truePos = feature_id > 2500) %>% 
  group_by(isPos, truePos) %>% 
  count() %>% 
  ungroup() %>% 
  full_join(res_labels) %>% 
  mutate(n = replace_na(n,0))

accuracy_frame_ssc %>% 
  filter(label %in% c("True Negative", "False Positive")) %>% 
  mutate(rate = n/sum(n)) %>% 
  filter(label == "False Positive")

#oracle

accuracy_frame_oracle <- run_models_to_test_oracle_roi[[2]] %>% 
  mutate(isPos = p.adjust( p.value_condition1, method = "fdr") < 0.05,
         truePos = feature_id > 2500) %>% 
  group_by(isPos, truePos) %>% 
  count()%>% 
  ungroup() %>% 
  full_join(res_labels) %>%
  mutate(n = replace_na(n,0))

accuracy_frame_oracle %>% 
  filter(label %in% c("True Negative", "False Positive")) %>% 
  mutate(rate = n/sum(n)) %>% 
  filter(label == "False Positive")

#skm
accuracy_frame_skm <- run_models_to_test_skm_roi_7c[[2]] %>% 
  mutate(isPos = p.adjust( p.value_condition1, method = "fdr") < 0.05,
         truePos = feature_id > 2500) %>% 
  group_by(isPos, truePos) %>% 
  count()%>% 
  ungroup() %>% 
  full_join(res_labels) %>%
  mutate(n = replace_na(n,0))

accuracy_frame_skm %>% 
  filter(label %in% c("True Negative", "False Positive")) %>% 
  mutate(rate = n/sum(n)) %>% 
  filter(label == "False Positive")


custom_colors <- c("#40B0A6", "#E1BE6A")  


bind_rows(run_models_to_test_oracle_roi[[2]] %>% 
            mutate(fdr.p.value_condition1 = p.adjust( p.value_condition1, method = "fdr"),
                   isPosAdj = fdr.p.value_condition1 < 0.05,
                   truePos = feature_id > 2500,
                   id = "oracle"),
          run_models_to_test_ssc_roi[[2]] %>% 
            mutate(fdr.p.value_condition1 = p.adjust( p.value_condition1, method = "fdr"),
                   isPosAdj = fdr.p.value_condition1 < 0.05,
                   truePos = feature_id > 2500,
                   id = "ssc"),
          run_models_to_test_skm_roi_7c[[2]] %>% 
            mutate(fdr.p.value_condition1 = p.adjust( p.value_condition1, method = "fdr"),
                   isPosAdj = fdr.p.value_condition1 < 0.05,
                   truePos = feature_id > 2500,
                   id = "skm")) %>% 
  mutate(fdr_finding = ifelse(isPosAdj,
                         ifelse(truePos,"True Positive", "False Positive"),
                         ifelse(truePos,"False Negative", "True Negative")),
         fdr_finding_pos_only = factor(ifelse(isPosAdj,
                              ifelse(truePos,"True Positive", "False Positive"),
                              "True/False Negative"),
                              levels = c("True Positive", "False Positive",
                              "True/False Negative")),
         id = factor(id, levels = c("oracle", "ssc", "skm"))) %>% 
  ggplot(aes(x = estimate_condition1, y = -log(fdr.p.value_condition1), color = fdr_finding_pos_only)) +
  geom_point(alpha = 0.65) +
  facet_grid(.~ id) +
  scale_color_manual(values = c("False Positive" = "#E1BE6A", "True Positive" = "#40B0A6", 
                                "True/False Negative" = "grey70"),
                     breaks = c("False Positive", "True Positive")) +
  labs(x = "Estimated difference of simulated condition means in ROI",
       y = "-Log(FDR Adjusted P Value)",
       color = "Finding") +
  theme_bw() +
  theme(legend.position = "bottom")
ggsave(file.path("Recreated_Figures","Figure-ROI-Simulation-1-Volcanos-by-Method.png"), height = 6, width = 10, units = "in")

##### other stuff ####



oracle_pos_features <- run_models_to_test_oracle_roi[[2]] %>% 
  mutate(isPos = p.value_condition1 < 0.05,
         truePos = feature_id > 2500,
         fdr_adj_pv = p.adjust(p.value_condition1, method = "fdr"),
         isPosAdj = fdr_adj_pv < 0.05) %>% 
  filter(isPosAdj) %>% pull(feature_id)

ssc_pos_features <- run_models_to_test_ssc_roi[[2]] %>% 
  mutate(isPos = p.value_condition1 < 0.05,
         truePos = feature_id > 2500,
         fdr_adj_pv = p.adjust(p.value_condition1, method = "fdr"),
         isPosAdj = fdr_adj_pv < 0.05) %>% 
  filter(isPosAdj) %>% filter(feature_id %in% as.integer(unique_top_20_all_samples)) %>% pull(feature_id)


library(VennDiagram)


library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")

ssc_control_charting_sample <- sscl[[which(samples$id == control_sample_id)]][[results_list[[1]][[which(results_list[[3]]$id == control_sample_id)]]$ssc_index]]

ssc_trt_charting_sample <- sscl[[which(samples$id == trt_sample_id)]][[results_list[[1]][[which(results_list[[3]]$id == trt_sample_id)]]$ssc_index]]


# Chart
venn.diagram(
  x = list(topFeatures(ssc_control_charting_sample, n = 500)$i, 
           topFeatures(ssc_trt_charting_sample, n = 500)$i, 2501:3000),
  category.names = c("Top 500 Features\nfrom a control image" , 
                     "Top 500 Features\nfrom a treatment image" , 
                     "True differentiating features"),
  filename = file.path("Recreated_Figures","Figure-ROI-Simulation-1-SSC-Venn-Shared-Features.png"),
  disable.logging = TRUE,
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 1000 , 
  width = 1000 , 
  resolution = 300,
  compression = "lzw",
  pointsize = 18,
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = brewer.pal(3, "Pastel2"),
  
  euler.d = T, sep.dist = 2,
  
  # Numbers
  # cex = .6,
  # fontface = "bold",
  fontfamily = "sans",
  
  # # Set names
  cat.cex = 0.6,
  # cat.fontface = "bold",
  # cat.default.pos = "outer",
  # cat.pos = c(-27, 27, 135),
  # cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  # rotation = 1,
  margin = 0.1
)

skm_control_charting_sample <- skml_7c[[which(samples$id == control_sample_id)]]

skm_trt_charting_sample <- skml_7c[[which(samples$id == trt_sample_id)]]

# Chart
venn.diagram(
  x = list(topFeatures(skm_control_charting_sample, n = 500)$i, 
           topFeatures(skm_trt_charting_sample, n = 500)$i, 2501:3000),
  category.names = c("Top 500 Features\nfrom a control image" , 
                     "Top 500 Features\nfrom a treatment image" , 
                     "True differentiating features"),
  filename = file.path("Recreated_Figures","Figure-ROI-Simulation-1-SKM-Venn-Shared-Features.png"),
  disable.logging = TRUE,
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 1000 , 
  width = 1000 , 
  resolution = 300,
  compression = "lzw",
  pointsize = 18,
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = brewer.pal(3, "Pastel2"),
  
  euler.d = T, sep.dist = 2,
  
  # Numbers
  # cex = .6,
  # fontface = "bold",
  fontfamily = "sans",
  
  # # Set names
  cat.cex = 0.6,
  # cat.fontface = "bold",
  # cat.default.pos = "outer",
  # cat.pos = c(-27, 27, 135),
  # cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  # rotation = 1,
  margin = 0.1
)


