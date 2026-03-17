library(matrixStats)
library(broom.mixed)
library(emmeans)
library(ggpubr)
library(lmerTest)
library(dplyr)
library(tidyr)
library(ggplot2)
library(Cardinal)
RNGkind(kind = "L'Ecuyer-CMRG")
source("Utilities.R")
source("SimulationUtilities.R")

##### Helper functions ##### 


##### Generate base images ##### 

baseline <- 50
pixel_sd <- 45
pdx <- PositionDataFrame(expand_grid(x = 1:100, y = 1:100))
pdx$all <- TRUE
fdx <- MassDataFrame(1:3000)
fdx$all <- baseline


samples <- data.frame(subject = c(1:8,1:8),
                      subtissue = c(rep(0,8),rep(1,8)),
                      condition = c(rep(0,4),rep(1,4),rep(0,4),rep(1,4))) %>%
  mutate(id = paste0(subject,"_",subtissue,"_",condition))

set.seed(123, kind = "L'Ecuyer-CMRG")

# sim_list_original <- list()
# for (i in 1:nrow(samples)) {
#   sim_list_original[[i]] <- simulateImage(pdx,
#                                           fdx, 
#                                           sdpixel = pixel_sd, 
#                                           sdrun = 0, 
#                                           sdmz = 0, 
#                                           sdmult = 0, 
#                                           sdpeaks = 0, 
#                                           sdnoise = 0.0, 
#                                           SAR = FALSE, 
#                                           baseline = 0,
#                                           spcorr = 0.5,
#                                           centroided = TRUE, 
#                                           verbose = TRUE,
#                                           BPPARAM=MulticoreParam())
# }
# 
# save(sim_list_original, file = "sim_list_original.RData")
sim_list_original_path <- "Simulation-Data/sim_list_original.RData"
stopifnot(file.exists(sim_list_original_path))
load(sim_list_original_path)
stopifnot(exists("sim_list_original"))

##### Add in effects and run models #####

delta <- 15
gamma <- sqrt(20)
resid <- sqrt(120)

effects <- tibble(mz = mz(sim_list_original[[1]]),
                  subtissue = c(rep(0,1500), rep(delta,1500)),
                  condition = c(rep(0,1500), rep(delta,1500)))

set.seed(2025)
gamma_mat_cad <- matrix(rnorm(4*nrow(effects), sd = gamma), nrow = 4)
gamma_mat_oa <- matrix(rnorm(4*nrow(effects), sd = gamma), nrow = 4)
gamma_mat_cad <- scale(gamma_mat_cad, center = TRUE, scale = FALSE)
gamma_mat_oa <- scale(gamma_mat_oa, center = TRUE, scale = FALSE)
gamma_mat <- rbind(gamma_mat_oa, gamma_mat_cad, gamma_mat_oa, gamma_mat_cad)

residual_mat <- matrix(rnorm(16*nrow(effects), sd = resid), nrow = 16)
residual_mat <- scale(residual_mat, center = TRUE, scale = FALSE)

effects_matrix <- sapply(1:16, function(i) {
  samples[[i, "subtissue"]]*effects$subtissue + 
    samples[[i, "condition"]]*effects$condition +
    gamma_mat[i,] +
    residual_mat[i,] +
    baseline
})

if (any(effects_matrix < 0)) {
  stop("Some effects have negative means. Increase baseline.")
}

sim_list_adj <- list()
for (i in 1:length(sim_list_original)) {
  print(paste0("Starting ", i))
  m <- sim_list_original[[i]]
  intensity(m) <- scale_row_sd_mean(as.matrix(intensity(m)), 
                                                         target_sd = pixel_sd, 
                                                         target_mean = effects_matrix[,i])
  sim_list_adj[[i]] <- summarizeFeatures(m, stat = c("mean", "sd"), na.rm = T, verbose = F)
}

mse_frame_sim_list <- extract_data_for_modeling_alt(sim_list_adj, samples)

# run_models_sim_list <- run_models(mse_frame_sim_list, 
#                                    formu = "mean ~ subtissue * condition + (1|subject)")
run_models_sim_list <- run_models_MEMORY_EFF(mse_frame_sim_list, 
                                              formu = "mean ~ subtissue * condition + (1|subject)")

##### Investigate and classify models by variation ##### 

feature_models_results_for_each_combo <- run_models_sim_list[[2]] %>% 
  mutate(sub_var = total_var * icc, 
         tech_var = total_var * (1-icc), 
         eff = as.numeric(feature_id) > 1500) %>% 
  dplyr::select(feature_id, icc, eff, total_var, sub_var, tech_var) %>%
  filter(icc > 0 & icc < 1 & total_var > 40 & total_var < 505) %>%
  group_by(eff, icc_group = icc < 0.5) %>%
  mutate(target_icc = if_else(icc_group, 0.25, 0.75),
         distance = abs(icc - target_icc)) %>%
  arrange(eff, icc_group, distance) %>%
  group_by(eff, icc_group) %>%
  slice_head(n = 150) %>%
  mutate(icc_label = if_else(icc_group, "low icc", "high icc")) %>%
  ungroup()

feature_models_results_for_each_combo %>% group_by(eff, icc_label) %>% summarise(median_icc = median(icc),
                                                                                 median_total = median(sub_var + tech_var),
                                                                                 median_tech = median(tech_var),
                                                                                 median_sub = median(sub_var))
feature_models_results_for_each_combo %>% ggplot(aes(x = total_var, fill = icc_label)) + geom_histogram(position = "identity", bins = 50, alpha = 0.4)

feature_models_results_for_each_combo %>% ggplot(aes(x = tech_var, fill = icc_label)) + geom_histogram(position = "identity", bins = 50, alpha = 0.4)
feature_models_results_for_each_combo %>% ggplot(aes(x = sub_var, fill = icc_label)) + geom_histogram(position = "identity", bins = 50, alpha = 0.4)



feature_models_results_for_each_combo %>%
  pivot_longer(total_var:tech_var, 
               names_to = "type", 
               names_pattern = "^(.*)_",
               values_to = "vari") %>%
  ggplot(aes(x = vari, fill = icc<0.5)) +
  geom_histogram(bins = 100, position = "identity", alpha = 0.4) +
  facet_wrap("type", ncol = 1) +
  xlim(-1,500)

##### Subset Tests #####

n_each <- 150

feature_ids_for_lowb_hight <- feature_models_results_for_each_combo %>% 
                                 filter(icc_label == "low icc") %>%
                                 arrange(feature_id) %>%
                                 pull(feature_id)

j_frame_lowb_hight <- mse_frame_sim_list %>% 
  filter(feature_id %in% feature_ids_for_lowb_hight)

# confirm there are exactly 150 features in each eff
j_frame_lowb_hight %>% 
  group_by(feature_id) %>% 
  summarize(n = n()) %>% 
  ungroup() %>% 
  mutate(eff = feature_id>1500) %>% 
  group_by(eff) %>% 
  summarize(n=n())

# this will essentilly reset the indicies - but thats ok because so long as there as 150 of each eff type and
# models are pulled in order of inc feature_id, the first 150 will be zero effect.
run_models_lowb_hight <- run_models_sim_list[[1]][feature_models_results_for_each_combo %>%
                                                    filter(icc_label == "low icc") %>% 
                                                    arrange(feature_id) %>%
                                                    pull(feature_id)]

feature_ids_for_highb_lowt <- feature_models_results_for_each_combo %>% 
  filter(icc_label == "high icc") %>% 
  arrange(feature_id) %>%
  pull(feature_id)


j_frame_highb_lowt <- mse_frame_sim_list %>% 
  filter(feature_id %in% feature_ids_for_highb_lowt)

j_frame_highb_lowt %>% 
  group_by(feature_id) %>% 
  summarize(n = n()) %>% 
  ungroup() %>% 
  mutate(eff = feature_id>1500) %>% 
  group_by(eff) %>% 
  summarize(n=n())

run_models_highb_lowt <- run_models_sim_list[[1]][feature_models_results_for_each_combo %>%
                                                    filter(icc_label == "high icc") %>% 
                                                    arrange(feature_id) %>%
                                                    pull(feature_id)]


tests_lowb_hight <- get_comparison_plotting_frame(j_frame_lowb_hight,run_models_lowb_hight, nmax = n_each)
tests_highb_lowt <- get_comparison_plotting_frame(j_frame_highb_lowt,run_models_highb_lowt, nmax = n_each)

##### Variation Fig for paper #####

top_labs <- c("Greater Technical Variation",
              "Greater Biological Variation")


add_zero_tests <- function(tests_df, nmax = 150) {
  if(nrow(tests_df %>% filter(n == nmax & finding %in% c("False Positive","True Negative")))>0) {
    tests_df <- bind_rows(tests_df,
                          tests_df %>% 
                            filter(n == nmax & finding %in% c("False Positive","True Negative")) %>% 
                            mutate(finding = ifelse(finding == "False Positive",
                                                    "True Negative",
                                                    "False Positive"),
                                   n = 0))
  }
  
  if(nrow(tests_df %>% filter(n == nmax & finding %in% c("True Positive","False Negative")))>0) {
    tests_df <- bind_rows(tests_df,
                          tests_df %>% 
                            filter(n == nmax & finding %in% c("True Positive","False Negative")) %>% 
                            mutate(finding = ifelse(finding == "True Positive",
                                                    "False Negative",
                                                    "True Positive"),
                                   n = 0))
  }
  
  return(tests_df)
  
}


df_test_vari <- bind_rows(add_zero_tests(tests_lowb_hight, nmax = n_each) %>% mutate(vari = top_labs[1]),
                          add_zero_tests(tests_highb_lowt, nmax = n_each) %>% mutate(vari = top_labs[2])) %>%
  mutate(vari = factor(vari, levels = top_labs)) %>%
  mutate(model_simple = factor(ifelse(test=="T Test of Conditional Means",
                                      "Mixed Effects Model",
                                      test),
                               levels = c("Mixed Effects Model",
                                          "Two Sample T Test",
                                          "Paired T Test"),
                               labels =c("Table 1 Model 3:\nMixed Effects Model\nwith both variables\n",
                                         "Table 1 Model 1:\nEquivalent to a\nTwo Sample T Test\nin this design",
                                         "Table 1 Model 2:\nEquivalent to a\nPaired T Test\nin this design")),
         model_family = factor(ifelse(test=="T Test of Conditional Means",
                                      "Mixed Effects Models",
                                      "Sample T Tests"),
                               levels = c("Mixed Effects Models",
                                          "Sample T Tests")),
         comparison_type = ifelse(term=="Subsection",
                                  "Within Subjects",
                                  "Between Subjects"),
         null = factor(null,
                       levels = c("In Medial:\nControl vs Osteoarthritis", "In Osteoarthritis:\nLateral vs Medial"),
                       labels =c("Hypothesis A, Between-subjects\nIn Medial: Osteoarthritis vs Control",
                                 "Hypothesis B, Within-subjects\nIn Osteoarthritis: Medial vs Lateral")))


itx_df <-bind_rows(j_frame_lowb_hight %>% mutate(vari = top_labs[1]), 
                   j_frame_highb_lowt %>% mutate(vari = top_labs[2])) %>%
  mutate(vari = factor(vari, levels = top_labs)) %>%
  dplyr::rename("Subsection"="subtissue", "Condition"="condition") %>%
  group_by(vari) %>%
  filter(feature_id == sort(unique(feature_id), decreasing = T)[5]) %>%
  ungroup() %>%
  mutate(Condition = factor(case_match(Condition, 
                                       "1"~"Osteoarthritis",
                                       "0"~"Control"),
                            levels = c("Osteoarthritis","Control")),
         Tissue = factor(case_match(Subsection, 
                                    "1"~"Medial",
                                    "0"~"Lateral"),
                         levels = c("Medial", "Lateral")))
custom_colors <- c("#40B0A6", "#E1BE6A")  


ggarrange(
  annotate_figure(
    ggarrange(
      annotate_figure(
        itx_df %>%
          filter(vari == top_labs[1]) %>%
          ggplot(aes(x = Tissue, y = mean, group = subject, color = Condition, fill = Condition)) +
          geom_point() +
          geom_line() +
          ylab("Mean Intensity") +
          xlab(NULL) +
          theme_classic() + 
          theme(legend.position = 'None'),
        top= top_labs[1]),
      annotate_figure(
        itx_df %>%
          filter(vari == top_labs[2]) %>%
          ggplot(aes(x = Tissue, y = mean, group = subject, color = Condition, fill = Condition)) +
          geom_point() +
          geom_line() +
          ylab(NULL) +
          xlab(NULL) +
          theme_classic() + 
          theme(legend.position = 'none'),
        top= top_labs[2]),
      ncol = 2, widths = c(4,4), common.legend = T, legend = "bottom", legend.grob = get_legend(itx_df %>%
                                                                                                  filter(vari == top_labs[1]) %>%
                                                                                                  ggplot(aes(x = Subsection, y = mean, group = subject, color = Condition, fill = Condition)) +
                                                                                                  geom_point() +
                                                                                                  geom_line() +
                                                                                                  ylab(NULL) +
                                                                                                  xlab(NULL) +
                                                                                                  theme_classic() +
                                                                                                  theme(legend.position = "bottom")))
  ),
  
  ggplot() + theme_minimal(),
  
  
  annotate_figure(
    ggarrange(
      annotate_figure(
        ggplot(df_test_vari %>% filter(vari == top_labs[1]) %>% 
                 filter(finding %in% c("False Positive","True Positive"))  %>% 
                 mutate(n = ifelse(finding == "False Positive",
                                   -1*n,
                                   n),
                        finding = factor(finding, levels= c("True Positive","False Positive"))), 
               aes(x = model_simple, y = n/n_each, fill = finding)) + 
          theme_classic() +
          geom_bar(stat = "identity", position = "stack") + 
          facet_grid(.~ null, space = "free", scales = "free_x") +
          theme(legend.position = 'none') +
          ylab(" ") +
          xlab(NULL) +
          geom_hline(yintercept = 0, color = "white") +
          geom_hline(yintercept = -0.05, color = "black", linetype = "dashed") +
          scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(-0.2,1)) + 
          scale_fill_manual(values = custom_colors),
        top= top_labs[1]),
      annotate_figure(
        ggplot(df_test_vari %>% filter(vari == top_labs[2]) %>% 
                 filter(finding %in% c("False Positive","True Positive"))  %>% 
                 mutate(n = ifelse(finding == "False Positive",
                                   -1*n,
                                   n),
                        finding = factor(finding, levels= c("True Positive","False Positive"))), 
               aes(x = model_simple, y = n/n_each, fill = finding)) + 
          theme_classic() +
          geom_bar(stat = "identity", position = "stack") + 
          facet_grid(.~ null, space = "free", scales = "free_x") +
          theme(legend.position = 'none') +
          ylab(NULL) +
          xlab(NULL) +
          geom_hline(yintercept = 0, color = "white") +
          geom_hline(yintercept = -0.05, color = "black", linetype = "dashed") +
          scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(-0.2,1)) + 
          scale_fill_manual(values = custom_colors),
        top= top_labs[2]),
      ncol = 2, widths = c(4,4), common.legend = T, legend = "bottom", legend.grob = get_legend(ggplot(df_test_vari %>% filter(vari == top_labs[2]) %>% 
                                                                                                         filter(finding %in% c("False Positive","True Positive"))  %>% 
                                                                                                         mutate(n = ifelse(finding == "False Positive",
                                                                                                                           -1*n,
                                                                                                                           n),
                                                                                                                finding = factor(finding, levels= c("True Positive","False Positive"))), 
                                                                                                       aes(x = model_simple, y = n/n_each, fill = finding)) + 
                                                                                                  theme_classic() +
                                                                                                  geom_bar(stat = "identity", position = "stack") + 
                                                                                                  facet_grid(.~ null, space = "free", scales = "free_x") +
                                                                                                  theme(legend.position = 'bottom') +
                                                                                                  ylab(NULL) +
                                                                                                  xlab(NULL) +
                                                                                                  labs(fill = "Rate") + 
                                                                                                  geom_hline(yintercept = 0, color = "white") +
                                                                                                  geom_hline(yintercept = -0.05, color = "black", linetype = "dashed") +
                                                                                                  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(-1,1)) + 
                                                                                                  scale_fill_manual(values = custom_colors)))
  ),
  nrow = 3, heights = c(1,0.05,1.2))

ggsave(file.path("Recreated_Figures","Figure-Variation-Simulation-2.png"), width = 11.5, height = 10, units = "in")

###### NEW Pix as rep figure #####

library(data.table)

results_list <- list()
for(i in 1:16) {
  m <- sim_list_adj[[i]]
  transposed_intensity_matrix <- t(as.matrix(intensity(m)))
  
  # Convert to data.table for better memory usage
  dt <- as.data.table(pData(sim_list_adj[[1]])[,1:2])
  
  # Add intensity columns
  for(col in 1:ncol(transposed_intensity_matrix)) {
    dt[[as.character(col)]] <- transposed_intensity_matrix[, col]
  }
  
  # Melt more efficiently than pivot_longer
  dt_long <- melt(dt, 
                  id.vars = c("x", "y"),
                  measure.vars = as.character(1:3000),
                  variable.name = "feature_id",
                  value.name = "intensity")
  
  # Add metadata
  dt_long[, subject := samples[i,1]]
  dt_long[, condition := samples[i,3]]
  dt_long[, subtissue := samples[i,2]]
  
  results_list[[i]] <- dt_long
  # Clean up memory
  rm(dt, dt_long, transposed_intensity_matrix)
  gc()
}

pix_as_rep_sim <- rbindlist(results_list)
rm(results_list)
# Convert factors if needed
cols_to_factor <- c("feature_id", "subject", "condition", "subtissue")
pix_as_rep_sim[, (cols_to_factor) := lapply(.SD, as.factor), .SDcols = cols_to_factor]

pix_as_rep_sim_lowb_hight <- pix_as_rep_sim[feature_id %in% as.factor(feature_ids_for_lowb_hight), ]


res_pixasrep_sim_lowb_hight <- run_models_MEMORY_EFF(pix_as_rep_sim_lowb_hight, 
                                          "intensity ~ condition * subtissue + (1|subject)")


pix_as_rep_sim_highb_lowt <- pix_as_rep_sim[feature_id %in% as.factor(feature_ids_for_highb_lowt), ]


#res_pixasrep_sim_highb_lowt  <- run_models(pix_as_rep_sim_highb_lowt, 
#                                           "intensity ~ condition * subtissue + (1|subject)")

res_pixasrep_sim_highb_lowt  <- run_models_MEMORY_EFF(pix_as_rep_sim_highb_lowt, 
                                           "intensity ~ condition * subtissue + (1|subject)")

get_comparison_plotting_frame_pix_as_rep <- function(mse_frame, model_res, nmax = 150) {
  ttest_frame_condition <- mse_frame %>% 
    filter(subtissue == "1") %>%
    group_by(feature_id) %>% 
    group_modify(~tidy(t.test(intensity ~ condition, data = .x, var.equal = T))) %>%
    ungroup() %>%
    mutate(SE = estimate / statistic,
           adj_p.value = p.adjust(p.value, method = "fdr"),
           old_feature_id = feature_id,
           feature_id = row_number())
  
  
  ttest_frame_subtissue <- mse_frame %>% 
    filter(condition == "1") %>%
    group_by(feature_id) %>% 
    arrange(subtissue, subject) %>%
    group_modify(~tidy(t.test(intensity ~ subtissue, data = .x, var.equal = T))) %>%
    ungroup() %>%
    mutate(SE = estimate / statistic,
           adj_p.value = p.adjust(p.value, method = "fdr"),
           old_feature_id = feature_id,
           feature_id = row_number())
  
  contrasts_means <- bind_rows(lapply(model_res[[1]], FUN = function(m){
    bind_rows(tidy(contrast(emmeans(m, ~ condition | subtissue, lmer.df = "satterthwaite"), "pairwise", adjust = "none")) %>% dplyr::rename("subgroup"="subtissue"),
              tidy(contrast(emmeans(m, ~ subtissue | condition, lmer.df = "satterthwaite"), "pairwise", adjust = "none")) %>% dplyr::rename("subgroup"="condition")) 
  })) %>% 
    group_by(term, subgroup) %>%
    mutate(feature_id = row_number(),
           adj_p.value = p.adjust(p.value, method = "fdr")) %>% ungroup()
  
  tests <- bind_rows(
    
    contrasts_means%>%
      filter(subgroup=="1") %>%
      mutate(significant=`p.value`< 0.05,
             truth = as.numeric(feature_id) > 150,
             finding = ifelse(significant,
                              ifelse(truth,
                                     "True Positive",
                                     "False Positive"),
                              ifelse(truth,
                                     "False Negative",
                                     "True Negative"))) %>% 
      group_by(finding, term, .drop = FALSE) %>% 
      summarise(n=n(), .groups = "keep") %>%
      ungroup() %>% 
      mutate(model = "Mixed Effects Model with All Variables",
             test = "T Test of Conditional Means",
             null = case_match(term,
                               "condition"~"In Medial:\nControl vs Osteoarthritis",
                               "subtissue"~"In Osteoarthritis:\nLateral vs Medial"),
             term = case_match(term,
                               "condition"~"Condition",
                               "condition:subtissue"~"Interaction",
                               "subtissue"~"Subsection")) %>%
      dplyr::select(finding, n, model, test, null, term),
    
    ttest_frame_condition %>% 
      mutate(significant=`p.value`< 0.05,
             truth = as.numeric(feature_id) > nmax,
             finding = ifelse(significant,
                              ifelse(truth,
                                     "True Positive",
                                     "False Positive"),
                              ifelse(truth,
                                     "False Negative",
                                     "True Negative"))) %>% 
      group_by(finding, .drop = FALSE) %>% 
      summarise(n=n(), .groups = "keep") %>%
      ungroup() %>% 
      mutate(model = "Linear Model with Condition",
             null = "In Medial:\nControl vs Osteoarthritis",
             test = "Two Sample T Test",
             term = "Condition") %>%
      dplyr::select(finding, n, model, test, null, term),
    
    ttest_frame_subtissue %>% 
      mutate(significant=`p.value`< 0.05,
             truth = as.numeric(feature_id) > nmax,
             finding = ifelse(significant,
                              ifelse(truth,
                                     "True Positive",
                                     "False Positive"),
                              ifelse(truth,
                                     "False Negative",
                                     "True Negative"))) %>% 
      group_by(finding,  .drop = FALSE) %>% 
      summarise(n=n(), .groups = "keep") %>%
      ungroup() %>%
      mutate(model = "Linear Model of Paired Differences\nBetween Subtissues",
             null = "In Osteoarthritis:\nLateral vs Medial",
             test = "Paired T Test",
             term = "Subsection") %>%
      dplyr::select(finding, n, model, test, null, term)
  ) %>% mutate(termc = case_match(term,
                                  "Condition"~"Control vs OA",
                                  "Interaction"~"Interaction",
                                  "Subsection"~"Lateral vs Medial")) %>%
    mutate(termc = ifelse(test == "F Test of Mean Squares",
                          term, 
                          termc)) %>%
    mutate(termc = factor(termc, 
                          levels = c("Condition",
                                     "Subsection",
                                     "Interaction",
                                     "Control vs OA",
                                     "Lateral vs Medial"),
                          labels = c("Condition\nAny Difference",
                                     "Subsection\nAny Difference",
                                     "Interaction\nAny Present",
                                     "Condition\nControl vs OA",
                                     "Subsection\nLateral vs Medial")))
  
  return(tests)
  
}

rm(sim_list_adj, pix_as_rep_sim)

tests_pix_as_reps_lowb_hight <- get_comparison_plotting_frame_pix_as_rep(pix_as_rep_sim_lowb_hight, res_pixasrep_sim_lowb_hight)
tests_pix_as_reps_highb_lowt <- get_comparison_plotting_frame_pix_as_rep(pix_as_rep_sim_highb_lowt, res_pixasrep_sim_highb_lowt)


##### Pix as reps Fig for paper #####

top_labs <- c("Greater Technical Variation",
              "Greater Biological Variation")




df_test_vari_pix_as_reps <- bind_rows(add_zero_tests(tests_pix_as_reps_lowb_hight) %>% mutate(vari = top_labs[1]),
                                      add_zero_tests(tests_pix_as_reps_highb_lowt) %>% mutate(vari = top_labs[2])) %>%
  mutate(vari = factor(vari, levels = top_labs)) %>%
  mutate(model_simple = factor(ifelse(test=="T Test of Conditional Means",
                                      "Mixed Effects Model",
                                      test),
                               levels = c("Mixed Effects Model",
                                          "Two Sample T Test",
                                          "Paired T Test"),
                               labels =c("Supplementary\nTable 1 Model 3:\nMixed Effects\nModel with\nboth variables",
                                         "Supplementary\nTable 1 Model 1:\nEquivalent to\nTwo Sample T Test\nin this design",
                                         "Supplementary\nTable 1 Model 2:\nEquivalent to\nTwo Sample T Test\nin this design")),
         model_family = factor(ifelse(test=="T Test of Conditional Means",
                                      "Mixed Effects Models",
                                      "Sample T Tests"),
                               levels = c("Mixed Effects Models",
                                          "Sample T Tests")),
         comparison_type = ifelse(term=="Subsection",
                                  "Within Subjects",
                                  "Between Subjects"),
         null = factor(null,
                       levels = c("In Medial:\nControl vs Osteoarthritis",
                                  "In Osteoarthritis:\nLateral vs Medial"),
                       labels =c("Hypothesis A, Between-subjects\nIn Medial: Osteoarthritis vs Control",
                                 "Hypothesis B, Within-subjects\nIn Osteoarthritis: Medial vs Lateral")))

annotate_figure(
  ggarrange(
    annotate_figure(
      ggplot(df_test_vari_pix_as_reps %>% filter(vari == top_labs[1]) %>% 
               filter(finding %in% c("False Positive","True Positive"))  %>% 
               mutate(n = ifelse(finding == "False Positive",
                                 -1*n,
                                 n),
                      finding = factor(finding, levels= c("True Positive","False Positive"))), 
             aes(x = model_simple, y = n/150, fill = finding)) + 
        theme_classic() +
        geom_bar(stat = "identity", position = "stack") + 
        facet_grid(.~ null, space = "free", scales = "free_x") +
        theme(legend.position = 'none') +
        ylab(" ") +
        xlab(NULL) +
        geom_hline(yintercept = 0, color = "white") +
        geom_hline(yintercept = -0.05, color = "black", linetype = "dashed") +
        scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(-1,1)) + 
        scale_fill_manual(values = custom_colors),
      top= top_labs[1]),
    annotate_figure(
      ggplot(df_test_vari_pix_as_reps %>% filter(vari == top_labs[2]) %>% 
               filter(finding %in% c("False Positive","True Positive"))  %>% 
               mutate(n = ifelse(finding == "False Positive",
                                 -1*n,
                                 n),
                      finding = factor(finding, levels= c("True Positive","False Positive"))), 
             aes(x = model_simple, y = n/150, fill = finding)) + 
        theme_classic() +
        geom_bar(stat = "identity", position = "stack") + 
        facet_grid(.~ null, space = "free", scales = "free_x") +
        theme(legend.position = 'none') +
        ylab(NULL) +
        xlab(NULL) +
        geom_hline(yintercept = 0, color = "white") +
        geom_hline(yintercept = -0.05, color = "black", linetype = "dashed") +
        scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(-1,1)) + 
        scale_fill_manual(values = custom_colors),
      top= top_labs[2]),
    ncol = 2, widths = c(4,4), common.legend = T, legend = "bottom", legend.grob = get_legend(ggplot(df_test_vari_pix_as_reps %>% filter(vari == top_labs[2]) %>% 
                                                                                                       filter(finding %in% c("False Positive","True Positive"))  %>% 
                                                                                                       mutate(n = ifelse(finding == "False Positive",
                                                                                                                         -1*n,
                                                                                                                         n),
                                                                                                              finding = factor(finding, levels= c("True Positive","False Positive"))), 
                                                                                                     aes(x = model_simple, y = n/150, fill = finding)) + 
                                                                                                theme_classic() +
                                                                                                geom_bar(stat = "identity", position = "stack") + 
                                                                                                facet_grid(.~ null, space = "free", scales = "free_x") +
                                                                                                theme(legend.position = 'bottom') +
                                                                                                ylab(NULL) +
                                                                                                xlab(NULL) +
                                                                                                labs(fill = "Rate") + 
                                                                                                geom_hline(yintercept = 0, color = "white") +
                                                                                                scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(-1,1)) + 
                                                                                                scale_fill_manual(values = custom_colors)))
)

ggsave(file.path("Recreated_Figures","Figure-Variation-Simulation-2-Pixels-as-Replicates.png"), width = 11.5, height = 6, units = "in")



##### Overestimation of effect sizes of significant values #####

with_effects <- run_models_sim_list[[1]][1501:3000]

contrasts_means_with_effects <- bind_rows(lapply(with_effects, FUN = function(m){
  bind_rows(tidy(contrast(emmeans(m, ~ condition | subtissue, lmer.df = "satterthwaite"), "pairwise", adjust = "none")) %>% dplyr::rename("subgroup"="subtissue"),
            tidy(contrast(emmeans(m, ~ subtissue | condition, lmer.df = "satterthwaite"), "pairwise", adjust = "none")) %>% dplyr::rename("subgroup"="condition")) 
})) %>% 
  group_by(term, subgroup) %>%
  mutate(feature_id = row_number(),
         adj_p.value = p.adjust(p.value, method = "fdr")) %>% ungroup()

means <- contrasts_means_with_effects %>%
  filter(subgroup == 1) %>% 
  mutate(group = p.value > 0.05) %>%
  group_by(group) %>%
  summarize(mean_estimate = -1*mean(abs(estimate), na.rm = TRUE))

# Plot with histogram and dashed vertical lines at group means
contrasts_means_with_effects %>%
  filter(subgroup == 1) %>% 
  mutate(group = p.value > 0.05) %>%
  ggplot(aes(x = -1*estimate, group = group, fill = group)) +
  theme_classic() +
  geom_histogram(bins = 100, position = "identity", alpha = 0.4) +
  geom_segment(x = 15, xend = 15, y = 0, yend = 64, size = .8) +
  geom_segment(data = means,
               aes(x = mean_estimate, xend = mean_estimate, y = 0, yend = 71, color = group),
               linetype = "31", size = 1) +
  labs(x = "Estimates of Effect",
       y = "Number of comparisons of from similations",
       group = "P < 0.05") +
  annotate("text", x = 15, y = 68, label = "True\nEffect") +
  annotate("text", x = means$mean_estimate, y = 74, label = c("Mean Effect\n when P>0.05","Mean Effect\nwhen P<0.05"))
