#Author: Lars Moret, MSc student @ FunGen
#Date: 29/10/2026
#Description:
#    R script for visualizing MSH6 G233A and within-clonal-group genetic diversity
#Inputs: clonal_group_pairwise_distances_untrimmed.tsv and clonal_vcf_MSH6_G233A_summary.tsv

## install packages
install.packages(c("readr", "dplyr", "ggplot2", "broom"))
library(readr)
library(dplyr)
library(ggplot2)
library(broom)

setwd("C:/Users/larsm/Downloads")
diversity_file <- "clonal_group_pairwise_distances_untrimmed.tsv"
msh6_file      <- "clonal_vcf_MSH6_G233A_summary.tsv"
if (!file.exists(diversity_file) || !file.exists(msh6_file)) {
  stop("One or both input files are missing. ")}

dist <- read_tsv(diversity_file, show_col_types = FALSE)
msh6 <- read_tsv(msh6_file, show_col_types = FALSE)

required_dist <- c("clonal_group", "group_size", "pair_count", "snp_sites", "observed_differences_per_pair", "pairwise_snp_difference_rate", "callability")
required_msh6 <- c("clonal_group", "group_status")

if (!all(required_dist %in% names(dist))) {
  stop("The diversity table is missing one or more required columns.")}

if (!all(required_msh6 %in% names(msh6))) {
  stop("The MSH6 summary is missing one or more required columns.")}

if (anyDuplicated(dist$clonal_group)) {
  stop("Duplicate clonal_group IDs found in the diversity table.")}

if (anyDuplicated(msh6$clonal_group)) {
  stop("Duplicate clonal_group IDs found in the MSH6 table.")}

## Join the diversity measrment
comparison_all <- dist %>%
  select(-any_of("MSH6_G233A_status")) %>%
  left_join(msh6 %>% select(clonal_group, group_status), by = "clonal_group") %>%
  mutate(MSH6_class = case_when(group_status == "Present" ~ "A233", group_status == "Absent"  ~ "G233", TRUE ~ NA_character_),
    MSH6_class = factor(MSH6_class, levels = c("A233", "G233")))

comparison <- comparison_all %>% filter(!is.na(MSH6_class)) %>% mutate(MSH6_class = droplevels(MSH6_class))

if (nlevels(comparison$MSH6_class) != 2) {
  stop("Both A233 and G233 groups must be present for the comparison.")}

write_tsv(comparison, "MSH6_clonal_diversity_joined.tsv")

## some quality checks
cat("\nMissing values in primary variables:\n")
print(
  comparison %>%
    summarise(missing_rate = sum(is.na(pairwise_snp_difference_rate)), missing_observed_differences = sum(is.na(observed_differences_per_pair)),
      missing_callability = sum(is.na(callability)), missing_group_size = sum(is.na(group_size))))

# Check for invalid rates
if (any(comparison$pairwise_snp_difference_rate < 0, na.rm = TRUE)) {
  stop("Negative pairwise SNP difference rates were found.")}

if (any(comparison$callability < 0 |
        comparison$callability > 1, na.rm = TRUE)) {
  stop("Callability values outside [0, 1] were found.")}

group_summary <- comparison %>%
  group_by(MSH6_class) %>%
  summarise(n_groups = n(),
    median_group_size = median(group_size, na.rm = TRUE),
    mean_group_size = mean(group_size, na.rm = TRUE),
    min_group_size = min(group_size, na.rm = TRUE),
    max_group_size = max(group_size, na.rm = TRUE),
    median_pairwise_rate = median(pairwise_snp_difference_rate, na.rm = TRUE),
    mean_pairwise_rate = mean(pairwise_snp_difference_rate, na.rm = TRUE),
    sd_pairwise_rate = sd(pairwise_snp_difference_rate, na.rm = TRUE),
    IQR_pairwise_rate = IQR(pairwise_snp_difference_rate, na.rm = TRUE),
    median_observed_differences = median(observed_differences_per_pair, na.rm = TRUE),
    mean_observed_differences = mean(observed_differences_per_pair, na.rm = TRUE),
    median_callability = median(callability, na.rm = TRUE),
    mean_callability = mean(callability, na.rm = TRUE),
    zero_rate_groups = sum(pairwise_snp_difference_rate == 0, na.rm = TRUE),
    .groups = "drop")

write_tsv(group_summary, "MSH6_clonal_diversity_group_summary.tsv") ## check this file for the summary

## unadjusted comparison using a two-sided wilcoxon rank-sum test
# Each clonal group contributes one diversity-rate observation
wilcox_primary <- wilcox.test(pairwise_snp_difference_rate ~ MSH6_class, data = comparison, alternative = "two.sided", exact = FALSE, conf.int = FALSE)
print(wilcox_primary)

#and the complementary comparison of observed differences per pair
wilcox_observed <- wilcox.test(observed_differences_per_pair ~ MSH6_class, data = comparison, alternative = "two.sided", exact = FALSE)
print(wilcox_observed)

# ------------------------------------------------------------
# 7. Adjusted model
# ------------------------------------------------------------


## Now the adjusted model, which excludes 0 rate groups since log(0) is undefined
# this is a log-linear model for positive pairwise rates, where A233 is the reference class
# This model estimates an association, not a causal effect

comparison_positive <- comparison %>% filter(pairwise_snp_difference_rate > 0, !is.na(group_size), group_size > 0, !is.na(callability)) %>% 
  mutate(MSH6_class = factor(MSH6_class, levels = c("A233", "G233")))

model_adjusted <- lm(log(pairwise_snp_difference_rate) ~ MSH6_class + log(group_size) + callability, data = comparison_positive)

print(summary(model_adjusted))

model_results <- tidy(model_adjusted, conf.int = TRUE)
write_tsv(model_results, "MSH6_clonal_diversity_adjusted_model.tsv")
print(model_results)

# Convert MSH6 effect from log scale to a ratio of rates, where G233 is relative to A233
msh6_effect <- model_results %>%
  filter(term == "MSH6_classG233") %>%
  mutate(rate_ratio_G233_vs_A233 = exp(estimate), ratio_ci_lower = exp(conf.low), ratio_ci_upper = exp(conf.high),
    rate_ratio_A233_vs_G233 = exp(-estimate), A233_vs_G233_ci_lower = exp(-conf.high), A233_vs_G233_ci_upper = exp(-conf.low))

print(msh6_effect)


## Sensivity analysis due to the large class imbalance between clonal groups 
## (i.e. a lot of smaller groups)

comparison_small <- comparison %>% filter(group_size >= 2, group_size <= 15) %>%
  mutate(MSH6_class = droplevels(MSH6_class))

if (nlevels(comparison_small$MSH6_class) == 2) {
  wilcox_small <- wilcox.test(pairwise_snp_difference_rate ~ MSH6_class,
    data = comparison_small, alternative = "two.sided", exact = FALSE)
  
  print(wilcox_small)
  
  comparison_small_positive <- comparison_small %>%
    filter(pairwise_snp_difference_rate > 0, !is.na(callability)) %>%
    mutate(MSH6_class = factor(MSH6_class, levels = c("A233", "G233")))
  
  if (nlevels(droplevels(comparison_small_positive$MSH6_class)) == 2 &&
      nrow(comparison_small_positive) > 3) {
    
    model_small <- lm(log(pairwise_snp_difference_rate) ~ MSH6_class + log(group_size) + callability, data = comparison_small_positive)
    print(summary(model_small))
    
    write_tsv(tidy(model_small, conf.int = TRUE), "MSH6_clonal_diversity_model_size2_15.tsv")}
  
} else {cat("Sensitivity test skipped")}

write_tsv(comparison_small, "MSH6_clonal_diversity_size2_15.tsv")


comparison_callable <- comparison %>% filter(callability >= 0.98) %>% mutate(MSH6_class = droplevels(MSH6_class))
print(table(comparison_callable$MSH6_class))

if (nlevels(comparison_callable$MSH6_class) == 2) {
  wilcox_callable <- wilcox.test(pairwise_snp_difference_rate ~ MSH6_class, 
    data = comparison_callable, alternative = "two.sided",exact = FALSE)
  
  print(wilcox_callable)
} else {cat("High-callability test skipped")}

write_tsv(comparison_callable, "MSH6_clonal_diversity_callability98.tsv")

## plotting for the diversity rate distributions, log10 axis for the positive rates
## 0-rate groups are, again, excluded

plot_data <- comparison %>% filter(pairwise_snp_difference_rate > 0, !is.na(pairwise_snp_difference_rate))

p <- ggplot(plot_data, aes(x = MSH6_class, y = pairwise_snp_difference_rate)) +
  geom_boxplot(outlier.shape = NA, width = 0.55) +
  geom_jitter(width = 0.16, height = 0, alpha = 0.45, size = 1.5) +
  scale_y_log10() + 
  labs(title = "Within-clonal-group SNP diversity by MSH6 class",
    subtitle = "Globally ascertained SNP panel; positive rates only",
    x = "MSH6 class", y = "Pairwise SNP difference rate (log10 scale)") + theme_classic(base_size = 12)

print(p)

ggsave("MSH6_clonal_diversity_boxplot.png", plot = p, width = 7, height = 5, dpi = 300)
