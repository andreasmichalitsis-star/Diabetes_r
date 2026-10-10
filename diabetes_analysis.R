################################################################################
# FINAL PROJECT: UNDERSTANDING DIABETES PROGRESSION
# Reproducible analysis script | Scenario 1
# Place this script beside diabetes.tab.txt (or select the data file when asked).
################################################################################

# Install these once if needed, then load them.
needed <- c("ggplot2", "dplyr", "broom")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing)
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(broom)})
options(digits = 4)

# Find the supplied data file in the project folder; otherwise open a chooser.
possible <- c("diabetes.tab.txt", "diabetes.tab(1).txt")
data_file <- possible[file.exists(possible)][1]
if (is.na(data_file)) {
  message("Choose diabetes.tab.txt")
  data_file <- file.choose()
}
diabetes <- read.delim(data_file, check.names = FALSE)
names(diabetes) <- tolower(names(diabetes))
if ("y" %in% names(diabetes)) names(diabetes)[names(diabetes) == "y"] <- "progression"
diabetes$sex <- factor(diabetes$sex)

# Create an outputs folder next to this script.
out_dir <- "project_outputs"
dir.create(out_dir, showWarnings = FALSE)

# Q1: data checks
cat("Dimensions:\n"); print(dim(diabetes))
cat("Variables and structure:\n"); print(str(diabetes))
cat("Missing values by variable:\n"); print(colSums(is.na(diabetes)))
cat("Duplicate rows:\n"); print(sum(duplicated(diabetes)))
cat("Summary:\n"); print(summary(diabetes))

# Q2: descriptive statistics
num <- diabetes |> select(where(is.numeric))
descriptive <- num |> summarise(across(everything(), list(
  mean = ~mean(.x), median = ~median(.x), sd = ~sd(.x),
  min = ~min(.x), q1 = ~quantile(.x, .25), q3 = ~quantile(.x, .75), max = ~max(.x)
)))
write.csv(descriptive, file.path(out_dir, "descriptive_statistics.csv"), row.names = FALSE)

p1 <- ggplot(diabetes, aes(bmi, progression)) + geom_point(alpha=.6, colour="#1D8A8A") +
  geom_smooth(method="lm", se=TRUE, colour="#15324B") +
  labs(title="BMI and one-year diabetes progression", x="BMI", y="Progression measure") + theme_minimal()
ggsave(file.path(out_dir,"bmi_progression.png"), p1, width=7, height=4, dpi=200)
p2 <- ggplot(diabetes, aes(progression)) + geom_histogram(bins=25, fill="#1D8A8A", colour="white") +
  labs(title="Distribution of diabetes progression", x="Progression measure", y="Patients") + theme_minimal()
ggsave(file.path(out_dir,"progression_distribution.png"), p2, width=7, height=4, dpi=200)

# Q3: correlations (sex excluded because it is categorical)
cor_matrix <- cor(num, use="complete.obs", method="pearson")
write.csv(round(cor_matrix, 3), file.path(out_dir,"correlation_matrix.csv"))
print(sort(cor_matrix[,"progression"], decreasing=TRUE))

# Q4: simple and multiple linear regression
simple_model <- lm(progression ~ bmi, data=diabetes)
multiple_model <- lm(progression ~ bmi + bp + s5, data=diabetes)
print(summary(simple_model)); print(tidy(simple_model, conf.int=TRUE)); print(glance(simple_model))
print(summary(multiple_model)); print(tidy(multiple_model, conf.int=TRUE)); print(glance(multiple_model))
model_comparison <- bind_rows(
  glance(simple_model) |> mutate(model="Simple: BMI"),
  glance(multiple_model) |> mutate(model="Multiple: BMI + BP + S5")
) |> select(model, r.squared, adj.r.squared, sigma)
print(model_comparison)
write.csv(tidy(multiple_model, conf.int=TRUE), file.path(out_dir,"multiple_model_coefficients.csv"), row.names=FALSE)

# Q5: residual diagnostics and VIF
aug <- broom::augment(multiple_model)
p3 <- ggplot(aug, aes(.fitted, .resid)) + geom_point(alpha=.6, colour="#1D8A8A") +
  geom_hline(yintercept=0, linetype=2) + geom_smooth(se=FALSE, colour="#D28A24") +
  labs(title="Residuals versus fitted values", x="Fitted progression", y="Residual") + theme_minimal()
ggsave(file.path(out_dir,"residuals_vs_fitted.png"), p3, width=7, height=4, dpi=200)
manual_vif <- function(model) {
  X <- model.matrix(model)[,-1,drop=FALSE]
  sapply(seq_len(ncol(X)), function(i) {
    others <- X[,-i,drop=FALSE]
    if (!ncol(others)) return(1)
    r2 <- summary(lm(X[,i] ~ others))$r.squared
    1/(1-r2)
  }) |> setNames(colnames(X))
}
print(manual_vif(multiple_model))
print(confint(multiple_model))

# Q6: observational design; interpret associations cautiously, not causally.
# All tables and figures above are generated from the supplied source data.
################################################################################
