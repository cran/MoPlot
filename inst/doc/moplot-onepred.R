## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message = FALSE,
  warning = TRUE,
  fig.width = 8,
  fig.height = 6,
  fig.retina = 1,
  fig.align = "center",
  out.width = "100%"
)
library(MoPlot)
set.seed(622)

## ----installation, eval = FALSE-----------------------------------------------
# install.packages("MoPlot")

## ----treatment-model----------------------------------------------------------
dat_iris <- datasets::iris
contrasts(dat_iris$Species) <- contr.treatment(3)
model_treatment <- lm(Sepal.Length ~ Species, data = dat_iris)

## ----treatment-table----------------------------------------------------------
model_contrast_table(model_treatment, type = "x")

## ----treatment-basic, fig.height = 4.3, fig.cap = "Observed sepal-length means under treatment coding."----
moplot_onepred(model_treatment, "Species", show_coef = FALSE, caption = FALSE)

## ----treatment-numbers--------------------------------------------------------
treatment_comparisons <- contrasts_onepred(model_treatment, "Species")
treatment_comparisons[, c("Contrast", "Levels", "Estimate", "CI_low", "CI_high",
                          "Significance", "Scale_factor")]

## ----treatment-complete, fig.cap = "Treatment contrasts with raw observations and the aligned coefficient panel."----
moplot_onepred(model_treatment, "Species", show_data = TRUE,
               style = "publication", caption = FALSE)

## ----treatment-base-----------------------------------------------------------
dat_base <- dat_iris
contrasts(dat_base$Species) <- contr.treatment(3, base = 3)
base_comparisons <- contrasts_onepred(lm(Sepal.Length ~ Species, data = dat_base), "Species")
unique(base_comparisons[, c("Contrast", "Levels")])

## ----sum-model, fig.cap = "Sum contrasts: each species against the Grand Mean."----
dat_sum <- dat_iris
contrasts(dat_sum$Species) <- contr.sum(3)
model_sum <- lm(Sepal.Length ~ Species, data = dat_sum)
moplot_onepred(model_sum, "Species", show_data = TRUE,
               style = "publication", caption = FALSE)

## ----helmert-model------------------------------------------------------------
set.seed(120)
dat_groups <- data.frame(Group = factor(rep(c("G1", "G2", "G3", "G4"), each = 30)))
dat_groups$Score <- rnorm(120, mean = rep(c(116, 136, 113, 119), each = 30), sd = 10)

dat_helmert <- dat_groups
contrasts(dat_helmert$Group) <- contr.helmert(4)
model_helmert <- lm(Score ~ Group, data = dat_helmert)

## ----helmert-plot, fig.cap = "Helmert contrasts. Each coefficient is a fraction of the displayed difference."----
moplot_onepred(model_helmert, "Group", show_data = TRUE,
               style = "publication", caption = FALSE)

## ----helmert-per-difference, fig.cap = "The same Helmert contrasts with per_difference = TRUE: the lower panel now equals the displayed differences."----
moplot_onepred(model_helmert, "Group", per_difference = TRUE,
               style = "publication", caption = FALSE)

## ----reverse-helmert, fig.cap = "Reverse Helmert contrasts."------------------
dat_reverse <- dat_groups
contrasts(dat_reverse$Group) <- reverse_helmert(4)
moplot_onepred(lm(Score ~ Group, data = dat_reverse), "Group",
               style = "publication", caption = FALSE)

## ----sliding, fig.cap = "Sliding differences between adjacent groups."--------
dat_sliding <- dat_groups
contrasts(dat_sliding$Group) <- MASS::contr.sdif(4)
moplot_onepred(lm(Score ~ Group, data = dat_sliding), "Group",
               style = "publication", caption = FALSE)

## ----custom, fig.cap = "Comparisons defined by an orthogonal custom coding."----
dat_custom <- dat_groups
contrasts(dat_custom$Group) <- cbind(
  first_vs_last = c(1, 1, -1, -1),
  within_first  = c(1, -1, 0, 0),
  within_last   = c(0, 0, 1, -1)
)
moplot_onepred(lm(Score ~ Group, data = dat_custom), "Group",
               style = "publication", caption = FALSE)

## ----polynomial, fig.height = 6.5, fig.cap = "Linear and quadratic components of an ordered factor."----
set.seed(7)
dat_dose <- data.frame(
  Dose = ordered(rep(c("Low", "Medium", "High"), each = 30),
                 levels = c("Low", "Medium", "High"))
)
dat_dose$Response <- rnorm(90, mean = c(10, 14, 15)[as.integer(dat_dose$Dose)], sd = 3)
moplot_onepred(lm(Response ~ Dose, data = dat_dose), "Dose", caption = FALSE)

## ----predicted, fig.cap = "Predicted means at the mean of the covariate."-----
set.seed(21)
dat_cov <- dat_groups
dat_cov$Age <- rnorm(120, mean = 40 + 5 * as.integer(dat_cov$Group), sd = 8)
dat_cov$Score <- dat_cov$Score + 0.5 * dat_cov$Age
contrasts(dat_cov$Group) <- MASS::contr.sdif(4)
model_cov <- lm(Score ~ Group + Age, data = dat_cov)

moplot_onepred(model_cov, "Group", style = "publication", caption = FALSE)

## ----predicted-check----------------------------------------------------------
predicted <- contrasts_onepred(model_cov, "Group")
attr(predicted, "estimate_type")
differences <- predicted$Estimate[c(1, 3, 5)] - predicted$Estimate[c(2, 4, 6)]
round(rbind(displayed = differences, coefficients = coef(model_cov)[2:4]), 3)

## ----binomial, fig.cap = "Predicted probabilities above and odds ratios below."----
set.seed(330)
dat_binomial <- data.frame(Group = factor(rep(c("G1", "G2", "G3"), each = 80)))
dat_binomial$Success <- rbinom(240, size = 1, prob = rep(c(0.25, 0.50, 0.70), each = 80))
model_binomial <- glm(Success ~ Group, family = binomial, data = dat_binomial)

moplot_onepred(model_binomial, "Group", estimate_type = "predicted",
               coefficient_scale = "interpreted", style = "publication",
               caption = FALSE)

## ----mixed, fig.cap = "Population-level predictions from a linear mixed-effects model."----
set.seed(440)
dat_mixed <- expand.grid(Subject = factor(seq_len(24)),
                         Group = factor(c("G1", "G2", "G3")),
                         Trial = seq_len(5))
subject_intercept <- rnorm(24, sd = 1.5)
dat_mixed$Score <- 10 + c(0, 1.5, 3)[as.integer(dat_mixed$Group)] +
  subject_intercept[as.integer(dat_mixed$Subject)] + rnorm(nrow(dat_mixed))
model_mixed <- lme4::lmer(Score ~ Group + (1 | Subject), data = dat_mixed)

moplot_onepred(model_mixed, "Group", style = "publication", caption = FALSE)

## ----glmm, fig.cap = "Population-level probabilities and odds ratios from a binomial mixed model."----
set.seed(550)
dat_glmm <- dat_mixed
eta <- -0.8 + c(0, 0.6, 1.2)[as.integer(dat_glmm$Group)] +
  rnorm(24, sd = 0.8)[as.integer(dat_glmm$Subject)]
dat_glmm$Success <- rbinom(nrow(dat_glmm), size = 1, prob = plogis(eta))
model_glmm <- lme4::glmer(Success ~ Group + (1 | Subject),
                          family = binomial, data = dat_glmm)

moplot_onepred(model_glmm, "Group", coefficient_scale = "interpreted",
               style = "publication", caption = FALSE)

## ----bounds, fig.height = 4.3, fig.cap = "Intervals of a 1-7 rating kept within the scale with bounds = c(1, 7)."----
set.seed(3)
dat_likert <- data.frame(Group = factor(rep(c("A", "B", "C"), each = 8)))
dat_likert$Rating <- c(sample(3:5, 8, TRUE), sample(4:6, 8, TRUE), c(7, 7, 7, 7, 7, 6, 7, 5))
contrasts(dat_likert$Group) <- contr.sum(3)
moplot_onepred(lm(Rating ~ Group, data = dat_likert), "Group",
               bounds = c(1, 7), show_coef = FALSE, caption = FALSE)

## ----caption, fig.show = "hide"-----------------------------------------------
figure <- moplot_onepred(model_helmert, "Group")

## ----extract-param------------------------------------------------------------
extract_param(model_binomial)[, c("Parameter", "Estimate", "Interpreted_Estimate",
                                  "Estimate_type", "p", "d", "d_method")]

## ----export, eval = FALSE-----------------------------------------------------
# moplot_export(figure, "helmert.pdf")

