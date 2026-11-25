## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message = FALSE,
  warning = TRUE,
  fig.width = 9,
  fig.height = 6,
  fig.retina = 1,
  fig.align = "center",
  out.width = "100%"
)
library(MoPlot)

## ----two-way-data-------------------------------------------------------------
set.seed(1)
dat <- expand.grid(Task = factor(c("A", "B", "C")),
                   Session = factor(c("First", "Second")),
                   replicate = 1:40)
dat$Score <- 3 + (dat$Task == "B") * 2 - (dat$Task == "C") +
  (dat$Session == "Second") + 2 * (dat$Task == "C") * (dat$Session == "Second") +
  rnorm(nrow(dat))

contrasts(dat$Task) <- contr.helmert(3)
contrasts(dat$Session) <- contr.sum(2)
model_two <- lm(Score ~ Task * Session, data = dat)

## ----two-way-table------------------------------------------------------------
decomposition <- contrasts_multpred(model_two, configuration = "semantic",
                                    methodological_note = FALSE)
decomposition[, c("coefficient", "split_side", "configuration_semantic",
                  "estimate", "ci_low", "ci_high", "scale_factor")]

## ----two-way-plot, fig.cap = "Two-way interaction with Session as the top-level grouping."----
moplot_multpred(model_two, interaction_info = FALSE, caption = FALSE)

## ----two-way-reversed, fig.cap = "The same interaction with Task as the top-level grouping."----
moplot_multpred(model_two, top_split = "Task", interaction_info = FALSE, caption = FALSE)

## ----two-way-per-difference, fig.cap = "Coefficients divided by their scale factors."----
moplot_multpred(model_two, per_difference = TRUE, interaction_info = FALSE,
                caption = FALSE)

## ----selection, fig.cap = "A single interaction coefficient."-----------------
moplot_multpred(model_two, coefficients = "Task2:Session1", caption = FALSE)

## ----custom, fig.width = 12, fig.height = 6.5, fig.cap = "Interaction between a sum-coded factor and a custom-coded factor."----
set.seed(2)
dat_custom <- expand.grid(Group = factor(c("A", "B", "C")),
                          Stimulus = factor(c("W", "X", "Y", "Z")),
                          replicate = 1:20)
dat_custom$y <- 5 + (dat_custom$Group == "B") * 2 +
  c(0, 1, 2, -2)[as.integer(dat_custom$Stimulus)] +
  3 * (dat_custom$Group == "C") * (dat_custom$Stimulus == "Y") +
  rnorm(nrow(dat_custom))

contrasts(dat_custom$Group) <- contr.sum(3)
contrasts(dat_custom$Stimulus) <- cbind(WX_vs_YZ = c(1, 1, -1, -1),
                                        W_vs_X = c(1, -1, 0, 0),
                                        Y_vs_Z = c(0, 0, 1, -1))
model_custom <- lm(y ~ Stimulus * Group, data = dat_custom)

moplot_multpred(model_custom, coefficients = c("StimulusWX_vs_YZ:Group1",
                                               "StimulusY_vs_Z:Group1"),
                interaction_info = FALSE, base_size = 11, caption = FALSE)

## ----three-way, fig.width = 10, fig.height = 7, fig.cap = "Three-way interaction in tree style."----
set.seed(3)
dat_three <- expand.grid(A = factor(c("a1", "a2", "a3")),
                         B = factor(c("b1", "b2")),
                         C = factor(c("c1", "c2")),
                         replicate = 1:20)
dat_three$y <- 5 + (dat_three$A == "a2") + 1.5 * (dat_three$A == "a3") * (dat_three$B == "b2") +
  (dat_three$C == "c2") * (dat_three$B == "b2") +
  2 * (dat_three$A == "a3") * (dat_three$B == "b2") * (dat_three$C == "c2") +
  rnorm(nrow(dat_three))

contrasts(dat_three$A) <- contr.helmert(3)
contrasts(dat_three$B) <- contr.sum(2)
contrasts(dat_three$C) <- contr.sum(2)
model_three <- lm(y ~ A * B * C, data = dat_three)

moplot_multpred(model_three, style = "tree", interaction_info = FALSE, caption = FALSE)

## ----glmm, fig.cap = "Predicted probabilities and odds ratios for a two-way interaction in a binomial mixed model."----
set.seed(4)
dat_glmm <- expand.grid(Subject = factor(1:40),
                        Condition = factor(c("Control", "Treatment")),
                        Phase = factor(c("Pre", "Post"), levels = c("Pre", "Post")),
                        trial = 1:6)
eta <- -0.5 + 0.8 * (dat_glmm$Condition == "Treatment") * (dat_glmm$Phase == "Post") +
  rnorm(40, sd = 0.7)[as.integer(dat_glmm$Subject)]
dat_glmm$Correct <- rbinom(nrow(dat_glmm), 1, plogis(eta))

contrasts(dat_glmm$Condition) <- contr.sum(2)
contrasts(dat_glmm$Phase) <- MASS::contr.sdif(2)
model_glmm <- lme4::glmer(Correct ~ Condition * Phase + (1 | Subject),
                          family = binomial, data = dat_glmm)

moplot_multpred(model_glmm, coefficient_scale = "interpreted",
                interaction_info = FALSE, caption = FALSE)

## ----caption, fig.show = "hide"-----------------------------------------------
figure <- moplot_multpred(model_two, interaction_info = FALSE)

