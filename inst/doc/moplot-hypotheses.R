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

## ----wrong-coding-------------------------------------------------------------
set.seed(2)
trial <- data.frame(group = factor(rep(c("Control", "T1", "T2"), each = 30),
                                   levels = c("Control", "T1", "T2")))
trial$score <- rnorm(90, mean = c(10, 12, 15)[trial$group], sd = 3)

intended <- rbind("Treat-Control" = c(Control = -1, T1 = 0.5, T2 = 0.5),
                  "T2-Control"    = c(Control = -1, T1 = 0,   T2 = 1))

trial_written <- trial
contrasts(trial_written$group) <- t(intended)
fit_written <- lm(score ~ group, data = trial_written)
round(coef(fit_written)[-1], 2)

## ----note, message = 1:2------------------------------------------------------
written <- contrasts_onepred(fit_written, "group")

## ----check-written, message = 1:2---------------------------------------------
written <- contrasts_onepred(fit_written, "group", hypotheses = intended)
attr(written, "hypothesis_check")

## ----plot-written, fig.cap = "The comparisons estimated by the coding written from the hypotheses."----
moplot_onepred(fit_written, "group", hypotheses = intended, style = "publication")

## ----right-coding, message = 1:2----------------------------------------------
coding <- MASS::ginv(intended)
dimnames(coding) <- list(levels(trial$group), rownames(intended))
round(coding, 2)

trial_derived <- trial
contrasts(trial_derived$group) <- coding
fit_derived <- lm(score ~ group, data = trial_derived)
derived <- contrasts_onepred(fit_derived, "group", hypotheses = intended)
round(coef(fit_derived)[-1], 2)

## ----plot-derived, fig.cap = "The comparisons estimated by the coding derived from the hypothesis matrix."----
moplot_onepred(fit_derived, "group", hypotheses = intended, style = "publication")

## ----hypr, eval = requireNamespace("hypr", quietly = TRUE), message = 1:2-----
h <- hypr::hypr(TreatControl = (T1 + T2) / 2 ~ Control,
                T2Control    = T2 ~ Control,
                levels = levels(trial$group))
trial_hypr <- trial
contrasts(trial_hypr$group) <- hypr::cmat(h)
fit_hypr <- lm(score ~ group, data = trial_hypr)
checked <- contrasts_onepred(fit_hypr, "group", hypotheses = hypr::hmat(h))

## ----helmert, message = 1:2---------------------------------------------------
trial_helmert <- trial
contrasts(trial_helmert$group) <- contr.helmert(3)
fit_helmert <- lm(score ~ group, data = trial_helmert)
plain_differences <- rbind("T1-Control"        = c(-1, 1, 0),
                           "T2-Mean(Control,T1)" = c(-0.5, -0.5, 1))
helmert <- contrasts_onepred(fit_helmert, "group", hypotheses = plain_differences)

## ----interaction, message = 1:3-----------------------------------------------
cells <- data.frame(
  treatment = factor(rep(c("Low", "Medium", "High"), 2), levels = c("Low", "Medium", "High")),
  context = factor(rep(c("Baseline", "Enhanced"), each = 3)),
  mu = c(5, 7, 3, 6, 10, 5))
dat_int <- cells[rep(seq_len(6), each = 5), ]
dat_int$outcome <- dat_int$mu + rep(c(-1, -0.5, 0, 0.5, 1), 6)

intended_treatment <- rbind("Medium-Low"  = c(Low = -1, Medium = 1, High = 0),
                            "High-Medium" = c(Low = 0, Medium = -1, High = 1))
intended_context <- rbind("Enhanced-Baseline" = c(Baseline = -1, Enhanced = 1))

contrasts(dat_int$treatment) <- MASS::ginv(intended_treatment)
contrasts(dat_int$context) <- MASS::contr.sdif(2)
model_int <- lm(outcome ~ treatment * context, data = dat_int)

decomposition <- contrasts_multpred(
  model_int,
  hypotheses = list(treatment = intended_treatment, context = intended_context),
  methodological_note = FALSE
)
attr(decomposition, "hypothesis_check")

