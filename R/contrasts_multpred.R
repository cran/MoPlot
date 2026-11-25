#' Decompose Interaction Contrasts into Interpretable Configurations
#'
#' @description
#' `contrasts_multpred()` decomposes model coefficients associated with
#' interactions among categorical predictors into interpretable
#' configurations of factor levels.
#'
#' For each selected interaction coefficient, the function identifies the
#' positive and negative level groups implied by the contrast matrices used
#' to fit the model. It then returns observed means or model-based
#' predictions for the resulting configurations.
#'
#' One variable involved in each interaction is selected as the
#' \emph{split variable}. The contrasts of the remaining target variable or
#' variables are then displayed within the levels or contrast-defined groups
#' of the split variable.
#'
#' The function supports two-way and higher-order interactions.
#'
#' @param model A fitted model object. Supported model classes include objects
#'   returned by [stats::lm()], [stats::glm()], [lme4::lmer()],
#'   [lme4::glmer()], [nlme::lme()], and [glmmTMB::glmmTMB()].
#'
#' @param data Optional data frame containing the variables used to fit
#'   `model`. If `NULL`, the data supplied in the model call are recovered,
#'   restricted to the rows used for fitting; the model frame is used only when
#'   they are unavailable. The original variables are needed for model-based
#'   estimates when the model contains transformed terms such as `log(x)`.
#'   Rows of a supplied data frame that were excluded from the fit, through
#'   `subset` or missing values, are dropped, matching rows by name.
#'
#' @param interaction Optional character string identifying the interaction
#'   term to decompose, for example `"A:B"` or `"A:B:C"`. When supplied,
#'   the interaction is resolved using the fitted model matrix and the
#'   coefficients associated with that interaction are selected automatically.
#'   If `coefficients` is also supplied, those coefficients are validated
#'   against the requested interaction. When both `interaction` and
#'   `coefficients` are `NULL`, all non-intercept coefficients containing
#'   `":"` are analyzed, preserving the default behavior of this function.
#'
#' @param coefficients Optional character vector containing the interaction
#'   coefficients to decompose. If `interaction` is `NULL` and `coefficients`
#'   is `NULL`, all non-intercept coefficients containing `":"` are analyzed.
#'   When `interaction` is supplied, `coefficients` may be used to select a
#'   subset of the coefficients belonging to that interaction.
#'
#' @param outcome Optional character string indicating the name of the
#'   dependent variable; by default, the dependent variable of `model`. It is
#'   used only for observed estimates, which summarize the response of the
#'   model evaluated in `data`: a transformed response such as `log(y)` is
#'   summarized on the scale of the model, and a logical response or a factor
#'   response of a binomial model is converted to 0/1 as in [stats::glm()].
#'
#' @param split_rule Character string indicating how the split variable is
#'   selected for each interaction coefficient. Available options are:
#'   \describe{
#'     \item{`"last_variable"`}{
#'       Selects the last variable appearing in the interaction coefficient.
#'     }
#'     \item{`"first_variable"`}{
#'       Selects the first variable appearing in the interaction coefficient.
#'     }
#'     \item{`"lowest_complexity"`}{
#'       Selects the variable with the lowest contrast-complexity score.
#'     }
#'     \item{`"highest_complexity"`}{
#'       Selects the variable with the highest contrast-complexity score.
#'     }
#'   }
#'
#' @param top_split Optional character string identifying a variable that must
#'   be used as the split variable. When supplied, `top_split`
#'   overrides `split_rule`.
#'
#' @param sig_level Numeric value indicating the significance threshold used
#'   to classify the original interaction coefficients as significant or
#'   non-significant. The default is `0.05`.
#'
#' @param conf_level Numeric value indicating the confidence level used for
#'   confidence intervals. The default is `0.95`.
#'
#' @param configuration Character string controlling which configuration
#'   labels are included in the output. Available options are:
#'   \describe{
#'     \item{`"none"`}{
#'       Does not return configuration-label columns.
#'     }
#'     \item{`"all"`}{
#'       Returns both explicit and semantic configuration labels.
#'     }
#'     \item{`"explicit"`}{
#'       Returns the complete mathematical composition of each configuration.
#'     }
#'     \item{`"semantic"`}{
#'       Returns a compact description of the levels or groups represented by
#'       each configuration.
#'     }
#'   }
#'
#' @param estimate_type Character string indicating which estimates are
#'   returned. Available options are:
#'   \describe{
#'     \item{`"auto"`}{
#'       Uses observed means when every analyzed interaction involves all the
#'       predictors of the fixed part, the model has no random effects and no
#'       unequal prior weights, and, for generalized models,
#'       `prediction_scale = "response"`. Model-based predictions are used
#'       otherwise, as in [contrasts_onepred()].
#'     }
#'     \item{`"observed"`}{
#'       Returns observed aggregated means computed directly from `data`.
#'     }
#'     \item{`"predicted"`}{
#'       Returns model-based predictions computed from the fixed-effect
#'       coefficients and their variance-covariance matrix.
#'     }
#'   }
#'
#' @param prediction_scale Character string indicating the scale on which
#'   model-based predictions are returned. Available options are:
#'   \describe{
#'     \item{`"response"`}{
#'       Returns predictions on the response scale. For generalized models,
#'       the inverse link function is applied.
#'     }
#'     \item{`"link"`}{
#'       Returns predictions on the linear-predictor scale.
#'     }
#'   }
#'   This argument is ignored when `estimate_type = "observed"`.
#'
#' @param methodological_note Logical. Should a methodological note describing
#'   the settings relevant to the current analysis be printed? The default is
#'   `TRUE`. The note is stored as an attribute of the returned object
#'   regardless of this setting.
#'
#' @param hypotheses Optional named list declaring, for one or more
#'   categorical predictors, the comparisons that their contrast coefficients
#'   are intended to estimate, for example
#'   `list(treatment = H_treatment, context = H_context)`. Each element is a
#'   numeric matrix with one row per contrast column of that predictor and one
#'   column per level, as the rows of a hypothesis matrix (see
#'   [contrasts_onepred()]). Each declared predictor is checked against the
#'   comparisons implied by its coding, and a message reports the result for
#'   each contrast column. Because the hypothesis of an interaction coefficient
#'   is the Kronecker product of the hypotheses of its components, each
#'   interaction coefficient estimates the product of the declared hypotheses
#'   when all its declared components do, multiplied by the product of their
#'   factors. Predictors not listed are not checked. The model is not
#'   modified. The result is stored in the `hypothesis_check` attribute, with
#'   one row per analyzed coefficient.
#'
#' @param bounds Optional numeric vector `c(lower, upper)` giving the natural
#'   range of a response analyzed with a Gaussian (identity link) model, such
#'   as `c(1, 7)` for a Likert score or `c(0, Inf)` for a reaction time. When
#'   supplied, confidence intervals of observed and predicted estimates are computed
#'   on a transformed scale (logit or logarithm) and transformed back, so that
#'   they remain within the bounds. Ignored, with a warning, for generalized
#'   models. See [contrasts_onepred()].
#'
#' @details
#' Each interaction coefficient is separated into its contrast-coded
#' components. The function then identifies:
#'
#' \itemize{
#'   \item the categorical predictor represented by each component;
#'   \item the corresponding contrast column;
#'   \item the factor levels receiving positive, negative, or zero weights.
#' }
#'
#' One predictor is selected as the split variable. The remaining
#' predictor or predictors are treated as target variables whose contrasts are
#' interpreted within the levels or contrast-defined groups of the
#' split variable.
#'
#' For example, suppose that `FactorA` has three levels and uses Helmert
#' contrasts, whereas `FactorB` has two levels and uses sliding-difference
#' contrasts. If `FactorB` is selected as the split variable, the
#' function displays the Helmert contrasts of `FactorA` separately within the
#' levels or contrast-defined groups of `FactorB`.
#'
#' Therefore, the split variable is not necessarily the variable of primary
#' interpretive interest. It defines the conditions under which the
#' contrasts of the remaining target variable or variables are displayed.
#'
#' The split variable is selected separately for each interaction coefficient.
#' Consequently, when all model interactions are analyzed simultaneously,
#' more than one value may appear in the `split_var` column. For example, in a
#' model containing `A:B`, `A:C`, `B:C`, and `A:B:C`, the variable selected by
#' `"highest_complexity"` may differ across coefficients because only the
#' variables involved in each coefficient are considered.
#'
#' Random-effect grouping variables are excluded when identifying the
#' categorical predictors involved in the fixed-effect interactions.
#'
#' @section Contrast complexity and split selection:
#' The options `"lowest_complexity"` and `"highest_complexity"` select the
#' split variable according to the contrast coding assigned to each
#' predictor.
#'
#' Contrast systems are ranked using the following approximate scores:
#'
#' \itemize{
#'   \item Treatment contrasts: `0`.
#'   \item Sum contrasts: `1`.
#'   \item Sliding-difference contrasts: `1`.
#'   \item Helmert contrasts: `2`.
#'   \item Reverse Helmert contrasts: `2`.
#'   \item Polynomial contrasts: `3`.
#'   \item Custom or unidentified contrast matrices: `4`.
#' }
#'
#' These scores affect only the choice of split variable. They do not
#' determine the order in which variables are displayed in the configuration
#' labels.
#'
#' If several variables have the same complexity score, the rightmost tied
#' variable in the interaction coefficient is selected. A warning is issued
#' to document the choice.
#'
#' The argument `top_split` bypasses the automatic rule and forces the
#' specified variable to be used as the split variable.
#'
#' @section Ordering of configuration labels:
#' Both explicit and semantic configuration labels use the following order:
#'
#' \preformatted{
#' target variable(s) -> split variable
#' }
#'
#' The target variables are displayed first, following their stable order
#' among the categorical predictors in the fixed-effects model. The
#' split variable selected by `split_rule` or `top_split` is always
#' displayed last.
#'
#' For a two-way interaction between `A` and `B`, if `A` is the split
#' variable, labels follow the order:
#'
#' \preformatted{
#' B_A
#' }
#'
#' If `B` is the split variable, labels follow the order:
#'
#' \preformatted{
#' A_B
#' }
#'
#' For a three-way interaction among `A`, `B`, and `C`, if `C` is the
#' split variable, labels follow the order:
#'
#' \preformatted{
#' A_B_C
#' }
#'
#' If `A` is the split variable, labels follow the order:
#'
#' \preformatted{
#' B_C_A
#' }
#'
#' This ordering rule affects only the presentation and interpretation of the
#' labels. It does not modify the contrast-complexity scores, the selected
#' split variable, or the numerical estimates.
#'
#' @section Explicit and semantic configurations:
#' `configuration_explicit` reports the individual cells entering an aggregate.
#' It is intended to show the exact mathematical composition of the displayed
#' estimate.
#'
#' For example:
#'
#' \preformatted{
#' (B2_A1 + B2_A2)/2
#' }
#'
#' indicates that the estimate is the arithmetic mean of cells `B2_A1` and
#' `B2_A2`.
#'
#' `configuration_semantic` reports the corresponding interpretation:
#'
#' \preformatted{
#' B2_Mean(A1,A2)
#' }
#'
#' Both labels use the same variable order. In this example, `B` is the target
#' variable and `A` is the split variable.
#'
#' With `B` as the split variable, the corresponding labels could be:
#'
#' \preformatted{
#' (A1_B2 + A2_B2)/2
#' Mean(A1,A2)_B2
#' }
#'
#' The explicit label answers the question, "Which cells were combined?",
#' whereas the semantic label answers the question, "What does the resulting
#' estimate represent?"
#'
#' @section Raw estimates:
#' For `estimate_type = "observed"`, estimates are computed directly from the
#' observed response values.
#'
#' At the individual-cell level, the function returns the observed mean,
#' standard error, and confidence interval.
#'
#' When several disjoint cells are combined, their means are aggregated using
#' normalized contrast weights. Standard errors are combined through classical
#' variance propagation.
#'
#' Confidence intervals for aggregated observed estimates use a
#' Satterthwaite-like approximation to the degrees of freedom when sufficient
#' information is available.
#'
#' Confidence intervals of observed means respect the natural range of the
#' response. For proportions (binary outcomes or binomial models), the Wilson
#' score interval is used; for averages of several cells it is
#' computed with the effective sample size \eqn{\hat{p}(1 - \hat{p}) /
#' \mathrm{SE}^2}, which reduces to the classical Wilson interval for a single
#' cell. For positive responses analyzed with count or positive-continuous
#' families (for example Poisson, negative binomial, or Gamma), the interval
#' is computed on the log scale, \eqn{\exp(\log m \pm z\,\mathrm{SE}/m)}; a
#' mean of zero receives the exact Poisson upper bound
#' \eqn{-\log(\alpha/2)/n}. Other responses use symmetric intervals.
#'
#' @section Model-based predictions:
#' For `estimate_type = "predicted"`, predictions are computed from the
#' fixed-effect model:
#'
#' \deqn{
#'   \hat{\eta} = X\hat{\beta}
#' }
#'
#' where \eqn{X} is the fixed-effect model matrix and \eqn{\hat{\beta}} is the
#' vector of fixed-effect coefficients.
#'
#' The covariance matrix of the linear predictors is obtained from:
#'
#' \deqn{
#'   \mathrm{Var}(\hat{\eta})
#'   =
#'   X V_{\hat{\beta}} X^\mathsf{T}
#' }
#'
#' where \eqn{V_{\hat{\beta}}} is the complete fixed-effect
#' variance-covariance matrix.
#'
#' In mixed-effects models, random effects are excluded. The returned
#' predictions are therefore population-level or fixed-effect predictions.
#'
#' Predictions are computed on a reference grid built for each coefficient.
#' The grid crosses the levels of the factors involved in the coefficient with
#' the levels of the other factors; numeric predictors are fixed at their
#' observed means. Each displayed cell averages over the other factors with
#' equal weights (a balanced average, independent of the observed cell
#' sizes). A factor that interacts with the variables of the coefficient, such
#' as `C` for a coefficient of `A:B` in a model containing `A:B:C`, is instead
#' held at the combination of its levels at which the coefficient is defined,
#' given by the first row of its hypothesis matrix: the balanced average for
#' sum, Helmert, and sliding-difference coding, and the reference level for
#' treatment coding.
#'
#' Averaging the predictions of the observed rows of each cell would instead
#' mix in the effects of covariates whose distribution differs between cells,
#' so that the displayed interaction contrast would no longer correspond to
#' the coefficient.
#'
#' For `glmmTMB` models, only the conditional component is analyzed.
#' Zero-inflation and dispersion components are not currently decomposed.
#'
#' @section Response-scale predictions:
#' For generalized models with `prediction_scale = "response"`, each
#' linear-predictor value is first transformed through the inverse link
#' function. The transformed predictions are then averaged.
#'
#' Therefore, the function computes:
#'
#' \deqn{
#'   \frac{1}{n}
#'   \sum_{i=1}^{n}
#'   g^{-1}(\hat{\eta}_i)
#' }
#'
#' rather than:
#'
#' \deqn{
#'   g^{-1}\left(
#'     \frac{1}{n}
#'     \sum_{i=1}^{n}
#'     \hat{\eta}_i
#'   \right)
#' }
#'
#' These two quantities are generally different when the inverse link is
#' nonlinear.
#'
#' @section Aggregation of predicted estimates:
#' Model-based predictions for different configurations are generally
#' correlated because they depend on the same fixed-effect coefficient vector.
#'
#' `contrasts_multpred()` therefore does not combine predicted standard errors
#' under an assumption of independence.
#'
#' Instead, each configuration is represented by its weights over the rows of
#' the reference grid, which are retained throughout the recursive
#' decomposition. When configurations are aggregated, their weights are
#' combined and the aggregate is calculated directly from the complete
#' fixed-effect variance-covariance matrix.
#'
#' On the link scale, an aggregated prediction with normalized weights
#' \eqn{a_i} is:
#'
#' \deqn{
#'   \hat{\theta}
#'   =
#'   \sum_i a_i X_i\hat{\beta}
#' }
#'
#' Defining:
#'
#' \deqn{
#'   L = \sum_i a_i X_i
#' }
#'
#' its variance is:
#'
#' \deqn{
#'   \mathrm{Var}(\hat{\theta})
#'   =
#'   L V_{\hat{\beta}} L^\mathsf{T}
#' }
#'
#' On the response scale, the aggregated prediction is:
#'
#' \deqn{
#'   \hat{\theta}
#'   =
#'   \sum_i
#'   a_i
#'   g^{-1}(X_i\hat{\beta})
#' }
#'
#' Its variance is approximated using the multivariate delta method. The
#' gradient is:
#'
#' \deqn{
#'   \nabla\hat{\theta}
#'   =
#'   \sum_i
#'   a_i
#'   g^{-1\prime}(\hat{\eta}_i)
#'   X_i
#' }
#'
#' and:
#'
#' \deqn{
#'   \mathrm{Var}(\hat{\theta})
#'   \approx
#'   \nabla\hat{\theta}
#'   V_{\hat{\beta}}
#'   \nabla\hat{\theta}^{\mathsf{T}}
#' }
#'
#' This procedure preserves the covariance among the predictions entering
#' each aggregate.
#'
#' On the response scale, the confidence interval of a prediction or of an
#' average of predictions \eqn{\hat{\theta}} is computed on the link scale and
#' transformed back,
#'
#' \deqn{g^{-1}\left(g(\hat{\theta}) \pm z \,
#'   \frac{\mathrm{SE}(\hat{\theta})}{\left|g^{-1\prime}(g(\hat{\theta}))\right|}\right),}
#'
#' so that its bounds remain within the range of the response (for example
#' between 0 and 1 for probabilities). For a single prediction this is the
#' interval computed on the linear predictor and transformed with the inverse
#' link. The interval is therefore generally asymmetric around the estimate,
#' whereas the reported standard error refers to the response scale.
#'
#' @section Supported contrast types:
#' The function supports treatment contrasts (with any reference level), sum
#' contrasts (with any omitted level), Helmert, reverse Helmert, and
#' sliding-difference contrasts, and custom contrast matrices.
#'
#' For every factor, the levels compared by each contrast column and their
#' weights are taken from the hypothesis matrix \eqn{H = [1, C]^{-1}} of its
#' coding matrix \eqn{C}, as described in [contrasts_onepred()]. Reading them
#' from the signs of \eqn{C} is correct only for orthogonal codings. Groups of
#' levels with unequal weights, which arise with non-orthogonal custom
#' codings, are labelled `W(a1 .20, a3 .80)`, and a message describes the
#' comparison actually estimated when it differs from the one suggested by
#' the signs of the coding column.
#'
#' Sum contrasts are interpreted as a focal level relative to the unweighted
#' mean across all levels of the corresponding factor.
#'
#' @section Scale factors:
#' The hypothesis of an interaction coefficient is the Kronecker product of
#' the hypothesis rows of its components. Each interaction coefficient
#' therefore equals the product of the scale factors of its components,
#' returned in `scale_factor`, multiplied by the displayed difference of
#' differences on the linear-predictor scale. For example, the interaction
#' between the second column of a Helmert-coded factor (\eqn{s = 1/3}) and a
#' sum-coded factor (\eqn{s = 1}) has \eqn{s = 1/3}.
#'
#' The function verifies this identity for each coefficient. For predicted
#' estimates it holds exactly unless a variable of the coefficient interacts
#' with a numeric predictor, in which case a warning is issued. For observed
#' estimates it holds only when the model contains the interacting factors
#' alone.
#'
#' Polynomial contrasts are not currently supported for interaction
#' decomposition because all levels generally contribute to each polynomial
#' component through nonzero weights and require a different interpretive
#' strategy.
#'
#' @section P-values:
#' P-values refer to the original interaction coefficients being decomposed.
#' They are not direct tests of the individual configuration-specific
#' estimates displayed in the output.
#'
#' The function recognizes coefficient-table columns named:
#'
#' \itemize{
#'   \item `"Pr(>|t|)"`;
#'   \item `"Pr(>|z|)"`;
#'   \item `"p-value"`.
#' }
#'
#' Models fitted with [lme4::lmer()] do not ordinarily provide p-values.
#' When `lmerTest` is installed, the model is converted internally using
#' [lmerTest::as_lmerModLmerTest()], and approximate p-values are extracted
#' from the resulting coefficient table.
#'
#' If no p-values are available, `p_value` is set to `NA` and
#' `significance` is set to `"NoPValue"`.
#'
#' @return
#' A data frame of class `"plot_ready_contrasts_df"` with one row per
#' decomposed configuration.
#'
#' The returned columns include:
#'
#' \describe{
#'   \item{`coefficient`}{
#'     Name of the original interaction coefficient being decomposed.
#'   }
#'   \item{`contrast_id`}{
#'     Sequential identifier assigned to the selected coefficient.
#'   }
#'   \item{`significance`}{
#'     Significance classification of the original coefficient:
#'     `"Sig"`, `"NotSig"`, or `"NoPValue"`.
#'   }
#'   \item{`p_value`}{
#'     P-value associated with the original interaction coefficient, when
#'     available.
#'   }
#'   \item{`scale_factor`}{
#'     Scale factor such that the coefficient equals `scale_factor`
#'     multiplied by the displayed difference of differences on the
#'     linear-predictor scale.
#'   }
#'   \item{`split_var`}{
#'     Split variable selected for the coefficient.
#'   }
#'   \item{factor columns}{
#'     One column for each categorical predictor involved in the analyzed
#'     interactions. Values identify specific levels or aggregated groups.
#'   }
#'   \item{`n`}{
#'     Number of observed rows underlying the displayed configuration.
#'   }
#'   \item{`estimate`}{
#'     Observed aggregate or model-based aggregate.
#'   }
#'   \item{`se`}{
#'     Standard error of the displayed estimate.
#'   }
#'   \item{`ci_low`}{
#'     Lower confidence-interval boundary.
#'   }
#'   \item{`ci_high`}{
#'     Upper confidence-interval boundary.
#'   }
#' }
#'
#' Depending on `configuration`, the output may also include:
#'
#' \describe{
#'   \item{`configuration_explicit`}{
#'     Mathematical description of the cells entering the aggregate.
#'   }
#'   \item{`configuration_semantic`}{
#'     Compact interpretation of the specific and averaged levels.
#'   }
#' }
#'
#' The returned object contains the following attributes:
#'
#' \describe{
#'   \item{`note`}{
#'     Methodological description of the estimation and aggregation
#'     procedures.
#'   }
#'   \item{`estimate_type`}{
#'     Whether the output contains `"observed"` or `"predicted"` estimates.
#'   }
#'   \item{`prediction_scale`}{
#'     Prediction scale used for model-based estimates. This is `NA` for observed
#'     estimates.
#'   }
#'   \item{`estimate_info`}{
#'     Human-readable description of the `estimate` column.
#'   }
#'   \item{`pvalue_source`}{
#'     Source of the p-values, such as `"model"`, `"lmerTest"`, or
#'     `"not_available"`.
#'   }
#'   \item{`split_rule`}{
#'     Rule used to select the split variable.
#'   }
#'   \item{`top_split`}{
#'     Split variable explicitly requested through `top_split`, when
#'     supplied.
#'   }
#'   \item{`scale_factors`}{
#'     Named vector of scale factors, one per coefficient.
#'   }
#'   \item{`identity_holds`}{
#'     Named logical vector indicating, for each coefficient, whether the
#'     displayed comparisons reconstruct the coefficient on the
#'     linear-predictor scale.
#'   }
#' }
#'
#' @examples
#' ###############################################################
#' ## Two-way interaction
#' ###############################################################
#' \donttest{
#'
#' library(lme4)
#' library(MASS)
#'
#' set.seed(123)
#'
#' dat <- expand.grid(
#'   ID = factor(paste0("S", 1:40)),
#'   A = factor(c("A1", "A2", "A3")),
#'   B = factor(c("B1", "B2", "B3", "B4")),
#'   Rep = 1:2
#' )
#'
#' contrasts(dat$A) <- stats::contr.helmert(3)
#' contrasts(dat$B) <- MASS::contr.sdif(4)
#'
#' random_intercept <- rnorm(
#'   nlevels(dat$ID),
#'   mean = 0,
#'   sd = 1
#' )
#'
#' effect_A <- c(
#'   A1 = 0,
#'   A2 = 1,
#'   A3 = 2
#' )
#'
#' effect_B <- c(
#'   B1 = 0,
#'   B2 = 0.5,
#'   B3 = 1.5,
#'   B4 = 2.5
#' )
#'
#' dat$Score <- 10 +
#'   effect_A[dat$A] +
#'   effect_B[dat$B] +
#'   0.7 * effect_A[dat$A] * effect_B[dat$B] +
#'   random_intercept[dat$ID] +
#'   rnorm(nrow(dat), mean = 0, sd = 1.5)
#'
#' mod <- lme4::lmer(
#'   Score ~ A * B + (1 | ID),
#'   data = dat
#' )
#'
#' # A is selected as the conditioning variable.
#' # Labels follow the order B -> A.
#' out_high <- contrasts_multpred(
#'   model = mod,
#'   outcome = "Score",
#'   estimate_type = "predicted",
#'   split_rule = "highest_complexity",
#'   configuration = "all"
#' )
#'
#' # B is selected as the conditioning variable.
#' # Labels follow the order A -> B.
#' out_low <- contrasts_multpred(
#'   model = mod,
#'   outcome = "Score",
#'   estimate_type = "predicted",
#'   split_rule = "lowest_complexity",
#'   configuration = "all"
#' )
#'
#'
#' ###############################################################
#' ## Three-way interaction
#' ###############################################################
#'
#' dat3 <- expand.grid(
#'   ID = factor(paste0("S", 1:30)),
#'   A = factor(c("A1", "A2", "A3")),
#'   B = factor(c("B1", "B2", "B3", "B4")),
#'   C = factor(c("C1", "C2")),
#'   Rep = 1:2
#' )
#'
#' contrasts(dat3$A) <- stats::contr.helmert(3)
#' contrasts(dat3$B) <- MASS::contr.sdif(4)
#' contrasts(dat3$C) <- stats::contr.treatment(2)
#'
#' dat3$Score <- 10 +
#'   effect_A[dat3$A] +
#'   effect_B[dat3$B] +
#'   1.2 * (dat3$C == "C2") +
#'   0.4 * effect_A[dat3$A] * effect_B[dat3$B] +
#'   rnorm(nrow(dat3), mean = 0, sd = 1.5)
#'
#' mod3 <- stats::lm(
#'   Score ~ A * B * C,
#'   data = dat3
#' )
#'
#' out3 <- contrasts_multpred(
#'   model = mod3,
#'   outcome = "Score",
#'   estimate_type = "predicted",
#'   split_rule = "highest_complexity",
#'   configuration = "all"
#' )
#'
#'
#' ###############################################################
#' ## glmmTMB model
#' ###############################################################
#'
#' eta_count <- 1.2 +
#'   0.2 * effect_A[dat$A] +
#'   0.15 * effect_B[dat$B]
#'
#' dat$count <- rpois(
#'   nrow(dat),
#'   lambda = exp(eta_count)
#' )
#'
#' mod_tmb <- glmmTMB::glmmTMB(
#'   count ~ A * B + (1 | ID),
#'   family = poisson,
#'   data = dat
#' )
#'
#' out_tmb <- contrasts_multpred(
#'   model = mod_tmb,
#'   outcome = "count",
#'   estimate_type = "predicted",
#'   prediction_scale = "response",
#'   split_rule = "highest_complexity",
#'   configuration = "semantic"
#' )
#' }
#'
#' @export


contrasts_multpred <- function(model,
                               data = NULL,
                               interaction = NULL,
                               coefficients = NULL,
                               outcome = NULL,
                               split_rule = c(
                                 "last_variable",
                                 "first_variable",
                                 "lowest_complexity",
                                 "highest_complexity"
                               ),
                               top_split = NULL,
                               sig_level = 0.05,
                               conf_level = 0.95,
                               configuration = c("none", "all", "explicit", "semantic"),
                               estimate_type = c("auto", "observed", "predicted"),
                             prediction_scale = c("response", "link"),
                             methodological_note = TRUE,
                             bounds = NULL,
                             hypotheses = NULL) {

  # ============================================================
  # 1. Argument matching and basic checks
  # ============================================================
  #
  # The function supports:
  #   - lm
  #   - glm
  #   - lmer/glmer models from lme4
  #
  # It is designed for models that include at least one interaction
  # involving categorical predictors.
  # ============================================================

  split_rule <- match.arg(split_rule)
  configuration <- match.arg(configuration)
  estimate_type <- match.arg(estimate_type)
  prediction_scale <- match.arg(prediction_scale)

  if (!.is_supported_model(model)) {

    stop(
      paste0(
        "Unsupported model class. Supported model classes include ",
        "lm, glm, lme4::lmer/glmer, lmerTest::lmer, nlme::lme, ",
        "and glmmTMB."
      ),
      call. = FALSE
    )
  }

  if (is.null(data)) {

    data <- .extract_analysis_data(model)

    if (is.null(data)) {
      stop(
        paste0(
          "`data` was not supplied and the analysis data could not be ",
          "recovered from `model`. Please provide the data frame used to fit ",
          "the model explicitly through `data`."
        ),
        call. = FALSE
      )
    }
  }

  if (!is.data.frame(data)) {
    stop("'data' must be a data frame.",call. = FALSE)
  }

  # Supplied data are restricted to the rows used for fitting and given
  # the contrasts used for fitting.
  data <- .align_data_with_model(model, data)

  if (is.null(outcome)) {
    outcome <- .outcome_name(model)
  }

  # Natural range declared by the user for Gaussian responses.
  bounds <- .validate_bounds(bounds, model)

  bound_critical <- stats::qnorm(1 - (1 - conf_level) / 2)

  # Model-based intervals within the bounds (identity-link models only).
  apply_bounds <- function(statistics) {

    if (!is.null(bounds)) {

      bounded <- .bounded_interval(statistics$estimate, statistics$se, bounds, bound_critical)

      statistics$ci_low <- bounded$ci_low
      statistics$ci_high <- bounded$ci_high
    }

    statistics
  }

  # ============================================================
  # 2. Extract model information
  # ============================================================
  #
  # This function extracts:
  #   - the model matrix;
  #   - coefficient names;
  #   - p-values, when available;
  #   - significance labels.
  #
  # For lmer models fitted with lme4, p-values are usually absent.
  # If lmerTest is installed, the function attempts to recover them.
  # ============================================================

  extract_model_info <- function(model,data, sig_level = 0.05) {

    fixed_effects <- .extract_fixef(model)

    coefficient_names <- setdiff(names(fixed_effects),"(Intercept)")

    p_information <- .extract_p_values(model = model,coefficient_names = coefficient_names, alpha = sig_level)

    list(
      coef_names = coefficient_names,
      sig_lookup = p_information$table,
      pvalue_source = p_information$source
    )
  }

  model_info <- extract_model_info(model, data, sig_level = sig_level)
  coef_names <- model_info$coef_names

  # Select the interaction coefficients to analyze.
  #
  # Backward compatibility is preserved:
  #   - interaction = NULL, coefficients = NULL:
  #       analyze all interaction coefficients among factors;
  #   - interaction = NULL, coefficients supplied:
  #       analyze exactly the requested coefficients.
  #
  # Interactions involving numeric variables are not decomposed; see
  # .extract_moplot_interaction_map().
  #
  # When 'interaction' is supplied, reuse the same model-matrix-based
  # selection logic used by moplot_multpred(). This maps an interaction
  # term (e.g., "A:B:C") to all coefficient columns belonging to that
  # term and, when 'coefficients' is also supplied, validates that the
  # requested coefficients belong to the selected interaction.
  if (!is.null(interaction)) {

    selection <- .resolve_moplot_interaction_selection(
      model = model,
      interaction = interaction,
      coefficients = coefficients,
      interaction_info = FALSE
    )

    interaction <- selection$interaction
    coefficients <- selection$coefficients

  } else {

    factor_coefficients <- .extract_moplot_interaction_map(model)$coefficient

    if (is.null(coefficients)) {

      coefficients <- factor_coefficients

    } else {

      unsupported_coefficients <- setdiff(
        intersect(coefficients, coef_names),
        factor_coefficients
      )

      if (length(unsupported_coefficients) > 0L) {
        stop(
          paste0(
            "The following coefficients do not belong to interactions among ",
            "factors: ",
            paste(unsupported_coefficients, collapse = ", "),
            ". Interactions involving numeric variables are not supported."
          ),
          call. = FALSE
        )
      }
    }
  }

  if (length(coefficients) == 0) {
    stop("No interaction coefficient to analyze.")
  }

  # Automatic estimate type: observed means only when every analyzed
  # interaction involves all the predictors of the fixed part, the model
  # has no random effects, and, for generalized models, the response
  # scale is requested. See .resolve_estimate_type().
  if (identical(estimate_type, "auto")) {

    interaction_map <- .extract_moplot_interaction_map(model)

    analyzed_interactions <- unique(interaction_map$interaction[interaction_map$coefficient %in% coefficients])

    resolved_types <- vapply(
      analyzed_interactions,
      function(term) {
        .resolve_estimate_type(
          model = model,
          predictor = all.vars(stats::as.formula(paste("~", term))),
          estimate_type = "auto",
          prediction_scale = prediction_scale
        )
      },
      character(1)
    )

    estimate_type <- if (length(resolved_types) > 0L && all(resolved_types == "observed")) "observed" else "predicted"
  }

  # Observed estimates summarize the response of the model evaluated in
  # 'data': a transformed response such as log(y) is summarized on the
  # scale of the model, and logical responses or factors of binomial
  # models are converted to 0/1. The values are stored in a reserved
  # column so that they follow every subset of the data. Model-based
  # estimates do not use the response values.
  response_column <- ".moplot_response"

  response_values <- if (estimate_type == "observed") {
    .response_values(
      model = model,
      data = data,
      outcome = outcome
    )
  } else {
    NULL
  }

  if (estimate_type == "observed") {
    data[[response_column]] <- response_values
    .check_values_within_bounds(response_values, bounds)
  }

  # Identify factor predictors involved in the model.
  fixed_formula <- .extract_fixed_formula(model)

  model_vars <- all.vars(fixed_formula)

  predictor_vars <- setdiff(model_vars,outcome)

  predictor_vars <- predictor_vars[predictor_vars %in% names(data)]

  factor_vars <- predictor_vars[vapply(data[predictor_vars],is.factor,logical(1))]

  if (length(factor_vars) == 0) {stop("There are no factor predictors in the model.")}

  # Natural range of the response, used for the intervals of raw
  # means (see .observed_interval()).
  response_range_type <- .response_range_type(
    family_object = .extract_family(model),
    values = response_values
  )

  # ============================================================
  # 3. Contrast hypotheses and classification
  # ============================================================
  #
  # The hypothesis matrix of each factor, H = [1, C]^-1, states what
  # each contrast coefficient compares. It is the source of truth for
  # the levels and weights on the positive and negative side of every
  # component; the coding matrix is used only to classify the contrast
  # type. See .contrast_hypotheses().
  # ============================================================

  hypotheses_cache <- list()

  get_variable_hypotheses <- function(var) {

    if (is.null(hypotheses_cache[[var]])) {

      if (!is.factor(data[[var]])) {
        stop(sprintf("The variable '%s' is not a factor in 'data'.", var), call. = FALSE)
      }

      contrast_matrix <- stats::contrasts(data[[var]])

      if (is.null(contrast_matrix)) {
        stop(sprintf("The variable '%s' has no available contrast matrix.", var), call. = FALSE)
      }

      hypotheses_cache[[var]] <<- .contrast_hypotheses(
        contrast_matrix = contrast_matrix,
        level_names = levels(data[[var]])
      )
    }

    hypotheses_cache[[var]]
  }

  variable_types <- vapply(
    factor_vars,
    function(var) {
      contrast_matrix <- stats::contrasts(data[[var]])
      if (is.null(contrast_matrix)) NA_character_ else .classify_contrast(contrast_matrix)$type
    },
    character(1)
  )

  # ============================================================
  # 3.1 Declared hypotheses
  # ============================================================
  #
  # The comparisons intended for one or more predictors are checked
  # against those implied by their coding. Interaction coefficients
  # inherit the result of their components (see the per-coefficient
  # check in decompose_one_coefficient()). The model is not modified.
  # ============================================================

  factor_checks <- list()

  if (!is.null(hypotheses)) {

    if (!is.list(hypotheses) ||
        is.data.frame(hypotheses) ||
        is.null(names(hypotheses)) ||
        any(!nzchar(names(hypotheses)))) {

      stop(
        paste0(
          "`hypotheses` must be a named list with one matrix per ",
          "categorical predictor, for example list(A = H_A, B = H_B)."
        ),
        call. = FALSE
      )
    }

    unknown_variables <- setdiff(names(hypotheses), factor_vars)

    if (length(unknown_variables) > 0L) {

      stop(
        paste0(
          "`hypotheses` refers to variables that are not categorical ",
          "predictors of the model: ", paste(unknown_variables, collapse = ", "), "."
        ),
        call. = FALSE
      )
    }

    for (var in names(hypotheses)) {

      factor_checks[[var]] <- .check_declared_hypotheses(
        declared = hypotheses[[var]],
        hypothesis_matrix = get_variable_hypotheses(var)$hypothesis_matrix,
        level_names = levels(data[[var]])
      )

      factor_checks[[var]]$messages <- sub(
        "^Contrast_([0-9]+)",
        paste0("Variable ", var, ", contrast \\1"),
        factor_checks[[var]]$messages
      )
    }
  }

  coefficient_checks <- list()

  # ============================================================
  # 4. Polynomial contrasts are not yet supported
  # ============================================================
  #
  # Interactions involving polynomial contrasts require a different
  # decomposition strategy because all levels contribute to each
  # contrast through polynomial weights.
  # ============================================================

  poly_vars <- factor_vars[variable_types %in% "Polynomial"]

  if (length(poly_vars) > 0) {
    stop(
      paste0(
        "Work in progress: this function does not yet support categorical variables ",
        "coded with polynomial contrasts. Variables involved: ",
        paste(poly_vars, collapse = ", ")
      )
    )
  }

  # ============================================================
  # 5. Contrast complexity score
  # ============================================================
  #
  # Used when split_rule is:
  #   - "lowest_complexity"
  #   - "highest_complexity"
  #
  # Lower values indicate simpler contrast codings.
  # Higher values indicate more complex or custom contrast codings.
  # ============================================================

  contrast_complexity <- function(var) {

    contrast_type <- variable_types[[var]]

    if (is.na(contrast_type)) { return(Inf) }

    switch(
      contrast_type,
      "Treatment" = 0,
      "Sum" = 1,
      "Sliding Differences" = 1,
      "Helmert" = 2,
      "Reverse Helmert" = 2,
      "Polynomial" = 3,
      # Custom or unidentified contrast matrices are assigned the
      # highest complexity score.
      4
    )
  }

  # ============================================================
  # 6. Methodological note
  # ============================================================

  disclaimer_note <- .build_multpred_methodological_note(
    model = model,
    estimate_type = estimate_type,
    prediction_scale = prediction_scale,
    split_rule = split_rule,
    top_split = top_split,
    configuration = configuration
  )

  if (isTRUE(methodological_note)) { message(disclaimer_note) }

  # ============================================================
  # 7. Build component lookup
  # ============================================================
  #
  # Each model coefficient is made of components.
  #
  # Example:
  #   A1:B2-1
  #
  # has two components:
  #   A1
  #   B2-1
  #
  # This lookup maps each component to:
  #   - the original variable;
  #   - the local contrast column.
  # ============================================================

  build_component_lookup <- function(var,
                                     data) {

    variable_formula <- stats::reformulate(var)

    terms_object <- stats::delete.response(stats::terms(variable_formula))

    variable_frame <- stats::model.frame(
      terms_object,
      data = data,
      drop.unused.levels = FALSE,
      na.action = stats::na.pass
    )

    # The contrast matrix is already stored in the factor contained
    # in 'data'. Therefore, it is not passed again through
    # contrasts.arg, avoiding duplicated contrast assignments and
    # the associated "contrasts dropped from factor" warning.
    model_matrix <- stats::model.matrix(terms_object,data = variable_frame)

    model_matrix <- model_matrix[ ,colnames(model_matrix) != "(Intercept)", drop = FALSE]

    if (ncol(model_matrix) == 0L) { return(NULL) }

    data.frame(
      component = colnames(model_matrix),
      variable = var,
      local_col = seq_len(ncol(model_matrix)),
      stringsAsFactors = FALSE
    )
  }

  lookup_list <- lapply(factor_vars, build_component_lookup, data = data)
  lookup <- do.call(rbind, lookup_list)

  if (is.null(lookup) || nrow(lookup) == 0) {
    stop("It was not possible to build the contrast-component lookup table.")
  }

  # ============================================================
  # 8. Data subsetting helper
  # ============================================================
  #
  # Subsets the original data according to a named list of
  # conditions, e.g. list(A = "A1", B = "B2").
  # ============================================================

  subset_by_conditions <- function(data, conditions) {

    out <- data

    if (length(conditions) == 0) return(out)

    for (nm in names(conditions)) {
      out <- out[out[[nm]] %in% conditions[[nm]], , drop = FALSE]
    }

    out
  }

  # ============================================================
  # 9. Label helpers
  # ============================================================
  #
  # These functions create readable labels describing which
  # configurations of levels are being averaged.
  # ============================================================

  collapse_source <- function(labels) {
    paste(unique(labels), collapse = "|||")
  }

  # Groups of levels with unequal hypothesis weights (non-orthogonal
  # custom codings) are labelled "W(a1 .20, a3 .80)". The labels are
  # registered while the current coefficient is decomposed; within one
  # coefficient each variable contributes a single contrast column,
  # so a set of levels identifies its group unambiguously.
  group_label_registry <- list()

  group_key <- function(var, levels) {
    paste0(var, "::", paste(sort(levels), collapse = "|"))
  }

  register_group_label <- function(var, weights) {

    if (length(weights) > 1L &&
        !isTRUE(all.equal(unname(weights), rep(1 / length(weights), length(weights))))) {

      group_label_registry[[group_key(var, names(weights))]] <<- .format_weighted_group(weights)
    }
  }

  group_token <- function(var, values) {

    label <- group_label_registry[[group_key(var, values)]]

    if (!is.null(label)) {
      return(label)
    }

    paste0("Mean(", paste(values, collapse = ","), ")")
  }

  make_configuration_explicit <- function(labels, weights = NULL) {

    if (!is.null(weights) &&
        length(weights) == length(labels) &&
        length(labels) > 1L) {

      weights <- weights / sum(weights)

      if (!isTRUE(all.equal(weights, rep(1 / length(weights), length(weights))))) {

        formatted_weights <- sub("^0\\.", ".", formatC(weights, format = "f", digits = 2))

        return(
          paste0("(", paste(paste0(formatted_weights, "*", labels), collapse = " + "), ")")
        )
      }
    }

    labs <- unique(labels)

    if (length(labs) == 1) return(labs)

    paste0("(", paste(labs, collapse = " + "), ")/", length(labs))
  }

  parse_named_label <- function(x) {

    parts <- strsplit(x, "\\|")[[1]]
    vars <- sub("=.*$", "", parts)
    vals <- sub("^.*=", "", parts)

    stats::setNames(vals, vars)
  }

  # named_to_display <- function(x) {
  #
  #   vals <- parse_named_label(x)
  #   vals <- vals[factor_vars[factor_vars %in% names(vals)]]
  #
  #   paste(unname(vals), collapse = "_")
  # }

  # ============================================================
  # Configuration-variable ordering
  # ============================================================
  #
  # This helper defines a stable display order for configuration
  # labels.
  #
  # All target contrast variables are displayed first, following
  # their order among the categorical predictors in the model.
  #
  # The split variable is the conditioning variable and is therefore
  # always displayed last.
  # ============================================================

  order_configuration_variables <- function(variables,
                                            split_var = NULL) {

    variables <- unique(variables)

    # Preserve the stable order of categorical predictors in the
    # fixed-effects model.
    ordered_variables <- factor_vars[factor_vars %in% variables]

    # Retain any unexpected variable that is not included in
    # factor_vars.
    remaining_variables <- setdiff(variables,ordered_variables)

    ordered_variables <- c(ordered_variables,remaining_variables)

    # The split variable is the conditioning variable and is placed
    # after all target contrast variables.
    if (!is.null(split_var) &&
        split_var %in% ordered_variables) {
          ordered_variables <- c(ordered_variables[ordered_variables != split_var],split_var)
    }

    ordered_variables
  }

  make_configuration_semantic <- function(named_labels,
                                          split_var = NULL) {

    labels <- unique(named_labels)

    parsed_labels <- lapply(labels,parse_named_label)

    variables_present <- unique(unlist(lapply(parsed_labels,names)))

    variable_order <- order_configuration_variables(
      variables = variables_present,
      split_var = split_var
    )

    same_variables <- all(
      vapply(
        parsed_labels,
        function(values) {
          all(
            variable_order %in% names(values)
          )
        },
        logical(1)
      )
    )

    if (!same_variables) {

      return(
        paste0(
          "Mean(",
          length(labels),
          " cells)"
        )
      )
    }

    value_matrix <- do.call(
      rbind,
      lapply(
        parsed_labels,
        function(values) {
          unname(
            values[variable_order]
          )
        }
      )
    )

    colnames(value_matrix) <- variable_order

    varying_variables <- apply(
      value_matrix,
      2L,
      function(values) {
        length(unique(values)) > 1L
      }
    )

    # No variable varies across the configurations.
    if (!any(varying_variables)) {

      return(
        paste(
          value_matrix[1L, ],
          collapse = "_"
        )
      )
    }

    semantic_tokens <- vapply(
      variable_order,
      function(variable_name) {

        values <- unique(
          value_matrix[, variable_name]
        )

        if (length(values) == 1L) {

          return(values)
        }

        group_token(variable_name, values)
      },
      character(1)
    )

    paste(
      semantic_tokens,
      collapse = "_"
    )
  }




  make_variable_columns <- function(named_labels, all_vars) {

    labs <- unique(named_labels)
    parsed <- lapply(labs, parse_named_label)

    out <- stats::setNames(as.list(rep("", length(all_vars))), all_vars)

    vars_present <- unique(unlist(lapply(parsed, names)))
    vars_present <- vars_present[vars_present %in% all_vars]

    for (v in vars_present) {

      vals <- unique(vapply(parsed, function(z) z[[v]], character(1)))

      if (length(vals) == 1) {
        out[[v]] <- vals
      } else {
        out[[v]] <- group_token(v, vals)
      }
    }

    as.data.frame(out, stringsAsFactors = FALSE)
  }




  make_ordered_label_info <- function(conditions,
                                      split_var = NULL) {

    if (length(conditions) == 0L) {

      return(
        list(
          display = "ALL",
          named = "ALL=ALL"
        )
      )
    }

    variables <- names(conditions)

    variable_order <- order_configuration_variables(
      variables = variables,
      split_var = split_var
    )

    values <- vapply(
      variable_order,
      function(variable_name) {
        as.character(
          conditions[[variable_name]]
        )
      },
      character(1)
    )

    list(
      display = paste(
        values,
        collapse = "_"
      ),
      named = paste(
        paste0(
          variable_order,
          "=",
          values
        ),
        collapse = "|"
      )
    )
  }

  # ============================================================
  # 10. Raw cell statistics
  # ============================================================
  #
  # Computes observed means and classical standard errors within
  # each cell/configuration.
  # ============================================================

  compute_raw_cell_stats <- function(sub) {

    x <- sub[[response_column]]
    x <- x[!is.na(x)]
    n <- length(x)

    if (n == 0) {
      return(data.frame(
        n = 0,
        estimate = NA_real_,
        sd = NA_real_,
        se = NA_real_,
        ci_low = NA_real_,
        ci_high = NA_real_
      ))
    }

    m <- mean(x)
    sdv <- if (n > 1) stats::sd(x) else NA_real_
    se <- if (n > 1) sdv / sqrt(n) else NA_real_

    alpha <- 1 - conf_level
    crit <- if (n > 1) stats::qt(1 - alpha / 2, df = n - 1) else NA_real_

    # Intervals respect the natural range of the response; for
    # unbounded responses they are the usual t intervals.
    interval <- .observed_interval(
      estimate = m,
      se = se,
      n = n,
      range_type = response_range_type,
      conf_level = conf_level,
      critical = crit,
      bounds = bounds
    )

    data.frame(
      n = n,
      estimate = m,
      sd = sdv,
      se = se,
      ci_low = if (n > 1 || response_range_type == "proportion") interval$ci_low else NA_real_,
      ci_high = if (n > 1 || response_range_type == "proportion") interval$ci_high else NA_real_
    )
  }


  # ============================================================
  # 13. Unified cell statistics
  # ============================================================
  #
  # Raw estimates are observed means of the rows of each cell.
  #
  # Predicted estimates are computed on the reference grid of the
  # current coefficient (see .reference_grid()): numeric predictors
  # are fixed at their mean and the other factors are averaged with
  # equal weights, or held at the levels at which the coefficient is
  # defined when they interact with its variables. Each predicted
  # cell is represented by its weights over the grid rows, which are
  # propagated through the recursive decomposition so that every
  # aggregate is computed directly from the fixed-effect covariance
  # matrix.
  # ============================================================

  prediction_context <- NULL

  compute_cell_stats <- function(conditions) {

    sub <- subset_by_conditions(data, conditions)

    if (estimate_type == "observed") {

      statistics <- compute_raw_cell_stats(sub)

      # Raw estimates do not require model-based prediction weights.
      statistics$.prediction_weights <- list(NULL)

      return(statistics)
    }

    weights <- .grid_cell_weights(prediction_context$grid, conditions)

    prediction <- apply_bounds(
      .summarize_grid_predictions(
        predictions = prediction_context$predictions,
        weights = weights,
        prediction_scale = prediction_scale,
        conf_level = conf_level
      )
    )

    statistics <- data.frame(
      n = nrow(sub),
      estimate = prediction$estimate,
      sd = NA_real_,
      se = prediction$se,
      ci_low = prediction$ci_low,
      ci_high = prediction$ci_high,
      stringsAsFactors = FALSE
    )

    statistics$.prediction_weights <- list(weights)

    statistics
  }

  # ============================================================
  # 15. Weighted aggregation of raw estimates and SEs
  # ============================================================
  #
  # This helper is used only when estimate_type = "observed".
  #
  # Raw cell means are based on disjoint groups and their standard
  # errors are combined through classical variance propagation.
  #
  # Predicted estimates are not aggregated with this helper because
  # model-based predictions share the same fixed-effect coefficients
  # and are generally correlated.
  # ============================================================

  combine_weighted_mean_stats <- function(estimates, ses, ns, weights) {

    ok <- !is.na(estimates) & !is.na(weights)

    estimates <- estimates[ok]
    ses <- ses[ok]
    ns <- ns[ok]
    weights <- weights[ok]

    if (length(estimates) == 0 || sum(weights) == 0) {
      return(data.frame(
        estimate = NA_real_,
        se = NA_real_,
        ci_low = NA_real_,
        ci_high = NA_real_,
        df = NA_real_
      ))
    }

    # Normalize weights so that they sum to 1.
    a <- weights / sum(weights)

    estimate_agg <- sum(a * estimates)

    # Variance propagation for a weighted mean.
    ses2 <- ifelse(is.na(ses), 0, ses)^2
    var_components <- (a^2) * ses2
    se_agg <- sqrt(sum(var_components, na.rm = TRUE))

    # Approximate degrees of freedom for raw observed estimates.
    valid_df <- !is.na(ns) & ns > 1 & !is.na(var_components)

    if (sum(valid_df) == 0 || sum(var_components[valid_df]) == 0) {
      df_agg <- NA_real_
    } else {
      num <- (sum(var_components[valid_df]))^2
      den <- sum((var_components[valid_df]^2) / (ns[valid_df] - 1))
      df_agg <- if (den > 0) num / den else NA_real_
    }

    alpha <- 1 - conf_level

    crit <- if (!is.na(df_agg) && estimate_type == "observed") {
      stats::qt(1 - alpha / 2, df = df_agg)
    } else {
      stats::qnorm(1 - alpha / 2)
    }

    interval <- .observed_interval(
      estimate = estimate_agg,
      se = se_agg,
      n = sum(ns, na.rm = TRUE),
      range_type = response_range_type,
      conf_level = conf_level,
      critical = crit,
      single_cell = FALSE,
      bounds = bounds
    )

    data.frame(
      estimate = estimate_agg,
      se = se_agg,
      ci_low = interval$ci_low,
      ci_high = interval$ci_high,
      df = df_agg
    )
  }


  # ============================================================
  # 16. Exact aggregation of model-based predictions
  # ============================================================
  #
  # Predicted cell estimates are correlated because they are obtained
  # from the same fitted coefficient vector, so their standard errors
  # must not be combined as if they were independent. Each child is
  # represented by its weights over the reference grid; the weights
  # are combined with the higher-level contrast weights and the
  # aggregate is computed once from the complete fixed-effect
  # covariance matrix.
  # ============================================================

  combine_prediction_bundles <- function(prediction_weights,
                                         weights) {

    valid_bundle <- vapply(
      prediction_weights,
      function(current_weights) {
        is.numeric(current_weights) && sum(current_weights) > 0
      },
      logical(1)
    ) &
      !is.na(weights) &
      weights > 0

    if (!any(valid_bundle)) {

      return(
        list(
          statistics = data.frame(
            estimate = NA_real_,
            se = NA_real_,
            ci_low = NA_real_,
            ci_high = NA_real_,
            stringsAsFactors = FALSE
          ),
          prediction_weights = NULL
        )
      )
    }

    outer_weights <- weights[valid_bundle] / sum(weights[valid_bundle])

    combined_weights <- Reduce(
      `+`,
      Map(
        function(current_weights, outer_weight) {
          outer_weight * current_weights / sum(current_weights)
        },
        prediction_weights[valid_bundle],
        outer_weights
      )
    )

    statistics <- apply_bounds(
      .summarize_grid_predictions(
        predictions = prediction_context$predictions,
        weights = combined_weights,
        prediction_scale = prediction_scale,
        conf_level = conf_level
      )
    )

    list(
      statistics = statistics,
      prediction_weights = combined_weights
    )
  }

  # ============================================================
  # 16. Choose split variable
  # ============================================================
  #
  # This function decides which variable should be used as the
  # main splitting variable when decomposing an interaction.
  #
  # Available rules:
  #
  #   first_variable:
  #     use the first variable appearing in the interaction term.
  #
  #   last_variable:
  #     use the last variable appearing in the interaction term.
  #
  #   lowest_complexity:
  #     choose the variable with the simplest contrast coding.
  #
  #   highest_complexity:
  #     choose the variable with the most complex contrast coding.
  #
  # If several variables have the same complexity score, the
  # rightmost tied variable is selected.
  #
  # top_split overrides all automatic choices.
  # ============================================================

  choose_split_var <- function(component_df,
                               rule = c(
                                 "last_variable",
                                 "first_variable",
                                 "lowest_complexity",
                                 "highest_complexity"
                               ),
                               forced = NULL) {

    rule <- match.arg(rule)

    vars <- component_df$variable
    vars <- vars[!duplicated(vars)]

    if (!is.null(forced)) {
      if (!forced %in% vars) {
        stop(sprintf("The variable '%s' is not involved in the coefficient.", forced))
      }
      return(forced)
    }

    if (length(vars) == 1) return(vars)

    if (rule == "first_variable") {
      return(vars[1])
    }

    if (rule == "last_variable") {
      return(vars[length(vars)])
    }

    scores <- vapply(vars, contrast_complexity, numeric(1))

    if (rule == "lowest_complexity") {
      target_score <- min(scores, na.rm = TRUE)
    } else {
      target_score <- max(scores, na.rm = TRUE)
    }

    candidate_vars <- vars[scores == target_score]

    if (length(candidate_vars) > 1) {
      warning(
        paste0(
          "Multiple variables have the same complexity score (",
          paste(candidate_vars, collapse = ", "),
          "). The rightmost tied variable in the interaction is selected: ",
          candidate_vars[length(candidate_vars)],
          "."
        )
      )
    }

    candidate_vars[length(candidate_vars)]
  }

  # ============================================================
  # 17. Extract contrast information for one variable
  # ============================================================
  #
  # Returns the levels and weights on the positive and negative side
  # of one contrast column of a factor, taken from its hypothesis
  # matrix. Sum contrasts compare one level with the Grand Mean of all
  # levels; all other codings use the disjoint positive and negative
  # groups of the hypothesis row.
  # ============================================================

  divergence_notes <- character()

  get_contrast_info <- function(var, local_col) {

    variable_hypotheses <- get_variable_hypotheses(var)

    if (local_col > length(variable_hypotheses$contrasts)) {
      stop(sprintf("Contrast index %d is out of range for variable '%s'.", local_col, var))
    }

    contrast <- variable_hypotheses$contrasts[[local_col]]

    positive_weights <- contrast$positive_weights[contrast$positive_weights > 0]
    negative_weights <- contrast$negative_weights[contrast$negative_weights > 0]

    register_group_label(var, positive_weights)
    register_group_label(var, negative_weights)

    if (!is.na(contrast$note) &&
        !var %in% names(factor_checks)) {

      divergence_notes <<- unique(
        c(
          divergence_notes,
          sub("^Contrast_[0-9]+", paste0("Variable ", var, ", contrast ", local_col), contrast$note)
        )
      )
    }

    list(
      positive_levels = names(positive_weights),
      positive_weights = unname(positive_weights),
      negative_levels = names(negative_weights),
      negative_weights = unname(negative_weights),
      scale = contrast$scale
    )
  }

  # ============================================================
  # 18. Aggregate across actual levels
  # ============================================================
  #
  # This is used when a contrast group contains multiple levels.
  #
  # For each actual level, the function recursively builds the
  # inner groups for the remaining variables, then combines the
  # resulting estimates using the current contrast weights.
  # ============================================================

  aggregate_across_actual_levels <- function(actual_levels,
                                             actual_weights,
                                             component_df_remaining,
                                             current_var,
                                             conditions = list(),
                                             split_rule = "last_variable",
                                             top_split_var = NULL) {

    nested_list <- lapply(seq_along(actual_levels), function(i) {

      lev <- actual_levels[i]

      tmp <- build_inner_groups_for_one_level(
        component_df = component_df_remaining,
        conditions = c(conditions, setNames(list(lev), current_var)),
        split_rule = split_rule,
        top_split_var = top_split_var
      )

      tmp$actual_level <- lev
      tmp$actual_weight <- actual_weights[i]
      tmp
    })

    nested_all <- do.call(rbind, nested_list)
    keys <- unique(nested_all$inner_side)

    do.call(rbind, lapply(keys, function(k) {

      sub <- nested_all[nested_all$inner_side == k, , drop = FALSE]

      if (estimate_type == "observed") {

        agg <- combine_weighted_mean_stats(
          estimates = sub$estimate,
          ses = sub$se,
          ns = sub$n,
          weights = sub$actual_weight
        )

        aggregated_prediction_weights <- NULL

      } else {

        prediction_aggregation <- combine_prediction_bundles(
          prediction_weights = sub$.prediction_weights,
          weights = sub$actual_weight
        )

        agg <- prediction_aggregation$statistics


        aggregated_prediction_weights <-
          prediction_aggregation$prediction_weights
      }

      labels_source <- unique(unlist(strsplit(sub$configuration_source, "\\|\\|\\|")))
      named_source <- unique(unlist(strsplit(sub$configuration_named_source, "\\|\\|\\|")))
      var_cols <- make_variable_columns(named_source, factor_vars)

      aggregated_result <- cbind(
        data.frame(
          inner_side = k,
          inner_variable_path = sub$inner_variable_path[1L],
          configuration_source = collapse_source(labels_source),
          configuration_named_source = collapse_source(named_source),
          configuration_explicit = make_configuration_explicit(
            sub$configuration_explicit,
            sub$actual_weight
          ),
          configuration_semantic = make_configuration_semantic(
            named_source,
            split_var = top_split_var
          ),
          n = sum(sub$n, na.rm = TRUE),
          estimate = agg$estimate,
          sd = NA_real_,
          se = agg$se,
          ci_low = agg$ci_low,
          ci_high = agg$ci_high,
          stringsAsFactors = FALSE
        ),
        var_cols
      )


      aggregated_result$.prediction_weights <-
        list(aggregated_prediction_weights)

      aggregated_result
    }))
  }

  # ============================================================
  # 19. Build inner groups for one level
  # ============================================================
  #
  # This is the recursive core of the function.
  #
  # It takes a set of remaining contrast components and:
  #   1. chooses the next variable to decompose;
  #   2. identifies positive and negative levels for its contrast;
  #   3. either computes cell statistics directly, or recursively
  #      decomposes the remaining variables.
  #
  # This allows the function to work with two-way and higher-order
  # interactions.
  # ============================================================

  build_inner_groups_for_one_level <- function(component_df,
                                               conditions = list(),
                                               split_rule = "last_variable",
                                               top_split_var = NULL) {

    # Base case: no more variables to decompose.
    # Compute the estimate for the current cell/configuration.
    if (nrow(component_df) == 0) {

      stats_df <- compute_cell_stats(conditions)
      lab <- make_ordered_label_info(conditions, split_var = top_split_var)
      var_cols <- make_variable_columns(lab$named, factor_vars)

      cell_result <- cbind(
        data.frame(
          inner_side = "cell",
          inner_variable_path = "",
          configuration_source = lab$display,
          configuration_named_source = lab$named,
          configuration_explicit = lab$display,
          configuration_semantic = make_configuration_semantic(
            lab$named,
            split_var = top_split_var
          ),
          n = stats_df$n,
          estimate = stats_df$estimate,
          sd = stats_df$sd,
          se = stats_df$se,
          ci_low = stats_df$ci_low,
          ci_high = stats_df$ci_high,
          stringsAsFactors = FALSE
        ),
        var_cols
      )

      # Internal prediction information used by higher recursive levels.

      cell_result$.prediction_weights <-
        stats_df$.prediction_weights

      return(cell_result)
    }

    # Choose the next variable to decompose.
    current_var <- choose_split_var(component_df, rule = split_rule)
    row_idx <- which(component_df$variable == current_var)[1]
    local_col <- component_df$local_col[row_idx]

    contrast_info <- get_contrast_info(current_var, local_col)
    remaining <- component_df[component_df$variable != current_var, , drop = FALSE]

    # Identify positive and negative levels of the current contrast.
    pos_levels <- contrast_info$positive_levels
    pos_weights <- contrast_info$positive_weights
    neg_levels <- contrast_info$negative_levels
    neg_weights <- contrast_info$negative_weights

    # If no variables remain, compute the positive and negative
    # contrast groups directly.
    if (nrow(remaining) == 0) {

      pos_cells <- lapply(
        seq_along(pos_levels),
        function(i) {

          lev <- pos_levels[i]

          conds <- c(conditions,setNames(list(lev),current_var))

          stats_df <- compute_cell_stats(conds)

          lab <- make_ordered_label_info(conds,split_var = top_split_var)

          positive_cell <- data.frame(
            label = lab$display,
            named_label = lab$named,
            n = stats_df$n,
            estimate = stats_df$estimate,
            se = stats_df$se,
            weight = pos_weights[i],
            stringsAsFactors = FALSE
          )


          positive_cell$.prediction_weights <- stats_df$.prediction_weights

          positive_cell
        }
      )

      pos_cells <- do.call(rbind,pos_cells)

      neg_cells <- lapply(
        seq_along(neg_levels),
        function(i) {

          lev <- neg_levels[i]

          conds <- c(
            conditions,
            setNames(
              list(lev),
              current_var
            )
          )

          stats_df <- compute_cell_stats(conds)

          lab <- make_ordered_label_info(conds,split_var = top_split_var)

          negative_cell <- data.frame(
            label = lab$display,
            named_label = lab$named,
            n = stats_df$n,
            estimate = stats_df$estimate,
            se = stats_df$se,
            weight = neg_weights[i],
            stringsAsFactors = FALSE
          )


          negative_cell$.prediction_weights <- stats_df$.prediction_weights

          negative_cell
        }
      )

      neg_cells <- do.call(rbind,neg_cells)

      if (estimate_type == "observed") {

        pos_stats <- combine_weighted_mean_stats(
          estimates = pos_cells$estimate,
          ses = pos_cells$se,
          ns = pos_cells$n,
          weights = pos_cells$weight
        )

        neg_stats <- combine_weighted_mean_stats(
          estimates = neg_cells$estimate,
          ses = neg_cells$se,
          ns = neg_cells$n,
          weights = neg_cells$weight
        )

        positive_prediction_weights <- NULL

        negative_prediction_weights <- NULL

      } else {

        positive_aggregation <- combine_prediction_bundles(
          prediction_weights = pos_cells$.prediction_weights,
          weights = pos_cells$weight
        )

        negative_aggregation <- combine_prediction_bundles(
          prediction_weights = neg_cells$.prediction_weights,
          weights = neg_cells$weight
        )

        pos_stats <- positive_aggregation$statistics
        neg_stats <- negative_aggregation$statistics


        positive_prediction_weights <- positive_aggregation$prediction_weights


        negative_prediction_weights <- negative_aggregation$prediction_weights
      }

      pos_var_cols <- make_variable_columns(pos_cells$named_label, factor_vars)
      neg_var_cols <- make_variable_columns(neg_cells$named_label, factor_vars)

      positive_result <- cbind(
        data.frame(
          inner_side = "positive",
          inner_variable_path = current_var,
          configuration_source = collapse_source(
            pos_cells$label
          ),
          configuration_named_source = collapse_source(
            pos_cells$named_label
          ),
          configuration_explicit = make_configuration_explicit(
            pos_cells$label,
            pos_cells$weight
          ),
          configuration_semantic = make_configuration_semantic(
            pos_cells$named_label,
            split_var = top_split_var
          ),
          n = sum(pos_cells$n, na.rm = TRUE),
          estimate = pos_stats$estimate,
          sd = NA_real_,
          se = pos_stats$se,
          ci_low = pos_stats$ci_low,
          ci_high = pos_stats$ci_high,
          stringsAsFactors = FALSE
        ),
        pos_var_cols
      )


      positive_result$.prediction_weights <- list(positive_prediction_weights)


      negative_result <- cbind(
        data.frame(
          inner_side = "negative",
          inner_variable_path = current_var,
          configuration_source = collapse_source(
            neg_cells$label
          ),
          configuration_named_source = collapse_source(
            neg_cells$named_label
          ),
          configuration_explicit = make_configuration_explicit(
            neg_cells$label,
            neg_cells$weight
          ),
          configuration_semantic = make_configuration_semantic(
            neg_cells$named_label,
            split_var = top_split_var
          ),
          n = sum(neg_cells$n, na.rm = TRUE),
          estimate = neg_stats$estimate,
          sd = NA_real_,
          se = neg_stats$se,
          ci_low = neg_stats$ci_low,
          ci_high = neg_stats$ci_high,
          stringsAsFactors = FALSE
        ),
        neg_var_cols
      )


      negative_result$.prediction_weights <-
        list(negative_prediction_weights)

      return(
        rbind(
          positive_result,
          negative_result
        )
      )
    }

    # Recursive case: there are still variables to decompose.
    pos_df <- aggregate_across_actual_levels(
      actual_levels = pos_levels,
      actual_weights = pos_weights,
      component_df_remaining = remaining,
      current_var = current_var,
      conditions = conditions,
      split_rule = split_rule,
      top_split_var = top_split_var
    )

    neg_df <- aggregate_across_actual_levels(
      actual_levels = neg_levels,
      actual_weights = neg_weights,
      component_df_remaining = remaining,
      current_var = current_var,
      conditions = conditions,
      split_rule = split_rule,
      top_split_var = top_split_var
    )

    pos_df$inner_side <- paste0(
      "positive",
      ifelse(pos_df$inner_side == "cell", "", paste0("_", pos_df$inner_side))
    )

    neg_df$inner_side <- paste0(
      "negative",
      ifelse(neg_df$inner_side == "cell", "", paste0("_", neg_df$inner_side))
    )

    pos_df$inner_variable_path <- ifelse(
      is.na(pos_df$inner_variable_path) |
        pos_df$inner_variable_path == "",
      current_var,
      paste0(
        current_var,
        .moplot_variable_separator,
        pos_df$inner_variable_path
      )
    )

    neg_df$inner_variable_path <- ifelse(
      is.na(neg_df$inner_variable_path) |
        neg_df$inner_variable_path == "",
      current_var,
      paste0(
        current_var,
        .moplot_variable_separator,
        neg_df$inner_variable_path
      )
    )

    rbind(pos_df, neg_df)
  }

  # ============================================================
  # 20. Build the full bundle for one side of the top split
  # ============================================================
  #
  # For the top split variable, this function builds all inner
  # configurations for each actual level and then aggregates them
  # according to the contrast weights.
  # ============================================================

  build_group_bundle <- function(actual_levels,
                                 actual_weights,
                                 coefficient,
                                 split_var,
                                 component_df_remaining,
                                 split_rule = "last_variable") {

    nested_list <- lapply(seq_along(actual_levels), function(i) {

      lev <- actual_levels[i]

      tmp <- build_inner_groups_for_one_level(
        component_df = component_df_remaining,
        conditions = setNames(list(lev), split_var),
        split_rule = split_rule,
        top_split_var = split_var
      )

      tmp$actual_level <- lev
      tmp$actual_weight <- actual_weights[i]
      tmp
    })

    nested_all <- do.call(rbind, nested_list)
    keys <- unique(nested_all$inner_side)

    do.call(rbind, lapply(keys, function(k) {

      sub <- nested_all[nested_all$inner_side == k, , drop = FALSE]

      if (estimate_type == "observed") {

        agg <- combine_weighted_mean_stats(
          estimates = sub$estimate,
          ses = sub$se,
          ns = sub$n,
          weights = sub$actual_weight
        )

        aggregated_prediction_weights <- NULL

      } else {

        prediction_aggregation <- combine_prediction_bundles(
          prediction_weights = sub$.prediction_weights,
          weights = sub$actual_weight
        )

        agg <- prediction_aggregation$statistics


        aggregated_prediction_weights <-
          prediction_aggregation$prediction_weights
      }

      labels_source <- unique(unlist(strsplit(sub$configuration_source, "\\|\\|\\|")))
      named_source <- unique(unlist(strsplit(sub$configuration_named_source, "\\|\\|\\|")))
      var_cols <- make_variable_columns(named_source, factor_vars)

      bundle_result <- cbind(
        data.frame(
          coefficient = coefficient,
          split_var = split_var,
          inner_side = k,
          inner_variable_path = sub$inner_variable_path[1L],
          configuration_source = collapse_source(labels_source),
          configuration_named_source = collapse_source(named_source),
          configuration_explicit = make_configuration_explicit(
            sub$configuration_explicit,
            sub$actual_weight
          ),
          configuration_semantic = make_configuration_semantic(
            named_source,
            split_var = split_var
          ),
          n = sum(sub$n, na.rm = TRUE),
          estimate = agg$estimate,
          se = agg$se,
          ci_low = agg$ci_low,
          ci_high = agg$ci_high,
          stringsAsFactors = FALSE
        ),
        var_cols
      )


      bundle_result$.prediction_weights <-
        list(aggregated_prediction_weights)

      bundle_result
    }))
  }

  # ============================================================
  # 20.1 Reference grid, scale factors, and coefficient identity
  # ============================================================
  #
  # The reference grid depends on the variables of the coefficient:
  # these are the focal factors, and any other factor interacting
  # with them is held at the levels at which the coefficient is
  # defined. Grids and predictions are cached by variable set.
  #
  # Each coefficient must equal its scale factor multiplied by the
  # displayed difference of differences on the linear-predictor
  # scale. For predicted estimates this holds exactly unless a
  # variable of the coefficient interacts with a numeric predictor;
  # for raw estimates it holds only for a model containing the
  # interacting factors alone.
  # ============================================================

  fixed_effects <- .extract_fixef(model)

  family_object <- .extract_family(model)

  prediction_cache <- list()
  coefficient_scales <- list()
  coefficient_identity <- list()

  get_prediction_context <- function(variables) {

    key <- paste(sort(unique(variables)), collapse = "|")

    if (is.null(prediction_cache[[key]])) {

      grid_information <- .reference_grid(
        model = model,
        data = data,
        focal_variables = unique(variables)
      )

      prediction_cache[[key]] <<- list(
        grid = grid_information,
        predictions = .compute_fixed_predictions(
          model = model,
          newdata = grid_information$grid,
          prediction_scale = prediction_scale,
          conf_level = conf_level
        )
      )
    }

    prediction_cache[[key]]
  }

  check_coefficient_identity <- function(coef_name,
                                         component_df,
                                         scale_factor) {

    signed_weights <- lapply(
      seq_len(nrow(component_df)),
      function(index) {
        contrast <- get_variable_hypotheses(component_df$variable[index])$contrasts[[component_df$local_col[index]]]
        contrast$positive_weights - contrast$negative_weights
      }
    )

    names(signed_weights) <- component_df$variable

    cells <- expand.grid(
      lapply(signed_weights, names),
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )

    cell_weights <- apply(
      cells,
      1L,
      function(cell) {
        prod(
          vapply(
            names(signed_weights),
            function(var) signed_weights[[var]][[cell[[var]]]],
            numeric(1)
          )
        )
      }
    )

    active_cells <- which(abs(cell_weights) > 0)

    link_values <- vapply(
      active_cells,
      function(cell_index) {

        conditions <- as.list(cells[cell_index, , drop = FALSE])

        if (estimate_type == "predicted") {

          return(
            sum(
              .grid_cell_weights(prediction_context$grid, conditions) *
                attr(prediction_context$predictions, "eta")
            )
          )
        }

        cell_mean <- mean(subset_by_conditions(data, conditions)[[response_column]], na.rm = TRUE)

        if (.is_generalized_model(model)) {
          cell_mean <- suppressWarnings(family_object$linkfun(cell_mean))
        }

        cell_mean
      },
      numeric(1)
    )

    if (!all(is.finite(link_values))) {
      return(FALSE)
    }

    implied <- scale_factor * sum(cell_weights[active_cells] * link_values)

    coefficient <- as.numeric(fixed_effects[[coef_name]])

    abs(implied - coefficient) <= 1e-6 * max(1, abs(coefficient))
  }

  # ============================================================
  # 21. Decompose one interaction coefficient
  # ============================================================
  #
  # This function takes a single model coefficient, such as:
  #
  #   A1:B2-1
  #
  # and decomposes it into interpretable positive/negative level
  # configurations according to the contrast matrices of the
  # variables involved.
  # ============================================================

  decompose_one_coefficient <- function(coef_name,
                                        split_rule = "last_variable",
                                        top_split = NULL) {

    parts <- strsplit(coef_name, ":", fixed = TRUE)[[1]]

    component_df <- lookup[lookup$component %in% parts, , drop = FALSE]
    component_df <- component_df[match(parts, component_df$component), , drop = FALSE]

    if (nrow(component_df) != length(parts)) {
      warning(sprintf("It was not possible to fully map coefficient '%s'. Skipped.", coef_name))
      return(NULL)
    }

    # Weighted-group labels are specific to the contrast columns of
    # the current coefficient.
    group_label_registry <<- list()

    # Model-based predictions are computed once per set of variables
    # on the corresponding reference grid.
    if (estimate_type == "predicted") {

      prediction_context <<- get_prediction_context(component_df$variable)
    }

    split_var <- choose_split_var(component_df, rule = split_rule, forced = top_split)

    split_row <- which(component_df$variable == split_var)[1]
    split_col <- component_df$local_col[split_row]

    contrast_info <- get_contrast_info(split_var, split_col)
    remaining <- component_df[component_df$variable != split_var, , drop = FALSE]

    pos_levels <- contrast_info$positive_levels
    pos_weights <- contrast_info$positive_weights
    neg_levels <- contrast_info$negative_levels
    neg_weights <- contrast_info$negative_weights

    pos_bundle <- build_group_bundle(
      actual_levels = pos_levels,
      actual_weights = pos_weights,
      coefficient = coef_name,
      split_var = split_var,
      component_df_remaining = remaining,
      split_rule = split_rule
    )

    neg_bundle <- build_group_bundle(
      actual_levels = neg_levels,
      actual_weights = neg_weights,
      coefficient = coef_name,
      split_var = split_var,
      component_df_remaining = remaining,
      split_rule = split_rule
    )

    pos_bundle$split_side <- "positive"
    neg_bundle$split_side <- "negative"

    # The interaction hypothesis is the Kronecker product of the
    # hypothesis rows of its components, so its scale factor is the
    # product of their scale factors.
    scale_factor <- prod(
      vapply(
        seq_len(nrow(component_df)),
        function(index) {
          get_variable_hypotheses(component_df$variable[index])$contrasts[[component_df$local_col[index]]]$scale
        },
        numeric(1)
      )
    )

    coefficient_scales[[coef_name]] <<- scale_factor

    if (length(factor_checks) > 0L) {

      declared <- component_df$variable %in% names(factor_checks)

      component_rows <- lapply(
        which(declared),
        function(index) {
          factor_checks[[component_df$variable[index]]]$table[component_df$local_col[index], ]
        }
      )

      statuses <- vapply(component_rows, `[[`, character(1), "Status")
      ratios <- vapply(component_rows, `[[`, numeric(1), "Ratio")

      status <- if (length(statuses) == 0L) {
        "not declared"
      } else if (any(statuses == "differs")) {
        "differs"
      } else if (!all(declared)) {
        "partially declared"
      } else if (prod(ratios) < 0) {
        "reversed"
      } else if (abs(prod(ratios) - 1) < 1e-6) {
        "matches"
      } else {
        "rescaled"
      }

      component_labels <- ifelse(
        declared,
        vapply(
          seq_len(nrow(component_df)),
          function(index) {
            if (!declared[index]) return(NA_character_)
            factor_checks[[component_df$variable[index]]]$table$Hypothesis[component_df$local_col[index]]
          },
          character(1)
        ),
        paste0(component_df$variable, " (not declared)")
      )

      coefficient_checks[[coef_name]] <<- data.frame(
        Coefficient = coef_name,
        Hypothesis = paste(component_labels, collapse = " x "),
        Status = status,
        Ratio = if (status %in% c("differs", "not declared")) NA_real_ else prod(ratios),
        stringsAsFactors = FALSE
      )
    }

    coefficient_identity[[coef_name]] <<- check_coefficient_identity(
      coef_name = coef_name,
      component_df = component_df,
      scale_factor = scale_factor
    )

    pos_bundle$scale_factor <- scale_factor
    neg_bundle$scale_factor <- scale_factor

    rbind(pos_bundle, neg_bundle)
  }

  # ============================================================
  # 22. Apply decomposition to all selected coefficients
  # ============================================================

  out_list <- lapply(
    coefficients,
    decompose_one_coefficient,
    split_rule = split_rule,
    top_split = top_split
  )

  out_list <- out_list[!sapply(out_list, is.null)]

  if (length(out_list) == 0) {
    stop("No coefficient was successfully processed.")
  }

  out <- do.call(rbind, out_list)

  if (estimate_type == "predicted") {
    .warn_estimates_outside_bounds(out$estimate, bounds)
  }

  coefficient_scales <- unlist(coefficient_scales)
  coefficient_identity <- unlist(coefficient_identity)

  failed_identity <- names(coefficient_identity)[!coefficient_identity]

  if (estimate_type == "predicted" &&
      length(failed_identity) > 0L) {

    warning(
      paste0(
        "The displayed comparisons do not reconstruct the following ",
        "coefficients: ",
        paste(failed_identity, collapse = ", "),
        ". This usually means that a variable of the coefficient interacts ",
        "with a numeric predictor: the coefficient refers to the value zero ",
        "of that predictor, whereas the predictions fix it at its mean."
      ),
      call. = FALSE
    )
  }

  # Notes on codings whose signs do not describe the estimated
  # comparisons point to the check through `hypotheses`.
  if (length(divergence_notes) > 0L) {
    divergence_notes <- paste(
      divergence_notes,
      "To check the intended comparisons, supply them through `hypotheses`."
    )
  }

  hypothesis_messages <- unlist(
    lapply(factor_checks, `[[`, "messages"),
    use.names = FALSE
  )

  hypothesis_check <- if (length(coefficient_checks) > 0L) {
    out_check <- do.call(rbind, coefficient_checks)
    rownames(out_check) <- NULL
    out_check
  } else {
    NULL
  }

  for (note in c(hypothesis_messages, divergence_notes)) {
    message(note)
  }

  # Internal prediction bundles are needed only during recursive
  # aggregation and are not part of the user-facing output.
  internal_prediction_columns <- ".prediction_weights"

  out <- out[ ,setdiff(names(out),internal_prediction_columns),drop = FALSE]

  # ============================================================
  # 23. Add significance information
  # ============================================================

  contrast_lookup <- data.frame(
    coefficient = coefficients,
    contrast_id = seq_along(coefficients),
    stringsAsFactors = FALSE
  )

  contrast_lookup <- merge(
    contrast_lookup,
    model_info$sig_lookup,
    by = "coefficient",
    all.x = TRUE,
    sort = FALSE
  )

  out <- merge(
    out,
    contrast_lookup,
    by = "coefficient",
    all.x = TRUE,
    sort = FALSE
  )

  # ============================================================
  # 24. Preserve hierarchical decomposition
  # ============================================================

  if (!"inner_side" %in% names(out)) {
    stop(
      "Internal hierarchical decomposition information was unavailable.",
      call. = FALSE
    )
  }

  decomposition_columns <- c(
    "coefficient",
    "contrast_id",
    "significance",
    "p_value",
    "scale_factor",
    "split_var",
    "split_side",
    "inner_side",
    "inner_variable_path",
    factor_vars,
    "n",
    "estimate",
    "se",
    "ci_low",
    "ci_high",
    "configuration_explicit",
    "configuration_semantic"
  )

  decomposition_columns <- intersect(decomposition_columns,names(out))

  decomposition <- out[ ,decomposition_columns,drop = FALSE]

  decomposition$decomposition_path <- ifelse(
    is.na(decomposition$inner_side) |
      decomposition$inner_side == "" |
      decomposition$inner_side == "cell",
    decomposition$split_side,
    paste0(
      decomposition$split_side,
      "_",
      decomposition$inner_side
    )
  )

  decomposition$variable_path <- ifelse(
    is.na(decomposition$inner_variable_path) |
      decomposition$inner_variable_path == "",
    decomposition$split_var,
    paste0(
      decomposition$split_var,
      .moplot_variable_separator,
      decomposition$inner_variable_path
    )
  )

  decomposition$decomposition_depth <- lengths(
    strsplit(
      decomposition$decomposition_path,
      "_",
      fixed = TRUE
    )
  )

  rownames(decomposition) <- NULL

  # ============================================================
  # 24.1 Select and order final columns
  # ============================================================

  base_cols <- c(
    "coefficient", "contrast_id", "significance", "p_value",
    "scale_factor", "split_var", "split_side",
    factor_vars,
    "n", "estimate",
    "se", "ci_low", "ci_high"
  )

  configuration_cols <- switch(
    configuration,
    none = character(0),
    all = c("configuration_explicit", "configuration_semantic"),
    explicit = "configuration_explicit",
    semantic = "configuration_semantic"
  )

  out <- out[, c(base_cols, configuration_cols)]
  rownames(out) <- NULL

  # ============================================================
  # 25. ESTIMATE INFORMATION
  # ============================================================

  if (estimate_type == "observed") {

    estimate_info <- paste(
      "Displayed estimate:",
      "observed aggregated means."
    )

  } else {

    family_object <- .extract_family(model)

    if (prediction_scale == "response") {

      if (!is.null(family_object) &&
          identical(family_object$family, "binomial")) {

        estimate_info <- paste(
          "Displayed estimate:",
          "average model-predicted probability on the response scale."
        )

      } else if (
        !is.null(family_object) &&
        family_object$family %in% c(
          "poisson",
          "Negative Binomial(1)",
          "Negative Binomial(2)"
        )
      ) {

        estimate_info <- paste(
          "Displayed estimate:",
          "average model-predicted expected count on the response scale."
        )

      } else {

        estimate_info <- paste(
          "Displayed estimate:",
          "average model-predicted value on the response scale."
        )
      }

    } else {

      estimate_info <- paste(
        "Displayed estimate:",
        "average model-predicted value on the linear-predictor/link scale."
      )
    }
  }

  # ============================================================
  # 26. Store metadata as attributes
  # ============================================================
  #
  # These attributes are useful for reproducibility and downstream
  # plotting functions.
  # ============================================================

  # Notes on non-orthogonal custom codings are part of the note.
  bounds_note <- if (is.null(bounds)) {
    NULL
  } else {
    paste0(
      "Confidence intervals were computed within bounds = c(",
      bounds[1L], ", ", bounds[2L], ") on a transformed scale, so that they ",
      "remain within the range of the response."
    )
  }

  attr(out, "note") <- paste(c(disclaimer_note, bounds_note, hypothesis_messages, divergence_notes), collapse = " ")
  attr(out, "hypothesis_check") <- hypothesis_check
  attr(out, "bounds") <- bounds
  attr(out, "estimate_type") <- estimate_type
  attr(out, "prediction_scale") <- if (estimate_type == "observed") NA_character_ else prediction_scale
  attr(out, "estimate_info") <- estimate_info
  attr(out, "pvalue_source") <- model_info$pvalue_source
  attr(out, "split_rule") <- split_rule
  attr(out, "top_split") <- top_split
  attr(out, "interaction") <- interaction
  attr(out, "decomposition") <- decomposition
  attr(out, "scale_factors") <- coefficient_scales
  attr(out, "identity_holds") <- coefficient_identity

  class(out) <- c("plot_ready_contrasts_df", class(out))

  return(out)
}


# Separator between variable names in hierarchical variable paths, such
# as "B > A". Variable names may contain underscores, so "_" cannot be
# used: "drug_type_time_point" could not be split unambiguously.
.moplot_variable_separator <- " > "
