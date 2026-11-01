# Multiple imputation sensitivity analysis for Cox Models 4 and 5 (PREREGISTRATION.md, Section 6.3).
# mice 3.19.0: m = 20, maxit = 20, seed 42; predictive mean matching for KPS, logistic regression for binary
# MGMT methylation status; imputation model: all model covariates, GPR153, MGMT expression, batch, the event
# indicator and the Nelson-Aalen cumulative hazard. Cox fits pooled by Rubin's rules.
# Usage: Rscript pipeline/a23_mice.R <mi_input.csv> <output.csv>
suppressPackageStartupMessages({library(mice); library(survival)})
a <- commandArgs(trailingOnly = TRUE)
d <- read.csv(a[1])
d$meth_status <- factor(d$meth_status)
meth <- make.method(d)
meth[] <- ""
meth["kps"] <- "pmm"
meth["meth_status"] <- "logreg"
pred <- make.predictorMatrix(d)
pred[, "time"] <- 0
imp <- mice(d, m = 20, maxit = 20, method = meth, predictorMatrix = pred, seed = 42, printFlag = FALSE)
out <- list()
fits <- list(
  "Model 4" = with(imp, coxph(Surv(time, event) ~ gpr153 + age + male + kps + meth_status, ties = "efron")),
  "Model 5" = with(imp, coxph(Surv(time, event) ~ gpr153 + age + male + kps + meth_status + mgmt, ties = "efron")))
for (mod in names(fits)) {
  s <- summary(pool(fits[[mod]]), conf.int = TRUE)
  out[[mod]] <- data.frame(model = mod, term = as.character(s$term), coef = s$estimate, se = s$std.error,
                           hr = exp(s$estimate), hr_low = exp(s[["2.5 %"]]), hr_high = exp(s[["97.5 %"]]),
                           p = s$p.value, df = s$df, m = 20, n = nrow(d))
}
write.csv(do.call(rbind, out), a[2], row.names = FALSE)
cat("mice: imputed KPS", sum(is.na(d$kps)), "and methylation status", sum(is.na(d$meth_status)), "of", nrow(d), "\n")
