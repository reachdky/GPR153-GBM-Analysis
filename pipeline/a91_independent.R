# Independent re-computation in R for a91_independent_check.py (PREREGISTRATION.md, Section 13.2).
# Usage: Rscript pipeline/a91_independent.R <workdir>
# Reads  <workdir>/m3.csv    time, event, gpr153, age, male          -> m3_out.csv    survival::coxph (Efron)
#        <workdir>/km.csv    time, event, batch, <gene>_log2tpm       -> km_logrank_out.csv, km_groups_out.csv,
#                                                                         km_at_risk_out.csv   survival::survdiff, survfit
#        <workdir>/meta.csv  analysis, yi, vi                         -> meta_out.csv  metafor::rma (DL), test "knha",
#                                                                         confint (Q-profile), predict (pi.type "Riley")
suppressPackageStartupMessages({library(survival); library(metafor)})
wd <- commandArgs(trailingOnly = TRUE)[1]
f <- function(x) file.path(wd, x)
DPM <- 365.25 / 12

d <- read.csv(f("m3.csv"))
m <- coxph(Surv(time, event) ~ gpr153 + age + male, data = d, ties = "efron",
           control = coxph.control(eps = 1e-12, iter.max = 100))
s <- summary(m)
write.csv(data.frame(term = rownames(s$coefficients), hr = s$conf.int[, "exp(coef)"], hr_low = s$conf.int[, "lower .95"],
                     hr_high = s$conf.int[, "upper .95"], p = s$coefficients[, "Pr(>|z|)"], n = m$n, events = m$nevent),
          f("m3_out.csv"), row.names = FALSE)

k <- read.csv(f("km.csv"))
months <- as.numeric(strsplit(readLines(f("km_months.txt")), ",")[[1]])
lr <- list(); gr <- list(); ar <- list()
for (gene in c("gpr153", "mgmt")) {
  x <- k[[paste0(gene, "_log2tpm")]]
  splits <- list("pooled median" = x > median(x), "within-batch median" = x > ave(x, k$batch, FUN = median))
  for (sp in names(splits)) {
    k$grp <- ifelse(splits[[sp]], "high", "low")
    sd <- survdiff(Surv(time, event) ~ grp, data = k)
    rev_all <- summary(survfit(Surv(time, 1 - event) ~ 1, data = k))$table["median"]
    lr[[length(lr) + 1]] <- data.frame(gene = toupper(gene), split = sp, n = sum(sd$n), chisq = sd$chisq,
                                       p = pchisq(sd$chisq, 1, lower.tail = FALSE), median_followup_days_reverse_km = rev_all)
    for (g in c("high", "low")) {
      kk <- k[k$grp == g, ]
      fit <- survfit(Surv(time, event) ~ 1, data = kk, conf.type = "log-log")
      tb <- summary(fit)$table
      rv <- summary(survfit(Surv(time, 1 - event) ~ 1, data = kk))$table["median"]
      gr[[length(gr) + 1]] <- data.frame(gene = toupper(gene), split = sp, group = g, n = nrow(kk), deaths = sum(kk$event),
                                         median_os_days = tb["median"], median_os_ci_low = tb["0.95LCL"],
                                         median_os_ci_high = tb["0.95UCL"], median_followup_days_reverse_km = rv)
      ar[[length(ar) + 1]] <- data.frame(gene = toupper(gene), split = sp, group = g, month = months,
                                         at_risk = summary(fit, times = months * DPM, extend = TRUE)$n.risk)
    }
  }
}
write.csv(do.call(rbind, lr), f("km_logrank_out.csv"), row.names = FALSE)
write.csv(do.call(rbind, gr), f("km_groups_out.csv"), row.names = FALSE)
write.csv(do.call(rbind, ar), f("km_at_risk_out.csv"), row.names = FALSE)

mt <- read.csv(f("meta.csv"), stringsAsFactors = FALSE)
out <- list()
for (a in unique(mt$analysis)) {
  e <- mt[mt$analysis == a, ]
  r <- rma(e$yi, e$vi, method = "DL")
  hk <- rma(e$yi, e$vi, method = "DL", test = "knha")
  ci <- confint(r, control = list(tol = 1e-12))$random
  pr <- if (r$k > 2) predict(r, pi.type = "Riley") else list(pi.lb = NA, pi.ub = NA)
  out[[length(out) + 1]] <- data.frame(analysis = a, k = r$k, mu = r$b[1], se = r$se, ci_low = r$ci.lb, ci_high = r$ci.ub,
                                       p = r$pval, Q = r$QE, Q_p = r$QEp, tau2 = r$tau2, I2 = r$I2 / 100,
                                       hk_ci_low = hk$ci.lb, hk_ci_high = hk$ci.ub, hk_p = hk$pval,
                                       tau2_ci_low = ci["tau^2", "ci.lb"], tau2_ci_high = ci["tau^2", "ci.ub"],
                                       pi_low = pr$pi.lb, pi_high = pr$pi.ub)
}
write.csv(do.call(rbind, out), f("meta_out.csv"), row.names = FALSE)
