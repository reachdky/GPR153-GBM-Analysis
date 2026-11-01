# Cox fits (survival::coxph, Efron ties) and proportional-hazards tests (survival::cox.zph, transform "km",
# per term and global) for a list of jobs.
# Usage: Rscript pipeline/a20_coxzph.R <jobs.csv>
# jobs.csv columns: input (CSV with columns time, event and the covariates), prefix, covariates
# (comma-separated), strata (empty or a column name). Writes <prefix>_coxph.csv and <prefix>_zph.csv per job.
suppressPackageStartupMessages(library(survival))
jobs <- read.csv(commandArgs(trailingOnly = TRUE)[1], stringsAsFactors = FALSE, na.strings = character(0))
for (i in seq_len(nrow(jobs))) {
  d <- read.csv(jobs$input[i])
  rhs <- paste(strsplit(jobs$covariates[i], ",")[[1]], collapse = " + ")
  if (!is.na(jobs$strata[i]) && nzchar(jobs$strata[i])) rhs <- paste0(rhs, " + strata(", jobs$strata[i], ")")
  fit <- coxph(as.formula(paste("Surv(time, event) ~", rhs)), data = d, ties = "efron",
               control = coxph.control(eps = 1e-12, iter.max = 100))
  z <- cox.zph(fit, transform = "km", terms = TRUE)
  write.csv(data.frame(term = rownames(z$table), chisq = z$table[, "chisq"], df = z$table[, "df"],
                       p = z$table[, "p"]), paste0(jobs$prefix[i], "_zph.csv"), row.names = FALSE)
  s <- summary(fit)
  ci <- s$conf.int
  write.csv(data.frame(term = rownames(s$coefficients), coef = s$coefficients[, "coef"], hr = ci[, "exp(coef)"],
                       hr_low = ci[, "lower .95"], hr_high = ci[, "upper .95"], se = s$coefficients[, "se(coef)"],
                       p = s$coefficients[, "Pr(>|z|)"], n = fit$n, events = fit$nevent, loglik = fit$loglik[2]),
            paste0(jobs$prefix[i], "_coxph.csv"), row.names = FALSE)
}
