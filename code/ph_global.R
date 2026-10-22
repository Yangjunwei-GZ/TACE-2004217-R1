args <- commandArgs(trailingOnly=TRUE)
input <- args[1]
out <- args[2]
suppressPackageStartupMessages(library(survival))
d <- read.csv(input, check.names=FALSE)
d$log_selected_volume <- log1p(d$selected_volume_cm3)
specs <- list(
 clinical=c('age','sex_female','child_b','log_afp','log_tumour_count','bclc_b'),
 primary_combined=c('age','sex_female','child_b','log_afp','log_tumour_count','bclc_b','selected_diameter_mm','mean_hu','radiomics_score')
)
rows <- list()
for (name in names(specs)) {
 vars <- specs[[name]]
 formula <- as.formula(paste('Surv(os_months, death_event) ~', paste(vars, collapse=' + ')))
 fit <- coxph(formula, data=d, ties='breslow', x=TRUE, y=TRUE)
 for (tr in c('rank','km')) {
  ph <- cox.zph(fit, transform=tr, global=TRUE)
  tab <- as.data.frame(ph$table)
  tab$variable <- rownames(tab)
  tab$model <- name
  tab$time_transform <- tr
  rows[[length(rows)+1]] <- tab
 }
}
result <- do.call(rbind, rows)
write.csv(result, file=out, row.names=FALSE)
age_signal <- any(result$variable=='age' & result$p<0.05)
age_path <- file.path(dirname(out), 'AGE_TIME_DEPENDENT_SENSITIVITY.csv')
if (age_signal) {
 d$radiomics_score_SD <- as.numeric(scale(d$radiomics_score, center=TRUE, scale=sqrt(mean((d$radiomics_score-mean(d$radiomics_score))^2))))
 fit_time <- coxph(Surv(os_months, death_event) ~ age + gender_woman + child_b + log_afp + log_tumour_count + bclc_b + selected_diameter_mm + mean_hu + radiomics_score_SD + tt(age), data=d, ties='breslow', tt=function(x,t,...) x*log(pmax(t,1e-8)))
 tab <- as.data.frame(summary(fit_time)$coefficients)
 ci <- as.data.frame(summary(fit_time)$conf.int)
 tab$variable <- rownames(tab)
 tab$CI_lower <- ci$`lower .95`
 tab$CI_upper <- ci$`upper .95`
 tab$directly_comparable_C_index_reported <- 'No'
 write.csv(tab, age_path, row.names=FALSE)
} else {
 write.csv(data.frame(status='No age PH signal at the frozen exploratory threshold; time-effect analysis not triggered',directly_comparable_C_index_reported='No'), age_path, row.names=FALSE)
}
cat('Global and per-variable PH diagnostics completed; exploratory unpenalised diagnostic fits only.\n')
