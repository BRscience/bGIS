util_cal_pinteractions <- function(
    data,
    survival_col,
    treatment_col,
    biomarker_col,
    ref = NULL,
    if_nacse = TRUE
) {
  if (!is.null(ref)) {
    data <-
      data %>%
      mutate(
        {{biomarker_col}} :=
          relevel(.data[[biomarker_col]], ref)
      )
  }
  if (if_nacse) {
    tb_ncase <- util_ncase_statistics(data, biomarker_col)
    freq <- tb_ncase[[1]] %>% pull(pct) %>% tail(1)
  } else {
    freq <- NA
  }
  formula <-
    as.formula(
      glue::glue("survival::Surv({survival_col},{survival_col}.status) ~
                 {treatment_col} + {biomarker_col} +
                 {treatment_col}:{biomarker_col}"))
  cox <- (survival::coxph(formula, data = data))
  p <- summary(cox)$coefficients[, 'Pr(>|z|)']
  HR <- summary(cox)$coefficients[, 'exp(coef)']
  z <- summary(cox)$coefficients[, 'z']
  df_res <- data.frame(
    treatment_col = treatment_col,
    biomarker_col = biomarker_col,
    freq = freq,
    survival = survival_col,
    HR = round(as.numeric(HR[3]), 3),
    subgroup_p = round(as.numeric(p[1]), 3),
    test_p = round(as.numeric(p[2]), 3),
    z_score = round(as.numeric(z[3]), 3),
    interaction_p = round(as.numeric(p[3]), 4)
  )
  lt_res <- list(
    cox_model = cox,
    res_table = df_res
  )
  return(lt_res)
}
