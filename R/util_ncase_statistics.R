util_ncase_statistics <- function(data, category_variates) {
  ncase <-
    purrr::map(
      .x = category_variates,
      .f = ~ util_ncase_statistics_basic(data, variate = .x)
    )
  names(ncase) <- category_variates
  return(ncase)
}

util_ncase_statistics_basic <- function(data, variate) {
  n_names <- data[[variate]] %>% unique() %>% length()
  if (n_names > 10L) {
    msg <-
      glue::glue("{variate}'s category class should be no more than 10!") %>%
      c(x = .) %>%
      rlang::format_error_bullets()
    rlang::abort(msg)
  }
  summary_res <-
    data %>%
    group_by(.data[[variate]]) %>%
    summarise(
      size = n()
    ) %>%
    transmute(
      group_name = .data[[variate]],
      size,
      label = glue::glue("{group_name} (n={size})") %>% as.character(),
      pct = size/sum(size)
    ) %>%
    ungroup()
  return(summary_res)
}
