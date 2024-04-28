#' Mark cols with specific types
#'
#' Mark cols with specific types
#'
#' @param report tibble, tables of whatever u gave.
#' @param col_to_mark character, colname that should be marked.
#' @param type_to_mark character, types that should be marked.
#'
#' @return tibble, marked tables.
#' @export
#'
#' @examples #TODO
util_mark_specific_col_types <- function(report,
                                         col_to_mark = "mutation_type",
                                         type_to_mark = "intron|synonymous"){

  report$dummy_order <- rownames(report)

  sample_to_mark_negative <-
    report %>%
    group_by(sample_id) %>%
    summarise(
      is_negative = all(stringr::str_detect(
        .data[[col_to_mark]], stringr::regex(type_to_mark, ignore_case = TRUE)
      )),
      .groups = "drop"
    ) %>%
    filter(is_negative) %>%
    pull(sample_id)

  if (length(sample_to_mark_negative) == 0) {
    report_to_filter_mark_content <-
      report %>%
      filter(!stringr::str_detect(
        .data[[col_to_mark]], stringr::regex(type_to_mark, ignore_case = TRUE)
      )) %>%
      select(-dummy_order)
    return(report_to_filter_mark_content)
  }

  report_to_mark_negative <-
    report %>%
    filter(sample_id %in% sample_to_mark_negative) %>%
    mutate(
      gene = "negative",
      mutation_type = "negative",
      description = "negative"
    ) %>%
    group_by(sample_id) %>%
    slice(1L) %>%
    ungroup() %>%
    select(sample_id, gene, mutation_type, description, dummy_order)

  report_to_filter_mark_content <-
    report %>%
    filter(!stringr::str_detect(
      .data[[col_to_mark]], stringr::regex(type_to_mark, ignore_case = TRUE)
    ))

  report_after_mark <-
    dplyr::bind_rows(report_to_filter_mark_content, report_to_mark_negative)

  report_after_mark <-
    report_after_mark %>%
    arrange(dummy_order) %>%
    select(-dummy_order)

  return(report_after_mark)
}
