#' Blood Based ITH calculation method in batch
#'
#' reference: http://10.10.11.111:10010/display/DS/bITH
#'
#' The method is based on ctDNA deep sequenced data.
#' Only SNV data is needed.
#'
#' Data requirement:
#'
#' sample_id : each sample should have a unique sample_id
#'
#' af : the original variant allele frequency.
#'
#' @param report report The report
#' @param maxaf_adjust default TRUE, for scale of variant allele frequency.
#' @param low_mutation_count_filter whether to filter low mutation counts data.
#' @param given_weights given weight vector, length should be the same as bins.
#' @param bin_width default 0.1, define the width of Shannon Index bins.
#' @param weighted_shannon_entropy default TRUE,whether to use weighted Shannon. Index
#' @param mutation_number_cutoff interger, numbers of mutations to filter
#'  samples.
#'
#' @export
#'
#' @examples
#'
#' util_bITH_calculation(
#'   report = data.frame(
#'     sample_id = rep("id1",5),
#'     mutation_type = rep("missense", 5),
#'     af_dbl = c(0.11,0.25,0.005,0.0075,0.15)
#'   )
#' )
#'
util_bITH_calculation <- function(report,
                                  maxaf_adjust = TRUE,
                                  low_mutation_count_filter = FALSE,
                                  mutation_number_cutoff = 5,
                                  weighted_shannon_entropy = TRUE,
                                  given_weights = NULL,
                                  bin_width = 0.1){
  report <-
    report %>%
    util_mark_specific_col_types(type_to_mark = "cn_|fusion|large")
  report <- util_trans_af(report)

  if(low_mutation_count_filter){
    report_count <-
      report %>%
      gsa(sample_id) %>%
      filter(n >= mutation_number_cutoff)
    report <-
      report %>%
      filter(sample_id %in% report_count$sample_id)
  }
  result <-
    purrr::map_dfr(
      .x = report %>%
        group_by(sample_id) %>%
        group_split(.keep = TRUE),
      .f = ~ cal_Shannon_Index(
        .x,
        maxaf_adjust = maxaf_adjust,
        weighted_shannon_entropy = weighted_shannon_entropy,
        given_weights = given_weights,
        bin_width = bin_width))
  return(result)
}

#' Blood Based ITH calculation method by sample
#'
#' This is only for single sample data, and is more customized.
#'
#'
#' @param tmp_report Reports only include one sample.
#' @param bin_width default 0.1, The width of Shannon Index bins.
#' @param start_pos default 0, The start position of  all bins.
#' @param end_pos default 1, The end position of  all bins.
#' @param maxaf_adjust default FALSE, for scale of variant allele frequency.
#' @param adjust_alpha default 1, rescale of bins, default no process.
#' @param weighted_shannon_entropy logiical, add weighted function or not.
#' @param given_weights vectors, length equals to number of bins.
cal_Shannon_Index <- function(tmp_report,
                              bin_width = 0.1,
                              start_pos = 0,
                              end_pos = 1,
                              maxaf_adjust = FALSE,
                              weighted_shannon_entropy = TRUE,
                              given_weights = NULL,
                              adjust_alpha = 1){
  if(weighted_shannon_entropy){
    if(!is.null(given_weights)){
      weights_vec <- given_weights
    }else{
      weights_vec <-
        seq(
          from=start_pos+bin_width/2,
          to=end_pos,
          by=bin_width
        )
    }
  }else{
    weights_vec <- rep(1,round(1/bin_width))
  }
  tmp_report <-
    tmp_report %>%
    dplyr::mutate(
      af_backup = af
    )
  max_af <- max(tmp_report$af,na.rm = TRUE)
  # min_af <- min(tmp_report$af,na.rm = TRUE)
  # af_dev <- max_af - min_af
  # mean_af <- mean(tmp_report$af,na.rm = TRUE)
  # fix_weight <- max_af/mean_af
  if(maxaf_adjust){
    tmp_report <-
      tmp_report %>%
      dplyr::mutate(
        af = af/max_af
      )
  }
  shannon_distribution <-
    data.frame(
      start = seq(from = start_pos,to = end_pos - bin_width,by=bin_width),
      end = seq(from = start_pos + bin_width,to = end_pos ,by=bin_width)
    )
  shannon_distribution <-
    shannon_distribution %>%
    dplyr::group_by(start,end) %>%
    dplyr::summarise(
      counts =
        tmp_report %>%
        filter(af > start,af <= end) %>%
        nrow(),
      .groups = "drop") %>%
    dplyr::mutate(pi = counts/nrow(tmp_report)) %>%
    dplyr::mutate(weight = weights_vec) %>%
    dplyr::mutate(
      pilnpi =
        case_when(
          pi != 0 ~ pi*log(pi,base=exp(1)),
          pi == 0 ~ 0
        )
    ) %>%
    dplyr::mutate(
      fipilnpi =
        pilnpi*weight
    )
  Hindex <- -sum(shannon_distribution$fipilnpi, na.rm = TRUE)
  shannon_index_dataframe <-
    data.frame(sample_id = unique(tmp_report$sample_id), bITH = Hindex)
  return(shannon_index_dataframe)
}

util_trans_af <- function(report) {
  if (!"af_dbl" %in% colnames(report)) {
    report <-
      report %>%
      mutate(
        af = case_when(
          stringr::str_detect(af, "^[\\d+\\.]+$") ~
            suppressWarnings(as.double(af)),
          stringr::str_detect(af, "^[\\d+\\.]+\\%+$") ~
            str_replace(af, "\\%+$", "") %>%
            {suppressWarnings(as.double(.)) / 100},
          TRUE ~ rlang::na_dbl
        )
      )
  } else {
    report <- report %>%
      mutate(af = af_dbl)
  }
  return(report)
}
