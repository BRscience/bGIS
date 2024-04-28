#' Standard Survival KM curve plot.
#'
#' Standard Survival KM curve plot.
#'
#' @param df_clicor tibble, clinical tables in clicor module
#' @param surfit_col character, colname of test variate like "TMB"
#' @param survival_col  character, colname of survival variate like "OS".
#'  BTW, if you give "OS", thene "OS.status" should also present in df_clicor.
#' @param time_unit character, like "days", to place in plots
#' @param remove_sep_HR_factorP logical, when surfit_col have more than 2 types,
#'  set to TURE too get comparison results for each two variates.
#' @param color_pal character vector, color hex in order.
#' @param legend_pos numeric vector of length 2, x pos and y pos of legend.
#' @param fontsize_of_legend numeric, fontsize_of_legend.
#' @param fontsize_of_axis numeric, fontsize_of_axis.
#' @param fontsize_of_labs numeric, fontsize_of_labs.
#' @param fontsize_of_annotation numeric, fontsize_of_annotation.
#' @param remove_overall_HR logical，to remove HR from pictures.
#' @param remove_overall_P logical，to remove P value from pictures.
#' @param fontfamily fontfamily for pictures
#' @param remove_legend_title logical，to remove legend_title from pictures.
#' @param anno_pos numeric vector of length 2, x pos and y pos of annotation.
#' @param ncase_in_legend logical，to annotate ncase in pictures or not.
#' @param fontsize_of_tabletitle numeric, fontsize_of_risk_table_title.
#' @param face_type character, default
#' @param tables.theme survminer theme, table themes for risk table.
#' @param surv.median.line character, survival median time lines, inherited
#'  from ggsurvplot.
#' @param linetype linetype vector or starta name inherited from ggsurvplot.
#' @param legend_labs inherited from legend.labs in ggsurvplot.
#' @param break_time breaks of time
#'
#' @return ggsurvplot object
#' @export
#'
#' @examples #TODO
plotutil_survival_standard <- function(
    df_clicor,
    surfit_col,
    survival_col,
    time_unit = "months",
    remove_overall_P = FALSE,
    remove_sep_HR_factorP = FALSE,
    remove_overall_HR = FALSE,
    remove_legend_title = FALSE,
    color_pal = NULL,
    legend_pos = c(0.75, 0.8),
    fontsize_of_legend = 14,
    fontsize_of_axis = 12,
    fontsize_of_labs = 14,
    fontsize_of_annotation = 5,
    fontsize_of_tabletitle = 13,
    anno_pos = NULL,
    ncase_in_legend = TRUE,
    fontfamily = "sans",
    face_type = c("bold", "plain", "italic")[1],
    tables.theme =
      survminer::theme_cleantable(),
    surv.median.line = c("none", "hv"),
    linetype = NULL,
    legend_labs = NULL,
    break_time = NULL,
    ...
) {
  surv.median.line <- match.arg(surv.median.line)
  # data filtering ----
  survival_col_status <- paste0(survival_col, ".status")

  df_clicor <-
    df_clicor %>%
    filter(
      !is.na(.data[[surfit_col]]),
      !is.na(.data[[survival_col]]),
      !is.na(.data[[survival_col_status]])
    ) %>%
    mutate(tmp_surfit := .data[[surfit_col]])
  # break time calculation ----
  if (is.null(break_time)) {
    break_time <-
      df_clicor %>%
      pull(.data[[survival_col]]) %>%
      max() %>%
      `/`(5) %>%
      ceiling()
  }
  # Labs modification ----
  time_lab_str <-
    paste0("Time in ", time_unit)
  if (survival_col == "PFS") {
    surv_lab_str <- "Progression-free survival ratio"
  } else if (survival_col == "OS") {
    surv_lab_str <- "Overall survival ratio"
  } else {
    surv_lab_str <- as.character(survival_col)
  }
  # Do Survfit & coxph test to prepare for plot objects ----
  formula <- as.formula(
    glue::glue(
      "survival::Surv({survival_col}, {survival_col_status}) ~ tmp_surfit"
    )
  )

  fit <-
    survminer::surv_fit(formula = formula, data = df_clicor)
  # names(fit$strata) <-
  #   stringr::str_replace(
  #     pattern = ".*=",
  #     replacement = "",
  #     string =  names(fit$strata)
  #   )
  res_cox <- survival::coxph(formula = formula, data = df_clicor)
  coxph_summary <- summary(res_cox)
  res_survdiff <- survival::survdiff(formula = formula, data = df_clicor)

  # plot ----
  ncase_summary <- util_ncase_statistics(df_clicor, surfit_col)
  n_of_surfit_col <- ncase_summary[[1]] %>% nrow()

  if (n_of_surfit_col > 2) {
    annotate_sep_HR_factorP <- TRUE
    annotate_overall_HR <- FALSE
  } else {
    annotate_overall_HR <- TRUE
    annotate_sep_HR_factorP <- FALSE
  }
  if (remove_overall_P) {
    annotate_overall_P <- FALSE
  } else {
    annotate_overall_P <- TRUE
  }
  if (remove_sep_HR_factorP) {
    annotate_sep_HR_factorP <- FALSE
  }
  if (remove_overall_HR) {
    annotate_overall_HR <- FALSE
  }

  pic <- plotutil_survival_with_survfit_obj(
    data = df_clicor,
    survfit_obj = fit,
    surfit_col = surfit_col,
    survival_col = survival_col,
    coxph_summary_obj = coxph_summary,
    survdiff_obj = res_survdiff,
    annotate_overall_P = annotate_overall_P,
    annotate_overall_HR = annotate_overall_HR,
    annotate_sep_HR_factorP = annotate_sep_HR_factorP,
    remove_legend_title = remove_legend_title,
    color_pal = color_pal,
    tables.theme = tables.theme,
    time_lab_str = time_lab_str,
    break_time = break_time,
    legend_pos = legend_pos,
    fontsize_of_legend = fontsize_of_legend,
    fontsize_of_axis = fontsize_of_axis,
    fontsize_of_labs = fontsize_of_labs,
    fontsize_of_annotation = fontsize_of_annotation,
    fontsize_of_tabletitle = fontsize_of_tabletitle,
    anno_pos = anno_pos,
    ncase_in_legend = ncase_in_legend,
    fontfamily = fontfamily,
    surv.median.line = surv.median.line,
    face_type = face_type,
    linetype = linetype,
    legend_labs = legend_labs,
    ...
  )

  return(pic)
}

plotutil_survival_with_survfit_obj <- function(
  data,
  surfit_col,
  survival_col,
  survfit_obj,
  coxph_summary_obj,
  survdiff_obj,
  annotate_overall_P = TRUE,
  annotate_overall_HR = TRUE,
  annotate_sep_HR_factorP = FALSE,
  ncase_in_legend = TRUE,
  break_time = 5,
  censor.size = 4.5,
  color_pal = NULL,
  legend_labs = NULL,
  time_lab_str = "Time in days",
  surv_lab_str = "Progression-free Survival (%)",
  risk.table = TRUE,
  risk.table.height = 0.25,
  tables.theme = survminer::theme_survminer(),
  legend_pos = c(0.75, 0.8),
  remove_legend_title = FALSE,
  fontsize_of_legend = 14,
  fontsize_of_axis = 12,
  fontsize_of_labs = 14,
  fontsize_of_annotation = 5,
  fontsize_of_tabletitle = 13,
  anno_pos = NULL,
  fontfamily = "sans",
  surv.median.line = "none",
  surv.scale = c("percent", "default")[1],
  face_type = c("bold", "plain", "italic")[1],
  linetype = linetype,
  ...
) {

  # check parameters ----

  if (annotate_overall_HR && annotate_sep_HR_factorP) {
    my_abort("Can't annotate overall HR and seperate HR in one pic. >_<|||")
  }

  # get color pal ----
  if (is.null(color_pal)) {
    # color_pal <- c("#e41a1c", "#377eb8", "#4daf4a",
    #                "#984ea3", "#ff7f00", "#ffff33")
    color_pal_default <- ggsci::pal_jama("default")(7)
  } else {
    color_pal_default <- color_pal
  }

  # show ncase

  if (is.null(legend_labs)) {
    ncase_summary <- util_ncase_statistics(data, surfit_col)
    nlabel_surfit_col <- ncase_summary[[surfit_col]] %>% pull(label)
    names(nlabel_surfit_col) <- ncase_summary[[surfit_col]] %>% pull(group_name)
    n_group_name_surfit_col <- ncase_summary[[surfit_col]] %>% pull(group_name)
    names(n_group_name_surfit_col) <- ncase_summary[[surfit_col]] %>% pull(group_name)
    if (ncase_in_legend) {
      legend_labs <- nlabel_surfit_col
    } else {
      legend_labs <- n_group_name_surfit_col
    }
  }

  # preprocess surv_lab_str ----

  if (survival_col == "PFS") {
    surv_lab_str <- "Progression-free Survival (%)"
  } else if (survival_col == "OS") {
    surv_lab_str <- "Overall Survival (%)"
  } else if (survival_col == "DFS") {
    surv_lab_str <- "Disease-free Survival (%)"
  } else {
    surv_lab_str <- survival_col
  }


  # generate ggsurv obj ----
  ggsurv <- survminer::ggsurvplot(
    fit = survfit_obj,
    data = data,
    risk.table = risk.table,
    risk.table.height = risk.table.height,
    risk.table.y.text = FALSE,
    surv.scale = surv.scale,
    surv.median.line = surv.median.line,
    conf.int = FALSE,
    censor.size = censor.size,
    tables.theme = tables.theme,
    fontsize = fontsize_of_axis/3,
    break.time.by = break_time,
    palette = color_pal,
    ggtheme = ggpubr::theme_pubr(),
    font.family = fontfamily,
    legend.labs = legend_labs,
    linetype = linetype,
    ...
  )

  # modify ggsurv plot ----

  if (!(!is.null(linetype) && !is.null(color_pal))) {
    ggsurv$plot <- ggsurv$plot +
      ggplot2::guides(
        color = ggplot2::guide_legend(
          title = ifelse(remove_legend_title, "", surfit_col),
          title.theme = ggplot2::element_text(
            size = fontsize_of_legend + 1,
            face = face_type,
            colour = "black"
          )
        )
      )
  }
  ggsurv$plot <- ggsurv$plot +
    #ggplot2::scale_y_continuous(labels = scales::percent_format(suffix = "")) +
    ggplot2::theme(
      legend.position = legend_pos,
      legend.margin = ggplot2::margin(),
      legend.text = ggplot2::element_text(size = fontsize_of_legend),
      axis.text = ggplot2::element_text(size = fontsize_of_axis),
      axis.title.x = ggplot2::element_text(
        face = face_type, color = "black", size = fontsize_of_labs
      ),
      axis.title.y = ggplot2::element_text(
        face = face_type, color = "black", size = fontsize_of_labs
      ),
    ) +
    ggplot2::labs(x = time_lab_str, y = surv_lab_str)

  # add annotation of p_value and HR ----
 if (is.null(anno_pos)) {
   ybase_for_overall_P <-
     ifelse(
       annotate_sep_HR_factorP,
       0.15,
       ifelse(annotate_overall_HR, 0.075, 0)
     )
   xbase <- 0.05
 } else {
   ybase_for_overall_P <- anno_pos[2]
   xbase <- anno_pos[1]
 }
  #TODO-solved: change oberall_p from sctest in coxph to survdiff logRank P
  #overall_P <- coxph_summary_obj$sctest["pvalue"]
  if (annotate_overall_P) {
    overall_P <- 1 - pchisq(survdiff_obj$chisq, length(survdiff_obj$n) - 1)
    ggsurv$plot <- ggsurv$plot +
      ggplot2::theme() +
      ggplot2::annotate(
        geom = "text",
        x =  xbase,
        hjust = 0,
        y = ybase_for_overall_P,
        size = fontsize_of_annotation,
        fontface = "italic",
        label = ifelse(
          overall_P <= 0.001,
          "p < 0.001",
          stringr::str_c("p = ", sprintf("%.3f", overall_P))
        )
      )
  }

  if (annotate_overall_HR) {
    overall_HR <- paste0(
      sprintf("%.2f", coxph_summary_obj$conf.int[, 1]), " [",
      sprintf("%.2f", coxph_summary_obj$conf.int[, 3]), "-",
      sprintf("%.2f", coxph_summary_obj$conf.int[, 4]), "]"
    )
    ggsurv$plot <- ggsurv$plot +
      ggplot2::annotate(
        geom = "text",
        x = xbase,
        y = ybase_for_overall_P - 0.075,
        hjust = 0,
        size = fontsize_of_annotation,
        fontface = "italic",
        label = paste0("HR = ", overall_HR)
      )
  }
  if (annotate_sep_HR_factorP) {
    sep_HR <- paste(
      sprintf("%.2f", coxph_summary_obj$conf.int[, 1]),
      collapse = ";"
    )
    sep_P <- paste(
      sprintf("%.2f", coxph_summary_obj$coefficients[, 5]),
      collapse = ";"
    )
    ggsurv$plot <- ggsurv$plot +
      ggplot2::annotate(
        geom = "text",
        x = xbase,
        y = ybase_for_overall_P - 0.075,
        hjust = 0,
        size = fontsize_of_annotation,
        fontface = "italic",
        label = paste0("HR = ", sep_HR)
      )
    ggsurv$plot <- ggsurv$plot +
      ggplot2::annotate(
        geom = "text",
        x = xbase,
        y = ybase_for_overall_P - 0.15,
        hjust = 0,
        size = fontsize_of_annotation,
        fontface = "italic",
        label = paste0("Factor P-value = ", sep_P)
      )
  }

  # modify ggsurv table ----
  if (risk.table) {
    ggsurv$table <-
      ggsurv$table +
      ggplot2::labs(
        x = ggplot2::element_blank(),
        y = ggplot2::element_blank()
      ) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = fontsize_of_labs)
      )
    ggsurv$table <-
      ggsurv$table +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = fontsize_of_tabletitle)
      )
  }

  return(ggsurv)
}
