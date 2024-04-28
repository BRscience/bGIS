util_cal_bGIS <- function(
    clinical_table,
    MSAF_cutoff = 0.005,
    bITH_cutoff = 0.65,
    bTMB_cutoff = 13,
    bCIN_cutoff = 0.3836
) {
  clinical_table <- clinical_table %>%
    mutate(
      BEP = case_when(
        MSAF >= MSAF_cutoff ~ "ctDNA+",
        TRUE ~ "ctDNA-"
      ),
      count_BEP = BEP == "ctDNA+",
      count_bITH = case_when(
        BEP == "ctDNA-" ~ FALSE,
        bITH <= bITH_cutoff ~ TRUE,
        TRUE ~ FALSE
      ),
      count_bTMB = case_when(
        BEP == "ctDNA-" ~ FALSE,
        bTMB >= bTMB_cutoff ~ TRUE,
        TRUE ~ FALSE
      ),
      count_CIN = case_when(
        BEP == "ctDNA-" ~ FALSE,
        bCIN <= bCIN_cutoff ~ TRUE,
        TRUE ~ FALSE
      )
    ) %>% mutate(
      pro_ici_flag = count_bITH + count_bTMB + count_CIN,
      group = case_when(
        BEP == "ctDNA-" ~ "bGIS-1",
        pro_ici_flag = 0 ~ "bGIS-3",
        pro_ici_flag > 0 ~ "bGIS-2",
        TRUE ~ NA_character_
      )
    )
  return(clinical_table)
}
