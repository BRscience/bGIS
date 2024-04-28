util_cal_bCIN <- function(tb_seg) {
  tb_seg <-
    tb_seg %>%
    mutate(length = end - start + 1) %>%
    group_by(sample, chr) %>%
    mutate(
      chr_length = sum(length),
      instable_length = case_when(
        event == "NEUT" ~ 0,
        TRUE ~ length
      ),
      chr_cin_length = sum(instable_length),
      chr_wgii = chr_cin_length/chr_length
    ) %>%
    ungroup()
  tb_bCIN <- tb_seg %>%
    group_by(sample, chr) %>%
    slice(1L) %>%
    ungroup() %>%
    group_by(sample) %>%
    summarise(
      wgii = mean(chr_wgii, na.rm = TRUE),
      .groups = "drop"
    )
  return(tb_bCIN)
}
