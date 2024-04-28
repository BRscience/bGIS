# bGIS

> bGIS subtyping strategy aims to distinguish patients benefit from the addition of PD-1 inhibitor to first-line chemotherapy. Use codes in this project to calculate bTMB, bITH and bCIN scores, and generate corresponding bGIS subtypes. And then calculate prediction efficacy with cox interaction analysis. Finally, visualize interaction effects with Kaplan-Meier Curves.

## 1. calculate bTMB

> bTMB was defined as numbers of non-intron/non-synonymous/non-hotspot SNVs/Indels divided by captured panel size (MB) .

```
library(dplyr)
library(stringr)
tb_bTMB <- muts_report %>% 
  filter(!str_detect(mutation_type, "synonymous|intron")) %>% 
  group_by(sample_id) %>% 
  summarise(
    bTMB_count = n(),
    bTMB = bTMB_count/1.003
  )
```

## 2. calculate bITH 

> bITH was defined as weighted Shannon index of SNVs/Indels allele frequency [1].

```
source("./util_mark_specific_col_types.R")
source("./util_bITH_calculation.R")
tb_bITH <- util_bITH_calculation(muts_report)
```

## 3. calculate bCIN

> bCIN was calculated with wgii index based on segments profiles [2].

segments format example, generated from ichorCNA [3]:
```
ID      chrom   start   end     num.mark        seg.median.logR copy.number     call    subclone.status logR_Copy_Number        Corrected_Copy_Number   Corrected_Call
s1   1       1000001 20000001       247     -0.107397104233856      2       NEUT    FALSE   0.967729472249722       2       NEUT
s1   1       20000001 50000001       189     0.188470794751652       4       AMP     FALSE   5.27198468677867        4       AMP
```

```
source("./util_cal_bCIN.R")
tb_bCIN <- util_cal_bCIN(tb_seg)
```

## 4. bGIS 

```
source("./util_cal_bGIS.R")
clinical_table <- util_cal_bGIS(clinical_table)
```

## 5. interaction effects

```
source("./util_ncase_statistics.R")
source("./util_cal_pinteractions.R")
util_cal_pinteractions(
  data = clinical_table,
  survival_col = "OS",
  treatment_col = "treatment",
  biomarker_col = "bGIS"
)
```

## 6. Visulization
```
source("plotutil_survival_standard.R")
source("util_ncase_statistics.R")
library(dplyr)
library(ggplot2)

x <- plotutil_survival_standard(
  df_clicor = clinical_table %>%
    mutate(
      group_comb = case_when(
        group == "bGIS-1" & treatment == "comb" ~ "T-bGIS-1",
        group == "bGIS-1" & treatment == "chemo" ~ "C-bGIS-1",
        group == "bGIS-2" & treatment == "comb" ~ "T-bGIS-2",
        group == "bGIS-2" & treatment == "chemo" ~ "C-bGIS-2",
        group == "bGIS-3" & treatment == "comb" ~ "T-bGIS-3",
        group == "bGIS-3" & treatment == "chemo" ~ "C-bGIS-3",
        TRUE ~ NA_character_
      ),
      group_comb = factor(
        group_comb,
        levels = c(
          "T-bGIS-1",
          "C-bGIS-1",
          "T-bGIS-2",
          "C-bGIS-2",
          "T-bGIS-3",
          "C-bGIS-3"
        )
      )
    ),
  surfit_col = "group",
  survival_col = "PFS",
  face_type = "plain",
  surv_scale = "percent",
  legend_pos = c(0.75, 0.8),
  remove_sep_HR_factorP = TRUE,
  remove_overall_P = TRUE,
  fontsize_of_legend = 8,
  linetype = c("solid", "twodash", "solid", "twodash", "solid", "twodash"),
  ncase_in_legend = FALSE,
  color_pal = c("#374E55FF", "#878787", "#BC3C29FF", "#fc9272", "#0072B5FF", "#9ecae1")
)

x$plot <- x$plot +
  guides(
    color = guide_legend(
      title = "chemo+anti-PD1\n       chemo"
    ),
    linetype = guide_legend(
      title = "chemo+anti-PD1\n       chemo"
    )
  ) +
  theme(
    legend.position = "top",
    legend.box="vertical",
    legend.text = element_text(family = "sans",size = 10,lineheight = 1.5),
    legend.title = element_text(family = "sans",size = 10,lineheight = 1.5)
  ) +
  ggplot2::scale_y_continuous(labels = scales::percent_format(suffix = "")) +
  xlab("Time since treatment start (months)")
```

### Reference
[1] Fan, Y., Liu, Y., Wang, L., Cai, Y., Cao, W., Sun, W., Zou, X., Li, B., Zhang, Z., Cai, S., et al. (2023). bITH, a blood-based metric of intratumor heterogeneity, is associated with clinical response to immune checkpoint blockade in non-small cell lung cancer. EBioMedicine 91, 104564. 10.1016/j.ebiom.2023.104564.

[2] Burrell, R.A., McClelland, S.E., Endesfelder, D., Groth, P., Weller, M.C., Shaikh, N., Domingo, E., Kanu, N., Dewhurst, S.M., Gronroos, E., et al. (2013). Replication stress links structural and numerical cancer chromosomal instability. Nature 494, 492-496. 10.1038/nature11935

[3] Adalsteinsson, V.A., Ha, G., Freeman, S.S., Choudhury, A.D., Stover, D.G., Parsons, H.A., Gydush, G., Reed, S.C., Rotem, D., Rhoades, J., et al. (2017). Scalable whole-exome sequencing of cell-free DNA reveals high concordance with metastatic tumors. Nat Commun 8, 1324. 10.1038/s41467-017-00965-y

