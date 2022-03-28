# Neocortical Localization and Thalamocortical Modulation of 
# Neuronal Hyperexcitability contribute to Fragile X Syndrome
# Communication Biology 
# Demographics
# Table 1, Supplemental Table 1, Figure 2A
# Author: E. Pedapati
# Version: 3/27/2022

pacman::p_load(tidyverse, labelled, compareGroups, lsmeans, BiocManager, filesstrings)
`%notin%` <- Negate(`%in%`)

source.demo <- "https://figshare.com/ndownloader/files/34517351"

demographics <- c("eegid", "group", "sex","mosaic","visitage")

selectvars <- c("iq_dev",
                "sbs_nvz", "sbs_vz", "adams_anxiety",
                "adams_ocd", "scq_total", "abc_irritable",
                "abc_hyperactivity", "abc_speech",
                "abc_lethargy", "abc_stereotypy", "wj3", "visitage"
)

colors_subgroup <- c(
  "FXS(F)" = "darkorange2",
  "FXS(M)" = "red",
  "Control(F)" = "darkviolet",
  "Control(M)" = "blue"
)

# exclude subjects with benzodiazapines 794, 1806, 2464
model.demo  <- 
  read_csv(source.demo) %>%
  dplyr::filter(eegid %notin% c("D0794_rest","D1806_rest","D2464_rest")) %>%
  dplyr::select(all_of(c(demographics, selectvars))) %>%
  mutate( NVIQ = (sbs_nvz * 15) + 100,
          VIQ  = (sbs_vz * 15) + 100) %>% 
  mutate(sex = factor(sex, levels=c("Male","Female"), labels=c("M","F")),
         subgroup = paste0(group,'_',sex),
         mgroup = paste0(subgroup,'_',mosaic)) %>%
  mutate(subgroup = factor(subgroup,
                           levels=c("TDC_M","FXS_M","TDC_F","FXS_F"),
                           labels=c("Control(M)","FXS(M)","Control(F)","FXS(F)")),
         mgroup=factor(mgroup, 
                       levels=c("TDC_M_4","FXS_F_3", "TDC_F_4", "FXS_M_1", "FXS_M_2"), 
                       labels=c("CM","FF","CF","FM","MM")),
         sex=factor(sex,
                    levels=c("M","F"), 
                    labels=c("Male","Female"))) 

model.demo %>% write_csv(file = "commbio_clinicalmeasures.csv")

# Model for short demographic table.
dat.demo_short <- model.demo %>%
  select(group, 
         visitage,
         iq_dev,
         VIQ,
         NVIQ,
         scq_total,
         adams_anxiety, 
         wj3) %>%
  set_variable_labels(visitage = "Age (y)", 
                      iq_dev = "FSIQ",
                      VIQ = "VIQ",
                      NVIQ = "NVIQ",
                      scq_total = "Social Score",
                      adams_anxiety = "Anxiety Score",
                      wj3="WJ-III")

# Model for long demographic table.
dat.demo_long <- model.demo %>%
  select(subgroup,
         visitage,
         iq_dev,
         VIQ,
         NVIQ,
         scq_total,
         abc_irritable,
         abc_hyperactivity, 
         abc_speech,
         abc_lethargy,
         abc_stereotypy,
         adams_ocd,
         adams_anxiety, 
         wj3) %>%
  set_variable_labels(visitage = "Age (Years)", subgroup="Group", iq_dev = "FSIQ",
                      scq_total = "SCQ", wj3 = "WJ-III", abc_irritable = "ABC-Irritability",
                      abc_lethargy = "ABC-Lethargy", abc_stereotypy = "ABC-Stereotypy",
                      abc_hyperactivity = "ABC-Hyperactivity", abc_speech = "ABC-Abnormal Speech",
                      adams_anxiety = "ADAMS-Anxiety", adams_ocd = "ADAMS-OCD")


# Model for NVIQ figure
dat.demo_iq <- model.demo %>% select(eegid, NVIQ, sbs_nvz, group, sex, subgroup) 

# Supplemental Table 1
# Caption: Short demographic table (Table)
result_table.demo_short <- descrTable(group ~ ., dat.demo_short, 
                                      hide.no = "no", 
                                      simplify = TRUE)

filename.supplTable1 = "commbio_supplTable1.docx"
export2word(result_table.demo_short, file = filename.supplTable1)
file.move(files = filename.supplTable1, "output/", overwrite = TRUE)
file.remove(str_replace(filename.supplTable1, "docx", "Rmd"))

# Table 1
# Caption: Long demographic table (Table)
result_table.demo_long <- descrTable(subgroup ~ ., dat.demo_long, 
                                     hide.no = "no", 
                                     simplify = TRUE)

filename.Table1 = "commbio_Table1.docx"
export2word(result_table.demo_long, file = filename.Table1)
file.move(files = filename.Table1, "output/", overwrite = TRUE)
file.remove(str_replace(filename.Table1, "docx", "Rmd"))


# Figure 2a
fit.nviq <- aov(NVIQ ~ group * sex, data=dat.demo_iq)
dat.table_nviq_estimates <- lsmeans::lsmeans(fit.nviq, ~ group:sex) %>% as_tibble()%>% 
  dplyr::select(-df) %>% mutate(across(where(is.numeric), round,2))

# Figure: Boxplot of NVIQ
p.demo_iq <- dat.demo_iq %>% ggplot(aes(x = subgroup, y = sbs_nvz, fill = subgroup)) +
  geom_boxplot(outlier.size = 0, alpha = .2) +
  scale_y_continuous(breaks = seq(-10, 2, 2)) +
  geom_point(aes(fill = subgroup, shape = sex),
             color = "black", size = 8,
             position = position_jitterdodge(jitter.width = 0.75), alpha = 0.8
  ) +
  scale_shape_manual(values = c(24, 21)) +
  scale_fill_manual(values = colors_subgroup) +
  scale_color_manual(values = colors_subgroup) +
  theme_bw(base_size = 30) +
  theme(
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black"), legend.position = "none"
  ) +
  xlab("Group") +
  ylab("Non-verbal IQ z-score")

write_csv(x = dat.demo_iq,file = "figshare/commbio_Figure2a_SourceData.csv")
ggsave("output/Figure2a.pdf",colormodel = "cmyk")


