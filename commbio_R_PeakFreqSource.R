# Peak Alpha Frequency
# Author: E. Pedapati
# Version: 3/27/2022
source("_Common.R")

pacman::p_load(
  R.matlab,
  tidyverse,
  here,
  weights,
  flextable,
  reshape2,
  rstatix,
  hdf5r,
  DescTools,
  nlme,
  emmeans,
  ggthemes,
  broom.mixed,
  dotwhisker,
  ggsci,
  ggsignif,
  units,
  ggsegDefaultExtra
)

pacman::p_load_current_gh("LCBC-UiO/ggsegDefaultExtra", "ggseg/ggseg")


# LOAD OTHER RDATA ============================================================#
# Group Assigments
source.demo <-
  "https://figshare.com/ndownloader/files/34517351"

source.group <-
  "https://figshare.com/ndownloader/files/34517366"
clininfo <- read_csv(source.group) %>%
  rename(subgroup2 = subgroup) %>%
  rename(sex2 = sex) %>%
  mutate(sex = str_extract(sex2, "^.{1}"),
         subgroup = paste0(group, "_", sex)) %>%
  relocate(eegid, group, sex, subgroup)

source.atlas <-
  "https://figshare.com/ndownloader/files/34517354"
nodeinfo <- read_csv(source.atlas)


# CUSTOM OUTPUTS   ============================================================#
target_file = "output/commbio_PAF.RData"
target_fig = "figshare/commbio_PAF.RData"
target_file.between_table <-
  str_replace(target_file, ".RData", "_between_table.docx")
target_file.within_table <-
  str_replace(target_file, ".RData", "_within_table.docx")
target_file.atlas_figure <-
  str_replace(target_fig, ".RData", "_atlas_figure.pdf")
target_file.intraregion_figure <-
  str_replace(target_file, ".RData", "_intraregion_figure.pdf")
target_file.group_table <-
  str_replace(target_file, ".RData", "_group_table.docx")
target_file.group_figure <-
  str_replace(target_fig, ".RData", "_group_figure.pdf")
target_file.nodetable <-
  str_replace(target_file, ".RData", "_node_table.docx")
target_file.node_figure <-
  str_replace(target_fig, ".RData", "_node_figure.pdf")
target_file.PafCorrelations <-
  str_replace(target_file, ".RData", "forCorr.RData")

source.mne.peakfreq <-
  "https://figshare.com/ndownloader/files/34517756"
import.peakfreq.source <- read_csv(source.mne.peakfreq)

# source level peak frequency
df.source <- import.peakfreq.source %>%
  left_join(nodeinfo, by = c("Electrode" = "labelclean")) %>%
  left_join(clininfo, by = c("eegid" = "eegid")) %>%
  rename(value = peakFreq) %>%
  drop_na(value) %>%
  mutate(
    side = ifelse(
      nchar(region) < 3,
      substring(region, 1, nchar(region) - 1),
      substring(region, 1, nchar(region) - 2)
    ),
    zone = substring(region, 2, nchar(region))
  ) %>%
  mutate(
    zone = factor(
      zone,
      levels = c("PF", "F", "L", "T", "P", "C", "O"),
      labels = c(
        "Prefrontal",
        "Frontal",
        "Limbic",
        "Temporal",
        "Parietal",
        "Central",
        "Occipital"
      )
    ),
    side = factor(
      side,
      levels = c("L", "R"),
      labels = c("Left", "Right")
    )
  ) %>%
  filter(zone %in% c("Occipital", "Parietal", "Central", "Frontal"))

unique(df.source$zone)

#==============================================================================#
# ANALYSIS #1      ============================================================#
#                  Between group analysis
#==============================================================================#

#==============================================================================#
# LME: Group, sex, and zone
#==============================================================================#

# Base LME model includes group, sex, lower band, and resting state network (RSN)
fit.1 <-
  nlme::lme(
    value ~ group * sex * zone,
    random = ~ 1 | eegid,
    method = "ML",
    correlation = corCompSymm(form =  ~ 1 |
                                eegid),
    data = df.source
  )
anova(fit.1)  # laterality did not matter

# numDF denDF  F-value p-value
# (Intercept)        1  4452 9761.976  <.0001
# group              1   140   22.539  <.0001
# sex                1   140    0.017  0.8967
# zone               3  4452   57.585  <.0001
# group:sex          1   140    0.018  0.8949
# group:zone         3  4452    7.246  0.0001
# sex:zone           3  4452    2.537  0.0549
# group:sex:zone     3  4452    3.481  0.0152

# selected model
fit.final <- fit.1
summary(fit.1)
anova.lme(fit.1)

# least-squared means
emc <- emmeans(fit.final, ~ group * sex * zone)
estimates <-
  broom::tidy(emc, conf.int = TRUE, conf.level = .95) %>%
  mutate(zone = factor(zone, levels = c(
    "Occipital", "Parietal", "Central", "Frontal"
  )))

estimates_sub <-
  broom::tidy(emc, conf.int = TRUE, conf.level = .95) %>%
  mutate(zone = factor(zone, levels = c(
    "Occipital", "Parietal", "Central", "Frontal"
  )))

df.sublevel <-
  df.source %>% dplyr::select(eegid, value, Electrode, group, sex, zone) %>%
  filter(zone %in%  c("Occipital", "Parietal", "Central", "Frontal")) %>%
  group_by(eegid, zone, group, sex) %>% dplyr::summarize(value = mean(value)) %>%
  ungroup() %>%
  mutate(subgroup = paste0(group, "_", sex),
         sex = factor(
           sex,
           levels = c("F", "M"),
           labels = c("Female", "Male")
         ))



# get paired significant contrasts
pairs_by_group <-
  broom::tidy(pairs(emc, by = c("zone", "sex"), adjust = 'none')) %>%
  group_by(contrast) %>% mutate(adj.p = p.adjust(p.value, method = "fdr", n = n())) %>%
  select(-term, -null.value) %>% ungroup() %>% select(-p.value) %>%
  mutate(
    estimate = weights::rd(estimate, 2),
    std.error = weights::rd(std.error, 2),
    statistic = weights::rd(statistic, 1),
    adj.p = scales::pvalue(adj.p)
  )

# get paired significant contrasts - across zones
pairs_by_group.zone <-
  broom::tidy(pairs(emc, by = c("group", "sex"), adjust = 'none')) %>%
  group_by(contrast) %>% mutate(adj.p = p.adjust(p.value, method = "fdr", n = n())) %>%
  select(-term, -null.value) %>% ungroup() %>% select(-p.value)

# get mean values for summary table
peakmeans <- estimates %>%
  mutate(mean = paste0(rd(estimate, 1), " \u00B1 ", rd(std.error, 2)),
         ploty = estimate) %>%
  select(group, sex, zone, mean, ploty) %>%
  pivot_wider(names_from = group, values_from = c(mean, ploty)) %>%
  left_join(pairs_by_group, by = c("zone", "sex")) %>%
  relocate(c(starts_with("FXS"), starts_with("TDC")), .before = estimate)

# Create significance stars for between group plots
peakmean.plot <-
  peakmeans %>% mutate(ploty = (ploty_FXS + ploty_TDC) / 2,
                       stars = add_sig_stars(adj.p)) %>%
  select(sex, zone, statistic, ploty,  stars, adj.p) %>%
  mutate(stars = ifelse(stars == "", "", stars))

#==============================================================================#
# LME 1: Group Table
#==============================================================================#
peakmeans <- estimates %>%
  mutate(mean = paste0(rd(estimate, 1), "\u00B1", rd(std.error, 2))) %>%
  select(group, sex, zone, mean) %>%
  pivot_wider(names_from = group, values_from = c(mean)) %>%
  left_join(pairs_by_group, by = c("zone", "sex")) %>% select(-contrast) %>%
  relocate(c(starts_with("FXS"), starts_with("TDC")), .before = estimate) %>%
  relocate(zone, .before = sex) %>%
  flextable() %>%
  set_header_labels(
    label = "Node",
    region = "Region",
    TDC = "Control",
    sex = "Sex",
    zone = "Region",
    estimate = "FXS-Control",
    std.error = "SE",
    df = "DF",
    statistic = "F",
    adj.p = "5% FDR"
  ) %>%
  flextable::merge_v(j = 1) %>%
  flextable::valign(j = 1, valign = "top") %>%
  autofit() %>% fix_border_issues() %>%
  save_as_docx(path = target_file.group_table)

out("Table", target_file.group_table)

peakmeans

#==============================================================================#
# LME 1: Plot
#==============================================================================#

# New facet label names for dose variable
sex.labs <- c(M = "Male", F = "Female")
subgroup.labs <- c(FXS_F = "FXS(F)", FXS_M = "FXS(M)")

#==============================================================================#
# FIGURE: WITHIN SUBJECT - BY REGION
#==============================================================================#

labels_zone = c("Occipital", "Parietal", "Central", "Frontal")

anno_line <- pairs_by_group.zone %>% rowwise() %>%
  mutate(
    subgroup = paste0(group, "_", sex),
    region1 = str_split(contrast, " - ", simplify = TRUE)[1],
    region2 = str_split(contrast, " - ", simplify = TRUE)[2],
    region1 = factor(region1, levels = labels_zone),
    region2 = factor(region2, levels = labels_zone),
    offset = ifelse(str_detect(group, "FXS"), -1, 1),
    ystart = ifelse(str_detect(contrast, "Frontal - Occipital"), .25, 0),
    ystart = ifelse(str_detect(contrast, "Central - Occipital"), .5, ystart),
    ystart = ifelse(str_detect(contrast, "Parietal - Occipital"), .75, ystart),
    ystart = ifelse(str_detect(contrast, "Parietal - Central"), 1, ystart),
    ystart = ifelse(str_detect(contrast, "Frontal - Parietal"), 1.25, ystart),
    y = ifelse(str_detect(group, "FXS"), 7.75 - ystart, 9.75 + ystart),
    vpos = ifelse(str_detect(group, "FXS"), -1.4, .4),
    label = add_sig_stars(adj.p),
    rand = rnorm(1)
  ) %>%
  relocate(region1, region2, .before = estimate) %>%
  select(subgroup,
         vpos,
         sex,
         rand,
         region1,
         statistic,
         region2,
         y,
         adj.p,
         label) %>%
  ungroup() %>%
  mutate(
    start = as.numeric(region1),
    end = as.numeric(region2),
    group = as.numeric(factor(subgroup))
  ) %>%
  filter(label != "")
anno_line

p.pafTopo <- estimates %>%
  mutate(subgroup = paste0(group, "_", sex)) %>%
  ggplot(aes(group = subgroup)) +
  geom_line(aes(
    group = subgroup,
    color = subgroup,
    y = estimate,
    x = zone
  ), size = 1.5) +
  geom_pointrange(
    aes(
      x = zone,
      fill = subgroup,
      y = estimate,
      ymin = estimate - std.error,
      ymax = estimate + std.error
    ),
    size = 1,
    shape = 21
  ) +
  scale_fill_manual(values = colors_subgroup2) +
  scale_color_manual(values = colors_subgroup2) +
  geom_signif(
    data = anno_line,
    aes(
      group = rand,
      xmin = start,
      xmax = end,
      annotations = label,
      y_position = y
    ),
    vjust = .6,
    textsize = 8,
    tip_length = 0,
    manual = TRUE
  ) +
  #geom_text(aes(x=zone, y=10, label=stars), nudge_x=0, angle=0, size=6, data=peakmean.plot) +
  facet_grid(cols = vars(sex), labeller = labeller(sex = sex.labs)) +
  theme_Publication() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 20),
    axis.text.x = element_text(
      size = 22,
      angle = 45,
      vjust = 1,
      hjust = 1
    ),
    axis.text.y = element_text(size = 22),
    axis.title.y = element_text(size = 24),
    axis.title.x = element_text(size = 24)
  ) +
  xlab("Cortical Region (Posterior to Anterior)") + ylab("Peak Alpha Frequency (Hz)") +
  coord_cartesian(ylim = c(6, 11.5)) +
  theme(aspect.ratio = 1)
ggsave(filename = target_file.intraregion_figure)
out("Figure", target_file.intraregion_figure)

ggsave(plot = p.pafTopo,
       filename = "figshare/Figure4b.pdf",
       colormodel = "cmyk")
write_csv(x = estimates, file = "figshare/Figure4b_SourceData.csv")

df.source.avg <- df.source %>%
  group_by(eegid, subgroup, zone) %>% dplyr::summarize(value = mean(value))
labels_subgroup <- c("FXS_F", "TDC_F", "FXS_M", "TDC_M")

anno_bar <- pairs_by_group %>%
  mutate(sex = factor(
    sex,
    levels = c("F", "M"),
    labels = c("Female", "Male")
  )) %>%
  mutate(zone = factor(zone, labels_zone)) %>%
  rowwise() %>%
  mutate(
    region1 = str_split(contrast, " - ", simplify = TRUE)[1],
    region2 = str_split(contrast, " - ", simplify = TRUE)[2],
    region1 = factor(region1, levels = labels_subgroup),
    region2 = factor(region2, levels = labels_subgroup),
    xoffset = .3,
    xcenter = as.numeric(zone),
    xstart = xcenter - xoffset,
    xend = xcenter + xoffset,
    ystart = ifelse(str_detect(contrast, "FXS - TDC"), 1.25, 1.25),
    y = 9.9,
    label = add_sig_stars(adj.p),
    rand = rnorm(1)
  ) %>%
  relocate(region1, region2, .before = estimate) %>% ungroup() %>%
  mutate(start = as.numeric(region1),
         end = as.numeric(region2))

#anno_bar %>% View()
as.numeric(anno_bar$zone)

anno_bar %>% select(zone, rand, sex, xoffset, xcenter, label, y)

p.topoRegion <- estimates %>%
  mutate(subgroup = paste0(group, "_", sex)) %>%
  mutate(sex = factor(
    sex,
    levels = c("F", "M"),
    labels = c("Female", "Male")
  )) %>%
  ggplot() +
  coord_cartesian(ylim = c(6.5, 15)) +
  geom_col(
    aes(
      x = zone,
      group = subgroup,
      fill = subgroup,
      y = estimate
    ),
    color = "black",
    position = position_dodge2(1),
    alpha = .3
  ) +
  geom_point(
    aes(
      x = zone,
      y = value,
      fill = subgroup,
      group = subgroup
    ),
    pch = 21,
    size = 1,
    color = "black",
    position = position_jitterdodge(jitter.width = .1),
    alpha = 0.8,
    data = df.sublevel
  ) +
  geom_errorbar(aes(
    x = zone,
    y = estimate,
    ymin = estimate,
    ymax = estimate + std.error
  ),
  position = position_dodge2()) +
  geom_signif(
    data = anno_bar,
    aes(
      group = zone,
      xmin = xstart,
      xmax = xend,
      annotations = label,
      y_position = 13
    ),
    vjust = .4,
    textsize = 6,
    tip_length = .02,
    manual = TRUE
  ) +
  scale_fill_manual(values = colors_subgroup2) +
  theme_Publication() +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(
          angle = 45,
          vjust = 1,
          hjust = 1
        )) +
  facet_grid(cols = vars(sex)) +
  xlab("Cortical Region") +
  ylab("Peak Alpha Frequency (Hz)") +
  theme(aspect.ratio = .75)

p.topoRegion
ggsave(filename = target_file.group_figure,
       width = 6,
       height = 4)
out("Figure", target_file.group_figure)

ggsave(plot = p.topoRegion,
       filename = "figshare/Figure4a.pdf",
       colormodel = "cmyk")
write_csv(x = estimates, file = "figshare/Figure4a_Estimates_SourceData.csv")
write_csv(x = df.sublevel, file = "figshare/Figure4a_SubjectLevel_SourceData.csv")

#==============================================================================#
# ANALYSIS #3      ============================================================#
#                  Node level view of Alpha Peak Frequency
#==============================================================================#

#==============================================================================#
# GGSEG ATLAS FOR SUPPLEMENT
#==============================================================================#
ggseg(atlas = dkextra, mapping = aes(fill = region)) +
  scale_fill_brain("dk")
ggsave(filename = target_file.atlas_figure)
out("Figure", target_file.atlas_figure)

# Higher resolution view
df.source2 <- import.peakfreq.source %>%
  left_join(nodeinfo, by = c("Electrode" = "labelclean")) %>%
  left_join(clininfo, by = c("eegid" = "eegid")) %>%
  rename(value = peakFreq) %>%
  drop_na()

## =============================================================================
## SAVE FOR CORRELATIONS
## =============================================================================
df.paf <- df.source2 %>% select(eegid, label, value, region, RSN)
save(df.paf, file = target_file.PafCorrelations)
out("Model", target_file.PafCorrelations)

#==============================================================================#
# LME 2: Group X Sex X Label
#==============================================================================#
# Base LME model includes group, sex, lower band, and resting state network (RSN)
fit.2 <-
  nlme::lme(
    value ~ group * sex * label,
    random = ~ 1 | eegid,
    method = "ML",
    correlation = corCompSymm(form =  ~ 1 |
                                eegid),
    data = df.source2
  )
anova(fit.2)  # laterality did not matter

# least-squared means
emc.2 <- emmeans(fit.2, ~ group * label)
estimates.2 <-
  broom::tidy(emc.2, conf.int = TRUE, conf.level = .95)

# get paired significant contrasts
pairs_by_group.2 <-
  broom::tidy(pairs(emc.2, by = c("label"), adjust = 'none')) %>%
  group_by(contrast) %>% mutate(adj.p = p.adjust(p.value, method = "fdr", n = n())) %>%
  select(-term, -null.value) %>% ungroup() %>% select(-p.value)

df.ggseg <- pairs_by_group.2 %>% rename(Name = label) %>%
  mutate(ggseglabel = paste0(
    str_split(Name, " ", simplify = TRUE)[, 2],
    "h_",
    str_split(Name, " ", simplify = TRUE)[, 1]
  ) %>% str_to_lower()) %>%
  rename(label = ggseglabel,
         value = statistic)
df.ggseg


#==============================================================================#
# LME 3: Node-level Plot
#==============================================================================#
# T-values

p.tvals <- df.ggseg %>%
  ggseg(atlas = dkextra,
        view = "lateral",
        mapping = aes(fill = value)) +
  scale_fill_gradientn(
    colours = c("blue", "white", "firebrick"),
    na.value = "gray",
    limits = c(-5, 5)
  ) +
  theme(legend.position = "bottom" , legend.text = element_text(size = 7))
ggsave(filename = target_file.node_figure)
out("Figure", target_file.node_figure)

#==============================================================================#
# LME 3: Node-Level Table
#==============================================================================#

est_for_merge <-
  estimates.2 %>% select(group, label, estimate, std.error) %>%
  mutate(meanse = paste0(rd(estimate, 2), "\u00B1", rd(std.error, 2))) %>%
  select(-estimate,-std.error) %>% pivot_wider(names_from = group, values_from = meanse)


ft.elecpeak <- pairs_by_group.2 %>% left_join(est_for_merge,
                                              by = c("label")) %>%
  left_join(nodeinfo %>% select(label, region), by = c("label")) %>%
  relocate(c(region, FXS, TDC), .before = estimate) %>% filter(adj.p <=  .05) %>%
  arrange(statistic) %>% select(-contrast) %>%
  mutate(
    estimate = weights::rd(estimate, 2),
    std.error = weights::rd(std.error, 2),
    statistic = weights::rd(statistic, 1),
    adj.p = scales::pvalue(adj.p),
    label = str_remove(label, " ")
  ) %>%
  left_join(nodeinfo %>% select(publish, labelclean),
            by = c("label" = "labelclean")) %>%
  relocate(publish) %>%
  arrange(publish) %>%
  select(-label) %>%
  rename(label = publish) %>%
  flextable() %>%
  set_header_labels(
    label = "Node",
    region = "Region",
    sex = "Sex",
    zone = "Region",
    estimate = "FXS-TDC",
    std.error = "SE",
    df = "DF",
    statistic = "F",
    adj.p = "5% FDR"
  ) %>%
  flextable::merge_v(j = 1) %>%
  flextable::valign(j = 1, valign = "top") %>%
  autofit() %>% fix_border_issues()

ft.elecpeak %>%
  save_as_docx(path = target_file.nodetable)

out("Table", target_file.nodetable)

#==============================================================================#
# EXPORT TABLE
# see above
#==============================================================================#
