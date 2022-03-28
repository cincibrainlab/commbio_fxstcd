

pacman::p_load(R.matlab,
               reshape2,
               tidyverse,
               dichromat,
               ggsci,
               viridis,
               nlme,
               emmeans,
               flextable)
source(file = "_Common.R")

target_file <- "output/commbio_powFig3_.RData"

# Load Datasets
csv.group <- "https://figshare.com/ndownloader/files/34517366"
dat.clin <- read_csv(csv.group) %>%
  rename(subgroup2 = subgroup) %>%
  rename(sex2 = sex) %>%
  mutate(sex = str_extract(sex2, "^.{1}"),
         subgroup = paste0(group, "_", sex)) %>%
  relocate(eegid, group, sex, subgroup) %>%
  mutate(subgroup = subgroup2)

csv.mne.relpow <- "https://figshare.com/ndownloader/files/34518185"
dat.pow <- read_csv(csv.mne.relpow)

# Models (computationally intensive can used precaluclated results)
modelByRegion <- "log(value)~group*sex*bandname*Region"
modelByRsn    <- "log(value)~group*sex*bandname*RSN"

# Function #1 LME (with Log link)
eegpower_lme <- function(modelformula) {
  fit <- lme(
    modelformula,
    random = ~ 1 | eegid,
    correlation = corCompSymm(form =  ~ 1 | eegid),
    data = dat.pow
  )
}

run_model = FALSE
if (run_model == TRUE) {
  # Region Model Results and EMMEANS
  fit.region <- eegpower_lme(as.formula(modelByRegion))
  emc.region <- emmeans(fit.region, ~ group * sex * bandname * region)
  # RSN Model Results and EMMEANS
  fit.rsn    <- eegpower_lme(as.formula(modelByRsn))
  emc.rsn <- emmeans(fit.rsn, ~ group * sex * bandname * RSN)
  save(fit.region, fit.rsn, emc.region, emc.rsn, file = "commbio_Fig3_modelResults.RData")
} else {
  # RSN Model Results and EMMEANS
  load(url("https://figshare.com/ndownloader/files/34532204"))
}

# all estimates
estimates.region <-
  broom::tidy(emc.region, conf.int = TRUE, conf.level = .95)
estimates.rsn    <-
  broom::tidy(emc.rsn, conf.int = TRUE, conf.level = .95)

fx_calcPairs <- function(emmeansObj, contrastString) {
  pairs_by_group <-
    broom::tidy(pairs(emmeansObj, by = contrastString, adjust = 'none')) %>%
    group_by(contrast) %>% mutate(adj.p = p.adjust(p.value, method = "fdr", n = n())) %>%
    select(-term, -null.value)
}

pairs.region.group  <-
  fx_calcPairs(emc.region, c("bandname", "region", "sex"))
pairs.rsn.group     <-
  fx_calcPairs(emc.rsn, c("bandname", "RSN", "sex"))
pairs.region.sex    <-
  fx_calcPairs(emc.region, c("bandname", "region", "group"))
pairs.rsn.sex       <-
  fx_calcPairs(emc.rsn, c("bandname", "RSN", "group"))

# Flextable for Group By Region
ft.pairs.tmp <-
  pairs.region.group %>% ungroup() %>%  select(sex, bandname, region, estimate, std.error, adj.p) %>%
  mutate(side = substr(region, 1, 1),
         region = substr(region, 2, nchar(region))) %>%
  mutate(
    estimate = weights::rd(estimate, 2),
    estimate = paste0(
      estimate,
      "\u00B1",
      weights::rd(std.error, 2),
      add_sig_stars(adj.p, cutoffs = c(0.01, 0.001, 0.0001))
    )
  ) %>%
  select(-adj.p, -std.error) %>% relocate(side, region, sex) %>%
  pivot_wider(names_from = bandname, values_from = estimate) %>%
  mutate(
    side = factor(
      side,
      levels = c("L", "R"),
      labels = c("Left", "Right")
    ),
    region = factor(
      region,
      levels = c("O", "L", "P", "T", "C",  "F",  "PF"),
      labels = c(
        "Occipital",
        "Limbic",
        "Parietal",
        "Temporal",
        "Central",
        "Frontal",
        "Prefrontal"
      )
    )
  ) %>%
  arrange(side, region)

ft.pairs <- ft.pairs.tmp %>% flextable()
ft.pairs <- set_header_labels(ft.pairs,
                              side = "L/R",
                              region = "Region",
                              sex = "Sex") %>% merge_v(j = c(1, 2, 3)) %>% valign(j =
                                                                                    c(1, 2, 3), valign = "top") %>%
  fix_border_issues() %>%
  autofit()

target_file.region.pairs <-
  str_replace(target_file, ".RData", "groupbyRegion_pairwise.docx")
ft.pairs %>% save_as_docx(path = target_file.region.pairs)


rm(ft.pairs.tmp, ft.pairs)

# Data for Group by Region Plot
plotdata.region <- pairs.region.group %>% mutate(side = substr(region,1,1),
                                             region = substr(region,2,nchar(region))) %>% 
  mutate(bandname = factor(bandname, levels=c("delta","theta","alpha1","alpha2","beta","gamma1","gamma2")),
         side = factor(side, levels=c("L","R"), labels = c("Left", "Right")),
         sex= factor(sex, levels=c("F","M"), labels=c("Females: FXS-Control", "Males: FXS-Control")),
         region = factor(region, levels=c("O","L","P", "T", "C",  "F",  "PF" ),
                         labels = c("Occipital","Limbic", "Parietal","Temporal","Central","Frontal","Prefrontal"))) %>% 
  rename(x=bandname, y=statistic, fill=region)

plotdata.region %>% select(y, adj.p) %>% filter(adj.p < .003) %>% arrange(-adj.p)
write_csv(x = plotdata.region, file = "figshare/Figure3a_SourceData.csv")

# Flextable for Group By Region

ft.pairs.tmp <- pairs.rsn.group %>% ungroup() %>%  
  select(sex, bandname,RSN, estimate, std.error, adj.p) %>% 
  mutate(bandname = factor(bandname, 
                           levels=c("delta","theta","alpha1",
                                    "alpha2","beta","gamma1","gamma2"))) %>% 
  arrange(RSN, bandname) %>%  
  mutate(estimate = weights::rd(estimate,2),
         estimate = paste0(estimate, "\u00B1", 
                           weights::rd(std.error,2), 
                           add_sig_stars(adj.p, cutoffs = c(0.01, 0.001, 0.0001)))) %>% 
  select(-adj.p,-std.error) %>% relocate(RSN, sex ) %>% 
  pivot_wider(names_from = bandname, values_from = estimate) %>% 
  mutate(sex= factor(sex, levels=c("M","F"), labels=c("Male","Female")))

ft.pairs.tmp$RSN
ft.pairs <- ft.pairs.tmp %>% flextable()
ft.pairs <- set_header_labels(ft.pairs,
                              RSN = "RSN",
                              sex = "Sex") %>% merge_v(j=c(1,2,3)) %>% 
  valign(j=c(1,2,3),valign = "top") %>% 
  fix_border_issues()%>% 
  autofit()
ft.pairs

target_file.rsn.pairs <-
  str_replace(target_file, ".RData", "groupbyRsn_pairwise.docx")
ft.pairs %>% save_as_docx(path = target_file.rsn.pairs)

rm(ft.pairs.tmp, ft.pairs)

# Create Plotting Data
plotdata.RSN <- pairs.rsn.group %>% 
  mutate(bandname = factor(bandname, levels=c("delta","theta","alpha1",
                                              "alpha2","beta","gamma1","gamma2")),
         RSN = factor(RSN, levels = c("DMN","DAN","SAN","AUD","VIS","other")),
         sex= factor(sex, levels=c("F","M"), 
                     labels=c("Females: FXS-Control", "Males: FXS-Control"))) %>% 
  rename(x=bandname, y=statistic, fill=RSN)

write_csv(x = plotdata.RSN, file = "figshare/Figure3b_SourceData.csv")

# Flextable for Sex by Region
ft.pairs.tmp <- pairs.region.sex %>% ungroup() %>%  
  select(group, bandname,region, estimate, std.error, adj.p) %>% 
  mutate(side = substr(region,1,1),
         region = substr(region,2, nchar(region))) %>% 
  mutate(estimate = weights::rd(estimate,2),
         estimate = paste0(estimate, "\u00B1", weights::rd(std.error,2), 
                           add_sig_stars(adj.p, cutoffs = c(0.01, 0.001, 0.0001)))) %>% 
  select(-adj.p,-std.error) %>% relocate(side, region, group ) %>% 
  pivot_wider(names_from = bandname, values_from = estimate) %>% 
  mutate(side = factor(side, levels=c("L","R"), labels = c("Left", "Right")),
         region = factor(region, levels=c("O","L","P", "T", "C",  "F",  "PF" ),
                         labels = c("Occipital","Limbic", "Parietal","Temporal",
                                    "Central","Frontal","Prefrontal"))) %>% 
  arrange(side,region)

ft.pairs <- ft.pairs.tmp %>% flextable()
ft.pairs <- set_header_labels(ft.pairs, 
                              side = "L/R",
                              region = "Region",
                              group = "Group") %>% 
  merge_v(j=c(1,2,3)) %>% valign(j=c(1,2,3),valign = "top") %>% 
  fix_border_issues()%>% 
  autofit()
ft.pairs

target_file.region.pairs.sex <-
  str_replace(target_file, ".RData", "sexbyRegion_pairwise.docx")
ft.pairs %>% save_as_docx(path = target_file.region.pairs.sex)
rm(ft.pairs.tmp, ft.pairs)

# Save data for plotting
plotdata.region.sex <- pairs.region.sex %>% mutate(side = substr(region,1,1),
                                             region = substr(region,2,nchar(region))) %>% 
  mutate(bandname = factor(bandname, levels=c("delta","theta","alpha1","alpha2","beta","gamma1","gamma2")),
         side = factor(side, levels=c("L","R"), labels = c("Left", "Right")),
         #sex= factor(sex, levels=c("Female","Male"), labels=c("Female: FXS-Control", "Males: FXS-Control")),
         region = factor(region, levels=c("O","L","P", "T", "C",  "F",  "PF" ),
                         labels = c("Occipital","Limbic", "Parietal","Temporal","Central","Frontal","Prefrontal")),
         group = factor(group, levels=c("TDC","FXS"), 
                        labels=c("Control: Females-Males","FXS: Females-Males"))) %>% 
  rename(x=bandname, y=statistic, fill=region)

plotdata.region.sex %>% select(y, adj.p) %>% filter(adj.p < .003) %>% arrange(-adj.p)
write_csv(x = plotdata.region.sex, file = "figshare/Figure3c_Region_SourceData.csv")

# Flextable for Sex by Rsn

ft.pairs.tmp <- pairs.rsn.sex %>% arrange(RSN) %>% ungroup() %>%  
  select(group, bandname,RSN, estimate, std.error, adj.p) %>% 
  mutate(estimate = weights::rd(estimate,2),
         estimate = paste0(estimate, "\u00B1", 
                           weights::rd(std.error,2), add_sig_stars(adj.p, cutoffs = c(0.01, 0.001, 0.0001)))) %>% 
  select(-adj.p,-std.error) %>% relocate(RSN, group ) %>% 
  pivot_wider(names_from = bandname, values_from = estimate)
ft.pairs.tmp$RSN
ft.pairs <- ft.pairs.tmp %>% flextable()
ft.pairs <- set_header_labels(ft.pairs,
                              RSN = "RSN",
                              sex = "Sex") %>% merge_v(j=c(1,2,3)) %>% 
  valign(j=c(1,2,3),valign = "top") %>% 
  fix_border_issues()%>% 
  autofit()
ft.pairs

target_file.rsn.pairs.sex <-
  str_replace(target_file, ".RData", "sexbyRsn_pairwise.docx")
ft.pairs %>% save_as_docx(path = target_file.rsn.pairs.sex)
rm(ft.pairs.tmp, ft.pairs)

# save data for plotting
plotdata.RSN <- pairs.rsn.sex %>% 
  mutate(bandname = factor(bandname, levels=c("delta","theta","alpha1","alpha2","beta","gamma1","gamma2")),
         RSN = factor(RSN, levels = c("DMN","DAN","SAN","AUD","VIS","other")),
         group = factor(group, levels=c("TDC","FXS"), 
                        labels=c("Control: Females-Males","FXS: Females-Males"))) %>% 
  rename(x=bandname, y=statistic, fill=RSN)
write_csv(x = plotdata.RSN, file = "figshare/Figure3c_Rsn_SourceData.csv")

# Create Figures
showtext_opts(dpi = 300)
showtext_auto(enable = TRUE)

fig3a_region <- read_csv("figshare/Figure3a_SourceData.csv") %>%
  mutate(plotname = "region") %>%
  rename(facetvar = sex) %>% 
  relocate(plotname)
fig3b_rsn <- read_csv("figshare/Figure3b_SourceData.csv") %>%
  mutate(plotname = "rsn") %>%
  rename(facetvar = sex) %>% 
  relocate(plotname)
fig3c_region <- read_csv("figshare/Figure3c_Region_SourceData.csv") %>%
  mutate(plotname = "regionsex") %>%
  rename(facetvar = group) %>% 
  relocate(plotname)
fig3c_rsn <- read_csv("figshare/Figure3c_Rsn_SourceData.csv") %>%
  mutate(plotname = "rsnsex") %>%
  rename(facetvar = group) %>% 
  relocate(plotname)

fig3df <- bind_rows(
  fig3a_region,
  fig3b_rsn,
  fig3c_region,
  fig3c_rsn) %>% 
  mutate(facetvar =factor(facetvar)) %>% 
  mutate(x = factor(x, 
                    levels=c("delta","theta","alpha1","alpha2","beta","gamma1","gamma2")),
         fill = factor(fill, 
                       levels = c("Occipital","Limbic", 
                                  "Parietal","Temporal",
                                  "Central","Frontal","Prefrontal", 
                                  "DMN","DAN","SAN","AUD","VIS","other"))) 
#  filter(facetvar != "Control: Females-Males")
fig3df %>% distinct(facetvar)
p.base = list()

for (i in unique(fig3df$plotname)) {
  print(i)
  p.base[[i]] <- fig3df %>%
    filter(plotname == i) %>%
    ggplot() +
    geom_col(aes(x = x, y = y, fill = fill),
             color = "black", size = 0.05,
             position = position_dodge2()
    ) +
    facet_wrap(~facetvar) +
    theme_Publication() +
    theme(aspect.ratio = .5) +
    coord_cartesian(ylim = c(-10, 10)) +
    theme(
      legend.direction = "horizontal",
      strip.background = element_blank(),
      strip.text = element_text(hjust = 0, size = 8),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
    )
}

plotspecific <- NA
plotspecific <- tribble(
  ~plotname,  ~ylab, ~filllab, ~p05cut, ~p001cut, ~colors,
  "region", "(FXS - Control)", "Region", 2, 3.38, "scale_fill_colorblind",
  "rsn", "(FXS - Control)", "RSN", 2, 3.54, "scale_fill_viridis_d",
  "regionsex", "(Female - Male)", "Region", 2, 3.36,  "scale_fill_colorblind",
  "rsnsex", "(Female - Male)", "RSN", 2.07, 3.3, "scale_fill_viridis_d"
)

p.final = list()
for(i in 1 : nrow(plotspecific)){
  ps <- plotspecific[i,]
  plotname <- ps$plotname
  p05 <- ps$p05cut
  p001 <- ps$p001cut
  p05l <- p05 + ((p001-p05)/2)
  p001l <- p001 +.2
  ylab <- paste0("Log Power T-scores\n", ps$ylab)
  xlab <- "Frequency Band"
  colorscale <- ps$colors
  filllab <- ps$filllab
  p.final[[i]] <- p.base[[i]] + 
    geom_hline(yintercept = c(p05, -p05), 
               linetype = "dotted", size=.25) +
    geom_hline(yintercept = c(p001,-p001), 
               linetype = "dotted",
               size=.25) +
    ylab(ylab) +
    xlab(xlab) +
    labs(fill = filllab) +
    eval(parse_expr(colorscale))() +
    theme(aspect.ratio = .5) +
    #coord_cartesian(ylim = c(-10, 10)) +
    theme(text = element_text(size=8),
          legend.direction = "horizontal", 
          strip.background = element_blank(),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
    ) +
    annotate(geom = "text", label="*", x = c(1,1), vjust=.5,
             y = c(p05l,-p05-.6), size=6/.pt) +
    annotate(geom = "text", label="***", x = c(1,1), vjust=.5,
             y = c(p001l,-p001-.6), size=6/.pt)
  savename <- paste0("figshare/Figure3_",plotname,".pdf")
  ggsave(
    plot =  p.final[[i]],
    filename = savename,
    colormodel = "cmyk", units = "in"
  )
}

p.final[[1]]
p.final[[2]]
p.final[[3]]
p.final[[4]]