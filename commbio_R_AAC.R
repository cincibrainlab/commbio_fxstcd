#==============================================================================#
# AAC: GLOBAL / LOCAL VERSION
#==============================================================================#
# remotes::install_github('m-clark/mixedup')
source("_Common.R")
pacman::p_load(R.matlab, tidyverse, here, weights, flextable, reshape2,
               rstatix,effects,ggsegDefaultExtra, glmmTMB,hdf5r, DescTools, nlme, emmeans,
               mixedup, ggthemes, ggseg, afex)
pacman::p_load(ggsci, ClinReport)


basename    <- "MneJunAAC" # Edit
prefix      <- paste0("model_", basename)
data_file   <- NA # any RData inputs (or NA)

target_file <- "output/commbio_AAC.RData"
target_fig<- "figshare/commbio_AAC.RData"

#==============================================================================#
# INPUT                                                                        #
#==============================================================================#
df.aac <- read_csv("https://figshare.com/ndownloader/files/34535618")

#==============================================================================#
# ASSIGNED OUTPUT FILES USING GENERIC target_file VARIABLE
#==============================================================================#

target_file_forCorr <- str_replace(target_file, ".RData", "_forCorr.RData") 
target_file_pairs <- str_replace(target_file, ".RData", "table_pairwise.docx")

#==============================================================================#
# COMPUTATIONS
#==============================================================================#

# factor order and labels
levels_subgroup  = c("TDC_M","FXS_M","TDC_F","FXS_F")
levels_RSN       = c("DMN","DAN","SAN","AUD","VIS","other")
levels_upperband = c("gamma1","gamma2","epsilon")
levels_lowerband = c("theta","alpha1","alpha2")

# # Convert to Fisher-Z scores
# df.aac.addZ <- import.aac %>%  mutate(z = DescTools::FisherZ(value))
# 
# # add atlas and subject assignments
# df.aac <- df.aac.addZ %>% left_join(clininfo, by=c("eegid"="eegid")) %>% 
#   left_join(nodeinfo %>% select(labelclean, RSN, region), 
#             by=c("label"="labelclean")) %>% 
#   mutate(subgroup = factor(subgroup, levels=levels_subgroup),
#          RSN = factor(RSN, levels=levels_RSN),
#          upperband = factor(upperband, levels=levels_upperband),
#          lowerband = factor(lowerband, levels=levels_lowerband)) %>% 
#   filter(upperband  %notin% c("epsilon","gamma2"))



# check number of rows  141 subjects * 3 lowerbands * 2 powertype  * 68 nodes
# = 57,528
stopifnot(nrow(df.aac) == 57528)

#==============================================================================#
# MODEL: EFFECT OF GROUP AND RSN ON AAC
#==============================================================================#

# Base LME model includes group, sex, lower band, and resting state network (RSN)
# fit.1 <- nlme::lme(value~group*sex*lowerband*RSN, random = ~1|eegid,
#              correlation=corCompSymm(form=~1|eegid), data=df.aac %>% 
#                filter(powertype == "relative"))
# anova(fit.1)



load(url("https://figshare.com/ndownloader/files/34535801"))
nodeinfo <- read_csv("https://figshare.com/ndownloader/files/34517354")

# ============================= LMES  ===========================================
target_file_fit_region <- str_replace(target_file, ".RData", "_fit_region.RData")

# Sex is non-contributory to model, so will remove
model.aac.sex <- value ~ group * lowerband * RSN * sex
model.aac.rsn <- value ~ group * lowerband * RSN
model.aac.region <- value ~ group * lowerband * label

run_model = FALSE
if(run_model == TRUE){
  fit.aac.region <- lme(model.aac.region, random = ~1|eegid,
                        correlation=corCompSymm(form=~1|eegid),
                        data = df.aac %>% filter(powertype == "relative"))
  save(fit.aac.region, file = target_file_fit_region)
}else{
load(target_file_fit_region)
}
anova(fit.aac.region)
#sjPlot:: tab_model(fit.aac.region)
#sjPlot::plot_model(fit.aac.region)
target_file_fit_rsn <- str_replace(target_file, ".RData", "_fit_rsn.RData")

if(run_model == TRUE){
  
  fit.aac.rsn <- lme(model.aac.rsn, random = ~1|eegid,
                     correlation=corCompSymm(form=~1|eegid), 
                     data = df.aac %>% filter(powertype == "relative"))
  save(fit.aac.rsn, file = target_file_fit_rsn)
}
load(target_file_fit_rsn)
anova(fit.aac.rsn)
# sjPlot:: tab_model(fit.aac.rsn)
# summary(fit.aac.rsn)


# =============================MODEL PLOT==================================

library(ggsci)
p.model <- sjPlot::plot_model(fit.aac.rsn, 
                              type = "eff",
                              terms = c("RSN","lowerband","group")) +
  theme_minimal() + 
  scale_color_nejm() +
  theme(legend.position = "bottom",
        axis.title = element_text(size=20),
        axis.text = element_text(size=16),
        strip.text = element_text(size=20),
        axis.text.x = element_text(angle = 90, vjust = .5),
        legend.text = element_text(size=18),
        legend.title = element_text(size=18),
        aspect.ratio = 4        ) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0)
p.model

tidy(fit.aac.rsn)

ggsave(filename = str_replace(target_file, ".RData", "_figure_region_model_plot.pdf"))
out("Figure", str_replace(target_fig, ".RData", "_figure_region_model_plot.pdf"))

ggsave(plot = p.model, filename = "figshare/Figure5c_model.pdf", colormodel = "cmyk")
write_csv(x = tidy(fit.aac.rsn), file = "figshare/Figure5c_model_SourceData.csv")

# ============================= EMMEANS ========================================

emmeans.aac.rsn <- ~ group * lowerband * RSN 
emmeans.aac.region <- ~ group * lowerband * label

emc.region <- emmeans(fit.aac.region, emmeans.aac.region) 
emc.rsn <- emmeans(fit.aac.rsn, emmeans.aac.rsn)

estimates.region <- broom::tidy(emc.region, conf.int = TRUE, conf.level = .95)
estimates.rsn <- broom::tidy(emc.rsn, conf.int = TRUE, conf.level = .95)

pairs.region <- tidy(pairs(emc.region, by = c("lowerband","label"), adjust='none')) %>% 
  group_by(contrast) %>% mutate(adj.p = p.adjust(p.value,method = "fdr", n = n())) %>% 
  select(-term,-null.value) %>% ungroup()

pairs.rsn <- tidy(pairs(emc.rsn, by = c("lowerband","RSN"), adjust='none')) %>% 
  group_by(contrast) %>% mutate(adj.p = p.adjust(p.value,method = "fdr", n = n())) %>% 
  select(-term,-null.value) %>% ungroup()


pairs.rsn2 <- tidy(pairs(emc.rsn, by = c("lowerband","group"), adjust='none')) %>% 
  group_by(contrast) %>% mutate(adj.p = p.adjust(p.value,method = "fdr", n = n())) %>% 
  select(-term,-null.value) %>% ungroup() %>% filter(adj.p < .05) %>% arrange(abs(estimate))
pairs.rsn2
# ============================= ESTIMATE TABLES =================================

table.pairs.region <- pairs.region  %>%  
  select(lowerband, contrast, label, estimate,std.error, adj.p) %>% 
  mutate(estimate = weights::rd(estimate,2),
         estimate = paste0(estimate, add_sig_stars(adj.p, 
                                                   cutoffs = c(0.05, 0.001, 0.0001)))) %>% 
  select(-adj.p,-std.error) %>% 
  pivot_wider(names_from = lowerband, values_from = estimate) %>% 
  mutate(side = str_extract(label,"R|L"),
         label = str_remove(label, "R|L")) %>% 
  pivot_wider(id_cols = label, names_from = side, values_from = c(theta,alpha1,alpha2)) %>% 
  group_by(label) %>% 
  summarise(theta_L = na.omit(theta_L), 
            theta_R  = na.omit(theta_R),
            alpha1_L = na.omit(alpha1_L),
            alpha1_R = na.omit(alpha1_R), 
            alpha2_L =na.omit(alpha2_L),
            alpha2_R = na.omit(alpha2_R)) %>% 
  left_join(nodeinfo %>% select(publish,labelclean) %>% 
              mutate(labelclean = str_remove(labelclean,"R|L"),
                     publish = str_remove(publish," R| L") ) %>% distinct(), by=c("label"="labelclean")) %>% 
  relocate(publish) %>% select(-label)

ft.pairs.region <- table.pairs.region %>% flextable() %>% set_header_labels("theta_L"="L",
                                                                            "theta_R"="R",
                                                                            "alpha1_L"="L",
                                                                            "alpha1_R"="R",
                                                                            "alpha2_L"="L",
                                                                            "alpha2_R"="R",
                                                                            "publish" = "Node") %>%
  merge_h(i=1, part="header") %>% 
  align(i=1, j=2:6, part="header", align = "center") %>% 
  add_header_row(values = c("","Theta","Alpha1","Alpha2"), colwidths = c(1,2,2,2)  ) %>% 
  add_header_lines(values = "Pairwise (FXS-Control) gamma1 power-power CFC comparisons by cortical node")


table.pairs.rsn <- pairs.rsn %>% 
  select(RSN, lowerband, estimate, std.error, statistic, df, p.value, adj.p) %>% 
  mutate(estimate = weights::rd(estimate,2),
         estimate = paste0(estimate, "\u00B1", weights::rd(std.error,2)),
         sig = add_sig_stars(adj.p, 
                             cutoffs = c(0.05, 0.01, 0.001)),
         statistic = weights::rd(statistic,2),
         p.value = formatC(p.value, format = "e", digits = 1),
         adj.p = formatC(adj.p, format = "e", digits = 1)) %>%  select(-std.error,-p.value)

ft.pairs.rsn <- table.pairs.rsn %>% flextable()  %>% 
  set_header_labels(lowerband = "Lower Band",
                    estimate = "Estimate",
                    statistic = "Statistic",
                    df = "DF",
                    p.value = "p",
                    adj.p = "5% FDR p",
                    sig = "Sig.")  %>% autofit()  %>% 
  merge_v(j=1) %>% valign(j=1, valign = 'top') %>% fix_border_issues()



target_file_table_est_region <- str_replace(target_file, ".RData", "_est_region.docx")
target_file_table_est_rsn <- str_replace(target_file, ".RData", "_est_rsn.docx")

ft.pairs.region %>% save_as_docx(path = target_file_table_est_region)
ft.pairs.rsn %>% save_as_docx(path = target_file_table_est_rsn)

out("Table", target_file_table_est_region)
out("Table", target_file_table_est_rsn)

# ============================= VISUALIZATION ==================================

# ============================= REGION ==================================

labels_lowerband <- c("theta-gamma1 AAC","alpha1-gamma1 AAC","alpha2-gamma1 AAC")

jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))


plotdata.region <- pairs.region %>% 
  left_join(nodeinfo %>% select(Name, labelclean), 
            by=c("label"="labelclean")) %>% ungroup() %>% 
  mutate(ggseglabel = paste0(str_split(Name, " ",simplify=TRUE)[,2],"h_",
                             str_split(Name, " ",simplify=TRUE)[,1]) %>% 
           str_to_lower()) %>% 
  select(-label,-Name) %>% 
  rename(label = ggseglabel)  
dkcustom <- dkextra
dkcustom$data <- dkextra$data %>% filter(side %in% c('superior','lateral'))

plotdata.region %>% 
  select(lowerband, statistic, label) %>% 
  mutate(lowerband = factor(lowerband, 
                            levels=levels_lowerband, 
                            labels=labels_lowerband),
         statistic = statistic * -1) %>% 
  rename(fill = statistic) %>% 
  group_by(lowerband) %>% 
  ggplot() + 
  geom_brain(atlas = dkcustom, 
             position = position_brain(hemi ~ side),
             aes(fill=fill), 
             color="black") +
  facet_grid(cols = vars(lowerband)) +
  scale_fill_gradientn(colors = jet.colors(7), na.value="lightgray",limits=c(-6,6)) +
  theme_Publication() +
  labs(fill = "T-score") + 
  theme(strip.background = element_blank(),
        strip.text = element_text(size=18))
ggsave(filename = str_replace(target_file, ".RData", "_figure_region.pdf"))
out("Figure", str_replace(target_file, ".RData", "_figure_region.pdf"))

# significant regions only
plotdata.region %>% 
  select(lowerband, statistic, adj.p, label) %>% 
  mutate(lowerband = factor(lowerband, 
                            levels=levels_lowerband, 
                            labels=labels_lowerband),
         statistic = statistic * -1,
         issig = ifelse(adj.p <= .05, .25,0),
         issig = ifelse(adj.p <= .01, .5, issig),
         issig = ifelse(adj.p <= .001, 1, issig),
         issig = factor(issig, levels=c(0,0.25,.5,1), labels=c("ns","<.05","<.01","<.001"))) %>% 
  rename(fill = issig) %>% 
  group_by(lowerband) %>% 
  ggplot() + 
  geom_brain(atlas = dkcustom, 
             position = position_brain(hemi ~ side),
             aes(fill=fill), 
             color="black") +
  facet_grid(cols = vars(lowerband)) +
  scale_fill_grey(start=1, end=0) +
  theme_Publication() +
  labs(fill = "Adj. p") + 
  theme(strip.background = element_blank(),
        strip.text = element_text(size=18))
ggsave(filename = str_replace(target_file, ".RData", "_figure_region_sig.pdf"))
out("Figure", str_replace(target_file, ".RData", "_figure_region_sig.pdf"))


# ============================= RSN ==================================

#==============================================================================#
# OUTPUT FOR CORRELATIONS
#==============================================================================#
target_file_region_corr <- str_replace(target_file, ".RData", "_forCorr.RData") 
save(df.aac, file = target_file_forCorr)
out("Model",target_file_forCorr)


plotdata.sub <- df.aac %>% filter(powertype == "relative") %>% 
  group_by(RSN, eegid, group, lowerband, upperband) %>% 
  dplyr::summarize( z = mean(z)) %>%  
  mutate(lowerband = factor(lowerband, 
                            levels=levels_lowerband, 
                            labels=labels_lowerband)) %T>% 
  save(file = target_file_region_corr) %>% 
  rename(x = RSN, fill = group, facet=lowerband, y=z ) %>% 
  arrange(facet)

plotdata.aac <- df.aac %>% filter(powertype == "relative") %>% 
  group_by(RSN, eegid, group, lowerband, upperband) %>% 
  dplyr::summarize( z = mean(z)) %>% group_by(RSN, group, lowerband) %>% 
  dplyr::summarize( mean = mean(z), 
                    sd = sd(z)/sqrt(n()),
                    n = n()) %>% 
  mutate(ymin = mean - sd, 
         ymax = mean + sd,
         lowerband = factor(lowerband, levels=levels_lowerband,
                            labels=labels_lowerband)) %>% 
  rename(x = RSN, fill = group, facet=lowerband, y=mean, ymin = ymin, ymax=ymax ) %>% 
  arrange(facet)

plotdata.sig <- pairs.rsn %>% ungroup() %>% select(lowerband, RSN, adj.p) %>% 
  mutate(sig = add_sig_stars(adj.p, cutoffs = c(0.05, 0.001, 0.0001))) %>% 
  mutate(lowerband = factor(lowerband, levels=levels_lowerband,
                            labels=labels_lowerband)) %>% 
  rename(x = RSN,facet=lowerband, label = sig)

ylab = "AAC Correlation Coefficient"
xlab = "RSN"
filllab = "Group"

p.aacregion<- ggplot() +
  geom_errorbar(aes(x=x, y=y, ymin = ymin, group=fill,
                    ymax = ymax), width =.5,
                position = position_dodge(1), plotdata.aac) +
  geom_col(aes(x=x, y=y,fill=fill, group=fill),  color="black",
           position = "dodge", data=plotdata.aac) +
  geom_jitter(aes(x=x, y=y,fill=fill, group=fill), alpha=.2,  color="black", shape=21,
              position = position_jitterdodge(.1), data=plotdata.sub) +
  geom_vline(aes(xintercept = .5), data = subset(plotdata.aac, facet != "theta-gamma1 AAC")) +
  geom_text(aes(x=x,label=label,y=.15 ), size=4, data=plotdata.sig) +
  facet_wrap(vars(facet)) +
  xlab(xlab) +
  ylab(ylab) +
  labs(fill = filllab) +
  scale_fill_manual(values = colors_group)+
  theme_Publication() +  
  theme(aspect.ratio = 2.5) +
  #coord_cartesian(ylim = c(-10, 10)) +
  theme(legend.direction = "horizontal", strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size=12),
        axis.text = element_text(size=18),
        axis.title = element_text(size=20),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20)) 
ggsave(filename = str_replace(target_file, ".RData", "_figure_rsn.pdf"))
out("Figure", str_replace(target_file, ".RData", "_figure_rsn.pdf"))

ggsave(plot = p.aacregion, filename = "figshare/Figure5_rsn.pdf", colormodel = "cmyk")
write_csv(x = plotdata.aac, file = "figshare/Figure5c_rsn__SourceData.csv")


# pairs_by_group %>% filter(RSN != "other") %>% 
#   ggplot() + 
#   geom_pointrange(aes(x=RSN, y=estimate, fill=lowerband,ymin=estimate-std.error,
#                       ymax=estimate+std.error),
#                   shape=21,size=1) + theme_Publication() +
#   ylab("AAC Contrast (FXS-TDC)")

