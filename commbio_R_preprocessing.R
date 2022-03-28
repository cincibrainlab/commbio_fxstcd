# Preprocessing Quality Assurance
# Supplemental Table 2
# Author: E. Pedapati
# Version: 3/27/2022

pacman::p_load(labelled, compareGroups, tidyverse, BiocManager, filesstrings)

# source(file = "_Common.R")

source.demo <- "https://figshare.com/ndownloader/files/34517351"
source.qi   <- "https://figshare.com/ndownloader/files/34517363"

df.clin <- read_csv(source.demo)

fx.add.group_assignments <- function(df, df.clin) {
  df %>% left_join( df.clin %>% dplyr::select(eegid, group, subgroup, sex, mgroup), 
                    by=c("eegid"="eegid"))}


model.qi <- read_csv(source.qi) %>% 
  fx.add.group_assignments(df.clin) %>% 
  select(-mgroup,-sex,-subgroup, -subj_subfolder)

# sanity check
model.qi %>% group_by(group) %>% summarize(n=n())



dat.qi <- model.qi %>%
  set_variable_labels(group   = "Group",
                      proc_xmax_raw    = "Total Duration(s)",
                      proc_xmax_epoch  = "Clean Duration(s)",
                      epoch_trials     = "Remaining Trials", 
                      proc_badchans    = "Bad Channels",
                      proc_removeComps = "Artifact Components")


dat.table_qi <- descrTable(group ~ .,dat.qi %>% 
                             select(-eegid,-id), hide.no = "no")

target_table = dat.table_qi
target_file = "commbio_supplTable2.docx"
export2word(target_table, file = target_file)
file.move(files = target_file, "output/", overwrite = TRUE)
file.remove(str_replace(target_file, "docx", "Rmd"))
