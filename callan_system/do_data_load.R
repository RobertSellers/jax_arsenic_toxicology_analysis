library(tidyverse)
library(janitor)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# contains harmony plate collection class builder
source("scripts/harmony_utils_v3.R")
# various add-on scripts and stats tools
source("scripts/custom_tools.R")

# Warning - this currently expects a temp_files data dir 
# Set overwrite to TRUE when underlying data dir has been changed
harmony_collection <- harmony_create_collection(
  dir = "./test",
  overwrite = FALSE
  )

## Data treatments

# consult https://cran.r-project.org/web/packages/janitor/vignettes/janitor.html
do_focus_collection$all_plates <- do_focus_collection$all_plates %>% 
  dplyr::rename(dose = concentration) %>%
  add_column(sex = NA, .after="individual") %>%
  filter(dir %notin% c("dC20244_A1__20200918T11_06_17Measurement_1_e6",
                         "dC20244_A2__20200914T11_13_52Measurement_1_e6",
                         "dC20244_B1__20200914T09_49_07Measurement_2_e6",
                         "dC20244_B2__20200918T09_22_05Measurement_1_e6",
                         "d20222_B2_Cal__20200818T11_40_06Measurement_3_e5"),
  dose != 0.5, individual != '11_20') %>%
  mutate(plate=ifelse((dir=='d20223_A1_Cal__20200817T17_03_14Measurement_2_e7'), 'A1', plate)) %>%
  mutate(plate=ifelse((dir=='dC20236_B2__20200828T14_45_19Measurement_1_e8'), 'B2', plate)) %>%
  mutate(plate=ifelse((dir=='dC20236_B1__20200828T13_12_27Measurement_1_e7'), 'B1', plate)) %>%
  mutate(plate=ifelse((dir=='d20223_B1_Cal__20200827T10_48_20Measurement_1_e7'), 'B1', plate)) %>%
  # remove last 3 columns from 
  mutate(keep = ifelse((dir == 'dC20255_B2__20200926T13_00_37Measurement_2_e6' & as.numeric(column) < 10),'keep',
                ifelse((dir == 'dC20264_B1__20200928T12_28_55Measurement_1_e5' & as.numeric(column) < 9), 'keep',
                ifelse((dir %notin% c('dC20264_B1__20200928T12_28_55Measurement_1_e5','dC20255_B2__20200926T13_00_37Measurement_2_e6')), 'keep','remove')))) %>%
  filter(keep == 'keep') %>%
  mutate(mouse = paste0(gsub('_','.',individual),"_", plate)) %>%
  dplyr::select(plate,row,column,dose,sex,individual, mouse, everything()) %>%
  # hardcode mapped replicates to paired IDs
   mutate(group_pair = LETTERS702[group_indices(., substr(plate,1,1), sub("\\_.*", "", dir))]) 