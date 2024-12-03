###
### Family outcomes and genetic explanations
### Basic script to set up data and analysis for discrete-time
### event history models (and related analyses) for fertility, union formation and
### dissolution outcomes from registers.
###
library(dplyr)
library(tidyr)
library(forcats)
library(tibble)
library(purrr)
library(readr)
library(stringr)
library(ggplot2)
library(data.table)
library(broom)
library(here)
library(huxtable)
###
### Data construction part
###
prefixit <- function(df,prefix, except = NULL) {
  orignames <- names(df)
  newnames <- str_c(prefix, orignames)
  newnames[orignames%in% except] <- orignames[orignames%in% except]
  
  setNames(df, nm=newnames)
}
myscale <- function(x) as.vector(scale.default(x))

### Get basic demography data
fasteoppl <- fread(here("../../../data/registers/SSB/01_data/data_v1.0/csv/POPULATION_FASTE_OPPL.csv")) %>%
  mutate(cohort = as.numeric(str_sub(foedsels_aar_mnd,1,4))) %>% 
  mutate(female = as.numeric(kjoenn)-1)
basic_info <- select(fasteoppl, w19_0634_lnr, cohort, female, fodeland, doeds_aar_mnd, mor_lnr, far_lnr)

### Get union history data
unions <- fread(here("../../../data/registers/SSB/01_data/data_v1.0/csv/POPULATION_FAMILIE_SAMLIV.csv"))


### Genetics ID numbers etc
genetics_info <- fread(here("../../../projects/openflux/data/prepared-genetic-data-individuals.csv"))  %>%
    filter(w19_0634_lnr!="")  %>%
    group_by(w19_0634_lnr)  %>%
    slice(1)  %>%
    ungroup  %>%
    mutate(across(contains("_PGS"), myscale))

nrow(genetics_info)

## Check up when first moba kids were born
mobafamilies <- genetics_info  %>%
    select(w19_0634_lnr, role, PREG_ID_2601, barn_nr)  %>%
    group_by(PREG_ID_2601,barn_nr)  %>%
    spread(role, w19_0634_lnr)  %>%
    rename(child_w19_0634_lnr = child_)  %>%
    rename(father_w19_0634_lnr = father)  %>%
    rename(mother_w19_0634_lnr = mother)  %>%
    ungroup

moba_birthyears <- mobafamilies %>%
    left_join(prefixit(basic_info, prefix="child_"))  %>%
    gather(key=parent, value=w19_0634_lnr, father_w19_0634_lnr, mother_w19_0634_lnr)  %>%
    arrange(w19_0634_lnr, child_cohort)  %>%
    group_by(w19_0634_lnr)  %>%
    slice(1)  %>%
    select(child_w19_0634_lnr, child_cohort, w19_0634_lnr)


### All relevant adults
relevant_adults <- genetics_info  %>%
    filter(role!="child_")  %>%
    filter(SENTRIXID!="")  %>%
    select(w19_0634_lnr) %>%
    distinct %>%
    filter(w19_0634_lnr!="")
nrow(relevant_adults)    


### Get personyears at risk of dissolution
our_unions <- relevant_adults  %>%
    left_join(unions,by="w19_0634_lnr")

personyears_unions <- our_unions %>%
    gather(var,info, -w19_0634_lnr) %>%
    mutate(year = as.numeric(str_sub(var, start=-4))) %>%
    mutate(var = str_sub(var,start=1,end=-6)) %>%
    arrange(w19_0634_lnr, year) %>%
    group_by(w19_0634_lnr, year) %>%
    spread(var,info)

### All personyears from age 18 to 2018
pys <- relevant_adults  %>%
    left_join(basic_info)  %>%
    left_join(personyears_unions)  %>%
    filter(year>(cohort+17))  %>%
    arrange(w19_0634_lnr, year)  %>%
    group_by(w19_0634_lnr)  %>%
    mutate(dissolved = case_when(
               ((lnr_sambo=="") & (lag(lnr_sambo,1)!="")) ~ 1,
               ((lnr_sambo!="") & (lag(lnr_sambo,1)!="") & (lnr_sambo!=lag(lnr_sambo,1))) ~ 1,
               TRUE ~ 0))  %>%
    mutate(atrisk_dissolution = if_else(((lnr_sambo!="") | (dissolved==1)),1,0))  %>%
    
    mutate(partner_w19_0634_lnr = lnr_sambo)  %>%
    mutate(partner_w19_0634_lnr = if_else(partner_w19_0634_lnr=="", lag(partner_w19_0634_lnr),partner_w19_0634_lnr))  %>%
    
    left_join(prefixit(basic_info, "partner_"))  %>%
    left_join(moba_birthyears)

personyears_atrisk <- pys  %>%
    filter(atrisk_dissolution==1)  %>%
    mutate(ego_w19_0634_lnr = w19_0634_lnr)  %>%
    left_join(prefixit(genetics_info, prefix="ego_"))  %>%
    left_join(prefixit(genetics_info, prefix="partner_"))

nrow(personyears_atrisk)
fwrite(personyears_atrisk, here("eventhistory/data/personyears-v5-all.csv"))


personyears_atrisk_prospective <- filter(personyears_atrisk, year>=child_cohort)
nrow(personyears_atrisk_prospective)
fwrite(personyears_atrisk_prospective, here("eventhistory/data/personyears-v5-prospective.csv"))
