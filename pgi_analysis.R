###
### Polygenic score analysis for genetics and partnership dissolution paper
###
library(tidyverse)
library(ggplot2)
library(tibble)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(readr)
library(forcats)
library(broom)
library(bife)
library(data.table)
library(broom)
library(here)
library(huxtable)
library(stargazer)

### Reset workspace
rm(list=ls())

### Scale function
myscale <- function(x) as.vector(scale.default(x))

# Get prepared data
if (!exists("personyears")) personyears <- fread(here("pgi_analysis/data/personyears-v5-prospective.csv"))

eventhistorydf <- personyears %>% 
  filter(female!=partner_female) %>% 
  rename(ego_sivilstand = sivilstand) %>% 
  rename(ego_mor_lnr = mor_lnr) %>% 
  rename(ego_far_lnr = far_lnr) %>% 
  rename(ego_doeds_aar_mnd = doeds_aar_mnd) %>% 
  rename(ego_lnr_ektefelle = lnr_ektefelle) %>% 
  rename(ego_lnr_familie = lnr_familie) %>% 
  rename(ego_lnr_sambo = lnr_sambo) %>% 
  rename(ego_dissolved = dissolved) %>% 
  rename(ego_atrisk_dissolution = atrisk_dissolution) %>% 
  rename(ego_female = female) %>% 
  rename(ego_fodeland = fodeland) %>% 
  rename(ego_cohort = cohort) %>% 
  rename(ego_year = year) %>% 
  select(-starts_with("child_"), -w19_0634_lnr) %>% 
  mutate(across(ends_with("_PGS"), myscale)) %>% 
  unite(col=parents_lnr, ego_far_lnr, ego_mor_lnr, sep="_", remove=TRUE, na.rm = FALSE)


pgi_names <- c("Educational attainment" = "ego_EA4_excl23andMe_mobaref_PGS",
               "Subjective well-being" = "ego_WellBeing_Baselmans_2019_PGS",
               "Age at 1st birth" = "ego_age_at_first_birth_pooled_mobaref_PGS",
               "Age at 1st sex" = "ego_age_at_first_sex_pooled_mobaref_PGS",
               "Cigarettes per day" = "ego_CigarettesPerDay_Liu_2019_PGS",
               "Drinks per week" = "ego_DrinksPerWeek_Liu_2019_PGS",
               "Number of sex partners" = "ego_NumSex_KarlssonLinnér_2019_mobaref_PGS",
               "Number of children ever born" = "ego_NEB_Mathieson_2023_mobaref_PGS",
               "Being a 'morning person'" = "ego_morningperson_jones_2019_PGS",
               "Loneliness" = "ego_loneliness_day_2018_mobaref_PGS",
               "Autism spectrum disorder" = "ego_ASD_Grove_2019_1kgreffrq_PGS",
               "ADHD" = "ego_ADHD_Demontis_2023_PGS",
               "Depression" = "ego_Depression_Howard_2019_PGS",
               "Neuroticism" = "ego_Neuroticism_Turley_2018_PGS",
               "BMI" = "ego_BMI_GIANT_2018_PGS",
               "Height" = "ego_GIANT_HEIGHT_YENGO_2022_PGS"
)
egopcs <- str_c("ego_PC",seq(1:10),collapse = " + ")
controls <- " + ego_genotyping_batch + ego_genotyping_batch "
multimod <- str_c("ego_dissolved ~ ",paste(pgi_names,collapse = " + ")," + ",egopcs,controls)


# Estimate models of PGIs for all, men and women
all_multimod <- glm(formula=multimod, family="binomial", data=eventhistorydf)
men_multimod <- glm(formula=multimod, family="binomial", data=filter(eventhistorydf, ego_female==0))
women_multimod <- glm(formula=multimod, family="binomial", data=filter(eventhistorydf, ego_female==1))

coefs_all <- tidy(all_multimod) %>%
  filter(grepl("_PGS$", term)) 
coefs_women <- tidy(women_multimod) %>%
  filter(grepl("_PGS$", term)) 
coefs_men <- tidy(men_multimod) %>%
  filter(grepl("_PGS$", term)) 

### Produce figure 1
fig1_data <- bind_rows(list("All" = coefs_all, "Women" = coefs_women, "Men" = coefs_men), .id="Sex") %>% 
  mutate(conf.low = estimate - 1.96 * std.error, 
         conf.high = estimate + 1.96 * std.error) %>% 
  mutate(or = exp(estimate), or_conf.low = exp(conf.low), or_conf.high = exp(conf.high)) %>% 
  mutate(is_significant = if_else(p.value<0.05, TRUE, FALSE)) %>% 
  left_join(data.frame(label = names(pgi_names), term = pgi_names), by="term") %>% 
  filter(Sex!="All") %>% 
  group_by(term) %>% 
  mutate(meanest = mean(estimate)) %>% 
  arrange(meanest) %>% 
  ungroup %>% 
  mutate(label= fct_reorder(label, -1*meanest))

fig1_colors <- c("Women" = "darkorchid1", "Men" = "#00CC99")
figure1 <- ggplot(fig1_data, aes(x = label, y = or, ymin = or_conf.low, ymax = or_conf.high, color = Sex, shape = Sex, alpha = is_significant)) +
  geom_pointrange(size = 0.7, position = position_dodge(width = 0.3)) + 
  scale_y_continuous(breaks = c(0.85, 0.9, 0.95, 1, 1.05, 1.1, 1.15)) +
  geom_hline(yintercept = 1, col="black", linewidth = 0.5, linetype = 2) +
  labs(x=NULL, y=NULL) +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.text=element_text(size=13,
                                 color = "#333333"),
        axis.text=element_text(size=13,
                               color = "black"),
        text = element_text(size=13),
        panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_line(color = "lightgrey",
                                  linetype = 3,
                                  linewidth = 0.3),
        axis.line = element_line(color = "gray", size = 0.3),
        legend.title = element_blank(),
        legend.key = element_rect(fill = "transparent", color = NA)) +
  coord_flip() + 
  scale_color_manual(values = fig1_colors) +
  scale_alpha_manual(values = c(0.3, 1), guide = "none") 

ggsave(filename = here("pgi_analysis/output/figure1.png"), figure1, device="png")
ggsave(filename = here("pgi_analysis/output/figure1.pdf"), figure1, device="pdf")


### Produce numerical results table for figure 1
tableS1 <- fig1_data %>%
  arrange(Sex, label) %>% 
  select(Sex = Sex, PGI = label, Beta = estimate, SE = std.error, P = p.value,LowerCI = conf.low, UpperCI = conf.high)
xlsx::write.xlsx(tableS1,file = here("pgi_analysis/output/tableS1.xlsx"))
knitr::kable(tableS1, format = "markdown", digits = 3) %>% 
  cat(file = here("pgi_analysis/output/table1S.md"))





###
### Within-family PGI analysis
###
siblingdf <- personyears  %>%
  group_by(w19_0634_lnr)  %>%
  arrange(w19_0634_lnr, year)  %>%
  mutate(everdissolved = max(dissolved))  %>%
  mutate(exposures = n())  %>%
  mutate(events = sum(dissolved))  %>%    
  mutate(ratedissolved = events/exposures)  %>%
  select(w19_0634_lnr, mor_lnr, far_lnr, events, exposures, ratedissolved, everdissolved, contains("_PGS"))  %>%
  slice(1)  %>%
  ungroup  %>%
  unite(far_lnr, mor_lnr, col="parents_lnr", remove=F, sep=":")  %>%
  group_by(parents_lnr) %>% 
  mutate(n_sibs = n()) %>% 
  filter(n_sibs>1) %>% 
  filter(n_sibs<10)
nrow(siblingdf)


### Esztimate family fixed-effects models 
withinest <- function(pgs)
{
  mod <- bife::bife(data=siblingdf, as.formula(str_c("everdissolved ~ ",pgs,"|parents_lnr")))
  return(mod)
}
single_within_models <- map(pgi_names, withinest) %>% setNames(nm = pgi_names)

fixbife <- function(m) {
  sm <- summary(m)
  dsm <-as.data.frame(sm$cm)  %>%
    setNames(nm=c("estimate","std.error","statistic","p.value"))  %>%
    rownames_to_column("term")  %>%
    mutate(conf.low = estimate-1.96*std.error, conf.high=estimate+1.96*std.error) %>% 
    mutate(or = exp(estimate), or_conf.low=exp(conf.low), or_conf.high = exp(conf.high))
  return(dsm)
}

single_within_dfs <- map(single_within_models, fixbife) %>% bind_rows %>% 
  mutate(Model = "Single") %>% left_join(tibble(term = pgi_names, PGI = names(pgi_names)), by="term")

### Within-family model with multiple PGIs
multiple_within_mod <- bife::bife(data=siblingdf, formula = as.formula(str_c("everdissolved~", str_c(pgi_names, collapse=" + "), " | parents_lnr")))
multiple_within_df <- multiple_within_mod %>%
  fixbife  %>%
  filter(str_detect(term, "_PGS")) %>% 
  arrange(desc(estimate)) %>% 
  left_join(tibble(term = pgi_names, PGI = names(pgi_names)), by="term") %>% 
  mutate(Model = "Multiple")

# Find order of within PGI effects
pgi_order <- multiple_within_df  %>%
  pull(PGI)

# Combined plot data
within_both_dfs <- bind_rows(multiple_within_df, single_within_dfs) %>% 
  arrange(desc(estimate)) %>% 
  mutate(is_significant = if_else(p.value<0.05,T,F)) %>% 
  mutate(PGI_factor = factor(PGI, level = pgi_order))

# Show lists of PGI names to make sure the ordering is correct.
select(within_both_dfs, contains("PGI")) %>% distinct

### Produce figure 2
within_colours <- c("goldenrod","deepskyblue2")
figure2 <- within_both_dfs %>% 
  ggplot(aes(x=PGI_factor, y = or, ymin=or_conf.low, ymax=or_conf.high, color=Model, alpha=is_significant)) +
  geom_pointrange(size=0.7, position=position_dodge2(width = 0.3)) + 
  geom_hline(yintercept = 1, col="black", linewidth = 0.5, linetype = 2) +
  scale_color_manual(values = within_colours) + 
  scale_alpha_manual(values = c(0.3,1), guide="none") +
  scale_y_continuous(breaks = c(0.5,0.75,1.0,1.25,1.5,1.75)) +
  xlab("") + ylab("") +
  coord_flip()  +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.text=element_text(size=13, color="#333333"),
        axis.text=element_text(size=13, color="black"),
        text = element_text(size=13),
        panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_line(color="lightgrey",linetype=3, size=0.3),
        axis.line = element_line(color="gray", size=0.3),
        legend.title=element_blank(),
        legend.key = element_rect(fill="transparent", color=NA)
  )

ggsave(here("pgi_analysis/output/figure2.png"), figure2, device="png")
ggsave(here("pgi_analysis/output/figure2.pdf"), figure2, device="pdf")


### Produce numerical results table for figure 2
tableS2 <- within_both_dfs %>%
  arrange(Model, PGI) %>% 
  select(Model, PGI, Beta = estimate, SE = std.error, P = p.value,LowerCI = conf.low, UpperCI = conf.high)
xlsx::write.xlsx(tableS1,file = here("pgi_analysis/output/tableS2.xlsx"))
knitr::kable(tableS1, format = "markdown", digits = 3) %>% 
  cat(file = here("pgi_analysis/output/table2S.md"))



###
### Supplementary materials
###

### Squared PGI difference predictors
diff_names <- gsub("^ego_","sqdiff_",pgi_names)
partner_pgis <- gsub("^ego_","partner_",pgi_names)

sqdiffs <- select(eventhistorydf, ego_w19_0634_lnr)
for (pgi in 1:length(pgi_names))  {
  pgi_ego <- pull(eventhistorydf, pgi_names[[pgi]])
  pgi_partner <- pull(eventhistorydf, partner_pgis[[pgi]])
  sqdiffs <- cbind(sqdiffs, (pgi_ego-pgi_partner)^2) %>% 
    setNames(nm=c("ego_w19_0634_lnr",diff_names[1:pgi]))
}
diff_df <-bind_cols(eventhistorydf, select(sqdiffs,-ego_w19_0634_lnr))
diffpcs <- str_c("ego_PC1","partner_PC1","ego_PC2","partner_PC2","ego_PC3","partner_PC3","ego_PC4","partner_PC4","ego_PC5","partner_PC5",sep=" + ",collapse = " + ")


# Estimate models of PGIs for all, single and multiple
single_diff_models <- map(diff_names, function(x) {
  model <- str_c("ego_dissolved ~ ", paste(x,collapse=" + "), " + ", diffpcs, controls)
})
est_single_sqdiff <- function(mod) {
  glm(formula = mod, family="binomial", data=diff_df)
}
single_sqdiff_models <- map(single_diff_models, est_single_sqdiff)
single_sqdiff_df <- map(single_sqdiff_models, tidy) %>% 
  bind_rows %>% 
  filter(str_detect(term, "sqdiff"))

sqdiff_multimod <- str_c("ego_dissolved ~ ",paste(diff_names,collapse = " + ")," + ",diffpcs,controls)
all_sqdiff_multimod <- glm(formula=sqdiff_multimod, family="binomial", data=diff_df)
sqdiff_coefs_multimod <- tidy(all_sqdiff_multimod) %>%
  filter(grepl("_PGS$", term)) 

sqdiff_coefs_all <- bind_rows(list("Single" = single_sqdiff_df, "Multiple" = sqdiff_coefs_multimod), .id="Model")

### Produce figure S5: squared difference models
figS5_data <- sqdiff_coefs_all %>% 
  mutate(conf.low = estimate - 1.96 * std.error, 
         conf.high = estimate + 1.96 * std.error) %>% 
  mutate(or = exp(estimate), or_conf.low = exp(conf.low), or_conf.high = exp(conf.high)) %>% 
  mutate(is_significant = if_else(p.value<0.05, TRUE, FALSE)) %>% 
  left_join(data.frame(label = names(pgi_names), term = diff_names), by="term") %>% 
  group_by(term) %>% 
  mutate(meanest = mean(estimate)) %>% 
  arrange(meanest) %>% 
  ungroup %>% 
  mutate(label= fct_reorder(label, -1*meanest))

figS5_colors <- c("Single" = "darkorchid1", "Multiple" = "#00CC99")
figureS5 <- ggplot(figS5_data, aes(x = label, y = or, ymin = or_conf.low, ymax = or_conf.high, color = Model, shape = Model, alpha = is_significant)) +
  geom_pointrange(size = 0.7, position = position_dodge(width = 0.3)) + 
  geom_hline(yintercept = 1, col="black", linewidth = 0.5, linetype = 2) +
  labs(x=NULL, y=NULL) +
  theme(legend.position = "none",
        legend.direction = "horizontal",
        legend.text=element_text(size=13,
                                 color = "#333333"),
        axis.text=element_text(size=13,
                               color = "black"),
        text = element_text(size=13),
        panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_line(color = "lightgrey",
                                  linetype = 3,
                                  size = 0.3),
        axis.line = element_line(color = "gray", size = 0.3),
        legend.title = element_blank(),
        legend.key = element_rect(fill = "transparent", color = NA)) +
  coord_flip() + 
  scale_y_continuous(breaks = c(0.925, 0.95, 0.975, 1.0, 1.025, 1.05, 1.075)) +
  scale_color_manual(values = figS5_colors, guide="none") +
  scale_alpha_manual(values = c(0.3, 1), guide = "none")

ggsave(filename = here("pgi_analysis/output/figureS5.png"), figureS5, device="png")
ggsave(filename = here("pgi_analysis/output/figureS5.pdf"), figureS5, device="pdf")


### Basic descriptive statistics for text and Table 1
num_men <- filter(eventhistorydf, ego_female==0) %>% pull(ego_w19_0634_lnr) %>% unique %>% length
num_women <- filter(eventhistorydf, ego_female==1) %>% pull(ego_w19_0634_lnr) %>% unique %>% length
num_persons <- num_men + num_women
prop_female_ind <- num_women/num_persons
numpy_men <- filter(eventhistorydf,ego_female==0, ego_IID!="") %>% nrow
numpy_women <- filter(eventhistorydf,ego_female==1, ego_IID!="") %>% nrow
numpy <- numpy_men+numpy_women
prop_female_py <- numpy_women/numpy

eventhistorydf$ego_cohort %>% mean

num_men + num_women
nrow(eventhistorydf)/num_persons

numdissolutions_py <- sum(eventhistorydf$ego_dissolved)
dissolutionrate_py <- numdissolutions_py/numpy
tstats <- eventhistorydf %>% group_by(ego_w19_0634_lnr) %>% summarize(n = n(),period = mean(ego_year), dissolved = max(ego_dissolved))
avgobservations <- tstats$n
meanperiod <- tstats$period
table(tstats$dissolved)/nrow(tstats)
