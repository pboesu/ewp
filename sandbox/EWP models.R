# Running Philipp's EWP models

library(tidyverse)

#wd<-"C:/Users/jerem/Dropbox/BTO/NRS_clutch_size/data/processed/dv/"

#wilwa<-readRDS(paste0(wd,"wilwa.rds"))

wilwa <- readRDS('sandbox/wilwa.rds')

# Creating a year factor for wilwa
wilwa$yearf<-as.factor(wilwa$year)

wilwa_trial<- wilwa %>% filter(max_num_eggs > 0 & max_num_eggs < 9, flag_error <1, year >=2000) #%>% droplevels()

# Re-defining wilwa2 as wilwa
wilwa2<- wilwa %>% filter(max_num_eggs > 0 & max_num_eggs < 9, flag_error <1, year >1960)


year_tab_wilwa<- count(wilwa_trial, yearf)

hist(wilwa_trial$max_num_eggs)

plot(wilwa_trial$max_num_eggs ~ wilwa_trial$min_first_egg)

library(ewp)

#### Comparable models with EWP ####
fit_null <- ewp_reg(max_num_eggs ~ 1, data = wilwa_trial)
summary(fit_null)

fit_min <- ewp_reg(max_num_eggs ~ min_first_egg, data = wilwa_trial)
summary(fit_min)

Model3_min<-glm(max_num_eggs~min_first_egg, family = quasipoisson, data = wilwa_trial)
summary(Model3_min)

fit_min2 <- ewp_reg(max_num_eggs ~ min_first_egg, data = wilwa2)
summary(fit_min2)

Model3_min2<-glm(max_num_eggs~min_first_egg, family = quasipoisson, data = wilwa2)
summary(Model3_min2)

fit_year<-ewp_reg(max_num_eggs~yearf, hessian = FALSE, data = wilwa_trial)
summary(fit_year)

Model3a_trial<-glm(max_num_eggs~yearf, family = quasipoisson, data = wilwa_trial)
summary(Model3a_trial)

fit_year2<-ewp_reg(max_num_eggs~yearf, hessian = FALSE, data = wilwa2)
summary(fit_year2)

Model3a<-glm(max_num_eggs~yearf, family = quasipoisson, data = wilwa2)
summary(Model3a)
