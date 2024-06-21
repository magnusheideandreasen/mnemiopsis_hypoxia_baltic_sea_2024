##### MODELS FOR MANUSCRIPT 1 UPDATED mid April 2024

library(readr)
library(glmmTMB)
library(DHARMa)
library(dbplyr)
library(ggplot2)
library(ggExtra)
library(sjPlot)
library(sjstats)
library(bbmle)
library(splines)
library(ggformula)
library(ggpubr)
library(dplyr)
library(parallel)
library(MuMIn)
library(patchwork)

## data prep youngs #########
setwd("C:/Users/mhean/OneDrive - Danmarks Tekniske Universitet/Dokumenter/Manuscripts/m1_balticjellycommunity/m1_data/m1_data_new/")
dat <- read.csv("al544_modeldata_april2024.csv")
View(dat)
dat <- dat[!dat$vlookup=="FL1multinetmidi1",] # lacks <10
dat <- dat[!dat$vlookup=="FL1multinetmidi2",] # lacks <10
dat <- dat[!dat$vlookup=="FL1multinetmidi3",] # lacks <10
dat <- dat[!dat$vlookup=="FL1multinetmidi4",] # lacks <10
dat <- dat[!dat$vlookup=="FL1multinetmidi5",] # lacks <10
dat <- dat[!dat$vlookup=="FL2multinetmidi1",] # lacks <10
dat <- dat[!dat$vlookup=="FL2multinetmidi2",] # lacks <10
dat <- dat[!dat$vlookup=="FL2multinetmidi3",] # lacks <10
dat <- dat[!dat$vlookup=="FL2multinetmidi4",] # lacks <10
dat <- dat[!dat$vlookup=="FL2multinetmidi5",] # lacks <10
dat <- dat[!dat$vlookup=="FL3multinetmidi1",] # lacks <10
dat <- dat[!dat$vlookup=="FL3multinetmidi2",] # lacks <10
dat <- dat[!dat$vlookup=="FL3multinetmidi3",] # lacks <10
dat <- dat[!dat$vlookup=="FL3multinetmidi4",] # lacks <10
dat <- dat[!dat$vlookup=="FL3multinetmidi5",] # lacks <10
dat <- dat[!dat$vlookup=="LB1Imultinetmidi2",] # outlier; low subsample vol
dat <- dat[!dat$vlookup=="BY32multinetmaxi1",] # only when scaling (currently misses ctd)
dat <- dat[!dat$vlookup=="SW8multinetmidi5",] # only when scaling (currently misses ctd)
dat <- dat[!dat$vlookup=="SW9multinetmidi5",] # no subsample vol
dat <- dat[!dat$vlookup=="SW2WP20",] # no subsample vol

dat$pos <- numFactor(dat$ycoord, dat$xcoord)
dat$region <- factor(dat$region)

dat$net_type[dat$net_type == "multinetmaxi"]<-"multinet" # calling both multinet types by the same name
dat$net_type[dat$net_type == "multinetmidi"]<-"multinet"
dat$surfsallevel <- as.character(dat$surfsallevel)

## keeping only relevant variables
keep <- c("surf_sal","surf_temp_c","oxy_mgl","bottom_oxy_mgl","avg_sample_depth","fluoro",
          "bel10_juv","bel10_tra","bel10_yadu",
          "subf_juv","subf_tra","subf_yadu","subf_adu","subf_aur",
          "nethaulsize","region","pos","vol_filt_m3",
          "n_juv","n_tra","n_yadu","n_adu","n_aur",
          "net_type","juv_dens","tra_dens","yadu_dens","adu_dens","aur_dens","surfsallevel","oxylevel_2.8_9.0")

# removing NAs
datnzs <- na.omit(dat[,keep])

# for later descaling of variables
m <- mean(datnzs$surf_sal)
sd <- sd(datnzs$surf_sal) 
m_oxy <- mean(datnzs$oxy_mgl)
sd_oxy <- sd(datnzs$oxy_mgl) 

#scaling pred variables
datnzs[,keep[1:6]] <- scale(datnzs[,keep[1:6]])

save(datnzs, m, sd, m_oxy, sd_oxy, datnzs, file="organized data.Rdata")

## data prep adult + aurelia ###################################################
dat <- read.csv("al544_modeldata_april2024.csv")
dat <- dat[!dat$vlookup=="LB1Imultinetmidi2",] # outlier; low subsample vol
dat <- dat[!dat$vlookup=="BY32multinetmaxi1",] # only when scaling (currently misses ctd)
dat <- dat[!dat$vlookup=="SW8multinetmidi5",] # only when scaling (currently misses ctd)

dat$pos <- numFactor(dat$ycoord, dat$xcoord)
dat$region <- factor(dat$region)
dat$net_type[dat$net_type == "multinetmaxi"]<-"multinet" # calling both multinet types by the same name
dat$net_type[dat$net_type == "multinetmidi"]<-"multinet"
dat <- dat[!dat$net_type=="bongo",] # omitting bongo net counts
dat$surfsallevel <- as.character(dat$surfsallevel)

## keeping only relevant variables
keep <- c("surf_sal","surf_temp_c","oxy_mgl","bottom_oxy_mgl","avg_sample_depth","fluoro",
          "bel10_juv","bel10_tra","bel10_yadu",
          "subf_juv","subf_tra","subf_yadu","subf_adu","subf_aur",
          "nethaulsize","region","pos","vol_filt_m3",
          "n_juv","n_tra","n_yadu","n_adu","n_aur",
          "net_type","juv_dens","tra_dens","yadu_dens","adu_dens","aur_dens","surfsallevel","oxylevel_2.8_9.0")

# removing NAs
datnzs2 <- na.omit(dat[,keep])

# for later descaling of variables
m2 <- mean(datnzs2$surf_sal)
sd2 <- sd(datnzs2$surf_sal) 
m_oxy2 <- mean(datnzs2$oxy_mgl)
sd_oxy2 <- sd(datnzs2$oxy_mgl) 

#scaling pred variables
datnzs2[,keep[1:6]] <- scale(datnzs2[,keep[1:6]])

save(datnzs2, m, sd, m_oxy, sd_oxy, datnzs2, file="organized data2.Rdata")

## custom theme (theme1) ############
saved_theme <- theme_bw() +   
  theme(panel.border = element_rect(colour = "black", linewidth = 1), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(plot.title = element_text(family = "Times New Roman",
                                  hjust=0.5, size = 12, face = "bold"),
        axis.title = element_text(family = "Times New Roman",size = 11),
        axis.text.x = element_text(family = "Times New Roman",size = 11),
        axis.text.y = element_text(family = "Times New Roman",size = 11),
        legend.title = element_text(family = "Times New Roman",size = 12),
        legend.text = element_text(family = "Times New Roman",size = 11),
        strip.background = element_rect(colour= "white", fill="white")) +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        plot.margin = unit(c(0.15, 0.15, 0.15, 0.15), "cm")) 

saved_theme %>% saveRDS('saved_theme.rds')
rm(saved_theme)
theme1 <- readRDS('saved_theme.rds')

## custom theme (theme2) ############ ARIAL
saved_theme2 <- theme_bw() +   
  theme(panel.border = element_rect(colour = "black", linewidth = 1), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(plot.title = element_text(family = "Arial",
                                  hjust=0.5, size = 12, face = "bold"),
        axis.title = element_text(family = "Arial",size = 11),
        axis.text.x = element_text(family = "Arial",size = 11),
        axis.text.y = element_text(family = "Arial",size = 11),
        legend.title = element_text(family = "Arial",size = 12),
        legend.text = element_text(family = "Arial",size = 11),
        strip.background = element_rect(colour= "white", fill="white")) +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        plot.margin = unit(c(0.15, 0.15, 0.15, 0.15), "cm")) 

saved_theme2 %>% saveRDS('saved_theme2.rds')
rm(saved_theme2)
theme2 <- readRDS('saved_theme2.rds')
