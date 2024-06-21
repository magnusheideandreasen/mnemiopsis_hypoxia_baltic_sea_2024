theme1 <- readRDS('saved_theme.rds')

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
library(gamm4)
library(ggeffects)

# for faster processing
ncores <- detectCores() # max available cores
clusterType <- ifelse(Sys.info()["sysname"] == "Windows", "PSOCK", "FORK")
mycl <- parallel::makeCluster(ncores, type=clusterType)  
doParallel::registerDoParallel(mycl)
invisible(clusterEvalQ(mycl, library(glmmTMB , logical = TRUE)))


## MODELS ###################

# JUVENILES ####################################################################
load("organized data.Rdata") #data for juv, tra, and yadu
clusterExport(mycl, "datnzs")

# full juv model
m_juv_1 <- glmmTMB(n_juv ~ s(surf_sal, k=3) + 
									 	surf_temp_c+
									 	oxy_mgl + 
									 	avg_sample_depth + 
									 	net_type+ 
									 	offset(log(subf_juv) + log(vol_filt_m3))+ 
									 	(1 | region) + (1| nethaulsize),
									 	# disp=~surf_sal,
									 	family = nbinom2(link = "log"),
									 data=datnzs, na.action = na.fail)

summary(m_juv_1)
dd1 <- dredge(m_juv_1, rank = "AICc", cluster = mycl, fixed = ~cond(offset(log(subf_juv) +log(vol_filt_m3) )))

ddb7 <- subset(dd1, delta<7)
pred_juv <- get.models(dd1, delta==0)[[1]] # best model
pred_juv <- update(pred_juv, .~.-(1|region))
summary(pred_juv)

## make oxy_mgl average per high and low##########
## oxygen juv prediction############
newdat <- datnzs[datnzs$surfsallevel=="HIGH"|datnzs$surfsallevel=="LOW",]
newdat <- newdat[newdat$net_type == "multinet",]
newdat[newdat$surfsallevel=="HIGH","surf_sal"] <- mean(newdat[newdat$surfsallevel=="HIGH","surf_sal"])
newdat[newdat$surfsallevel=="LOW","surf_sal"] <- mean(newdat[newdat$surfsallevel=="LOW","surf_sal"])
newdat <- transform(newdat, 
										n_juv_scaled = n_juv/subf_juv/vol_filt_m3,
										avg_sample_depth = 0,
										fluoro=0,
										region = NA,
										nethaulsize = NA,
										vol_filt_m3 = 1,
										subf_juv = 1)

newdat$predjuv <- predict(pred_juv, newdat, type = "response") # juv mnemiopsis

ggplot(aes(x = oxy_mgl*sd_oxy+m_oxy, color=surfsallevel), data = newdat, size = 0.9)+
	geom_point(aes(y = n_juv_scaled))+
	geom_line(aes(y = predjuv))+	
	labs(x = "Bottom oxygen mg L-1", y = "Juvenile m-3")+
	theme1+
	theme(legend.spacing.y = unit(0.2, "cm"),
				legend.spacing.x = unit(0.0,"cm"),
				legend.title.align = 0.8)+
	scale_color_manual(labels = c("High (>15.0)",
																"Low (<9.0)"),
										 name = "Surface salinity",
										 values = c("grey10","grey70"))+
	labs(x = "", 
			 y = expression(paste(italic("M. leidyi"),juvenile~m^"-3", sep="")))

## salinity juv prediction########## 
newdat <- datnzs[datnzs$oxylevel_2.8_9.0=="HIGH"|datnzs$oxylevel_2.8_9.0=="LOW",]
newdat <- newdat[newdat$net_type == "multinet",]
newdat[newdat$oxylevel_2.8_9.0=="HIGH","oxy_mgl"] <- mean(newdat[newdat$oxylevel_2.8_9.0=="HIGH","oxy_mgl"])
newdat[newdat$oxylevel_2.8_9.0=="LOW","oxy_mgl"] <- mean(newdat[newdat$oxylevel_2.8_9.0=="LOW","oxy_mgl"])
newdat <- transform(newdat, 
										n_juv_scaled = n_juv/subf_juv/vol_filt_m3,
										avg_sample_depth = 0,
										region = NA,
										nethaulsize = NA,
										vol_filt_m3 = 1,
										subf_juv = 1)

newdat$predjuv <- predict(pred_juv, newdat, type = "response",re.form = NULL) # juv mnemiopsis

a_reg <- 
	ggplot(aes(x = surf_sal*sd+m, color=oxylevel_2.8_9.0), data = newdat, size = 0.9)+
	geom_point(aes(y = n_juv_scaled))+
	geom_line(aes(y = predjuv))+
	theme1+
	theme(legend.spacing.y = unit(0.2, "cm"),
				legend.spacing.x = unit(0.0,"cm"),
				legend.title.align = 0.8)+
	scale_color_manual(labels = c("High (>9.3)",
																"Low (<2.8)"),
										 name = expression(paste("Dissolved oxygen",~(mg~L^"-1"), sep="")),
										 values = c("grey10","grey70"))+
	labs(x = "", 
			 y = expression(paste(italic("M. leidyi"),~juvenile~m^"-3", sep="")))+
  ylim(0,75)
a_reg



# TRANSITIONALS ################################################################

## full tra models
temp =subset(datnzs, bel10_tra>0)
clusterExport(mycl, "temp")

m_tra_1 <- glmmTMB(n_tra ~ s(surf_sal, k=3) + 
									 	surf_temp_c+
									 	oxy_mgl + 
									 	avg_sample_depth + 
									 	net_type+ 
									 	offset(log(subf_tra) + log(vol_filt_m3) + log(bel10_tra))+ 
									 	(1 | region) + (1| nethaulsize),
									 									 	disp=~surf_sal,
									 family = nbinom2(link = "log"),
									 data=temp, na.action = na.fail)

summary(m_tra_1)

dd1 <- dredge(m_tra_1, rank = "AICc", cluster = mycl, fixed = ~cond(offset(log(subf_tra) +log(vol_filt_m3) + log(bel10_tra) )))
summary(get.models(dd1, delta==0)[[1]]) # best model 

ddb7 <- subset(dd1, delta<7)
View(ddb7)
pred_tra <- get.models(dd1, delta==0)[[1]] # best model
pred_tra <- update(pred_tra, .~.-(1|region))


## salinity tra prediction#############
newdat <- datnzs[datnzs$oxylevel_2.8_9.0=="HIGH"|datnzs$oxylevel_2.8_9.0=="LOW",]
newdat <- newdat[newdat$net_type == "multinet",]
newdat[newdat$oxylevel_2.8_9.0=="HIGH","oxy_mgl"] <- mean(newdat[newdat$oxylevel_2.8_9.0=="HIGH","oxy_mgl"])
newdat[newdat$oxylevel_2.8_9.0=="LOW","oxy_mgl"] <- mean(newdat[newdat$oxylevel_2.8_9.0=="LOW","oxy_mgl"])
newdat <- transform(newdat, 
										n_tra_scaled = n_tra/subf_tra/vol_filt_m3/bel10_tra,
										avg_sample_depth = 0,
										region = NA,
										nethaulsize = NA,
										vol_filt_m3 = 1,
										subf_tra = 1,
										bel10_tra=1)

newdat$predtra <- predict(pred_tra, newdat, type = "response", re.form = NULL) # tra mnemiopsis

b <- 
	ggplot(aes(x = surf_sal*sd+m, color=oxylevel_2.8_9.0), data = newdat, size = 0.9)+
	geom_point(aes(y = n_tra_scaled))+
	geom_line(aes(y = predtra))+
	theme1+
	theme(legend.spacing.y = unit(0.2, "cm"),
				legend.spacing.x = unit(0.0,"cm"),
				legend.title.align = 0.8)+
	scale_color_manual(labels = c("High (>9.3)",
																"Low (<2.8)"),
										 name = expression(paste("Dissolved oxygen",~(mg~L^"-1"), sep="")),
										 values = c("grey10","grey70"))+
	labs(x = "", 
			 y = expression(paste(italic("M. leidyi"),~transitional~m^"-3", sep="")))+
	ylim(0,40)
b 

# YOUNG ADULTS #################################################################

## full yadu model
m_yadu_1 <- glmmTMB(n_yadu ~ s(surf_sal, k = 3) + 
									 	surf_temp_c+
									 	oxy_mgl +
									 	avg_sample_depth + 
									 	net_type+ 
									 	offset(log(subf_yadu) + log(vol_filt_m3) + log(bel10_yadu))+ 
									 	(1 | region) + (1| nethaulsize),
									 	disp=~surf_sal, 
									 family = nbinom2(link = "log"),
									 data=datnzs, na.action = na.fail)

summary(m_yadu_1)

dd1 <- dredge(m_yadu_1, rank = "AICc", cluster = mycl, fixed = ~cond(offset(log(subf_yadu) +log(vol_filt_m3) + log(bel10_yadu) )))

ddb7 <- subset(dd1, delta<7)
pred_yadu <- get.models(dd1, delta==0)[[1]] # best model
pred_yadu <- update(pred_yadu, .~.-(1|region))
summary(pred_yadu)

## salinity yadu prediction######### 
newdat <- datnzs[datnzs$oxylevel_2.8_9.0=="HIGH"|datnzs$oxylevel_2.8_9.0=="LOW",]
newdat <- newdat[newdat$net_type == "multinet",]
newdat[newdat$oxylevel_2.8_9.0=="HIGH","oxy_mgl"] <- mean(newdat[newdat$oxylevel_2.8_9.0=="HIGH","oxy_mgl"])
newdat[newdat$oxylevel_2.8_9.0=="LOW","oxy_mgl"] <- mean(newdat[newdat$oxylevel_2.8_9.0=="LOW","oxy_mgl"])
newdat <- transform(newdat, 
                    n_yadu_scaled = n_yadu/subf_yadu/vol_filt_m3/bel10_yadu,
                    avg_sample_depth = 0,
                    surf_temp_c=0,
                    region = NA,
                    nethaulsize = NA,
                    vol_filt_m3 = 1,
                    subf_yadu = 1,
                    bel10_yadu=1)

newdat$predyadu <- predict(pred_yadu, newdat, type = "response", re.form = NULL) # yadu mnemiopsis

c_reg <- ggplot(aes(x = surf_sal*sd+m, color=oxylevel_2.8_9.0), data = newdat, size = 0.9)+
	geom_point(aes(y = n_yadu_scaled))+
	geom_line(aes(y = predyadu))+
	theme1+
	theme(legend.spacing.y = unit(0.2, "cm"),
				legend.spacing.x = unit(0.0,"cm"),
				legend.title = element_text(hjust=0.8))+
	scale_color_manual(labels = c("High (>9.3)",
																"Low (<2.8)"),
										 name = expression(paste("Dissolved oxygen",~(mg~L^"-1"), sep="")),
										 values = c("grey10","grey70"))+
	labs(x = "Salinity", 
			 y = expression(paste(italic("M. leidyi"),~young~adult~m^"-3", sep="")))+
  ylim(0,7)
c_reg


# ADULTS #######################################################################
load("organized data2.Rdata") #data for adults and aurelia

# full adu model
m_adu_1 <- glmmTMB(n_adu ~ s(surf_sal, k=3) + 
											surf_temp_c+
											s(oxy_mgl, k=3) + 
											avg_sample_depth + 
											net_type+ 
											offset(log(subf_adu) + log(vol_filt_m3))+ 
											(1 | region) + (1| nethaulsize),
											# disp=~surf_sal,
										    family = nbinom2(link = "log"),
										      data=datnzs, na.action = na.fail)

summary(m_adu_1)

dd1 <- dredge(m_adu_1, rank = "AICc", cluster = mycl, fixed = ~cond(offset(log(subf_adu) +log(vol_filt_m3) )))

# ddb7 <- subset(dd1, delta<7)
# View(ddb7)
# pred_adu <- get.models(dd1, delta==0)[[1]] # best model
# summary(pred_adu)

# second best adu model deltaAICc = 0.15 incl. sal spline
pred_adu <- glmmTMB(n_adu ~ s(surf_sal, k=3) +
                     s(oxy_mgl, k=3) + 
                     net_type+ 
                     offset(log(subf_adu) + log(vol_filt_m3))+ 
                     (1 | region) + (1| nethaulsize),
                   # disp=~surf_sal,
                   family = nbinom2(link = "log"),
                   data=datnzs, na.action = na.fail)

pred_adu <- update(pred_adu, .~.-(1|region))
summary(pred_adu)

## oxygen adu prediction######################
newdat <- datnzs[datnzs$surfsallevel=="HIGH"|datnzs$surfsallevel=="LOW",]
newdat <- newdat[newdat$net_type == "multinet",]
newdat[newdat$surfsallevel=="HIGH","surf_sal"] <- mean(newdat[newdat$surfsallevel=="HIGH","surf_sal"])
newdat[newdat$surfsallevel=="LOW","surf_sal"] <- mean(newdat[newdat$surfsallevel=="LOW","surf_sal"])
newdat <- transform(newdat, 
										n_adu_scaled = n_adu/subf_adu/vol_filt_m3,
										avg_sample_depth = 0,
										fluoro=0,
										region = NA,
										nethaulsize = NA,
										vol_filt_m3 = 1,
										subf_adu = 1)

newdat$predadu <- predict(pred_adu, newdat, type = "response") # adu mnemiopsis

ggplot(aes(x = oxy_mgl*sd_oxy+m_oxy, color=surfsallevel), data = newdat, size = 0.9)+
	geom_point(aes(y = n_adu_scaled))+
	geom_line(aes(y = predadu))+
	labs(x = "Bottom oxygen mg L-1", y = "adult m-3")+
	theme1+
	theme(legend.spacing.y = unit(0.2, "cm"),
				legend.spacing.x = unit(0.0,"cm"),
				legend.title.align = 0.8)+
	scale_color_manual(labels = c("High (>15.0)",
																"Low (<9.0)"),
										 name = "Surface salinity",
										 values = c("grey10","grey70"))+
	labs(x = "", 
			 y = expression(paste(italic("M. leidyi"),aduenile~m^"-3", sep="")))

## salinity adu prediction##############
newdat <- datnzs[datnzs$oxylevel_2.8_9.0=="HIGH"|datnzs$oxylevel_2.8_9.0=="LOW",]
newdat <- newdat[newdat$net_type == "multinet",]
newdat[newdat$oxylevel_2.8_9.0=="HIGH","oxy_mgl"] <- mean(newdat[newdat$oxylevel_2.8_9.0=="HIGH","oxy_mgl"])
newdat[newdat$oxylevel_2.8_9.0=="LOW","oxy_mgl"] <- mean(newdat[newdat$oxylevel_2.8_9.0=="LOW","oxy_mgl"])
newdat <- transform(newdat, 
										n_adu_scaled = n_adu/subf_adu/vol_filt_m3,
										avg_sample_depth = 0,
										fluoro=0,
										region = NA,
										nethaulsize = NA,
										vol_filt_m3 = 1,
										subf_adu = 1)

newdat$predadu <- predict(pred_adu, newdat, type = "response")#  adu mnemiopsis

d_reg <- ggplot(aes(x = surf_sal*sd+m, color=oxylevel_2.8_9.0), data = newdat, size = 0.9)+
	geom_point(aes(y = n_adu_scaled))+
	geom_line(aes(y = predadu))+
	theme1+
	theme(legend.spacing.y = unit(0.2, "cm"),
				legend.spacing.x = unit(0.0,"cm"),
				legend.title.align = 0.8)+
	scale_color_manual(labels = c("High (>9.3)",
																"Low (<2.8)"),
										 name = expression(paste("Dissolved oxygen",~(mg~L^"-1"), sep="")),
										 values = c("grey10","grey70"))+
	labs(x = "", 
			 y = expression(paste(italic("M. leidyi"),~adult~m^"-3", sep="")))+
	ylim(0,4)
d_reg 

# AURELIA ######################################################################

## full aur model
m_aur_1 <- glmmTMB(n_aur ~ surf_sal +
									 	surf_temp_c +
									 	oxy_mgl + 
									 	avg_sample_depth + 
									 	net_type+ 
									 	offset(log(subf_aur) + log(vol_filt_m3))+ 
									 	(1 | region) + 
									 	(1| nethaulsize),
									 #									 	disp=~surf_sal,
									 family = nbinom2(link = "log"),
									 data=datnzs, na.action = na.fail)

summary(m_aur_1)

dd1 <- dredge(m_aur_1, rank = "AICc", cluster = mycl, fixed = ~cond(offset(log(subf_aur) +log(vol_filt_m3) )))
# summary(get.models(dd1, delta==0)[[1]]) # best model
# ddb7 <- subset(dd1, delta<7)
# View(ddb7)
pred_aur <- get.models(dd1, delta==0)[[1]] # best model
pred_aur <- update(pred_aur, .~.-(1|region))

## salinity aur prediction###############
newdat <- datnzs[datnzs$oxylevel_2.8_9.0=="HIGH"|datnzs$oxylevel_2.8_9.0=="LOW",]
newdat <- newdat[newdat$net_type == "multinet",]
newdat[newdat$oxylevel_2.8_9.0=="HIGH","oxy_mgl"] <- mean(newdat[newdat$oxylevel_2.8_9.0=="HIGH","oxy_mgl"])
newdat[newdat$oxylevel_2.8_9.0=="LOW","oxy_mgl"] <- mean(newdat[newdat$oxylevel_2.8_9.0=="LOW","oxy_mgl"])

newdat <- transform(newdat, 
									 n_aur_scaled = n_aur/subf_aur/vol_filt_m3,
									 avg_sample_depth = 0,
										region = NA,
										nethaulsize = NA,
										vol_filt_m3 = 1,
									 subf_aur = 1)

newdat$predaur <- predict(pred_aur, newdat, type = "response") # aur 

e <- ggplot(aes(x = surf_sal*sd+m, color=oxylevel_2.8_9.0), data = newdat, size = 0.9)+
	geom_point(aes(y = n_aur_scaled))+
	geom_line(aes(y = predaur))+
	theme1+
	theme(legend.spacing.y = unit(0.2, "cm"),
				legend.spacing.x = unit(0.0,"cm"),
				legend.title.align = 0.8)+
	scale_color_manual(labels = c("High (>9.3)",
																"Low (<2.8)"),
										 name = expression(paste("Dissolved oxygen",~(mg~L^"-1"), sep="")),
										 values = c("grey10","grey70"))+
	labs(x = "Salinity", 
			 y = expression(paste(italic("A. aurita"),~m^"-3", sep="")))
e

## figures ############# 
#assign a, b, c etc before using

## validation figure
library(patchwork)
a+b+c+d+e+guide_area()+
	plot_layout(nrow = 3, byrow = FALSE, guides = "collect")+
	plot_annotation(tag_levels = 'A')
ggsave("junemodels_sal_al544.png", dpi=1200,dev='png',height=20,width = 16.4,units = 'cm')

## validation figure with region (requires redefining models, after adding update term, to 'a_reg' etc.)
a_reg+b_reg+c_reg+d_reg+e_reg+guide_area()+
  plot_layout(nrow = 3, byrow = FALSE, guides = "collect")+
  plot_annotation(tag_levels = 'A')
ggsave("junemodels_sal_al544_region.png", dpi=1200,dev='png',height=20,width = 16.4,units = 'cm')


## effect plots
a <- plot_model(m_juv_1,
                vline.color = "darkred",
                sort.est = F,
                show.values = TRUE,
                value.offset = .2,
                title = "juv ML full")
b <- plot_model(m_tra_1,
                vline.color = "darkred",
                sort.est = F,
                show.values = TRUE,
                value.offset = .2,
                title = "tra ML full")
c <- plot_model(m_yadu_1,
                vline.color = "darkred",
                sort.est = F,
                show.values = TRUE,
                value.offset = .2,
                title = "yadu ML full")
d <- plot_model(m_adu_1,
                vline.color = "darkred",
                sort.est = F,
                show.values = TRUE,
                value.offset = .2,
                title = "adu ML full")
e <- plot_model(m_aur_1,
                vline.color = "darkred",
                sort.est = F,
                show.values = TRUE,
                value.offset = .2,
                title = "aur ML full")

a+b+c+d+e+
  plot_layout(nrow = 1, byrow = FALSE, guides = "collect")
ggsave("effectplots_2024_june.png", dpi=300,dev='png',height=20,width = 38,units = 'cm')
