setwd("~/Desktop/Garibaldi_Analysis-1")

library(ggplot2)
#library(ggvegan)
library(plyr)
library(dplyr)
library(tidyr)
library(forcats)
library(gridExtra)
library(MCMCglmm)
library(readxl)
library(vegan)
library(cowplot)

##Climate data
clim<-read.table("./Raw data/climate_Iceland.txt", h=T, row.names= NULL)%>% 
          filter(month ==7)%>% 
          mutate(type=case_when(type=="temperature"~"July temperature (°C)",
                                type=="precipitation"~"July precipitation (cm)"))%>% 
          mutate(site=case_when(site=="Thingsvellir"~"Thingvellir",
                       TRUE~site))


clim_plot<-ggplot(data=clim, aes(x=year, y=value, col = site)) +
          geom_line() +
          ylab("") +
          xlab("Year")+
          labs(linetype="Sites",col="Sites")+
          scale_colour_brewer("", palette="Dark2")+
          facet_wrap(scale='free_y',~ type, ncol = 2)+
          theme_bw()+
          theme(strip.background = element_rect(fill = 'white'))


ggsave("clim_plot.tiff", clim_plot, dpi = 300)

##################################################

#Calculate changes in vegetation
dataset <-read_xlsx("/Users/lizmccleary/Desktop/Garibaldi_Analysis-1/merged_2022_2024.xlsx")%>% 
  group_by(SITE,DATE,TRTMT,PLOT,X,Y,SPECIES)%>% 
  mutate(position=case_when(HitOrder=="Bottom" ~ "understory",
                            HitOrder=="Top"~"canopy",
                            is.na(HitOrder)~"middle"))%>% 
  ungroup()

                    
allDataCanopy <- dataset %>%
                  filter(position=="canopy")

allDataUnder<- dataset %>%
  filter(position=="understory")

allDataTopBottom <- dataset %>%
  filter(!position=="middle")

##Top and bottom total species cover
coverTopBottom <- allDataTopBottom%>% 
  select(SITE, TRTMT, PLOT, DATE, SPECIES, HitOrder)%>% 
  mutate(SPECIES=as.factor(SPECIES))%>% 
  group_by(DATE, SITE, TRTMT, PLOT, SPECIES) %>% 
  tally(name = "cover") %>%
  ungroup()

coverCanopy <- allDataCanopy%>% 
  select(SITE, TRTMT, PLOT, DATE, SPECIES, HitOrder)%>% 
  mutate(SPECIES=as.factor(SPECIES))%>% 
  group_by(DATE,SITE,TRTMT,PLOT, SPECIES)%>% 
  tally(name = "cover")%>% 
  mutate(cover = cover/100)%>%
  ungroup()

coverUnder <- allDataUnder%>% 
  select(SITE, TRTMT, PLOT, DATE, SPECIES, HitOrder)%>% 
  mutate(SPECIES=as.factor(SPECIES))%>% 
  group_by(DATE,SITE, TRTMT,PLOT, SPECIES)%>% 
  tally(name = "cover")%>%
  mutate(cover = cover/100)%>%
  ungroup()

#####################NMDS changes in species cover
#Make wide format

SITES<- coverTopBottom %>% 
        filter(SITE %in% c("Salix", "Cassiope","Meadow"))%>% 
        select(DATE,SITE,TRTMT,PLOT, SPECIES, cover)%>% 
        pivot_wider(names_from = SPECIES, values_from = cover)

#Remove columns with all 0 or NA values 
#Swap all NA for zeros
SITES[is.na(SITES)] <- 0


###NMDS
#Do NMDS for all sites
year <- SITES [,1]
SITE <- SITES[,2]
TRTMT <- SITES [,3]
veg<- SITES [,5:49]

veg.mds <- metaMDS(veg, distance = "bray",autotransform = F,k=3,trymax=300)

site.scrs <- as.data.frame(scores(veg.mds, display = "sites"))
site.scrsSITES <- cbind(site.scrs, year, SITE, TRTMT)
site.scrsSITES <- site.scrsSITES[-28,]


#site.scrsSITES$site<- c("Salix", "Cassiope", "Meadow") #figure generated without this

#fit environmental variables with envfit
#env<- cbind(year, TRTMT)
#fit<-envfit(veg.mds, env)

en_coord_cont = as.data.frame(scores(fit, "vectors")) * ordiArrowMul(fit)
en_coord_cat = as.data.frame(scores(fit, "factors")) * ordiArrowMul(fit)

#Plot SITES with envfit

ggplot()+ 
  geom_point(data=site.scrsSITES, aes(NMDS1, NMDS2, colour = as.factor(SITE), shape=factor(TRTMT)),
             size=3) 
#+
  #geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               #data = en_coord_cont, size =1, alpha = 0.5, colour = "grey30") +
  #geom_point(data = en_coord_cat, aes(x = NMDS1, y = NMDS2), 
             #shape = "diamond", size = 4, alpha = 0.6, colour = "navy") +
  #geom_text(data = en_coord_cat, aes(x = NMDS1, y = NMDS2+0.04), 
            #label = row.names(en_coord_cat), colour = "navy", fontface = "bold") + 
  #geom_text(data = en_coord_cont, aes(x = NMDS1, y = NMDS2), colour = "grey30", 
            #fontface = "bold", label = row.names(en_coord_cont))

#Run permanova
veg.dist <- vegdist(veg, method="bray")
adonis2(veg ~ TRTMT*DATE, data = env, permutations = 999, method="bray", by = NULL)

data(dune.env)


#NMDS for CASS

year <- CASS [,1]
TRTMT <- CASS[,3]
veg<- CASS[,5:60]

veg.mds <- metaMDS(veg, distance = "bray",autotransform = F,k=3,trymax=300)

env<- cbind(year, TRTMT)
fit<-envfit(veg.mds, env)

en_coord_cont = as.data.frame(scores(fit, "vectors")) * ordiArrowMul(fit)
en_coord_cat = as.data.frame(scores(fit, "factors")) * ordiArrowMul(fit)

site.scrs <- as.data.frame(scores(veg.mds, display = "sites"))
site.scrsCASS <- cbind(site.scrs, year,TRTMT)
site.scrsCASS$site<-"Cassiope"

#NMDS for MEAD

year <- MEAD [,1]
TRTMT <- MEAD[,3]
veg<- MEAD[,5:60]

veg.mds <- metaMDS(veg, distance = "bray",autotransform = F,k=3,trymax=300)

env<- cbind(year, TRTMT)
fit<-envfit(veg.mds, env)

en_coord_cont = as.data.frame(scores(fit, "vectors")) * ordiArrowMul(fit)
en_coord_cat = as.data.frame(scores(fit, "factors")) * ordiArrowMul(fit)

site.scrs <- as.data.frame(scores(veg.mds, display = "sites"))
site.scrsMEAD <- cbind(site.scrs, year,TRTMT)
site.scrsMEAD$site<-"Meadow"


ggplot()+ 
  geom_point(data=site.scrsAUD, aes(NMDS1, NMDS2, colour = as.factor(YEAR), shape=factor(TRTMT)),
             size=3)+
geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
             data = en_coord_cont, size =1, alpha = 0.5, colour = "grey30") +
  geom_point(data = en_coord_cat, aes(x = NMDS1, y = NMDS2), 
             shape = "diamond", size = 4, alpha = 0.6, colour = "navy") +
  geom_text(data = en_coord_cat, aes(x = NMDS1, y = NMDS2+0.04), 
            label = row.names(en_coord_cat), colour = "navy", fontface = "bold") + 
  geom_text(data = en_coord_cont, aes(x = NMDS1, y = NMDS2), colour = "grey30", 
            fontface = "bold", label = row.names(en_coord_cont))



site.scrsTotal<- rbind(site.scrsSALIX, site.scrsCASS, site.scrsMEAD)%>% 
                  mutate(TRTMT = factor(TRTMT, levels = c("control", "warming"))) 


ggplot()+ 
  geom_point(data=site.scrsTotal, aes(NMDS1, NMDS2, colour = as.factor(YEAR), shape=factor(TRTMT)),
             size=3)+
  facet_grid(rows = vars(site) ,scales="free_y")+
  scale_colour_brewer("", palette="Dark2")+
  theme_classic()+ 
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+
  labs(colour = factor("Year"), shape="Treatment",title ="" )

ggsave("NMDS_topBottom.tiff", last_plot(), dpi = 500)
##############################################
###Changes in GFBROAD
library(RColorBrewer)
library(forcats)

GFBROAD<- allDataTopBottom%>% 
  select(SITE, TRTMT, PLOT, YEAR, GFBROAD, HitOrder)%>% 
  filter(!GFBROAD=="ABIOTIC",
    #     !GFBROAD=="SOIL",
         !GFBROAD=="FUNGI")%>% 
  mutate(GFBROAD=as.factor(GFBROAD),
         YEAR=as.factor(YEAR))%>% 
  group_by(YEAR,SITE, TRTMT,PLOT)%>% 
  count(GFBROAD, .drop = FALSE)%>%
  mutate(cover=n/sum(n)*100)%>%
  ungroup()%>%
  group_by(YEAR,SITE, TRTMT,GFBROAD)%>%
  summarise(coverSITE=mean(cover))%>%
mutate(TRTMT = factor(TRTMT, levels = c("control", "warming")))%>%
  mutate(YEAR = factor(YEAR, levels = c("2022", "2024")))%>%
  mutate(GFBROAD = factor(GFBROAD, levels = c("FORBSV", "GRAMINOID", "ESHRUB", 
                                            "DSHRUB", "MOSS", "LICHEN", 
                                            "LITTER", "SOIL")))%>%
  filter_at(vars(GFBROAD), all_vars(!is.na(.)))


ggsave("Cover_topBottom.tiff", last_plot(), dpi = 300)

GFBROADtop<- allDataCanopy%>% 
  select(SITE, TRTMT, PLOT, YEAR, GFBROAD, HitOrder)%>% 
  filter(!GFBROAD=="ABIOTIC",
         #     !GFBROAD=="SOIL",
         !GFBROAD=="FUNGI")%>% 
  mutate(GFBROAD=as.factor(GFBROAD),
         YEAR=as.factor(YEAR))%>% 
  group_by(YEAR,SITE, TRTMT,PLOT)%>% 
  count(GFBROAD, .drop = FALSE)%>%
  mutate(cover=n/sum(n)*100)%>%
  ungroup()%>%
  group_by(YEAR,SITE, TRTMT,GFBROAD)%>%
  summarise(coverSITE=mean(cover))%>%
  mutate(TRTMT = factor(TRTMT, levels = c("control", "warming")))%>%
  mutate(YEAR = factor(YEAR, levels = c("2022", "2024")))%>%
  mutate(GFBROAD = factor(GFBROAD, levels = c("FORBSV", "GRAMINOID", "ESHRUB", 
                                              "DSHRUB", "MOSS", "LICHEN", 
                                              "LITTER", "SOIL")))
ggsave("Cover_topBottom.tiff", last_plot(), dpi = 500)


#Plot Change in GFBROAD
ggplot(data=GFBROAD,aes(x=YEAR, y=coverSITE, fill=GFBROAD)) + 
  geom_bar(position="stack",stat="identity") +
  facet_grid(TRTMT~SITE ,scales="free_y")+
  scale_fill_brewer("", palette="Dark2")+
  xlab("Sampling year")+ 
  ylab("Cover (%)")+
  theme_bw()+ 
  theme(
    strip.background = element_blank(),
    axis.text.x=element_text(angle = 45, vjust = 0.5),
    axis.title.x=element_blank(),
    axis.ticks.x=element_blank())

ggplot(data=GFBROADtop,aes(x=YEAR, y=coverSITE, fill=GFBROAD)) + 
  geom_bar(position="stack",stat="identity") +
  facet_grid(TRTMT~SITE ,scales="free_y")+
  scale_fill_brewer("", palette="Dark2")+
  xlab("Sampling year")+ 
  ylab("Cover (%)")+
  theme_bw()+ 
  theme(
    strip.background = element_blank(),
    axis.text.x=element_text(angle = 45, vjust = 0.5),
    axis.title.x=element_blank(),
    axis.ticks.x=element_blank())


###########################################
####CHANGE IN HEIGHT

dataset <-read_xlsx("merged_2022_2024.xlsx")

Height_SAL<-dataset%>%
              filter(SITE=="SALIX")%>%
              filter(Hit==1)%>%
              select(YEAR,TRTMT, PLOT, CanopyHeight.mm.)%>%
              mutate(CanopyHeight.mm.=as.numeric(CanopyHeight.mm.))%>%
              na.omit()%>%
              group_by(YEAR,TRTMT, PLOT)%>%
              summarise(mean_height=mean(CanopyHeight.mm.))%>%
              mutate(TRTMT = factor(TRTMT, levels = c("control", "warming")))

prior2 <- list(R = list(V = 1, nu = 0.002), 
               G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.v = 10000), 
                        G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.v = 10000)))


SALIX_Height_CTL <- MCMCglmm(mean_height ~ I(YEAR-2022), 
                             random = ~ YEAR+PLOT, data = Height_SAL[Height_SAL$TRTMT == "control",], 
                             family = "gaussian", pr = TRUE, nitt = 100000, 
                             burnin = 20000, prior = prior2)

SALIX_Height_OTC <- MCMCglmm(mean_height ~ I(YEAR-2022), 
                           random = ~ YEAR+PLOT, data = Height_SAL[Height_SAL$TRTMT == "warming",], 
                           family = "gaussian", pr = TRUE, nitt = 100000, 
                           burnin = 20000, prior = prior2)

#Test difference between treatment
canopy_m <- MCMCglmm(mean_height ~ I(YEAR-2022)+TRTMT-1, random = ~ YEAR + PLOT,
                            data = Height_SAL, 
                            family = "gaussian", pr = TRUE, nitt = 100000, 
                            burnin = 20000, prior = prior2)



# Calculating model predictions

#CTL
nyears <- 23
niter <- length(SALIX_Height_CTL$Sol[,"(Intercept)"])

SALIX_Height_CTLpreds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    SALIX_Height_CTLpreds[i,j] <- SALIX_Height_CTL$Sol[i,"(Intercept)"] + SALIX_Height_CTL$Sol[i,"I(YEAR - 2022)"]*j
  }
}


SALIX_Height_CTLpreds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  SALIX_Height_CTLpreds_df [i,] <- quantile(SALIX_Height_CTLpreds[,i], c(0.025, 0.5, 0.975))
}

SALIX_Height_CTLpreds_df <- cbind.data.frame(lower = SALIX_Height_CTLpreds_df[,1], 
                                               mean = SALIX_Height_CTLpreds_df[,2], 
                                               upper = SALIX_Height_CTLpreds_df[,3], year = seq(1:23))
#OTC
nyears <- 23
niter <- length(SALIX_Height_OTC$Sol[,"(Intercept)"])

SALIX_Height_OTCpreds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    SALIX_Height_OTCpreds[i,j] <- SALIX_Height_OTC$Sol[i,"(Intercept)"] + SALIX_Height_OTC$Sol[i,"I(YEAR - 2022)"]*j
  }
}


SALIX_Height_OTCpreds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  SALIX_Height_OTCpreds_df [i,] <- quantile(SALIX_Height_OTCpreds[,i], c(0.025, 0.5, 0.975))
}

SALIX_Height_OTCpreds_df <- cbind.data.frame(lower = SALIX_Height_OTCpreds_df[,1], 
                                               mean = SALIX_Height_OTCpreds_df[,2], 
                                               upper = SALIX_Height_OTCpreds_df[,3], year = seq(1:23))


###Plot change in height
ggplot() +
  geom_point(data= Height_SAL, 
             aes(x = YEAR, y = mean_height, color=TRTMT),
             alpha = 0.8, size = 2)+
  scale_color_manual(values = c("#1B9E77", "#7570B3", "#E6AB02"), labels = c("control","warming")) +
  scale_fill_manual(values = c("#1B9E77", "#7570B3", "#E6AB02"),  labels =c("control","warming")) +
  scale_x_continuous(breaks = c(1997, 2000, 2007, 2014, 2019)) +
  geom_ribbon(data = SALIX_Height_CTLpreds_df, aes(x = year + 1996, ymin = lower, ymax = upper), fill = "#1B9E77", alpha = 0.2) +
  geom_line(data = SALIX_Height_CTLpreds_df, aes(x = year + 1996, y = mean), colour = "#1B9E77") +
  geom_ribbon(data = SALIX_Height_OTCpreds_df, aes(x = year + 1996, ymin = lower, ymax = upper), fill = "#7570B3", alpha = 0.2) +
  geom_line(data = SALIX_Height_OTCpreds_df, aes(x = year + 1996, y = mean), colour = "#7570B3",  linetype="dashed") +
  theme_bw()+
  xlab("Year")+
  ylab("Mean height (cm)")+
  theme(
    #axis.text.x=element_blank(),
     #   axis.title.x=element_blank(),
      #  axis.title.y=element_blank(),
       # axis.ticks.x =element_blank(),
        text = element_text(size=9),
        #legend.position="none",
        plot.margin = unit(c(0,0.5,0,0.5), "cm"))

ggsave("Height_regression.tiff", last_plot(), dpi = 500)

####Alternative plot for change in height: BOXPLOT
ggplot(data=Height_SAL, aes(x=factor(YEAR), y=mean_height,fill = factor(TRTMT))) + 
  geom_boxplot() +
  scale_fill_brewer("Treatment", palette="Dark2")+
  xlab("Sampling year")+ 
  ylab("Height (mm)")+
  theme_bw()+
  theme(
    strip.background = element_blank(),
    axis.title.x=element_blank(),
    axis.ticks.x=element_blank(),
    legend.position = "right",
    axis.text.x = element_text(angle = 45, vjust = 0.5))

ggsave("Height_boxplot.tiff", last_plot(), dpi = 300)

################################################
#Species specific changes
###Get cover data from line 34-73

prior2 <- list(R = list(V = 1, nu = 0.002), 
               G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.v = 10000), 
                        G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.v = 10000)))

### Linear model for all species AUD Top and Bottom
#CTL

SALIX <- coverTopBottom%>%
  filter(SITE=="Salix")%>%
#  mutate(cover=cover*100)%>%
  na.omit()


AUD_plot_all_CTL <- MCMCglmm(cover ~ I(YEAR-1996) * SPP, 
                             random = ~ YEAR+PLOT, data = AUD[AUD$TRTMT == "CTL",], 
                             family = "gaussian", pr = TRUE, nitt = 100000, 
                             burnin = 20000, prior = prior2)
summary(AUD_plot_all_CTL)

#OTC
SAL_plot_all_OTC <- MCMCglmm(cover ~ I(YEAR-1996) * SPP, 
                             random = ~ YEAR+PLOT, data = SALIX[SALIX$TRTMT == "warming",], 
                             family = "gaussian", pr = TRUE, nitt = 100000, 
                             burnin = 20000, prior = prior2)

summary(SAL_plot_all_OTC)

#Linear model for main species to check for difference in treatment
#IS THIS THE SPOT I WANT TO GROUP SPECIES BY FAMILY??

SAL_plot_CAR <- MCMCglmm(cover ~ I(YEAR-2022)+TRTMT-1, 
                            random = ~ YEAR + PLOT,
                            data = SALIX[SALIX$SPP == "car",], 
                            family = "gaussian", pr = TRUE, nitt = 100000, 
                            burnin = 20000, prior = prior2)
summary(SAL_plot_CAR)

SAL_plot_JUN <- MCMCglmm(cover ~ I(YEAR-2022)+TRTMT-1, 
                            random = ~ YEAR + PLOT,
                            data = SALIX[SALIX$SPP == "jun",], 
                            family = "gaussian", pr = TRUE, nitt = 100000, 
                            burnin = 20000, prior = prior2)
summary(SAL_plot_JUN)

AUD_plot_CETISL <- MCMCglmm(cover ~ I(YEAR-1996)+TRTMT-1, 
                            random = ~ YEAR + PLOT,
                            data = AUD[AUD$SPP == "CETISL",], 
                            family = "gaussian", pr = TRUE, nitt = 100000, 
                            burnin = 20000, prior = prior2)
summary(AUD_plot_CETISL)

AUD_plot_BISVIV <- MCMCglmm(cover ~ I(YEAR-1996)+TRTMT-1, 
                            random = ~ YEAR + PLOT,
                            data = AUD[AUD$SPP == "BISVIV",], 
                            family = "gaussian", pr = TRUE, nitt = 100000, 
                            burnin = 20000, prior = prior2)
summary(AUD_plot_BISVIV)

##Linear model for main species to plot model predictions
AUD_res1 <- subset(AUD, SPP =="BETNAN")
AUD_res2 <- subset(AUD, SPP =="RACLAN")
AUD_res3 <- subset(AUD, SPP =="CETISL")
AUD_res4 <- subset(AUD, SPP =="BISVIV")
AUD_res5 <- subset(AUD, SPP =="LITTER")
AUD_res6 <- subset(AUD, SPP =="SANUNC")
#BETNAN

AUD_plot_BETNAN_CTL <- MCMCglmm(cover ~ I(YEAR-1996), random = ~ YEAR + PLOT,
                                data = AUD_res1[AUD_res1$TRTMT == "CTL",], 
                                family = "gaussian", pr = TRUE, nitt = 100000, 
                                burnin = 20000, prior = prior2)
summary(AUD_plot_BETNAN_CTL)

AUD_plot_BETNAN_OTC <- MCMCglmm(cover ~ I(YEAR-1996), random = ~ YEAR + PLOT,
                                data = AUD_res1[AUD_res1$TRTMT == "OTC",], 
                                family = "gaussian", pr = TRUE, nitt = 100000, 
                                burnin = 20000, prior = prior2)
summary(AUD_plot_BETNAN_OTC)

AUD_plot_BETNAN_GRA <- MCMCglmm(cover ~ I(YEAR-1996), random = ~ YEAR + PLOT,
                                data = AUD_res1[AUD_res1$TRTMT == "GRA",], 
                                family = "gaussian", pr = TRUE, nitt = 100000, 
                                burnin = 20000, prior = prior2)
summary(AUD_plot_BETNAN_GRA)

#RACLAN
AUD_plot_RACLAN_CTL <- MCMCglmm(cover ~ I(YEAR-1996), random = ~ YEAR + PLOT,
                                data = AUD_res2[AUD_res2$TRTMT == "CTL",], 
                                family = "gaussian", pr = TRUE, nitt = 100000, 
                                burnin = 20000, prior = prior2)
summary(AUD_plot_RACLAN_CTL)

AUD_plot_RACLAN_OTC <- MCMCglmm(cover ~ I(YEAR-1996), random = ~ YEAR + PLOT,
                                data = AUD_res2[AUD_res2$TRTMT == "OTC",], 
                                family = "gaussian", pr = TRUE, nitt = 100000, 
                                burnin = 20000, prior = prior2)
summary(AUD_plot_RACLAN_OTC)

AUD_plot_RACLAN_GRA <- MCMCglmm(cover ~ I(YEAR-1996), random = ~ YEAR + PLOT,
                                data = AUD_res2[AUD_res2$TRTMT == "GRA",], 
                                family = "gaussian", pr = TRUE, nitt = 100000, 
                                burnin = 20000, prior = prior2)
summary(AUD_plot_RACLAN_GRA)

#CETISL
AUD_plot_CETISL_CTL <- MCMCglmm(cover ~ I(YEAR-1996), random = ~ YEAR + PLOT,
                                data = AUD_res3[AUD_res3$TRTMT == "CTL",], 
                                family = "gaussian", pr = TRUE, nitt = 100000, 
                                burnin = 20000, prior = prior2)
summary(AUD_plot_CETISL_CTL)

AUD_plot_CETISL_OTC <- MCMCglmm(cover ~ I(YEAR-1996), random = ~ YEAR + PLOT,
                                data = AUD_res3[AUD_res3$TRTMT == "OTC",], 
                                family = "gaussian", pr = TRUE, nitt = 100000, 
                                burnin = 20000, prior = prior2)
summary(AUD_plot_CETISL_OTC)

AUD_plot_CETISL_GRA <- MCMCglmm(cover ~ I(YEAR-1996), random = ~ YEAR + PLOT,
                                data = AUD_res3[AUD_res3$TRTMT == "GRA",], 
                                family = "gaussian", pr = TRUE, nitt = 100000, 
                                burnin = 20000, prior = prior2)
summary(AUD_plot_CETISL_GRA)

#BISVIV
AUD_plot_BISVIV_CTL <- MCMCglmm(cover ~ I(YEAR-1996), random = ~ YEAR + PLOT,
                                data = AUD_res4[AUD_res4$TRTMT == "CTL",], 
                                family = "gaussian", pr = TRUE, nitt = 100000, 
                                burnin = 20000, prior = prior2)
summary(AUD_plot_BISVIV_CTL)

AUD_plot_BISVIV_OTC <- MCMCglmm(cover ~ I(YEAR-1996), random = ~ YEAR + PLOT,
                                data = AUD_res4[AUD_res4$TRTMT == "OTC",], 
                                family = "gaussian", pr = TRUE, nitt = 100000, 
                                burnin = 20000, prior = prior2)
summary(AUD_plot_BISVIV_OTC)

AUD_plot_BISVIV_GRA <- MCMCglmm(cover ~ I(YEAR-1996), random = ~ YEAR + PLOT,
                                data = AUD_res4[AUD_res4$TRTMT == "GRA",], 
                                family = "gaussian", pr = TRUE, nitt = 100000, 
                                burnin = 20000, prior = prior2)
summary(AUD_plot_BISVIV_GRA)

#LITTER
AUD_plot_LITTER_CTL <- MCMCglmm(cover ~ I(YEAR-1996), random = ~ YEAR + PLOT,
                                data = AUD_res5[AUD_res5$TRTMT == "CTL",], 
                                family = "gaussian", pr = TRUE, nitt = 100000, 
                                burnin = 20000, prior = prior2)
summary(AUD_plot_LITTER_CTL)

AUD_plot_LITTER_OTC <- MCMCglmm(cover ~ I(YEAR-1996), random = ~ YEAR + PLOT,
                                data = AUD_res5[AUD_res5$TRTMT == "OTC",], 
                                family = "gaussian", pr = TRUE, nitt = 100000, 
                                burnin = 20000, prior = prior2)
summary(AUD_plot_LITTER_OTC)

AUD_plot_LITTER_GRA <- MCMCglmm(cover ~ I(YEAR-1996), random = ~ YEAR + PLOT,
                                data = AUD_res5[AUD_res5$TRTMT == "GRA",], 
                                family = "gaussian", pr = TRUE, nitt = 100000, 
                                burnin = 20000, prior = prior2)
summary(AUD_plot_LITTER_GRA)

#SANUNC
AUD_plot_SANUNC_CTL <- MCMCglmm(cover ~ I(YEAR-1996), random = ~ YEAR + PLOT,
                                data = AUD_res6[AUD_res6$TRTMT == "CTL",], 
                                family = "gaussian", pr = TRUE, nitt = 100000, 
                                burnin = 20000, prior = prior2)
summary(AUD_plot_SANUNC_CTL)

AUD_plot_SANUNC_OTC <- MCMCglmm(cover ~ I(YEAR-1996), random = ~ YEAR + PLOT,
                                data = AUD_res6[AUD_res6$TRTMT == "OTC",], 
                                family = "gaussian", pr = TRUE, nitt = 100000, 
                                burnin = 20000, prior = prior2)
summary(AUD_plot_SANUNC_OTC)

AUD_plot_SANUNC_GRA <- MCMCglmm(cover ~ I(YEAR-1996), random = ~ YEAR + PLOT,
                                data = AUD_res6[AUD_res6$TRTMT == "GRA",], 
                                family = "gaussian", pr = TRUE, nitt = 100000, 
                                burnin = 20000, prior = prior2)
summary(AUD_plot_SANUNC_GRA)

# Calculating model predictions
#BETNAN

#CTL
nyears <- 23
niter <- length(AUD_plot_BETNAN_CTL$Sol[,"(Intercept)"])

AUD_plot_BETNANCTLpreds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    AUD_plot_BETNANCTLpreds[i,j] <- AUD_plot_BETNAN_CTL$Sol[i,"(Intercept)"] + AUD_plot_BETNAN_CTL$Sol[i,"I(YEAR - 1996)"]*j
  }
}

AUD_plot_BETNANCTLpreds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  AUD_plot_BETNANCTLpreds_df [i,] <- quantile(AUD_plot_BETNANCTLpreds[,i], c(0.025, 0.5, 0.975))
}

AUD_plot_BETNANCTLpreds_df <- cbind.data.frame(lower = AUD_plot_BETNANCTLpreds_df[,1], 
                                               mean = AUD_plot_BETNANCTLpreds_df[,2], 
                                               upper = AUD_plot_BETNANCTLpreds_df[,3], year = seq(1:23))
#OTC
nyears <- 23
niter <- length(AUD_plot_BETNAN_OTC$Sol[,"(Intercept)"])

AUD_plot_BETNANOTCpreds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    AUD_plot_BETNANOTCpreds[i,j] <- AUD_plot_BETNAN_OTC$Sol[i,"(Intercept)"] + AUD_plot_BETNAN_OTC$Sol[i,"I(YEAR - 1996)"]*j
  }
}

AUD_plot_BETNANOTCpreds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  AUD_plot_BETNANOTCpreds_df [i,] <- quantile(AUD_plot_BETNANOTCpreds[,i], c(0.025, 0.5, 0.975))
}

AUD_plot_BETNANOTCpreds_df <- cbind.data.frame(lower = AUD_plot_BETNANOTCpreds_df[,1], 
                                               mean = AUD_plot_BETNANOTCpreds_df[,2], 
                                               upper = AUD_plot_BETNANOTCpreds_df[,3], year = seq(1:23))

#GRA
nyears <- 23
niter <- length(AUD_plot_BETNAN_GRA$Sol[,"(Intercept)"])

AUD_plot_BETNANGRApreds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    AUD_plot_BETNANGRApreds[i,j] <- AUD_plot_BETNAN_GRA$Sol[i,"(Intercept)"] + AUD_plot_BETNAN_GRA$Sol[i,"I(YEAR - 1996)"]*j
  }
}


AUD_plot_BETNANGRApreds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  AUD_plot_BETNANGRApreds_df [i,] <- quantile(AUD_plot_BETNANGRApreds[,i], c(0.025, 0.5, 0.975))
}

AUD_plot_BETNANGRApreds_df <- cbind.data.frame(lower = AUD_plot_BETNANGRApreds_df[,1], 
                                               mean = AUD_plot_BETNANGRApreds_df[,2], 
                                               upper = AUD_plot_BETNANGRApreds_df[,3], year = seq(1:23))

#Making graph
AUD_BETNAN<- 
  ggplot() +
  geom_point(data= AUD_res1, 
             aes(x = YEAR, y = cover, color=TRTMT),
             alpha = 0.8, size = 2)+
  scale_color_manual(values = c("#1B9E77", "#7570B3", "#E6AB02"), labels = c("CTL","GRA", "OTC")) +
  scale_fill_manual(values = c("#1B9E77", "#7570B3", "#E6AB02"),  labels =c("CTL","GRA", "OTC")) +
  scale_x_continuous(breaks = c(1997, 2000, 2007, 2014, 2019)) +
 # geom_ribbon(data = AUD_plot_BETNANCTLpreds_df, aes(x = year + 1996, ymin = lower, ymax = upper), fill = "#1B9E77", alpha = 0.2) +
 # geom_line(data = AUD_plot_BETNANCTLpreds_df, aes(x = year + 1996, y = mean), colour = "#1B9E77") +
  geom_ribbon(data = AUD_plot_BETNANGRApreds_df, aes(x = year + 1996, ymin = lower, ymax = upper), fill = "#7570B3", alpha = 0.2) +
  geom_line(data = AUD_plot_BETNANGRApreds_df, aes(x = year + 1996, y = mean), colour = "#7570B3") +
  geom_ribbon(data = AUD_plot_BETNANOTCpreds_df, aes(x = year + 1996, ymin = lower, ymax = upper), fill = "#E6AB02", alpha = 0.2) +
  geom_line(data = AUD_plot_BETNANOTCpreds_df, aes(x = year + 1996, y = mean), colour = "#E6AB02",  linetype="dashed") +
  theme_bw()+
  ylab("Number of hits")+
  labs(title="",subtitle=expression(italic("Betula nana")), size=9)+
  theme(axis.text.x=element_blank(),
        axis.title.x=element_blank(),
      #  axis.title.y=element_blank(),
        axis.ticks.x =element_blank(),
        text = element_text(size=9),
        legend.position="none",
        plot.margin = unit(c(0,0.5,0,0.5), "cm"))

#RACLAN
#CTL
nyears <- 23
niter <- length(AUD_plot_RACLAN_CTL$Sol[,"(Intercept)"])

AUD_plot_RACLANCTLpreds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    AUD_plot_RACLANCTLpreds[i,j] <- AUD_plot_RACLAN_CTL$Sol[i,"(Intercept)"] + 
      AUD_plot_RACLAN_CTL$Sol[i,"I(YEAR - 1996)"]*j
  }
}


AUD_plot_RACLANCTLpreds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  AUD_plot_RACLANCTLpreds_df [i,] <- quantile(AUD_plot_RACLANCTLpreds[,i], c(0.025, 0.5, 0.975))
}

AUD_plot_RACLANCTLpreds_df <- cbind.data.frame(lower = AUD_plot_RACLANCTLpreds_df[,1], 
                                               mean = AUD_plot_RACLANCTLpreds_df[,2], 
                                               upper = AUD_plot_RACLANCTLpreds_df[,3], year = seq(1:23))
#OTC
nyears <- 23
niter <- length(AUD_plot_RACLAN_OTC$Sol[,"(Intercept)"])

AUD_plot_RACLANOTCpreds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    AUD_plot_RACLANOTCpreds[i,j] <- AUD_plot_RACLAN_OTC$Sol[i,"(Intercept)"] + 
      AUD_plot_RACLAN_OTC$Sol[i,"I(YEAR - 1996)"]*j
  }
}


AUD_plot_RACLANOTCpreds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  AUD_plot_RACLANOTCpreds_df [i,] <- quantile(AUD_plot_RACLANOTCpreds[,i], c(0.025, 0.5, 0.975))
}

AUD_plot_RACLANOTCpreds_df <- cbind.data.frame(lower = AUD_plot_RACLANOTCpreds_df[,1], 
                                               mean = AUD_plot_RACLANOTCpreds_df[,2], 
                                               upper = AUD_plot_RACLANOTCpreds_df[,3], year = seq(1:23))

#GRA
nyears <- 23
niter <- length(AUD_plot_RACLAN_GRA$Sol[,"(Intercept)"])

AUD_plot_RACLANGRApreds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    AUD_plot_RACLANGRApreds[i,j] <- AUD_plot_RACLAN_GRA$Sol[i,"(Intercept)"] + 
      AUD_plot_RACLAN_GRA$Sol[i,"I(YEAR - 1996)"]*j
  }
}

AUD_plot_RACLANGRApreds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  AUD_plot_RACLANGRApreds_df [i,] <- quantile(AUD_plot_RACLANGRApreds[,i], c(0.025, 0.5, 0.975))
}

AUD_plot_RACLANGRApreds_df <- cbind.data.frame(lower = AUD_plot_RACLANGRApreds_df[,1], 
                                               mean = AUD_plot_RACLANGRApreds_df[,2], 
                                               upper = AUD_plot_RACLANGRApreds_df[,3], year = seq(1:23))
#Making graph
AUD_RACLAN<- 
  ggplot() +
  geom_point(data= AUD_res2, 
             aes(x = YEAR, y = cover, color=TRTMT),
             alpha = 0.8, size = 2)+
  scale_color_manual(values = c("#1B9E77", "#7570B3", "#E6AB02"), labels = c("CTL","GRA", "OTC")) +
  scale_fill_manual(values = c("#1B9E77", "#7570B3", "#E6AB02"),  labels =c("CTL","GRA", "OTC")) +
  scale_x_continuous(breaks = c(1997, 2000, 2007, 2014, 2019)) +
  geom_ribbon(data = AUD_plot_RACLANCTLpreds_df, aes(x = year + 1996, ymin = lower, ymax = upper), fill = "#1B9E77", alpha = 0.2) +
  geom_line(data = AUD_plot_RACLANCTLpreds_df, aes(x = year + 1996, y = mean), colour = "#1B9E77") +
  geom_ribbon(data = AUD_plot_RACLANGRApreds_df, aes(x = year + 1996, ymin = lower, ymax = upper), fill = "#7570B3", alpha = 0.2) +
  geom_line(data = AUD_plot_RACLANGRApreds_df, aes(x = year + 1996, y = mean), colour = "#7570B3") +
  geom_ribbon(data = AUD_plot_RACLANOTCpreds_df, aes(x = year + 1996, ymin = lower, ymax = upper), fill = "#E6AB02", alpha = 0.2) +
  geom_line(data = AUD_plot_RACLANOTCpreds_df, aes(x = year + 1996, y = mean), colour = "#E6AB02") +
  ylab("Number of hits")+
  theme_bw()+
  labs(title= "Audkuluheidi", subtitle=expression(italic("Racomitrium lanuginosum")))+
  theme(axis.text.x=element_blank(),
        axis.title.x=element_blank(),
     #   axis.title.y=element_blank(),
        axis.ticks.x =element_blank(),
        legend.position = "none",
        text = element_text(size=9),
        plot.title=element_text(size=12),
        plot.margin = unit(c(0.5,0.5,0,0.5), "cm"))

#CETISL
#CTL
nyears <- 23
niter <- length(AUD_plot_CETISL_CTL$Sol[,"(Intercept)"])

AUD_plot_CETISLCTLpreds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    AUD_plot_CETISLCTLpreds[i,j] <- AUD_plot_CETISL_CTL$Sol[i,"(Intercept)"] + 
      AUD_plot_CETISL_CTL$Sol[i,"I(YEAR - 1996)"]*j
  }
}


AUD_plot_CETISLCTLpreds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  AUD_plot_CETISLCTLpreds_df [i,] <- quantile(AUD_plot_CETISLCTLpreds[,i], c(0.025, 0.5, 0.975))
}

AUD_plot_CETISLCTLpreds_df <- cbind.data.frame(lower = AUD_plot_CETISLCTLpreds_df[,1], 
                                               mean = AUD_plot_CETISLCTLpreds_df[,2], 
                                               upper = AUD_plot_CETISLCTLpreds_df[,3], year = seq(1:23))
#OTC
nyears <- 23
niter <- length(AUD_plot_CETISL_OTC$Sol[,"(Intercept)"])

AUD_plot_CETISLOTCpreds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    AUD_plot_CETISLOTCpreds[i,j] <- AUD_plot_CETISL_OTC$Sol[i,"(Intercept)"] + 
      AUD_plot_CETISL_OTC$Sol[i,"I(YEAR - 1996)"]*j
  }
}


AUD_plot_CETISLOTCpreds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  AUD_plot_CETISLOTCpreds_df [i,] <- quantile(AUD_plot_CETISLOTCpreds[,i], c(0.025, 0.5, 0.975))
}

AUD_plot_CETISLOTCpreds_df <- cbind.data.frame(lower = AUD_plot_CETISLOTCpreds_df[,1], 
                                               mean = AUD_plot_CETISLOTCpreds_df[,2], 
                                               upper = AUD_plot_CETISLOTCpreds_df[,3], year = seq(1:23))

#GRA
nyears <- 23
niter <- length(AUD_plot_CETISL_GRA$Sol[,"(Intercept)"])

AUD_plot_CETISLGRApreds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    AUD_plot_CETISLGRApreds[i,j] <- AUD_plot_CETISL_GRA$Sol[i,"(Intercept)"] + 
      AUD_plot_CETISL_GRA$Sol[i,"I(YEAR - 1996)"]*j
  }
}

AUD_plot_CETISLGRApreds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  AUD_plot_CETISLGRApreds_df [i,] <- quantile(AUD_plot_CETISLGRApreds[,i], c(0.025, 0.5, 0.975))
}

AUD_plot_CETISLGRApreds_df <- cbind.data.frame(lower = AUD_plot_CETISLGRApreds_df[,1], 
                                               mean = AUD_plot_CETISLGRApreds_df[,2], 
                                               upper = AUD_plot_CETISLGRApreds_df[,3], year = seq(1:23))
#Making graph
AUD_CETISL<- 
  ggplot() +
  geom_point(data= AUD_res3, 
             aes(x = YEAR, y = cover, color=TRTMT),
             alpha = 0.8, size = 2)+
  scale_color_manual(values = c("#1B9E77", "#7570B3", "#E6AB02"), labels = c("CTL","GRA", "OTC")) +
  scale_fill_manual(values = c("#1B9E77", "#7570B3", "#E6AB02"),  labels =c("CTL","GRA", "OTC")) +
  scale_x_continuous(breaks = c(1997, 2000, 2007, 2014, 2019)) +
#  geom_ribbon(data = AUD_plot_CETISLCTLpreds_df, aes(x = year + 1996, ymin = lower, ymax = upper), fill = "#1B9E77", alpha = 0.2) +
#  geom_line(data = AUD_plot_CETISLCTLpreds_df, aes(x = year + 1996, y = mean), colour = "#1B9E77") +
 # geom_ribbon(data = AUD_plot_CETISLGRApreds_df, aes(x = year + 1996, ymin = lower, ymax = upper), fill = "#7570B3", alpha = 0.2) +
#  geom_line(data = AUD_plot_CETISLGRApreds_df, aes(x = year + 1996, y = mean), colour = "#7570B3") +
  geom_ribbon(data = AUD_plot_CETISLOTCpreds_df, aes(x = year + 1996, ymin = lower, ymax = upper), fill = "#E6AB02", alpha = 0.2) +
  geom_line(data = AUD_plot_CETISLOTCpreds_df, aes(x = year + 1996, y = mean), colour = "#E6AB02") +
  ylab("Number of hits")+
  theme_bw()+
  labs(title= "", subtitle=expression(italic("CETISL")))+
  theme(axis.text.x=element_blank(),
        axis.title.x=element_blank(),
   #     axis.title.y=element_blank(),
        axis.ticks.x =element_blank(),
        legend.position = "none",
        text = element_text(size=9),
        plot.title=element_text(size=12),
        plot.margin = unit(c(0.5,0.5,0,0.5), "cm"))

#BISVIV
#CTL
nyears <- 23
niter <- length(AUD_plot_BISVIV_CTL$Sol[,"(Intercept)"])

AUD_plot_BISVIVCTLpreds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    AUD_plot_BISVIVCTLpreds[i,j] <- AUD_plot_BISVIV_CTL$Sol[i,"(Intercept)"] + 
      AUD_plot_BISVIV_CTL$Sol[i,"I(YEAR - 1996)"]*j
  }
}


AUD_plot_BISVIVCTLpreds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  AUD_plot_BISVIVCTLpreds_df [i,] <- quantile(AUD_plot_BISVIVCTLpreds[,i], c(0.025, 0.5, 0.975))
}

AUD_plot_BISVIVCTLpreds_df <- cbind.data.frame(lower = AUD_plot_BISVIVCTLpreds_df[,1], 
                                               mean = AUD_plot_BISVIVCTLpreds_df[,2], 
                                               upper = AUD_plot_BISVIVCTLpreds_df[,3], year = seq(1:23))
#OTC
nyears <- 23
niter <- length(AUD_plot_BISVIV_OTC$Sol[,"(Intercept)"])

AUD_plot_BISVIVOTCpreds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    AUD_plot_BISVIVOTCpreds[i,j] <- AUD_plot_BISVIV_OTC$Sol[i,"(Intercept)"] + 
      AUD_plot_BISVIV_OTC$Sol[i,"I(YEAR - 1996)"]*j
  }
}


AUD_plot_BISVIVOTCpreds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  AUD_plot_BISVIVOTCpreds_df [i,] <- quantile(AUD_plot_BISVIVOTCpreds[,i], c(0.025, 0.5, 0.975))
}

AUD_plot_BISVIVOTCpreds_df <- cbind.data.frame(lower = AUD_plot_BISVIVOTCpreds_df[,1], 
                                               mean = AUD_plot_BISVIVOTCpreds_df[,2], 
                                               upper = AUD_plot_BISVIVOTCpreds_df[,3], year = seq(1:23))

#GRA
nyears <- 23
niter <- length(AUD_plot_BISVIV_GRA$Sol[,"(Intercept)"])

AUD_plot_BISVIVGRApreds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    AUD_plot_BISVIVGRApreds[i,j] <- AUD_plot_BISVIV_GRA$Sol[i,"(Intercept)"] + 
      AUD_plot_BISVIV_GRA$Sol[i,"I(YEAR - 1996)"]*j
  }
}

AUD_plot_BISVIVGRApreds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  AUD_plot_BISVIVGRApreds_df [i,] <- quantile(AUD_plot_BISVIVGRApreds[,i], c(0.025, 0.5, 0.975))
}

AUD_plot_BISVIVGRApreds_df <- cbind.data.frame(lower = AUD_plot_BISVIVGRApreds_df[,1], 
                                               mean = AUD_plot_BISVIVGRApreds_df[,2], 
                                               upper = AUD_plot_BISVIVGRApreds_df[,3], year = seq(1:23))
#Making graph
AUD_BISVIV<- 
  ggplot() +
  geom_point(data= AUD_res3, 
             aes(x = YEAR, y = cover, color=TRTMT),
             alpha = 0.8, size = 2)+
  scale_color_manual(values = c("#1B9E77", "#7570B3", "#E6AB02"), labels = c("CTL","GRA", "OTC")) +
  scale_fill_manual(values = c("#1B9E77", "#7570B3", "#E6AB02"),  labels =c("CTL","GRA", "OTC")) +
  scale_x_continuous(breaks = c(1997, 2000, 2007, 2014, 2019)) +
  #  geom_ribbon(data = AUD_plot_BISVIVCTLpreds_df, aes(x = year + 1996, ymin = lower, ymax = upper), fill = "#1B9E77", alpha = 0.2) +
  #  geom_line(data = AUD_plot_BISVIVCTLpreds_df, aes(x = year + 1996, y = mean), colour = "#1B9E77") +
  # geom_ribbon(data = AUD_plot_BISVIVGRApreds_df, aes(x = year + 1996, ymin = lower, ymax = upper), fill = "#7570B3", alpha = 0.2) +
  #  geom_line(data = AUD_plot_BISVIVGRApreds_df, aes(x = year + 1996, y = mean), colour = "#7570B3") +
  geom_ribbon(data = AUD_plot_BISVIVOTCpreds_df, aes(x = year + 1996, ymin = lower, ymax = upper), fill = "#E6AB02", alpha = 0.2) +
  geom_line(data = AUD_plot_BISVIVOTCpreds_df, aes(x = year + 1996, y = mean), colour = "#E6AB02") +
  ylab("Number of hits")+
  theme_bw()+
  labs(title= "", subtitle=expression(italic("BISVIV")))+
  theme(axis.text.x=element_blank(),
        axis.title.x=element_blank(),
     #   axis.title.y=element_blank(),
        axis.ticks.x =element_blank(),
        legend.position = "none",
        text = element_text(size=9),
        plot.title=element_text(size=12),
        plot.margin = unit(c(0.5,0.5,0,0.5), "cm"))

#LITTER
#CTL
nyears <- 23
niter <- length(AUD_plot_LITTER_CTL$Sol[,"(Intercept)"])

AUD_plot_LITTERCTLpreds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    AUD_plot_LITTERCTLpreds[i,j] <- AUD_plot_LITTER_CTL$Sol[i,"(Intercept)"] + 
      AUD_plot_LITTER_CTL$Sol[i,"I(YEAR - 1996)"]*j
  }
}


AUD_plot_LITTERCTLpreds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  AUD_plot_LITTERCTLpreds_df [i,] <- quantile(AUD_plot_LITTERCTLpreds[,i], c(0.025, 0.5, 0.975))
}

AUD_plot_LITTERCTLpreds_df <- cbind.data.frame(lower = AUD_plot_LITTERCTLpreds_df[,1], 
                                               mean = AUD_plot_LITTERCTLpreds_df[,2], 
                                               upper = AUD_plot_LITTERCTLpreds_df[,3], year = seq(1:23))
#OTC
nyears <- 23
niter <- length(AUD_plot_LITTER_OTC$Sol[,"(Intercept)"])

AUD_plot_LITTEROTCpreds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    AUD_plot_LITTEROTCpreds[i,j] <- AUD_plot_LITTER_OTC$Sol[i,"(Intercept)"] + 
      AUD_plot_LITTER_OTC$Sol[i,"I(YEAR - 1996)"]*j
  }
}


AUD_plot_LITTEROTCpreds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  AUD_plot_LITTEROTCpreds_df [i,] <- quantile(AUD_plot_LITTEROTCpreds[,i], c(0.025, 0.5, 0.975))
}

AUD_plot_LITTEROTCpreds_df <- cbind.data.frame(lower = AUD_plot_LITTEROTCpreds_df[,1], 
                                               mean = AUD_plot_LITTEROTCpreds_df[,2], 
                                               upper = AUD_plot_LITTEROTCpreds_df[,3], year = seq(1:23))

#GRA
nyears <- 23
niter <- length(AUD_plot_LITTER_GRA$Sol[,"(Intercept)"])

AUD_plot_LITTERGRApreds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    AUD_plot_LITTERGRApreds[i,j] <- AUD_plot_LITTER_GRA$Sol[i,"(Intercept)"] + 
      AUD_plot_LITTER_GRA$Sol[i,"I(YEAR - 1996)"]*j
  }
}

AUD_plot_LITTERGRApreds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  AUD_plot_LITTERGRApreds_df [i,] <- quantile(AUD_plot_LITTERGRApreds[,i], c(0.025, 0.5, 0.975))
}

AUD_plot_LITTERGRApreds_df <- cbind.data.frame(lower = AUD_plot_LITTERGRApreds_df[,1], 
                                               mean = AUD_plot_LITTERGRApreds_df[,2], 
                                               upper = AUD_plot_LITTERGRApreds_df[,3], year = seq(1:23))
#Making graph
AUD_LITTER<- 
  ggplot() +
  geom_point(data= AUD_res3, 
             aes(x = YEAR, y = cover, color=TRTMT),
             alpha = 0.8, size = 2)+
  scale_color_manual(values = c("#1B9E77", "#7570B3", "#E6AB02"), labels = c("CTL","GRA", "OTC")) +
  scale_fill_manual(values = c("#1B9E77", "#7570B3", "#E6AB02"),  labels =c("CTL","GRA", "OTC")) +
  scale_x_continuous(breaks = c(1997, 2000, 2007, 2014, 2019)) +
  #  geom_ribbon(data = AUD_plot_LITTERCTLpreds_df, aes(x = year + 1996, ymin = lower, ymax = upper), fill = "#1B9E77", alpha = 0.2) +
  #  geom_line(data = AUD_plot_LITTERCTLpreds_df, aes(x = year + 1996, y = mean), colour = "#1B9E77") +
  # geom_ribbon(data = AUD_plot_LITTERGRApreds_df, aes(x = year + 1996, ymin = lower, ymax = upper), fill = "#7570B3", alpha = 0.2) +
  #  geom_line(data = AUD_plot_LITTERGRApreds_df, aes(x = year + 1996, y = mean), colour = "#7570B3") +
  geom_ribbon(data = AUD_plot_LITTEROTCpreds_df, aes(x = year + 1996, ymin = lower, ymax = upper), fill = "#E6AB02", alpha = 0.2) +
  geom_line(data = AUD_plot_LITTEROTCpreds_df, aes(x = year + 1996, y = mean), colour = "#E6AB02") +
  ylab("Number of hits")+
  theme_bw()+
  labs(title= "", subtitle=expression(italic("LITTER")))+
  theme(axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        #   axis.title.y=element_blank(),
        axis.ticks.x =element_blank(),
        legend.position = "none",
        text = element_text(size=9),
        plot.title=element_text(size=12),
        plot.margin = unit(c(0.5,0.5,0,0.5), "cm"))

#SANUNC
#CTL
nyears <- 23
niter <- length(AUD_plot_SANUNC_CTL$Sol[,"(Intercept)"])

AUD_plot_SANUNCCTLpreds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    AUD_plot_SANUNCCTLpreds[i,j] <- AUD_plot_SANUNC_CTL$Sol[i,"(Intercept)"] + 
      AUD_plot_SANUNC_CTL$Sol[i,"I(YEAR - 1996)"]*j
  }
}


AUD_plot_SANUNCCTLpreds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  AUD_plot_SANUNCCTLpreds_df [i,] <- quantile(AUD_plot_SANUNCCTLpreds[,i], c(0.025, 0.5, 0.975))
}

AUD_plot_SANUNCCTLpreds_df <- cbind.data.frame(lower = AUD_plot_SANUNCCTLpreds_df[,1], 
                                               mean = AUD_plot_SANUNCCTLpreds_df[,2], 
                                               upper = AUD_plot_SANUNCCTLpreds_df[,3], year = seq(1:23))
#OTC
nyears <- 23
niter <- length(AUD_plot_SANUNC_OTC$Sol[,"(Intercept)"])

AUD_plot_SANUNCOTCpreds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    AUD_plot_SANUNCOTCpreds[i,j] <- AUD_plot_SANUNC_OTC$Sol[i,"(Intercept)"] + 
      AUD_plot_SANUNC_OTC$Sol[i,"I(YEAR - 1996)"]*j
  }
}


AUD_plot_SANUNCOTCpreds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  AUD_plot_SANUNCOTCpreds_df [i,] <- quantile(AUD_plot_SANUNCOTCpreds[,i], c(0.025, 0.5, 0.975))
}

AUD_plot_SANUNCOTCpreds_df <- cbind.data.frame(lower = AUD_plot_SANUNCOTCpreds_df[,1], 
                                               mean = AUD_plot_SANUNCOTCpreds_df[,2], 
                                               upper = AUD_plot_SANUNCOTCpreds_df[,3], year = seq(1:23))

#GRA
nyears <- 23
niter <- length(AUD_plot_SANUNC_GRA$Sol[,"(Intercept)"])

AUD_plot_SANUNCGRApreds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    AUD_plot_SANUNCGRApreds[i,j] <- AUD_plot_SANUNC_GRA$Sol[i,"(Intercept)"] + 
      AUD_plot_SANUNC_GRA$Sol[i,"I(YEAR - 1996)"]*j
  }
}

AUD_plot_SANUNCGRApreds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  AUD_plot_SANUNCGRApreds_df [i,] <- quantile(AUD_plot_SANUNCGRApreds[,i], c(0.025, 0.5, 0.975))
}

AUD_plot_SANUNCGRApreds_df <- cbind.data.frame(lower = AUD_plot_SANUNCGRApreds_df[,1], 
                                               mean = AUD_plot_SANUNCGRApreds_df[,2], 
                                               upper = AUD_plot_SANUNCGRApreds_df[,3], year = seq(1:23))
#Making graph
AUD_SANUNC<- 
  ggplot() +
  geom_point(data= AUD_res3, 
             aes(x = YEAR, y = cover, color=TRTMT),
             alpha = 0.8, size = 2)+
  scale_color_manual(values = c("#1B9E77", "#7570B3", "#E6AB02"), labels = c("CTL","GRA", "OTC")) +
  scale_fill_manual(values = c("#1B9E77", "#7570B3", "#E6AB02"),  labels =c("CTL","GRA", "OTC")) +
  scale_x_continuous(breaks = c(1997, 2000, 2007, 2014, 2019)) +
  #  geom_ribbon(data = AUD_plot_SANUNCCTLpreds_df, aes(x = year + 1996, ymin = lower, ymax = upper), fill = "#1B9E77", alpha = 0.2) +
  #  geom_line(data = AUD_plot_SANUNCCTLpreds_df, aes(x = year + 1996, y = mean), colour = "#1B9E77") +
  # geom_ribbon(data = AUD_plot_SANUNCGRApreds_df, aes(x = year + 1996, ymin = lower, ymax = upper), fill = "#7570B3", alpha = 0.2) +
  #  geom_line(data = AUD_plot_SANUNCGRApreds_df, aes(x = year + 1996, y = mean), colour = "#7570B3") +
  geom_ribbon(data = AUD_plot_SANUNCOTCpreds_df, aes(x = year + 1996, ymin = lower, ymax = upper), fill = "#E6AB02", alpha = 0.2) +
  geom_line(data = AUD_plot_SANUNCOTCpreds_df, aes(x = year + 1996, y = mean), colour = "#E6AB02") +
  ylab("Number of hits")+
  theme_bw()+
  labs(title= "", subtitle=expression(italic("SANUNC")))+
  theme(axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        #   axis.title.y=element_blank(),
        axis.ticks.x =element_blank(),
        legend.position = "none",
        text = element_text(size=9),
        plot.title=element_text(size=12),
        plot.margin = unit(c(0.5,0.5,0,0.5), "cm"))

###################Thingsvellir
THIN <- coverTopBottom%>%
  filter(SITE=="THINGVELLIR")%>%
#  mutate(cover=cover*100)%>%
  na.omit()

# Linear model for all species THIN per treatment
#CTL

THIN_plot_all_CTL <- MCMCglmm(cover ~ I(YEAR-1994) * SPP, 
                              random = ~ YEAR+PLOT, 
                              data = THIN[THIN$TRTMT == "CTL",], 
                              family = "gaussian", pr = TRUE, nitt = 100000, 
                              burnin = 20000, prior = prior2)
summary(THIN_plot_all_CTL)


#OTC
THIN_plot_all_OTC <- MCMCglmm(cover ~ I(YEAR-1994) * SPP, 
                              random = ~ YEAR+PLOT, 
                              data = THIN[THIN$TRTMT == "OTC",], 
                              family = "gaussian", pr = TRUE, nitt = 100000, 
                              burnin = 20000, prior = prior2)
summary(THIN_plot_all_OTC)

#Linear model for main species to check for difference in treatment
#RACLAN
THIN_plot_RACLAN <- MCMCglmm(cover ~ I(YEAR-1994)+TRTMT-1,
                             random = ~ YEAR + PLOT,
                             data = THIN[THIN$SPP == "RACLAN",], 
                             family = "gaussian", pr = TRUE, nitt = 100000, 
                             burnin = 20000, prior = prior2)
summary(THIN_plot_RACLAN)

#GALBOR
THIN_plot_GALBOR<- MCMCglmm(cover ~ I(YEAR-1994)+TRTMT-1,
                             random = ~ YEAR + PLOT,
                             data = THIN[THIN$SPP == "GALBOR",], 
                             family = "gaussian", pr = TRUE, nitt = 100000, 
                             burnin = 20000, prior = prior2)
summary(THIN_plot_GALBOR)

#FESRUB
THIN_plot_FESRUB<- MCMCglmm(cover ~ I(YEAR-1994)+TRTMT-1,
                            random = ~ YEAR + PLOT,
                            data = THIN[THIN$SPP == "FESRUB",], 
                            family = "gaussian", pr = TRUE, nitt = 100000, 
                            burnin = 20000, prior = prior2)
summary(THIN_plot_FESRUB)

##Linear model for main species to plot model predictions
THIN_res1 <- subset(THIN, SPP =="RACLAN")
THIN_res2 <- subset(THIN, SPP =="GALBOR")
THIN_res3 <- subset(THIN, SPP =="FESRUB")

#RACLAN
THIN_plot_RACLAN_CTL <- MCMCglmm(cover ~ I(YEAR-1994), random = ~ YEAR + PLOT,
                                 data = THIN_res1[THIN_res1$TRTMT == "CTL",], 
                                 family = "gaussian", pr = TRUE, nitt = 100000, 
                                 burnin = 20000, prior = prior2)
summary(THIN_plot_RACLAN_CTL)

THIN_plot_RACLAN_OTC <- MCMCglmm(cover ~ I(YEAR-1994), random = ~ YEAR + PLOT,
                                 data = THIN_res1[THIN_res1$TRTMT == "OTC",], 
                                 family = "gaussian", pr = TRUE, nitt = 100000, 
                                 burnin = 20000, prior = prior2)
summary(THIN_plot_RACLAN_OTC)

#GALBOR
THIN_plot_GALBOR_CTL <- MCMCglmm(cover ~ I(YEAR-1994), random = ~ YEAR + PLOT,
                                 data = THIN_res2[THIN_res2$TRTMT == "CTL",], 
                                 family = "gaussian", pr = TRUE, nitt = 100000, 
                                 burnin = 20000, prior = prior2)
summary(THIN_plot_GALBOR_CTL)

THIN_plot_GALBOR_OTC <- MCMCglmm(cover ~ I(YEAR-1994), random = ~ YEAR + PLOT,
                                 data = THIN_res2[THIN_res2$TRTMT == "OTC",], 
                                 family = "gaussian", pr = TRUE, nitt = 100000, 
                                 burnin = 20000, prior = prior2)
summary(THIN_plot_GALBOR_OTC)


#FESRUB
THIN_plot_FESRUB_CTL <- MCMCglmm(cover ~ I(YEAR-1994), random = ~ YEAR + PLOT,
                                 data = THIN_res3[THIN_res3$TRTMT == "CTL",], 
                                 family = "gaussian", pr = TRUE, nitt = 100000, 
                                 burnin = 20000, prior = prior2)
summary(THIN_plot_FESRUB_CTL)

THIN_plot_FESRUB_OTC <- MCMCglmm(cover ~ I(YEAR-1994), random = ~ YEAR + PLOT,
                                 data = THIN_res3[THIN_res3$TRTMT == "OTC",], 
                                 family = "gaussian", pr = TRUE, nitt = 100000, 
                                 burnin = 20000, prior = prior2)
summary(THIN_plot_FESRUB_OTC)

#RACLAN
#CTL
nyears <- 25
niter <- length(THIN_plot_RACLAN_CTL$Sol[,"(Intercept)"])

THIN_plot_RACLANCTLpreds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    THIN_plot_RACLANCTLpreds[i,j] <- THIN_plot_RACLAN_CTL$Sol[i,"(Intercept)"] + 
      THIN_plot_RACLAN_CTL$Sol[i,"I(YEAR - 1994)"]*j
  }
}


THIN_plot_RACLANCTLpreds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  THIN_plot_RACLANCTLpreds_df [i,] <- quantile(THIN_plot_RACLANCTLpreds[,i], c(0.025, 0.5, 0.975))
}

THIN_plot_RACLANCTLpreds_df <- cbind.data.frame(lower = THIN_plot_RACLANCTLpreds_df[,1], 
                                                mean = THIN_plot_RACLANCTLpreds_df[,2], 
                                                upper = THIN_plot_RACLANCTLpreds_df[,3], year = seq(1:25))
#OTC
nyears <- 25
niter <- length(THIN_plot_RACLAN_OTC$Sol[,"(Intercept)"])

THIN_plot_RACLANOTCpreds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    THIN_plot_RACLANOTCpreds[i,j] <- THIN_plot_RACLAN_OTC$Sol[i,"(Intercept)"] + 
      THIN_plot_RACLAN_OTC$Sol[i,"I(YEAR - 1994)"]*j
  }
}

THIN_plot_RACLANOTCpreds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  THIN_plot_RACLANOTCpreds_df [i,] <- quantile(THIN_plot_RACLANOTCpreds[,i], c(0.025, 0.5, 0.975))
}

THIN_plot_RACLANOTCpreds_df <- cbind.data.frame(lower = THIN_plot_RACLANOTCpreds_df[,1], 
                                                mean = THIN_plot_RACLANOTCpreds_df[,2], 
                                                upper = THIN_plot_RACLANOTCpreds_df[,3], year = seq(1:25))

# Making GRAPH
THIN_RACLAN<- 
  ggplot() +
  geom_point(data= THIN_res1, 
             aes(x = YEAR, y = cover, color=TRTMT),
             alpha = 0.8, size = 2)+
  scale_color_manual(values = c("#1B9E77",  "#E6AB02"), labels = c("CTL", "OTC")) +
  scale_fill_manual(values = c("#1B9E77", "#E6AB02"),  labels =c("CTL", "OTC")) +
  scale_x_continuous(breaks = c(1995, 1998,2000, 2007, 2014, 2019)) +
  geom_ribbon(data = THIN_plot_RACLANCTLpreds_df, aes(x = year + 1994, ymin = lower, ymax = upper), fill = "#1B9E77", alpha = 0.2) +
  geom_line(data = THIN_plot_RACLANCTLpreds_df, aes(x = year + 1994, y = mean), colour = "#1B9E77") +
  geom_ribbon(data = THIN_plot_RACLANOTCpreds_df, aes(x = year + 1994, ymin = lower, ymax = upper), fill = "#E6AB02", alpha = 0.2) +
  geom_line(data = THIN_plot_RACLANOTCpreds_df, aes(x = year + 1994, y = mean), colour = "#E6AB02") +
  theme_bw()+
  ylab("Cover (%)")+
  xlab("Year")+
  labs(title= "Thingvellir", subtitle=expression(italic("Racomitrium lanuginosum")), size=9)+
  theme(legend.position="none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(size=9),
        plot.title=element_text(size=12),
        plot.margin = unit(c(0.5,0.5,0,0.5), "cm"))

#GALBOR
#CTL
nyears <- 25
niter <- length(THIN_plot_GALBOR_CTL$Sol[,"(Intercept)"])

THIN_plot_GALBORCTLpreds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    THIN_plot_GALBORCTLpreds[i,j] <- THIN_plot_GALBOR_CTL$Sol[i,"(Intercept)"] + 
      THIN_plot_GALBOR_CTL$Sol[i,"I(YEAR - 1994)"]*j
  }
}


THIN_plot_GALBORCTLpreds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  THIN_plot_GALBORCTLpreds_df [i,] <- quantile(THIN_plot_GALBORCTLpreds[,i], c(0.025, 0.5, 0.975))
}

THIN_plot_GALBORCTLpreds_df <- cbind.data.frame(lower = THIN_plot_GALBORCTLpreds_df[,1], 
                                                mean = THIN_plot_GALBORCTLpreds_df[,2], 
                                                upper = THIN_plot_GALBORCTLpreds_df[,3], year = seq(1:25))
#OTC
nyears <- 25
niter <- length(THIN_plot_GALBOR_OTC$Sol[,"(Intercept)"])

THIN_plot_GALBOROTCpreds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    THIN_plot_GALBOROTCpreds[i,j] <- THIN_plot_GALBOR_OTC$Sol[i,"(Intercept)"] + 
      THIN_plot_GALBOR_OTC$Sol[i,"I(YEAR - 1994)"]*j
  }
}

THIN_plot_GALBOROTCpreds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  THIN_plot_GALBOROTCpreds_df [i,] <- quantile(THIN_plot_GALBOROTCpreds[,i], c(0.025, 0.5, 0.975))
}

THIN_plot_GALBOROTCpreds_df <- cbind.data.frame(lower = THIN_plot_GALBOROTCpreds_df[,1], 
                                                mean = THIN_plot_GALBOROTCpreds_df[,2], 
                                                upper = THIN_plot_GALBOROTCpreds_df[,3], year = seq(1:25))

# Making GRAPH
THIN_GALBOR<- 
  ggplot() +
  geom_point(data= THIN_res2, 
             aes(x = YEAR, y = cover, color=TRTMT),
             alpha = 0.8, size = 2)+
  scale_color_manual(values = c("#1B9E77",  "#E6AB02"), labels = c("CTL", "OTC")) +
  scale_fill_manual(values = c("#1B9E77", "#E6AB02"),  labels =c("CTL", "OTC")) +
  scale_x_continuous(breaks = c(1995, 1998,2000, 2007, 2014, 2019)) +
 # geom_ribbon(data = THIN_plot_GALBORCTLpreds_df, aes(x = year + 1994, ymin = lower, ymax = upper), fill = "#1B9E77", alpha = 0.2) +
#  geom_line(data = THIN_plot_GALBORCTLpreds_df, aes(x = year + 1994, y = mean), colour = "#1B9E77") +
  geom_ribbon(data = THIN_plot_GALBOROTCpreds_df, aes(x = year + 1994, ymin = lower, ymax = upper), fill = "#E6AB02", alpha = 0.2) +
  geom_line(data = THIN_plot_GALBOROTCpreds_df, aes(x = year + 1994, y = mean), colour = "#E6AB02") +
  theme_bw()+
  ylab("Cover (%)")+
  xlab("Year")+
  labs(title= "", subtitle=expression(italic("GALBOR")), size=9)+
  theme(legend.position="none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(size=9),
        plot.title=element_text(size=12),
        plot.margin = unit(c(0.5,0.5,0,0.5), "cm"))

##
#FESRUB
#CTL
nyears <- 25
niter <- length(THIN_plot_FESRUB_CTL$Sol[,"(Intercept)"])

THIN_plot_FESRUBCTLpreds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    THIN_plot_FESRUBCTLpreds[i,j] <- THIN_plot_FESRUB_CTL$Sol[i,"(Intercept)"] + 
      THIN_plot_FESRUB_CTL$Sol[i,"I(YEAR - 1994)"]*j
  }
}


THIN_plot_FESRUBCTLpreds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  THIN_plot_FESRUBCTLpreds_df [i,] <- quantile(THIN_plot_FESRUBCTLpreds[,i], c(0.025, 0.5, 0.975))
}

THIN_plot_FESRUBCTLpreds_df <- cbind.data.frame(lower = THIN_plot_FESRUBCTLpreds_df[,1], 
                                                mean = THIN_plot_FESRUBCTLpreds_df[,2], 
                                                upper = THIN_plot_FESRUBCTLpreds_df[,3], year = seq(1:25))
#OTC
nyears <- 25
niter <- length(THIN_plot_FESRUB_OTC$Sol[,"(Intercept)"])

THIN_plot_FESRUBOTCpreds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    THIN_plot_FESRUBOTCpreds[i,j] <- THIN_plot_FESRUB_OTC$Sol[i,"(Intercept)"] + 
      THIN_plot_FESRUB_OTC$Sol[i,"I(YEAR - 1994)"]*j
  }
}

THIN_plot_FESRUBOTCpreds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  THIN_plot_FESRUBOTCpreds_df [i,] <- quantile(THIN_plot_FESRUBOTCpreds[,i], c(0.025, 0.5, 0.975))
}

THIN_plot_FESRUBOTCpreds_df <- cbind.data.frame(lower = THIN_plot_FESRUBOTCpreds_df[,1], 
                                                mean = THIN_plot_FESRUBOTCpreds_df[,2], 
                                                upper = THIN_plot_FESRUBOTCpreds_df[,3], year = seq(1:25))

# Making GRAPH
THIN_FESRUB<- 
  ggplot() +
  geom_point(data= THIN_res3, 
             aes(x = YEAR, y = cover, color=TRTMT),
             alpha = 0.8, size = 2)+
  scale_color_manual(values = c("#1B9E77",  "#E6AB02"), labels = c("CTL", "OTC")) +
  scale_fill_manual(values = c("#1B9E77", "#E6AB02"),  labels =c("CTL", "OTC")) +
  scale_x_continuous(breaks = c(1995, 1998,2000, 2007, 2014, 2019)) +
  geom_ribbon(data = THIN_plot_FESRUBCTLpreds_df, aes(x = year + 1994, ymin = lower, ymax = upper), fill = "#1B9E77", alpha = 0.2) +
  geom_line(data = THIN_plot_FESRUBCTLpreds_df, aes(x = year + 1994, y = mean), colour = "#1B9E77") +
  #geom_ribbon(data = THIN_plot_FESRUBOTCpreds_df, aes(x = year + 1994, ymin = lower, ymax = upper), fill = "#E6AB02", alpha = 0.2) +
  #geom_line(data = THIN_plot_FESRUBOTCpreds_df, aes(x = year + 1994, y = mean), colour = "#E6AB02") +
  theme_bw()+
  ylab("Cover (%)")+
  xlab("Year")+
  labs(title= "", subtitle=expression(italic("FESRUB")), size=9)+
  theme(legend.position="none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(size=9),
        plot.title=element_text(size=12),
        plot.margin = unit(c(0.5,0.5,0,0.5), "cm"))

prow <- plot_grid(AUD_RACLAN, THIN_RACLAN + theme(legend.position="none"), 
                  AUD_BETNAN, THIN_GALBOR ,
                  AUD_CETISL, THIN_FESRUB ,
                  align = 'v',
                  hjust = -1,
                  ncol=2,
                  nrow = 3)

legend_b <- get_legend(AUD_RACLAN + theme(legend.position="bottom"))
p<-plot_grid( prow, legend_b, ncol = 1, rel_heights = c(1, 0.1))

##############################Species change with top only################

### Linear model for all species AUD Top
#CTL

AUD <- coverCanopy%>%
  filter(SITE=="AUDKULUHEIDI")%>%
  mutate(cover=cover*100)%>%
  na.omit()


AUD_plot_all_CTL <- MCMCglmm(cover ~ I(YEAR-1996) * SPP, 
                             random = ~ YEAR+PLOT, data = AUD[AUD$TRTMT == "CTL",], 
                             family = "gaussian", pr = TRUE, nitt = 100000, 
                             burnin = 20000, prior = prior2)
summary(AUD_plot_all_CTL)

#OTC
AUD_plot_all_OTC <- MCMCglmm(cover ~ I(YEAR-1996) * SPP, 
                             random = ~ YEAR+PLOT, data = AUD[AUD$TRTMT == "OTC",], 
                             family = "gaussian", pr = TRUE, nitt = 100000, 
                             burnin = 20000, prior = prior2)


summary(AUD_plot_all_OTC)

#GRA
AUD_plot_all_GRA <- MCMCglmm(cover ~ I(YEAR-1996) * SPP, 
                             random = ~ YEAR+PLOT, data = AUD[AUD$TRTMT == "GRA",], 
                             family = "gaussian", pr = TRUE, nitt = 100000, 
                             burnin = 20000, prior = prior2)
summary(AUD_plot_all_GRA)

#Linear model for main species to check for difference in treatment

AUD_plot_BETNAN <- MCMCglmm(cover ~ I(YEAR-1996)+TRTMT-1, 
                            random = ~ YEAR + PLOT,
                            data = AUD[AUD$SPP == "BETNAN",], 
                            family = "gaussian", pr = TRUE, nitt = 100000, 
                            burnin = 20000, prior = prior2)
summary(AUD_plot_BETNAN)

AUD_plot_RACLAN <- MCMCglmm(cover ~ I(YEAR-1996)+TRTMT-1, 
                            random = ~ YEAR + PLOT,
                            data = AUD[AUD$SPP == "RACLAN",], 
                            family = "gaussian", pr = TRUE, nitt = 100000, 
                            burnin = 20000, prior = prior2)
summary(AUD_plot_RACLAN)

AUD_plot_BISVIV <- MCMCglmm(cover ~ I(YEAR-1996)+TRTMT-1, 
                            random = ~ YEAR + PLOT,
                            data = AUD[AUD$SPP == "BISVIV",], 
                            family = "gaussian", pr = TRUE, nitt = 100000, 
                            burnin = 20000, prior = prior2)
summary(AUD_plot_BISVIV)

AUD_plot_EMPNIG <- MCMCglmm(cover ~ I(YEAR-1996)+TRTMT-1, 
                            random = ~ YEAR + PLOT,
                            data = AUD[AUD$SPP == "EMPNIG",], 
                            family = "gaussian", pr = TRUE, nitt = 100000, 
                            burnin = 20000, prior = prior2)
summary(AUD_plot_EMPNIG)

##Linear model for main species to plot model predictions
AUD_res1 <- subset(AUD, SPP =="BETNAN")
AUD_res2 <- subset(AUD, SPP =="RACLAN")
AUD_res3 <- subset(AUD, SPP =="BISVIV")
AUD_res4 <- subset(AUD, SPP =="EMPNIG")

#BETNAN
AUD_plot_BETNAN_CTL <- MCMCglmm(cover ~ I(YEAR-1996), random = ~ YEAR + PLOT,
                                data = AUD_res1[AUD_res1$TRTMT == "CTL",], 
                                family = "gaussian", pr = TRUE, nitt = 100000, 
                                burnin = 20000, prior = prior2)
summary(AUD_plot_BETNAN_CTL)

AUD_plot_BETNAN_OTC <- MCMCglmm(cover ~ I(YEAR-1996), random = ~ YEAR + PLOT,
                                data = AUD_res1[AUD_res1$TRTMT == "OTC",], 
                                family = "gaussian", pr = TRUE, nitt = 100000, 
                                burnin = 20000, prior = prior2)
summary(AUD_plot_BETNAN_OTC)

AUD_plot_BETNAN_GRA <- MCMCglmm(cover ~ I(YEAR-1996), random = ~ YEAR + PLOT,
                                data = AUD_res1[AUD_res1$TRTMT == "GRA",], 
                                family = "gaussian", pr = TRUE, nitt = 100000, 
                                burnin = 20000, prior = prior2)
summary(AUD_plot_BETNAN_GRA)

#RACLAN
AUD_plot_RACLAN_CTL <- MCMCglmm(cover ~ I(YEAR-1996), random = ~ YEAR + PLOT,
                                data = AUD_res2[AUD_res2$TRTMT == "CTL",], 
                                family = "gaussian", pr = TRUE, nitt = 100000, 
                                burnin = 20000, prior = prior2)
summary(AUD_plot_RACLAN_CTL)

AUD_plot_RACLAN_OTC <- MCMCglmm(cover ~ I(YEAR-1996), random = ~ YEAR + PLOT,
                                data = AUD_res2[AUD_res2$TRTMT == "OTC",], 
                                family = "gaussian", pr = TRUE, nitt = 100000, 
                                burnin = 20000, prior = prior2)
summary(AUD_plot_RACLAN_OTC)

AUD_plot_RACLAN_GRA <- MCMCglmm(cover ~ I(YEAR-1996), random = ~ YEAR + PLOT,
                                data = AUD_res2[AUD_res2$TRTMT == "GRA",], 
                                family = "gaussian", pr = TRUE, nitt = 100000, 
                                burnin = 20000, prior = prior2)
summary(AUD_plot_RACLAN_GRA)

#BISVIV
AUD_plot_BISVIV_CTL <- MCMCglmm(cover ~ I(YEAR-1996), random = ~ YEAR + PLOT,
                                data = AUD_res3[AUD_res3$TRTMT == "CTL",], 
                                family = "gaussian", pr = TRUE, nitt = 100000, 
                                burnin = 20000, prior = prior2)
summary(AUD_plot_BISVIV_CTL)

AUD_plot_BISVIV_OTC <- MCMCglmm(cover ~ I(YEAR-1996), random = ~ YEAR + PLOT,
                                data = AUD_res3[AUD_res3$TRTMT == "OTC",], 
                                family = "gaussian", pr = TRUE, nitt = 100000, 
                                burnin = 20000, prior = prior2)
summary(AUD_plot_BISVIV_OTC)

AUD_plot_BISVIV_GRA <- MCMCglmm(cover ~ I(YEAR-1996), random = ~ YEAR + PLOT,
                                data = AUD_res3[AUD_res3$TRTMT == "GRA",], 
                                family = "gaussian", pr = TRUE, nitt = 100000, 
                                burnin = 20000, prior = prior2)
summary(AUD_plot_BISVIV_GRA)

#EMPNIG
AUD_plot_EMPNIG_CTL <- MCMCglmm(cover ~ I(YEAR-1996), random = ~ YEAR + PLOT,
                                data = AUD_res4[AUD_res4$TRTMT == "CTL",], 
                                family = "gaussian", pr = TRUE, nitt = 100000, 
                                burnin = 20000, prior = prior2)
summary(AUD_plot_EMPNIG_CTL)

AUD_plot_EMPNIG_OTC <- MCMCglmm(cover ~ I(YEAR-1996), random = ~ YEAR + PLOT,
                                data = AUD_res4[AUD_res4$TRTMT == "OTC",], 
                                family = "gaussian", pr = TRUE, nitt = 100000, 
                                burnin = 20000, prior = prior2)
summary(AUD_plot_EMPNIG_OTC)

AUD_plot_EMPNIG_GRA <- MCMCglmm(cover ~ I(YEAR-1996), random = ~ YEAR + PLOT,
                                data = AUD_res4[AUD_res4$TRTMT == "GRA",], 
                                family = "gaussian", pr = TRUE, nitt = 100000, 
                                burnin = 20000, prior = prior2)
summary(AUD_plot_EMPNIG_GRA)


# Calculating model predictions
#BETNAN
#CTL
nyears <- 23
niter <- length(AUD_plot_BETNAN_CTL$Sol[,"(Intercept)"])

AUD_plot_BETNANCTLpreds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    AUD_plot_BETNANCTLpreds[i,j] <- AUD_plot_BETNAN_CTL$Sol[i,"(Intercept)"] + AUD_plot_BETNAN_CTL$Sol[i,"I(YEAR - 1996)"]*j
  }
}

AUD_plot_BETNANCTLpreds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  AUD_plot_BETNANCTLpreds_df [i,] <- quantile(AUD_plot_BETNANCTLpreds[,i], c(0.025, 0.5, 0.975))
}

AUD_plot_BETNANCTLpreds_df <- cbind.data.frame(lower = AUD_plot_BETNANCTLpreds_df[,1], 
                                               mean = AUD_plot_BETNANCTLpreds_df[,2], 
                                               upper = AUD_plot_BETNANCTLpreds_df[,3], year = seq(1:23))
#OTC
nyears <- 23
niter <- length(AUD_plot_BETNAN_OTC$Sol[,"(Intercept)"])

AUD_plot_BETNANOTCpreds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    AUD_plot_BETNANOTCpreds[i,j] <- AUD_plot_BETNAN_OTC$Sol[i,"(Intercept)"] + AUD_plot_BETNAN_OTC$Sol[i,"I(YEAR - 1996)"]*j
  }
}

AUD_plot_BETNANOTCpreds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  AUD_plot_BETNANOTCpreds_df [i,] <- quantile(AUD_plot_BETNANOTCpreds[,i], c(0.025, 0.5, 0.975))
}

AUD_plot_BETNANOTCpreds_df <- cbind.data.frame(lower = AUD_plot_BETNANOTCpreds_df[,1], 
                                               mean = AUD_plot_BETNANOTCpreds_df[,2], 
                                               upper = AUD_plot_BETNANOTCpreds_df[,3], year = seq(1:23))

#GRA
nyears <- 23
niter <- length(AUD_plot_BETNAN_GRA$Sol[,"(Intercept)"])

AUD_plot_BETNANGRApreds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    AUD_plot_BETNANGRApreds[i,j] <- AUD_plot_BETNAN_GRA$Sol[i,"(Intercept)"] + AUD_plot_BETNAN_GRA$Sol[i,"I(YEAR - 1996)"]*j
  }
}

AUD_plot_BETNANGRApreds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  AUD_plot_BETNANGRApreds_df [i,] <- quantile(AUD_plot_BETNANGRApreds[,i], c(0.025, 0.5, 0.975))
}

AUD_plot_BETNANGRApreds_df <- cbind.data.frame(lower = AUD_plot_BETNANGRApreds_df[,1], 
                                               mean = AUD_plot_BETNANGRApreds_df[,2], 
                                               upper = AUD_plot_BETNANGRApreds_df[,3], year = seq(1:23))

#Making graph
AUD_BETNAN<- 
  ggplot() +
  geom_point(data= AUD_res1, 
             aes(x = YEAR, y = cover, color=TRTMT),
             alpha = 0.8, size = 2)+
  scale_color_manual(values = c("#1B9E77", "#7570B3", "#E6AB02"), labels = c("CTL","GRA", "OTC")) +
  scale_fill_manual(values = c("#1B9E77", "#7570B3", "#E6AB02"),  labels =c("CTL","GRA", "OTC")) +
  scale_x_continuous(breaks = c(1997, 2000, 2007, 2014, 2019)) +
 #   geom_ribbon(data = AUD_plot_BETNANCTLpreds_df, aes(x = year + 1996, ymin = lower, ymax = upper), fill = "#1B9E77", alpha = 0.2) +
 #   geom_line(data = AUD_plot_BETNANCTLpreds_df, aes(x = year + 1996, y = mean), colour = "#1B9E77") +
 # geom_ribbon(data = AUD_plot_BETNANGRApreds_df, aes(x = year + 1996, ymin = lower, ymax = upper), fill = "#7570B3", alpha = 0.2) +
 # geom_line(data = AUD_plot_BETNANGRApreds_df, aes(x = year + 1996, y = mean), colour = "#7570B3") +
  geom_ribbon(data = AUD_plot_BETNANOTCpreds_df, aes(x = year + 1996, ymin = lower, ymax = upper), fill = "#E6AB02", alpha = 0.2) +
  geom_line(data = AUD_plot_BETNANOTCpreds_df, aes(x = year + 1996, y = mean), colour = "#E6AB02") +
  theme_bw()+
  labs(title="",subtitle=expression(italic("Betula nana")), size=9)+
  theme(axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.x =element_blank(),
        text = element_text(size=9),
        legend.position="none",
        plot.margin = unit(c(0,0.5,0,0.5), "cm"))

#RACLAN
#CTL
nyears <- 23
niter <- length(AUD_plot_RACLAN_CTL$Sol[,"(Intercept)"])

AUD_plot_RACLANCTLpreds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    AUD_plot_RACLANCTLpreds[i,j] <- AUD_plot_RACLAN_CTL$Sol[i,"(Intercept)"] + 
      AUD_plot_RACLAN_CTL$Sol[i,"I(YEAR - 1996)"]*j
  }
}

AUD_plot_RACLANCTLpreds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  AUD_plot_RACLANCTLpreds_df [i,] <- quantile(AUD_plot_RACLANCTLpreds[,i], c(0.025, 0.5, 0.975))
}

AUD_plot_RACLANCTLpreds_df <- cbind.data.frame(lower = AUD_plot_RACLANCTLpreds_df[,1], 
                                               mean = AUD_plot_RACLANCTLpreds_df[,2], 
                                               upper = AUD_plot_RACLANCTLpreds_df[,3], year = seq(1:23))
#OTC
nyears <- 23
niter <- length(AUD_plot_RACLAN_OTC$Sol[,"(Intercept)"])

AUD_plot_RACLANOTCpreds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    AUD_plot_RACLANOTCpreds[i,j] <- AUD_plot_RACLAN_OTC$Sol[i,"(Intercept)"] + 
      AUD_plot_RACLAN_OTC$Sol[i,"I(YEAR - 1996)"]*j
  }
}

AUD_plot_RACLANOTCpreds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  AUD_plot_RACLANOTCpreds_df [i,] <- quantile(AUD_plot_RACLANOTCpreds[,i], c(0.025, 0.5, 0.975))
}

AUD_plot_RACLANOTCpreds_df <- cbind.data.frame(lower = AUD_plot_RACLANOTCpreds_df[,1], 
                                               mean = AUD_plot_RACLANOTCpreds_df[,2], 
                                               upper = AUD_plot_RACLANOTCpreds_df[,3], year = seq(1:23))

#GRA
nyears <- 23
niter <- length(AUD_plot_RACLAN_GRA$Sol[,"(Intercept)"])

AUD_plot_RACLANGRApreds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    AUD_plot_RACLANGRApreds[i,j] <- AUD_plot_RACLAN_GRA$Sol[i,"(Intercept)"] + 
      AUD_plot_RACLAN_GRA$Sol[i,"I(YEAR - 1996)"]*j
  }
}

AUD_plot_RACLANGRApreds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  AUD_plot_RACLANGRApreds_df [i,] <- quantile(AUD_plot_RACLANGRApreds[,i], c(0.025, 0.5, 0.975))
}

AUD_plot_RACLANGRApreds_df <- cbind.data.frame(lower = AUD_plot_RACLANGRApreds_df[,1], 
                                               mean = AUD_plot_RACLANGRApreds_df[,2], 
                                               upper = AUD_plot_RACLANGRApreds_df[,3], year = seq(1:23))
#Making graph
AUD_RACLAN<- 
  ggplot() +
  geom_point(data= AUD_res2, 
             aes(x = YEAR, y = cover, color=TRTMT),
             alpha = 0.8, size = 2)+
  scale_color_manual(values = c("#1B9E77", "#7570B3", "#E6AB02"), labels = c("CTL","GRA", "OTC")) +
  scale_fill_manual(values = c("#1B9E77", "#7570B3", "#E6AB02"),  labels =c("CTL","GRA", "OTC")) +
  scale_x_continuous(breaks = c(1997, 2000, 2007, 2014, 2019)) +
  geom_ribbon(data = AUD_plot_RACLANCTLpreds_df, aes(x = year + 1996, ymin = lower, ymax = upper), fill = "#1B9E77", alpha = 0.2) +
  geom_line(data = AUD_plot_RACLANCTLpreds_df, aes(x = year + 1996, y = mean), colour = "#1B9E77") +
  geom_ribbon(data = AUD_plot_RACLANGRApreds_df, aes(x = year + 1996, ymin = lower, ymax = upper), fill = "#7570B3", alpha = 0.2) +
  geom_line(data = AUD_plot_RACLANGRApreds_df, aes(x = year + 1996, y = mean), colour = "#7570B3") +
  geom_ribbon(data = AUD_plot_RACLANOTCpreds_df, aes(x = year + 1996, ymin = lower, ymax = upper), fill = "#E6AB02", alpha = 0.2) +
  geom_line(data = AUD_plot_RACLANOTCpreds_df, aes(x = year + 1996, y = mean), colour = "#E6AB02") +
  ylab("Cover (%)")+
  theme_bw()+
  labs(title= "Audkuluheidi", subtitle=expression(italic("Racomitrium lanuginosum")))+
  theme(axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.x =element_blank(),
        legend.position = "none",
        text = element_text(size=9),
        plot.title=element_text(size=12),
        plot.margin = unit(c(0.5,0.5,0,0.5), "cm"))

#BISVIV
#CTL
nyears <- 23
niter <- length(AUD_plot_BISVIV_CTL$Sol[,"(Intercept)"])

AUD_plot_BISVIVCTLpreds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    AUD_plot_BISVIVCTLpreds[i,j] <- AUD_plot_BISVIV_CTL$Sol[i,"(Intercept)"] + 
      AUD_plot_BISVIV_CTL$Sol[i,"I(YEAR - 1996)"]*j
  }
}

AUD_plot_BISVIVCTLpreds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  AUD_plot_BISVIVCTLpreds_df [i,] <- quantile(AUD_plot_BISVIVCTLpreds[,i], c(0.025, 0.5, 0.975))
}

AUD_plot_BISVIVCTLpreds_df <- cbind.data.frame(lower = AUD_plot_BISVIVCTLpreds_df[,1], 
                                               mean = AUD_plot_BISVIVCTLpreds_df[,2], 
                                               upper = AUD_plot_BISVIVCTLpreds_df[,3], year = seq(1:23))
#OTC
nyears <- 23
niter <- length(AUD_plot_BISVIV_OTC$Sol[,"(Intercept)"])

AUD_plot_BISVIVOTCpreds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    AUD_plot_BISVIVOTCpreds[i,j] <- AUD_plot_BISVIV_OTC$Sol[i,"(Intercept)"] + 
      AUD_plot_BISVIV_OTC$Sol[i,"I(YEAR - 1996)"]*j
  }
}

AUD_plot_BISVIVOTCpreds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  AUD_plot_BISVIVOTCpreds_df [i,] <- quantile(AUD_plot_BISVIVOTCpreds[,i], c(0.025, 0.5, 0.975))
}

AUD_plot_BISVIVOTCpreds_df <- cbind.data.frame(lower = AUD_plot_BISVIVOTCpreds_df[,1], 
                                               mean = AUD_plot_BISVIVOTCpreds_df[,2], 
                                               upper = AUD_plot_BISVIVOTCpreds_df[,3], year = seq(1:23))

#GRA
nyears <- 23
niter <- length(AUD_plot_BISVIV_GRA$Sol[,"(Intercept)"])

AUD_plot_BISVIVGRApreds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    AUD_plot_BISVIVGRApreds[i,j] <- AUD_plot_BISVIV_GRA$Sol[i,"(Intercept)"] + 
      AUD_plot_BISVIV_GRA$Sol[i,"I(YEAR - 1996)"]*j
  }
}

AUD_plot_BISVIVGRApreds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  AUD_plot_BISVIVGRApreds_df [i,] <- quantile(AUD_plot_BISVIVGRApreds[,i], c(0.025, 0.5, 0.975))
}

AUD_plot_BISVIVGRApreds_df <- cbind.data.frame(lower = AUD_plot_BISVIVGRApreds_df[,1], 
                                               mean = AUD_plot_BISVIVGRApreds_df[,2], 
                                               upper = AUD_plot_BISVIVGRApreds_df[,3], year = seq(1:23))
#Making graph
AUD_BISVIV<- 
  ggplot() +
  geom_point(data= AUD_res3, 
             aes(x = YEAR, y = cover, color=TRTMT),
             alpha = 0.8, size = 2)+
  scale_color_manual(values = c("#1B9E77", "#7570B3", "#E6AB02"), labels = c("CTL","GRA", "OTC")) +
  scale_fill_manual(values = c("#1B9E77", "#7570B3", "#E6AB02"),  labels =c("CTL","GRA", "OTC")) +
  scale_x_continuous(breaks = c(1997, 2000, 2007, 2014, 2019)) +
  geom_ribbon(data = AUD_plot_BISVIVCTLpreds_df, aes(x = year + 1996, ymin = lower, ymax = upper), fill = "#1B9E77", alpha = 0.2) +
  geom_line(data = AUD_plot_BISVIVCTLpreds_df, aes(x = year + 1996, y = mean), colour = "#1B9E77") +
  geom_ribbon(data = AUD_plot_BISVIVGRApreds_df, aes(x = year + 1996, ymin = lower, ymax = upper), fill = "#7570B3", alpha = 0.2) +
  geom_line(data = AUD_plot_BISVIVGRApreds_df, aes(x = year + 1996, y = mean), colour = "#7570B3") +
  geom_ribbon(data = AUD_plot_BISVIVOTCpreds_df, aes(x = year + 1996, ymin = lower, ymax = upper), fill = "#E6AB02", alpha = 0.2) +
  geom_line(data = AUD_plot_BISVIVOTCpreds_df, aes(x = year + 1996, y = mean), colour = "#E6AB02") +
  ylab("Cover (%)")+
  theme_bw()+
  labs(title= "", subtitle=expression(italic("BISVIV")))+
  theme(axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.x =element_blank(),
        legend.position = "none",
        text = element_text(size=9),
        plot.title=element_text(size=12),
        plot.margin = unit(c(0.5,0.5,0,0.5), "cm"))

#EMPNIG
#CTL
nyears <- 23
niter <- length(AUD_plot_EMPNIG_CTL$Sol[,"(Intercept)"])

AUD_plot_EMPNIGCTLpreds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    AUD_plot_EMPNIGCTLpreds[i,j] <- AUD_plot_EMPNIG_CTL$Sol[i,"(Intercept)"] + 
      AUD_plot_EMPNIG_CTL$Sol[i,"I(YEAR - 1996)"]*j
  }
}

AUD_plot_EMPNIGCTLpreds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  AUD_plot_EMPNIGCTLpreds_df [i,] <- quantile(AUD_plot_EMPNIGCTLpreds[,i], c(0.025, 0.5, 0.975))
}

AUD_plot_EMPNIGCTLpreds_df <- cbind.data.frame(lower = AUD_plot_EMPNIGCTLpreds_df[,1], 
                                               mean = AUD_plot_EMPNIGCTLpreds_df[,2], 
                                               upper = AUD_plot_EMPNIGCTLpreds_df[,3], year = seq(1:23))
#OTC
nyears <- 23
niter <- length(AUD_plot_EMPNIG_OTC$Sol[,"(Intercept)"])

AUD_plot_EMPNIGOTCpreds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    AUD_plot_EMPNIGOTCpreds[i,j] <- AUD_plot_EMPNIG_OTC$Sol[i,"(Intercept)"] + 
      AUD_plot_EMPNIG_OTC$Sol[i,"I(YEAR - 1996)"]*j
  }
}

AUD_plot_EMPNIGOTCpreds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  AUD_plot_EMPNIGOTCpreds_df [i,] <- quantile(AUD_plot_EMPNIGOTCpreds[,i], c(0.025, 0.5, 0.975))
}

AUD_plot_EMPNIGOTCpreds_df <- cbind.data.frame(lower = AUD_plot_EMPNIGOTCpreds_df[,1], 
                                               mean = AUD_plot_EMPNIGOTCpreds_df[,2], 
                                               upper = AUD_plot_EMPNIGOTCpreds_df[,3], year = seq(1:23))

#GRA
nyears <- 23
niter <- length(AUD_plot_EMPNIG_GRA$Sol[,"(Intercept)"])

AUD_plot_EMPNIGGRApreds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    AUD_plot_EMPNIGGRApreds[i,j] <- AUD_plot_EMPNIG_GRA$Sol[i,"(Intercept)"] + 
      AUD_plot_EMPNIG_GRA$Sol[i,"I(YEAR - 1996)"]*j
  }
}

AUD_plot_EMPNIGGRApreds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  AUD_plot_EMPNIGGRApreds_df [i,] <- quantile(AUD_plot_EMPNIGGRApreds[,i], c(0.025, 0.5, 0.975))
}

AUD_plot_EMPNIGGRApreds_df <- cbind.data.frame(lower = AUD_plot_EMPNIGGRApreds_df[,1], 
                                               mean = AUD_plot_EMPNIGGRApreds_df[,2], 
                                               upper = AUD_plot_EMPNIGGRApreds_df[,3], year = seq(1:23))
#Making graph
AUD_EMPNIG<- 
  ggplot() +
  geom_point(data= AUD_res3, 
             aes(x = YEAR, y = cover, color=TRTMT),
             alpha = 0.8, size = 2)+
  scale_color_manual(values = c("#1B9E77", "#7570B3", "#E6AB02"), labels = c("CTL","GRA", "OTC")) +
  scale_fill_manual(values = c("#1B9E77", "#7570B3", "#E6AB02"),  labels =c("CTL","GRA", "OTC")) +
  scale_x_continuous(breaks = c(1997, 2000, 2007, 2014, 2019)) +
#  geom_ribbon(data = AUD_plot_EMPNIGCTLpreds_df, aes(x = year + 1996, ymin = lower, ymax = upper), fill = "#1B9E77", alpha = 0.2) +
#  geom_line(data = AUD_plot_EMPNIGCTLpreds_df, aes(x = year + 1996, y = mean), colour = "#1B9E77") +
   geom_ribbon(data = AUD_plot_EMPNIGGRApreds_df, aes(x = year + 1996, ymin = lower, ymax = upper), fill = "#7570B3", alpha = 0.2) +
   geom_line(data = AUD_plot_EMPNIGGRApreds_df, aes(x = year + 1996, y = mean), colour = "#7570B3") +
   geom_ribbon(data = AUD_plot_EMPNIGOTCpreds_df, aes(x = year + 1996, ymin = lower, ymax = upper), fill = "#E6AB02", alpha = 0.2) +
   geom_line(data = AUD_plot_EMPNIGOTCpreds_df, aes(x = year + 1996, y = mean), colour = "#E6AB02") +
  ylab("Cover (%)")+
  theme_bw()+
  labs(title= "", subtitle=expression(italic("EMPNIG")))+
  theme(axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.x =element_blank(),
        legend.position = "none",
        text = element_text(size=9),
        plot.title=element_text(size=12),
        plot.margin = unit(c(0.5,0.5,0,0.5), "cm"))

###################Thingsvellir
THIN <- coverCanopy%>%
  filter(SITE=="THINGVELLIR")%>%
  mutate(cover=cover*100)%>%
  na.omit()

# Linear model for all species THIN per treatment
#CTL
THIN_plot_all_CTL <- MCMCglmm(cover ~ I(YEAR-1994) * SPP, 
                              random = ~ YEAR+PLOT, 
                              data = THIN[THIN$TRTMT == "CTL",], 
                              family = "gaussian", pr = TRUE, nitt = 100000, 
                              burnin = 20000, prior = prior2)
summary(THIN_plot_all_CTL)


#OTC
THIN_plot_all_OTC <- MCMCglmm(cover ~ I(YEAR-1994) * SPP, 
                              random = ~ YEAR+PLOT, 
                              data = THIN[THIN$TRTMT == "OTC",], 
                              family = "gaussian", pr = TRUE, nitt = 100000, 
                              burnin = 20000, prior = prior2)
summary(THIN_plot_all_OTC)

#Linear model for main species to check for difference in treatment
#GALBOR
THIN_plot_GALBOR <- MCMCglmm(cover ~ I(YEAR-1994)+TRTMT-1,
                             random = ~ YEAR + PLOT,
                             data = THIN[THIN$SPP == "GALBOR",], 
                             family = "gaussian", pr = TRUE, nitt = 100000, 
                             burnin = 20000, prior = prior2)
summary(THIN_plot_GALBOR)

#FESRUB
THIN_plot_FESRUB <- MCMCglmm(cover ~ I(YEAR-1994)+TRTMT-1,
                             random = ~ YEAR + PLOT,
                             data = THIN[THIN$SPP == "FESRUB",], 
                             family = "gaussian", pr = TRUE, nitt = 100000, 
                             burnin = 20000, prior = prior2)
summary(THIN_plot_FESRUB)


#CARBIG
THIN_plot_CARBIG <- MCMCglmm(cover ~ I(YEAR-1994)+TRTMT-1,
                             random = ~ YEAR + PLOT,
                             data = THIN[THIN$SPP == "CARBIG",], 
                             family = "gaussian", pr = TRUE, nitt = 100000, 
                             burnin = 20000, prior = prior2)
summary(THIN_plot_CARBIG)

#AVEFLE
THIN_plot_AVEFLE <- MCMCglmm(cover ~ I(YEAR-1994)+TRTMT-1,
                             random = ~ YEAR + PLOT,
                             data = THIN[THIN$SPP == "AVEFLE",], 
                             family = "gaussian", pr = TRUE, nitt = 100000, 
                             burnin = 20000, prior = prior2)
summary(THIN_plot_AVEFLE)

##Linear model for main species to plot model predictions
THIN_res1 <- subset(THIN, SPP =="GALBOR")
THIN_res2 <- subset(THIN, SPP =="FESRUB")
THIN_res3 <- subset(THIN, SPP =="CARBIG")
THIN_res4 <- subset(THIN, SPP =="AVEFLE")

#GALBOR
THIN_plot_GALBOR_CTL <- MCMCglmm(cover ~ I(YEAR-1994), random = ~ YEAR + PLOT,
                                 data = THIN_res1[THIN_res1$TRTMT == "CTL",], 
                                 family = "gaussian", pr = TRUE, nitt = 100000, 
                                 burnin = 20000, prior = prior2)
summary(THIN_plot_GALBOR_CTL)

THIN_plot_GALBOR_OTC <- MCMCglmm(cover ~ I(YEAR-1994), random = ~ YEAR + PLOT,
                                 data = THIN_res1[THIN_res1$TRTMT == "OTC",], 
                                 family = "gaussian", pr = TRUE, nitt = 100000, 
                                 burnin = 20000, prior = prior2)
summary(THIN_plot_GALBOR_OTC)

#FESRUB
THIN_plot_FESRUB_CTL <- MCMCglmm(cover ~ I(YEAR-1994), random = ~ YEAR + PLOT,
                                 data = THIN_res2[THIN_res2$TRTMT == "CTL",], 
                                 family = "gaussian", pr = TRUE, nitt = 100000, 
                                 burnin = 20000, prior = prior2)
summary(THIN_plot_FESRUB_CTL)

THIN_plot_FESRUB_OTC <- MCMCglmm(cover ~ I(YEAR-1994), random = ~ YEAR + PLOT,
                                 data = THIN_res2[THIN_res2$TRTMT == "OTC",], 
                                 family = "gaussian", pr = TRUE, nitt = 100000, 
                                 burnin = 20000, prior = prior2)
summary(THIN_plot_FESRUB_OTC)

#CARBIG
THIN_plot_CARBIG_CTL <- MCMCglmm(cover ~ I(YEAR-1994), random = ~ YEAR + PLOT,
                                 data = THIN_res3[THIN_res3$TRTMT == "CTL",], 
                                 family = "gaussian", pr = TRUE, nitt = 100000, 
                                 burnin = 20000, prior = prior2)
summary(THIN_plot_CARBIG_CTL)

THIN_plot_CARBIG_OTC <- MCMCglmm(cover ~ I(YEAR-1994), random = ~ YEAR + PLOT,
                                 data = THIN_res3[THIN_res3$TRTMT == "OTC",], 
                                 family = "gaussian", pr = TRUE, nitt = 100000, 
                                 burnin = 20000, prior = prior2)
summary(THIN_plot_CARBIG_OTC)

#AVEFLE
THIN_plot_AVEFLE_CTL <- MCMCglmm(cover ~ I(YEAR-1994), random = ~ YEAR + PLOT,
                                 data = THIN_res4[THIN_res4$TRTMT == "CTL",], 
                                 family = "gaussian", pr = TRUE, nitt = 100000, 
                                 burnin = 20000, prior = prior2)
summary(THIN_plot_AVEFLE_CTL)

THIN_plot_AVEFLE_OTC <- MCMCglmm(cover ~ I(YEAR-1994), random = ~ YEAR + PLOT,
                                 data = THIN_res4[THIN_res4$TRTMT == "OTC",], 
                                 family = "gaussian", pr = TRUE, nitt = 100000, 
                                 burnin = 20000, prior = prior2)
summary(THIN_plot_AVEFLE_OTC)

#GALBOR
#CTL
nyears <- 25
niter <- length(THIN_plot_GALBOR_CTL$Sol[,"(Intercept)"])

THIN_plot_GALBORCTLpreds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    THIN_plot_GALBORCTLpreds[i,j] <- THIN_plot_GALBOR_CTL$Sol[i,"(Intercept)"] + 
      THIN_plot_GALBOR_CTL$Sol[i,"I(YEAR - 1994)"]*j
  }
}

THIN_plot_GALBORCTLpreds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  THIN_plot_GALBORCTLpreds_df [i,] <- quantile(THIN_plot_GALBORCTLpreds[,i], c(0.025, 0.5, 0.975))
}

THIN_plot_GALBORCTLpreds_df <- cbind.data.frame(lower = THIN_plot_GALBORCTLpreds_df[,1], 
                                                mean = THIN_plot_GALBORCTLpreds_df[,2], 
                                                upper = THIN_plot_GALBORCTLpreds_df[,3], year = seq(1:25))
#OTC
nyears <- 25
niter <- length(THIN_plot_GALBOR_OTC$Sol[,"(Intercept)"])

THIN_plot_GALBOROTCpreds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    THIN_plot_GALBOROTCpreds[i,j] <- THIN_plot_GALBOR_OTC$Sol[i,"(Intercept)"] + 
      THIN_plot_GALBOR_OTC$Sol[i,"I(YEAR - 1994)"]*j
  }
}

THIN_plot_GALBOROTCpreds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  THIN_plot_GALBOROTCpreds_df [i,] <- quantile(THIN_plot_GALBOROTCpreds[,i], c(0.025, 0.5, 0.975))
}

THIN_plot_GALBOROTCpreds_df <- cbind.data.frame(lower = THIN_plot_GALBOROTCpreds_df[,1], 
                                                mean = THIN_plot_GALBOROTCpreds_df[,2], 
                                                upper = THIN_plot_GALBOROTCpreds_df[,3], year = seq(1:25))

# Making GRAPH
THIN_GALBOR<- 
  ggplot() +
  geom_point(data= THIN_res1, 
             aes(x = YEAR, y = cover, color=TRTMT),
             alpha = 0.8, size = 2)+
  scale_color_manual(values = c("#1B9E77",  "#E6AB02"), labels = c("CTL", "OTC")) +
  scale_fill_manual(values = c("#1B9E77", "#E6AB02"),  labels =c("CTL", "OTC")) +
  scale_x_continuous(breaks = c(1995, 1998,2000, 2007, 2014, 2019)) +
  geom_ribbon(data = THIN_plot_GALBORCTLpreds_df, aes(x = year + 1994, ymin = lower, ymax = upper), fill = "#1B9E77", alpha = 0.2) +
  geom_line(data = THIN_plot_GALBORCTLpreds_df, aes(x = year + 1994, y = mean), colour = "#1B9E77") +
  geom_ribbon(data = THIN_plot_GALBOROTCpreds_df, aes(x = year + 1994, ymin = lower, ymax = upper), fill = "#E6AB02", alpha = 0.2) +
  geom_line(data = THIN_plot_GALBOROTCpreds_df, aes(x = year + 1994, y = mean), colour = "#E6AB02") +
  theme_bw()+
  ylab("Cover (%)")+
  xlab("Year")+
  labs(title= "Thingvellir", subtitle=expression(italic("GALBOR")), size=9)+
  theme(legend.position="none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(size=9),
        plot.title=element_text(size=12),
        plot.margin = unit(c(0.5,0.5,0,0.5), "cm"))

#FESRUB
#CTL
nyears <- 25
niter <- length(THIN_plot_FESRUB_CTL$Sol[,"(Intercept)"])

THIN_plot_FESRUBCTLpreds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    THIN_plot_FESRUBCTLpreds[i,j] <- THIN_plot_FESRUB_CTL$Sol[i,"(Intercept)"] + 
      THIN_plot_FESRUB_CTL$Sol[i,"I(YEAR - 1994)"]*j
  }
}

THIN_plot_FESRUBCTLpreds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  THIN_plot_FESRUBCTLpreds_df [i,] <- quantile(THIN_plot_FESRUBCTLpreds[,i], c(0.025, 0.5, 0.975))
}

THIN_plot_FESRUBCTLpreds_df <- cbind.data.frame(lower = THIN_plot_FESRUBCTLpreds_df[,1], 
                                                mean = THIN_plot_FESRUBCTLpreds_df[,2], 
                                                upper = THIN_plot_FESRUBCTLpreds_df[,3], year = seq(1:25))
#OTC
nyears <- 25
niter <- length(THIN_plot_FESRUB_OTC$Sol[,"(Intercept)"])

THIN_plot_FESRUBOTCpreds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    THIN_plot_FESRUBOTCpreds[i,j] <- THIN_plot_FESRUB_OTC$Sol[i,"(Intercept)"] + 
      THIN_plot_FESRUB_OTC$Sol[i,"I(YEAR - 1994)"]*j
  }
}

THIN_plot_FESRUBOTCpreds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  THIN_plot_FESRUBOTCpreds_df [i,] <- quantile(THIN_plot_FESRUBOTCpreds[,i], c(0.025, 0.5, 0.975))
}

THIN_plot_FESRUBOTCpreds_df <- cbind.data.frame(lower = THIN_plot_FESRUBOTCpreds_df[,1], 
                                                mean = THIN_plot_FESRUBOTCpreds_df[,2], 
                                                upper = THIN_plot_FESRUBOTCpreds_df[,3], year = seq(1:25))

# Making GRAPH
THIN_FESRUB<- 
  ggplot() +
  geom_point(data= THIN_res2, 
             aes(x = YEAR, y = cover, color=TRTMT),
             alpha = 0.8, size = 2)+
  scale_color_manual(values = c("#1B9E77",  "#E6AB02"), labels = c("CTL", "OTC")) +
  scale_fill_manual(values = c("#1B9E77", "#E6AB02"),  labels =c("CTL", "OTC")) +
  scale_x_continuous(breaks = c(1995, 1998,2000, 2007, 2014, 2019)) +
 # geom_ribbon(data = THIN_plot_FESRUBCTLpreds_df, aes(x = year + 1994, ymin = lower, ymax = upper), fill = "#1B9E77", alpha = 0.2) +
#  geom_line(data = THIN_plot_FESRUBCTLpreds_df, aes(x = year + 1994, y = mean), colour = "#1B9E77") +
  geom_ribbon(data = THIN_plot_FESRUBOTCpreds_df, aes(x = year + 1994, ymin = lower, ymax = upper), fill = "#E6AB02", alpha = 0.2) +
  geom_line(data = THIN_plot_FESRUBOTCpreds_df, aes(x = year + 1994, y = mean), colour = "#E6AB02") +
  theme_bw()+
  ylab("Cover (%)")+
  xlab("Year")+
  labs(title= "", subtitle=expression(italic("FESRUB")), size=9)+
  theme(legend.position="none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(size=9),
        plot.title=element_text(size=12),
        plot.margin = unit(c(0.5,0.5,0,0.5), "cm"))

#CARBIG
#CTL
nyears <- 25
niter <- length(THIN_plot_CARBIG_CTL$Sol[,"(Intercept)"])

THIN_plot_CARBIGCTLpreds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    THIN_plot_CARBIGCTLpreds[i,j] <- THIN_plot_CARBIG_CTL$Sol[i,"(Intercept)"] + 
      THIN_plot_CARBIG_CTL$Sol[i,"I(YEAR - 1994)"]*j
  }
}

THIN_plot_CARBIGCTLpreds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  THIN_plot_CARBIGCTLpreds_df [i,] <- quantile(THIN_plot_CARBIGCTLpreds[,i], c(0.025, 0.5, 0.975))
}

THIN_plot_CARBIGCTLpreds_df <- cbind.data.frame(lower = THIN_plot_CARBIGCTLpreds_df[,1], 
                                                mean = THIN_plot_CARBIGCTLpreds_df[,2], 
                                                upper = THIN_plot_CARBIGCTLpreds_df[,3], year = seq(1:25))
#OTC
nyears <- 25
niter <- length(THIN_plot_CARBIG_OTC$Sol[,"(Intercept)"])

THIN_plot_CARBIGOTCpreds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    THIN_plot_CARBIGOTCpreds[i,j] <- THIN_plot_CARBIG_OTC$Sol[i,"(Intercept)"] + 
      THIN_plot_CARBIG_OTC$Sol[i,"I(YEAR - 1994)"]*j
  }
}

THIN_plot_CARBIGOTCpreds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  THIN_plot_CARBIGOTCpreds_df [i,] <- quantile(THIN_plot_CARBIGOTCpreds[,i], c(0.025, 0.5, 0.975))
}

THIN_plot_CARBIGOTCpreds_df <- cbind.data.frame(lower = THIN_plot_CARBIGOTCpreds_df[,1], 
                                                mean = THIN_plot_CARBIGOTCpreds_df[,2], 
                                                upper = THIN_plot_CARBIGOTCpreds_df[,3], year = seq(1:25))

# Making GRAPH
THIN_CARBIG<- 
  ggplot() +
  geom_point(data= THIN_res3, 
             aes(x = YEAR, y = cover, color=TRTMT),
             alpha = 0.8, size = 2)+
  scale_color_manual(values = c("#1B9E77",  "#E6AB02"), labels = c("CTL", "OTC")) +
  scale_fill_manual(values = c("#1B9E77", "#E6AB02"),  labels =c("CTL", "OTC")) +
  scale_x_continuous(breaks = c(1995, 1998,2000, 2007, 2014, 2019)) +
  geom_ribbon(data = THIN_plot_CARBIGCTLpreds_df, aes(x = year + 1994, ymin = lower, ymax = upper), fill = "#1B9E77", alpha = 0.2) +
  geom_line(data = THIN_plot_CARBIGCTLpreds_df, aes(x = year + 1994, y = mean), colour = "#1B9E77") +
  geom_ribbon(data = THIN_plot_CARBIGOTCpreds_df, aes(x = year + 1994, ymin = lower, ymax = upper), fill = "#E6AB02", alpha = 0.2) +
  geom_line(data = THIN_plot_CARBIGOTCpreds_df, aes(x = year + 1994, y = mean), colour = "#E6AB02") +
  theme_bw()+
  ylab("Cover (%)")+
  xlab("Year")+
  labs(title= "", subtitle=expression(italic("CARBIG")), size=9)+
  theme(legend.position="none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(size=9),
        plot.title=element_text(size=12),
        plot.margin = unit(c(0.5,0.5,0,0.5), "cm"))

#AVEFLE
#CTL
nyears <- 25
niter <- length(THIN_plot_AVEFLE_CTL$Sol[,"(Intercept)"])

THIN_plot_AVEFLECTLpreds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    THIN_plot_AVEFLECTLpreds[i,j] <- THIN_plot_AVEFLE_CTL$Sol[i,"(Intercept)"] + 
      THIN_plot_AVEFLE_CTL$Sol[i,"I(YEAR - 1994)"]*j
  }
}

THIN_plot_AVEFLECTLpreds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  THIN_plot_AVEFLECTLpreds_df [i,] <- quantile(THIN_plot_AVEFLECTLpreds[,i], c(0.025, 0.5, 0.975))
}

THIN_plot_AVEFLECTLpreds_df <- cbind.data.frame(lower = THIN_plot_AVEFLECTLpreds_df[,1], 
                                                mean = THIN_plot_AVEFLECTLpreds_df[,2], 
                                                upper = THIN_plot_AVEFLECTLpreds_df[,3], year = seq(1:25))
#OTC
nyears <- 25
niter <- length(THIN_plot_AVEFLE_OTC$Sol[,"(Intercept)"])

THIN_plot_AVEFLEOTCpreds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    THIN_plot_AVEFLEOTCpreds[i,j] <- THIN_plot_AVEFLE_OTC$Sol[i,"(Intercept)"] + 
      THIN_plot_AVEFLE_OTC$Sol[i,"I(YEAR - 1994)"]*j
  }
}

THIN_plot_AVEFLEOTCpreds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  THIN_plot_AVEFLEOTCpreds_df [i,] <- quantile(THIN_plot_AVEFLEOTCpreds[,i], c(0.025, 0.5, 0.975))
}

THIN_plot_AVEFLEOTCpreds_df <- cbind.data.frame(lower = THIN_plot_AVEFLEOTCpreds_df[,1], 
                                                mean = THIN_plot_AVEFLEOTCpreds_df[,2], 
                                                upper = THIN_plot_AVEFLEOTCpreds_df[,3], year = seq(1:25))

# Making GRAPH
THIN_AVEFLE<- 
  ggplot() +
  geom_point(data= THIN_res3, 
             aes(x = YEAR, y = cover, color=TRTMT),
             alpha = 0.8, size = 2)+
  scale_color_manual(values = c("#1B9E77",  "#E6AB02"), labels = c("CTL", "OTC")) +
  scale_fill_manual(values = c("#1B9E77", "#E6AB02"),  labels =c("CTL", "OTC")) +
  scale_x_continuous(breaks = c(1995, 1998,2000, 2007, 2014, 2019)) +
  #  geom_ribbon(data = THIN_plot_AVEFLECTLpreds_df, aes(x = year + 1994, ymin = lower, ymax = upper), fill = "#1B9E77", alpha = 0.2) +
  #  geom_line(data = THIN_plot_AVEFLECTLpreds_df, aes(x = year + 1994, y = mean), colour = "#1B9E77") +
  geom_ribbon(data = THIN_plot_AVEFLEOTCpreds_df, aes(x = year + 1994, ymin = lower, ymax = upper), fill = "#E6AB02", alpha = 0.2) +
  geom_line(data = THIN_plot_AVEFLEOTCpreds_df, aes(x = year + 1994, y = mean), colour = "#E6AB02") +
  theme_bw()+
  ylab("Cover (%)")+
  xlab("Year")+
  labs(title= "", subtitle=expression(italic("AVEFLE")), size=9)+
  theme(legend.position="none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(size=9),
        plot.title=element_text(size=12),
        plot.margin = unit(c(0.5,0.5,0,0.5), "cm"))



prow <- plot_grid(AUD_RACLAN, THIN_GALBOR, 
                  AUD_BETNAN, THIN_FESRUB ,
                  AUD_BISVIV,THIN_CARBIG,
                  AUD_EMPNIG,THIN_AVEFLE,
                  align = 'v',
                  hjust = -1,
                  ncol=2,
                  nrow = 4)
legend_b <- get_legend(AUD_RACLAN + theme(legend.position="bottom"))
p<-plot_grid( prow, legend_b, ncol = 1, rel_heights = c(1, 0.1))


##############################Species change with bottom only################

AUD <- coverUnder%>%
  filter(SITE=="AUDKULUHEIDI")%>%
  mutate(cover=cover*100)%>%
  na.omit()

#CTL
AUD_plot_all_CTL <- MCMCglmm(cover ~ I(YEAR-1996) * SPP, 
                             random = ~ YEAR+PLOT, data = AUD[AUD$TRTMT == "CTL",], 
                             family = "gaussian", pr = TRUE, nitt = 100000, 
                             burnin = 20000, prior = prior2)
summary(AUD_plot_all_CTL)

#OTC
AUD_plot_all_OTC <- MCMCglmm(cover ~ I(YEAR-1996) * SPP, 
                             random = ~ YEAR+PLOT, data = AUD[AUD$TRTMT == "OTC",], 
                             family = "gaussian", pr = TRUE, nitt = 100000, 
                             burnin = 20000, prior = prior2)


summary(AUD_plot_all_OTC)

#GRA
AUD_plot_all_GRA <- MCMCglmm(cover ~ I(YEAR-1996) * SPP, 
                             random = ~ YEAR+PLOT, data = AUD[AUD$TRTMT == "GRA",], 
                             family = "gaussian", pr = TRUE, nitt = 100000, 
                             burnin = 20000, prior = prior2)
summary(AUD_plot_all_GRA)

#Linear model for main species to check for difference in treatment
AUD_plot_RACLAN <- MCMCglmm(cover ~ I(YEAR-1996)+TRTMT-1, 
                            random = ~ YEAR + PLOT,
                            data = AUD[AUD$SPP == "RACLAN",], 
                            family = "gaussian", pr = TRUE, nitt = 100000, 
                            burnin = 20000, prior = prior2)
summary(AUD_plot_RACLAN)

AUD_plot_CETISL <- MCMCglmm(cover ~ I(YEAR-1996)+TRTMT-1, 
                            random = ~ YEAR + PLOT,
                            data = AUD[AUD$SPP == "CETISL",], 
                            family = "gaussian", pr = TRUE, nitt = 100000, 
                            burnin = 20000, prior = prior2)
summary(AUD_plot_CETISL)

AUD_plot_CLAARB <- MCMCglmm(cover ~ I(YEAR-1996)+TRTMT-1, 
                            random = ~ YEAR + PLOT,
                            data = AUD[AUD$SPP == "CLAARB",], 
                            family = "gaussian", pr = TRUE, nitt = 100000, 
                            burnin = 20000, prior = prior2)
summary(AUD_plot_CLAARB)

##Linear model for main species to plot model predictions
AUD_res1 <- subset(AUD, SPP =="RACLAN")
AUD_res2 <- subset(AUD, SPP =="CETISL")
AUD_res3 <- subset(AUD, SPP =="CLAARB")

#RACLAN
AUD_plot_RACLAN_CTL <- MCMCglmm(cover ~ I(YEAR-1996), random = ~ YEAR + PLOT,
                                data = AUD_res1[AUD_res1$TRTMT == "CTL",], 
                                family = "gaussian", pr = TRUE, nitt = 100000, 
                                burnin = 20000, prior = prior2)
summary(AUD_plot_RACLAN_CTL)

AUD_plot_RACLAN_OTC <- MCMCglmm(cover ~ I(YEAR-1996), random = ~ YEAR + PLOT,
                                data = AUD_res1[AUD_res1$TRTMT == "OTC",], 
                                family = "gaussian", pr = TRUE, nitt = 100000, 
                                burnin = 20000, prior = prior2)
summary(AUD_plot_RACLAN_OTC)

AUD_plot_RACLAN_GRA <- MCMCglmm(cover ~ I(YEAR-1996), random = ~ YEAR + PLOT,
                                data = AUD_res1[AUD_res1$TRTMT == "GRA",], 
                                family = "gaussian", pr = TRUE, nitt = 100000, 
                                burnin = 20000, prior = prior2)
summary(AUD_plot_RACLAN_GRA)

#CETISL
AUD_plot_CETISL_CTL <- MCMCglmm(cover ~ I(YEAR-1996), random = ~ YEAR + PLOT,
                                data = AUD_res2[AUD_res2$TRTMT == "CTL",], 
                                family = "gaussian", pr = TRUE, nitt = 100000, 
                                burnin = 20000, prior = prior2)
summary(AUD_plot_CETISL_CTL)

AUD_plot_CETISL_OTC <- MCMCglmm(cover ~ I(YEAR-1996), random = ~ YEAR + PLOT,
                                data = AUD_res2[AUD_res2$TRTMT == "OTC",], 
                                family = "gaussian", pr = TRUE, nitt = 100000, 
                                burnin = 20000, prior = prior2)
summary(AUD_plot_CETISL_OTC)

AUD_plot_CETISL_GRA <- MCMCglmm(cover ~ I(YEAR-1996), random = ~ YEAR + PLOT,
                                data = AUD_res2[AUD_res2$TRTMT == "GRA",], 
                                family = "gaussian", pr = TRUE, nitt = 100000, 
                                burnin = 20000, prior = prior2)
summary(AUD_plot_CETISL_GRA)

#CLAARB
AUD_plot_CLAARB_CTL <- MCMCglmm(cover ~ I(YEAR-1996), random = ~ YEAR + PLOT,
                                data = AUD_res3[AUD_res3$TRTMT == "CTL",], 
                                family = "gaussian", pr = TRUE, nitt = 100000, 
                                burnin = 20000, prior = prior2)
summary(AUD_plot_CLAARB_CTL)

AUD_plot_CLAARB_OTC <- MCMCglmm(cover ~ I(YEAR-1996), random = ~ YEAR + PLOT,
                                data = AUD_res3[AUD_res3$TRTMT == "OTC",], 
                                family = "gaussian", pr = TRUE, nitt = 100000, 
                                burnin = 20000, prior = prior2)
summary(AUD_plot_CLAARB_OTC)

AUD_plot_CLAARB_GRA <- MCMCglmm(cover ~ I(YEAR-1996), random = ~ YEAR + PLOT,
                                data = AUD_res3[AUD_res3$TRTMT == "GRA",], 
                                family = "gaussian", pr = TRUE, nitt = 100000, 
                                burnin = 20000, prior = prior2)
summary(AUD_plot_CLAARB_GRA)


# Calculating model predictions
#RACLAN
#CTL
nyears <- 23
niter <- length(AUD_plot_RACLAN_CTL$Sol[,"(Intercept)"])

AUD_plot_RACLANCTLpreds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    AUD_plot_RACLANCTLpreds[i,j] <- AUD_plot_RACLAN_CTL$Sol[i,"(Intercept)"] + AUD_plot_RACLAN_CTL$Sol[i,"I(YEAR - 1996)"]*j
  }
}

AUD_plot_RACLANCTLpreds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  AUD_plot_RACLANCTLpreds_df [i,] <- quantile(AUD_plot_RACLANCTLpreds[,i], c(0.025, 0.5, 0.975))
}

AUD_plot_RACLANCTLpreds_df <- cbind.data.frame(lower = AUD_plot_RACLANCTLpreds_df[,1], 
                                               mean = AUD_plot_RACLANCTLpreds_df[,2], 
                                               upper = AUD_plot_RACLANCTLpreds_df[,3], year = seq(1:23))
#OTC
nyears <- 23
niter <- length(AUD_plot_RACLAN_OTC$Sol[,"(Intercept)"])

AUD_plot_RACLANOTCpreds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    AUD_plot_RACLANOTCpreds[i,j] <- AUD_plot_RACLAN_OTC$Sol[i,"(Intercept)"] + AUD_plot_RACLAN_OTC$Sol[i,"I(YEAR - 1996)"]*j
  }
}

AUD_plot_RACLANOTCpreds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  AUD_plot_RACLANOTCpreds_df [i,] <- quantile(AUD_plot_RACLANOTCpreds[,i], c(0.025, 0.5, 0.975))
}

AUD_plot_RACLANOTCpreds_df <- cbind.data.frame(lower = AUD_plot_RACLANOTCpreds_df[,1], 
                                               mean = AUD_plot_RACLANOTCpreds_df[,2], 
                                               upper = AUD_plot_RACLANOTCpreds_df[,3], year = seq(1:23))

#GRA
nyears <- 23
niter <- length(AUD_plot_RACLAN_GRA$Sol[,"(Intercept)"])

AUD_plot_RACLANGRApreds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    AUD_plot_RACLANGRApreds[i,j] <- AUD_plot_RACLAN_GRA$Sol[i,"(Intercept)"] + AUD_plot_RACLAN_GRA$Sol[i,"I(YEAR - 1996)"]*j
  }
}

AUD_plot_RACLANGRApreds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  AUD_plot_RACLANGRApreds_df [i,] <- quantile(AUD_plot_RACLANGRApreds[,i], c(0.025, 0.5, 0.975))
}

AUD_plot_RACLANGRApreds_df <- cbind.data.frame(lower = AUD_plot_RACLANGRApreds_df[,1], 
                                               mean = AUD_plot_RACLANGRApreds_df[,2], 
                                               upper = AUD_plot_RACLANGRApreds_df[,3], year = seq(1:23))

#Making graph
AUD_RACLAN<- 
  ggplot() +
  geom_point(data= AUD_res1, 
             aes(x = YEAR, y = cover, color=TRTMT),
             alpha = 0.8, size = 2)+
  scale_color_manual(values = c("#1B9E77", "#7570B3", "#E6AB02"), labels = c("CTL","GRA", "OTC")) +
  scale_fill_manual(values = c("#1B9E77", "#7570B3", "#E6AB02"),  labels =c("CTL","GRA", "OTC")) +
  scale_x_continuous(breaks = c(1997, 2000, 2007, 2014, 2019)) +
 # geom_ribbon(data = AUD_plot_RACLANCTLpreds_df, aes(x = year + 1996, ymin = lower, ymax = upper), fill = "#1B9E77", alpha = 0.2) +
#  geom_line(data = AUD_plot_RACLANCTLpreds_df, aes(x = year + 1996, y = mean), colour = "#1B9E77") +
 # geom_ribbon(data = AUD_plot_RACLANGRApreds_df, aes(x = year + 1996, ymin = lower, ymax = upper), fill = "#7570B3", alpha = 0.2) +
#  geom_line(data = AUD_plot_RACLANGRApreds_df, aes(x = year + 1996, y = mean), colour = "#7570B3") +
  geom_ribbon(data = AUD_plot_RACLANOTCpreds_df, aes(x = year + 1996, ymin = lower, ymax = upper), fill = "#E6AB02", alpha = 0.2) +
  geom_line(data = AUD_plot_RACLANOTCpreds_df, aes(x = year + 1996, y = mean), colour = "#E6AB02") +
  theme_bw()+
  ylab("Cover (%)")+
  labs(title="Audkuluheidi",subtitle=expression(italic("Racomitrium lanuginosum")), size=9)+
  theme(axis.text.x=element_blank(),
        axis.title.x=element_blank(),
      #  axis.title.y=element_blank(),
        axis.ticks.x =element_blank(),
        text = element_text(size=9),
        legend.position="none",
        plot.margin = unit(c(0,0.5,0,0.5), "cm"))

#CETISL
#CTL
nyears <- 23
niter <- length(AUD_plot_CETISL_CTL$Sol[,"(Intercept)"])

AUD_plot_CETISLCTLpreds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    AUD_plot_CETISLCTLpreds[i,j] <- AUD_plot_CETISL_CTL$Sol[i,"(Intercept)"] + 
      AUD_plot_CETISL_CTL$Sol[i,"I(YEAR - 1996)"]*j
  }
}

AUD_plot_CETISLCTLpreds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  AUD_plot_CETISLCTLpreds_df [i,] <- quantile(AUD_plot_CETISLCTLpreds[,i], c(0.025, 0.5, 0.975))
}

AUD_plot_CETISLCTLpreds_df <- cbind.data.frame(lower = AUD_plot_CETISLCTLpreds_df[,1], 
                                               mean = AUD_plot_CETISLCTLpreds_df[,2], 
                                               upper = AUD_plot_CETISLCTLpreds_df[,3], year = seq(1:23))
#OTC
nyears <- 23
niter <- length(AUD_plot_CETISL_OTC$Sol[,"(Intercept)"])

AUD_plot_CETISLOTCpreds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    AUD_plot_CETISLOTCpreds[i,j] <- AUD_plot_CETISL_OTC$Sol[i,"(Intercept)"] + 
      AUD_plot_CETISL_OTC$Sol[i,"I(YEAR - 1996)"]*j
  }
}

AUD_plot_CETISLOTCpreds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  AUD_plot_CETISLOTCpreds_df [i,] <- quantile(AUD_plot_CETISLOTCpreds[,i], c(0.025, 0.5, 0.975))
}

AUD_plot_CETISLOTCpreds_df <- cbind.data.frame(lower = AUD_plot_CETISLOTCpreds_df[,1], 
                                               mean = AUD_plot_CETISLOTCpreds_df[,2], 
                                               upper = AUD_plot_CETISLOTCpreds_df[,3], year = seq(1:23))

#GRA
nyears <- 23
niter <- length(AUD_plot_CETISL_GRA$Sol[,"(Intercept)"])

AUD_plot_CETISLGRApreds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    AUD_plot_CETISLGRApreds[i,j] <- AUD_plot_CETISL_GRA$Sol[i,"(Intercept)"] + 
      AUD_plot_CETISL_GRA$Sol[i,"I(YEAR - 1996)"]*j
  }
}

AUD_plot_CETISLGRApreds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  AUD_plot_CETISLGRApreds_df [i,] <- quantile(AUD_plot_CETISLGRApreds[,i], c(0.025, 0.5, 0.975))
}

AUD_plot_CETISLGRApreds_df <- cbind.data.frame(lower = AUD_plot_CETISLGRApreds_df[,1], 
                                               mean = AUD_plot_CETISLGRApreds_df[,2], 
                                               upper = AUD_plot_CETISLGRApreds_df[,3], year = seq(1:23))
#Making graph
AUD_CETISL<- 
  ggplot() +
  geom_point(data= AUD_res2, 
             aes(x = YEAR, y = cover, color=TRTMT),
             alpha = 0.8, size = 2)+
  scale_color_manual(values = c("#1B9E77", "#7570B3", "#E6AB02"), labels = c("CTL","GRA", "OTC")) +
  scale_fill_manual(values = c("#1B9E77", "#7570B3", "#E6AB02"),  labels =c("CTL","GRA", "OTC")) +
  scale_x_continuous(breaks = c(1997, 2000, 2007, 2014, 2019)) +
  geom_ribbon(data = AUD_plot_CETISLCTLpreds_df, aes(x = year + 1996, ymin = lower, ymax = upper), fill = "#1B9E77", alpha = 0.2) +
  geom_line(data = AUD_plot_CETISLCTLpreds_df, aes(x = year + 1996, y = mean), colour = "#1B9E77") +
#  geom_ribbon(data = AUD_plot_CETISLGRApreds_df, aes(x = year + 1996, ymin = lower, ymax = upper), fill = "#7570B3", alpha = 0.2) +
#  geom_line(data = AUD_plot_CETISLGRApreds_df, aes(x = year + 1996, y = mean), colour = "#7570B3") +
#  geom_ribbon(data = AUD_plot_CETISLOTCpreds_df, aes(x = year + 1996, ymin = lower, ymax = upper), fill = "#E6AB02", alpha = 0.2) +
#  geom_line(data = AUD_plot_CETISLOTCpreds_df, aes(x = year + 1996, y = mean), colour = "#E6AB02",  linetype="dashed") +
  ylab("Cover (%)")+
  theme_bw()+
  labs(title= "", subtitle=expression(italic("CETISL")))+
  theme(axis.text.x=element_blank(),
        axis.title.x=element_blank(),
      #  axis.title.y=element_blank(),
        axis.ticks.x =element_blank(),
        legend.position = "none",
        text = element_text(size=9),
        plot.title=element_text(size=12),
        plot.margin = unit(c(0.5,0.5,0,0.5), "cm"))

#CLAARB
#CTL
nyears <- 23
niter <- length(AUD_plot_CLAARB_CTL$Sol[,"(Intercept)"])

AUD_plot_CLAARBCTLpreds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    AUD_plot_CLAARBCTLpreds[i,j] <- AUD_plot_CLAARB_CTL$Sol[i,"(Intercept)"] + 
      AUD_plot_CLAARB_CTL$Sol[i,"I(YEAR - 1996)"]*j
  }
}

AUD_plot_CLAARBCTLpreds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  AUD_plot_CLAARBCTLpreds_df [i,] <- quantile(AUD_plot_CLAARBCTLpreds[,i], c(0.025, 0.5, 0.975))
}

AUD_plot_CLAARBCTLpreds_df <- cbind.data.frame(lower = AUD_plot_CLAARBCTLpreds_df[,1], 
                                               mean = AUD_plot_CLAARBCTLpreds_df[,2], 
                                               upper = AUD_plot_CLAARBCTLpreds_df[,3], year = seq(1:23))
#OTC
nyears <- 23
niter <- length(AUD_plot_CLAARB_OTC$Sol[,"(Intercept)"])

AUD_plot_CLAARBOTCpreds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    AUD_plot_CLAARBOTCpreds[i,j] <- AUD_plot_CLAARB_OTC$Sol[i,"(Intercept)"] + 
      AUD_plot_CLAARB_OTC$Sol[i,"I(YEAR - 1996)"]*j
  }
}

AUD_plot_CLAARBOTCpreds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  AUD_plot_CLAARBOTCpreds_df [i,] <- quantile(AUD_plot_CLAARBOTCpreds[,i], c(0.025, 0.5, 0.975))
}

AUD_plot_CLAARBOTCpreds_df <- cbind.data.frame(lower = AUD_plot_CLAARBOTCpreds_df[,1], 
                                               mean = AUD_plot_CLAARBOTCpreds_df[,2], 
                                               upper = AUD_plot_CLAARBOTCpreds_df[,3], year = seq(1:23))

#GRA
nyears <- 23
niter <- length(AUD_plot_CLAARB_GRA$Sol[,"(Intercept)"])

AUD_plot_CLAARBGRApreds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    AUD_plot_CLAARBGRApreds[i,j] <- AUD_plot_CLAARB_GRA$Sol[i,"(Intercept)"] + 
      AUD_plot_CLAARB_GRA$Sol[i,"I(YEAR - 1996)"]*j
  }
}

AUD_plot_CLAARBGRApreds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  AUD_plot_CLAARBGRApreds_df [i,] <- quantile(AUD_plot_CLAARBGRApreds[,i], c(0.025, 0.5, 0.975))
}

AUD_plot_CLAARBGRApreds_df <- cbind.data.frame(lower = AUD_plot_CLAARBGRApreds_df[,1], 
                                               mean = AUD_plot_CLAARBGRApreds_df[,2], 
                                               upper = AUD_plot_CLAARBGRApreds_df[,3], year = seq(1:23))
#Making graph
AUD_CLAARB<- 
  ggplot() +
  geom_point(data= AUD_res3, 
             aes(x = YEAR, y = cover, color=TRTMT),
             alpha = 0.8, size = 2)+
  scale_color_manual(values = c("#1B9E77", "#7570B3", "#E6AB02"), labels = c("CTL","GRA", "OTC")) +
  scale_fill_manual(values = c("#1B9E77", "#7570B3", "#E6AB02"),  labels =c("CTL","GRA", "OTC")) +
  scale_x_continuous(breaks = c(1997, 2000, 2007, 2014, 2019)) +
 # geom_ribbon(data = AUD_plot_CLAARBCTLpreds_df, aes(x = year + 1996, ymin = lower, ymax = upper), fill = "#1B9E77", alpha = 0.2) +
 # geom_line(data = AUD_plot_CLAARBCTLpreds_df, aes(x = year + 1996, y = mean), colour = "#1B9E77") +
  geom_ribbon(data = AUD_plot_CLAARBGRApreds_df, aes(x = year + 1996, ymin = lower, ymax = upper), fill = "#7570B3", alpha = 0.2) +
  geom_line(data = AUD_plot_CLAARBGRApreds_df, aes(x = year + 1996, y = mean), colour = "#7570B3") +
  geom_ribbon(data = AUD_plot_CLAARBOTCpreds_df, aes(x = year + 1996, ymin = lower, ymax = upper), fill = "#E6AB02", alpha = 0.2) +
  geom_line(data = AUD_plot_CLAARBOTCpreds_df, aes(x = year + 1996, y = mean), colour = "#E6AB02") +
  ylab("Cover (%)")+
  theme_bw()+
  labs(title= "", subtitle=expression(italic("CLAARB")))+
  theme(legend.position = "none",
        text = element_text(size=9),
        plot.title=element_text(size=12),
        plot.margin = unit(c(0.5,0.5,0,0.5), "cm"))



###################Thingsvellir
THIN <- coverUnder%>%
  filter(SITE=="THINGVELLIR")%>%
  mutate(cover=cover*100)%>%
  na.omit()

# Linear model for all species THIN per treatment
#CTL
THIN_plot_all_CTL <- MCMCglmm(cover ~ I(YEAR-1994) * SPP, 
                              random = ~ YEAR+PLOT, 
                              data = THIN[THIN$TRTMT == "CTL",], 
                              family = "gaussian", pr = TRUE, nitt = 100000, 
                              burnin = 20000, prior = prior2)
summary(THIN_plot_all_CTL)

#OTC
THIN_plot_all_OTC <- MCMCglmm(cover ~ I(YEAR-1994) * SPP, 
                              random = ~ YEAR+PLOT, 
                              data = THIN[THIN$TRTMT == "OTC",], 
                              family = "gaussian", pr = TRUE, nitt = 100000, 
                              burnin = 20000, prior = prior2)
summary(THIN_plot_all_OTC)

#Linear model for main species to check for difference in treatment
#RACLAN
THIN_plot_RACLAN <- MCMCglmm(cover ~ I(YEAR-1994)+TRTMT-1,
                             random = ~ YEAR + PLOT,
                             data = THIN[THIN$SPP == "RACLAN",], 
                             family = "gaussian", pr = TRUE, nitt = 100000, 
                             burnin = 20000, prior = prior2)
summary(THIN_plot_RACLAN)

#CETISL
THIN_plot_CETISL <- MCMCglmm(cover ~ I(YEAR-1994)+TRTMT-1,
                             random = ~ YEAR + PLOT,
                             data = THIN[THIN$SPP == "CETISL",], 
                             family = "gaussian", pr = TRUE, nitt = 100000, 
                             burnin = 20000, prior = prior2)
summary(THIN_plot_CETISL)


##Linear model for main species to plot model predictions
THIN_res1 <- subset(THIN, SPP =="RACLAN")
THIN_res2 <- subset(THIN, SPP =="CETISL")

#RACLAN
THIN_plot_RACLAN_CTL <- MCMCglmm(cover ~ I(YEAR-1994), random = ~ YEAR + PLOT,
                                 data = THIN_res1[THIN_res1$TRTMT == "CTL",], 
                                 family = "gaussian", pr = TRUE, nitt = 100000, 
                                 burnin = 20000, prior = prior2)
summary(THIN_plot_RACLAN_CTL)

THIN_plot_RACLAN_OTC <- MCMCglmm(cover ~ I(YEAR-1994), random = ~ YEAR + PLOT,
                                 data = THIN_res1[THIN_res1$TRTMT == "OTC",], 
                                 family = "gaussian", pr = TRUE, nitt = 100000, 
                                 burnin = 20000, prior = prior2)
summary(THIN_plot_RACLAN_OTC)

#CETISL
THIN_plot_CETISL_CTL <- MCMCglmm(cover ~ I(YEAR-1994), random = ~ YEAR + PLOT,
                                 data = THIN_res2[THIN_res2$TRTMT == "CTL",], 
                                 family = "gaussian", pr = TRUE, nitt = 100000, 
                                 burnin = 20000, prior = prior2)
summary(THIN_plot_CETISL_CTL)

THIN_plot_CETISL_OTC <- MCMCglmm(cover ~ I(YEAR-1994), random = ~ YEAR + PLOT,
                                 data = THIN_res2[THIN_res2$TRTMT == "OTC",], 
                                 family = "gaussian", pr = TRUE, nitt = 100000, 
                                 burnin = 20000, prior = prior2)
summary(THIN_plot_CETISL_OTC)


#RACLAN
#CTL
nyears <- 25
niter <- length(THIN_plot_RACLAN_CTL$Sol[,"(Intercept)"])

THIN_plot_RACLANCTLpreds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    THIN_plot_RACLANCTLpreds[i,j] <- THIN_plot_RACLAN_CTL$Sol[i,"(Intercept)"] + 
      THIN_plot_RACLAN_CTL$Sol[i,"I(YEAR - 1994)"]*j
  }
}

THIN_plot_RACLANCTLpreds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  THIN_plot_RACLANCTLpreds_df [i,] <- quantile(THIN_plot_RACLANCTLpreds[,i], c(0.025, 0.5, 0.975))
}

THIN_plot_RACLANCTLpreds_df <- cbind.data.frame(lower = THIN_plot_RACLANCTLpreds_df[,1], 
                                                mean = THIN_plot_RACLANCTLpreds_df[,2], 
                                                upper = THIN_plot_RACLANCTLpreds_df[,3], year = seq(1:25))
#OTC
nyears <- 25
niter <- length(THIN_plot_RACLAN_OTC$Sol[,"(Intercept)"])

THIN_plot_RACLANOTCpreds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    THIN_plot_RACLANOTCpreds[i,j] <- THIN_plot_RACLAN_OTC$Sol[i,"(Intercept)"] + 
      THIN_plot_RACLAN_OTC$Sol[i,"I(YEAR - 1994)"]*j
  }
}

THIN_plot_RACLANOTCpreds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  THIN_plot_RACLANOTCpreds_df [i,] <- quantile(THIN_plot_RACLANOTCpreds[,i], c(0.025, 0.5, 0.975))
}

THIN_plot_RACLANOTCpreds_df <- cbind.data.frame(lower = THIN_plot_RACLANOTCpreds_df[,1], 
                                                mean = THIN_plot_RACLANOTCpreds_df[,2], 
                                                upper = THIN_plot_RACLANOTCpreds_df[,3], year = seq(1:25))

# Making GRAPH
THIN_RACLAN<- 
  ggplot() +
  geom_point(data= THIN_res1, 
             aes(x = YEAR, y = cover, color=TRTMT),
             alpha = 0.8, size = 2)+
  scale_color_manual(values = c("#1B9E77",  "#E6AB02"), labels = c("CTL", "OTC")) +
  scale_fill_manual(values = c("#1B9E77", "#E6AB02"),  labels =c("CTL", "OTC")) +
  scale_x_continuous(breaks = c(1995, 1998,2000, 2007, 2014, 2019)) +
 # geom_ribbon(data = THIN_plot_RACLANCTLpreds_df, aes(x = year + 1994, ymin = lower, ymax = upper), fill = "#1B9E77", alpha = 0.2) +
 # geom_line(data = THIN_plot_RACLANCTLpreds_df, aes(x = year + 1994, y = mean), colour = "#1B9E77") +
  geom_ribbon(data = THIN_plot_RACLANOTCpreds_df, aes(x = year + 1994, ymin = lower, ymax = upper), fill = "#E6AB02", alpha = 0.2) +
  geom_line(data = THIN_plot_RACLANOTCpreds_df, aes(x = year + 1994, y = mean), colour = "#E6AB02") +
  theme_bw()+
  ylab("Cover (%)")+
  xlab("Year")+
  labs(title= "Thingsvellir", subtitle=expression(italic("Racomitrium lanuginosum")), size=9)+
  theme(legend.position="none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(size=9),
        plot.title=element_text(size=12),
        plot.margin = unit(c(0.5,0.5,0,0.5), "cm"))

#CETISL
#CTL
nyears <- 25
niter <- length(THIN_plot_CETISL_CTL$Sol[,"(Intercept)"])

THIN_plot_CETISLCTLpreds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    THIN_plot_CETISLCTLpreds[i,j] <- THIN_plot_CETISL_CTL$Sol[i,"(Intercept)"] + 
      THIN_plot_CETISL_CTL$Sol[i,"I(YEAR - 1994)"]*j
  }
}

THIN_plot_CETISLCTLpreds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  THIN_plot_CETISLCTLpreds_df [i,] <- quantile(THIN_plot_CETISLCTLpreds[,i], c(0.025, 0.5, 0.975))
}

THIN_plot_CETISLCTLpreds_df <- cbind.data.frame(lower = THIN_plot_CETISLCTLpreds_df[,1], 
                                                mean = THIN_plot_CETISLCTLpreds_df[,2], 
                                                upper = THIN_plot_CETISLCTLpreds_df[,3], year = seq(1:25))
#OTC
nyears <- 25
niter <- length(THIN_plot_CETISL_OTC$Sol[,"(Intercept)"])

THIN_plot_CETISLOTCpreds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    THIN_plot_CETISLOTCpreds[i,j] <- THIN_plot_CETISL_OTC$Sol[i,"(Intercept)"] + 
      THIN_plot_CETISL_OTC$Sol[i,"I(YEAR - 1994)"]*j
  }
}

THIN_plot_CETISLOTCpreds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  THIN_plot_CETISLOTCpreds_df [i,] <- quantile(THIN_plot_CETISLOTCpreds[,i], c(0.025, 0.5, 0.975))
}

THIN_plot_CETISLOTCpreds_df <- cbind.data.frame(lower = THIN_plot_CETISLOTCpreds_df[,1], 
                                                mean = THIN_plot_CETISLOTCpreds_df[,2], 
                                                upper = THIN_plot_CETISLOTCpreds_df[,3], year = seq(1:25))

# Making GRAPH
THIN_CETISL<- 
  ggplot() +
  geom_point(data= THIN_res2, 
             aes(x = YEAR, y = cover, color=TRTMT),
             alpha = 0.8, size = 2)+
  scale_color_manual(values = c("#1B9E77",  "#E6AB02"), labels = c("CTL", "OTC")) +
  scale_fill_manual(values = c("#1B9E77", "#E6AB02"),  labels =c("CTL", "OTC")) +
  scale_x_continuous(breaks = c(1995, 1998,2000, 2007, 2014, 2019)) +
 # geom_ribbon(data = THIN_plot_CETISLCTLpreds_df, aes(x = year + 1994, ymin = lower, ymax = upper), fill = "#1B9E77", alpha = 0.2) +
#  geom_line(data = THIN_plot_CETISLCTLpreds_df, aes(x = year + 1994, y = mean), colour = "#1B9E77") +
  # geom_ribbon(data = THIN_plot_CETISLOTCpreds_df, aes(x = year + 1994, ymin = lower, ymax = upper), fill = "#E6AB02", alpha = 0.2) +
  #  geom_line(data = THIN_plot_CETISLOTCpreds_df, aes(x = year + 1994, y = mean), colour = "#E6AB02") +
  theme_bw()+
  ylab("Cover (%)")+
  xlab("Year")+
  labs(title= "", subtitle=expression(italic("CETISL")), size=9)+
  theme(legend.position="none",
        axis.ticks.x = element_blank(),
        text = element_text(size=9),
        plot.title=element_text(size=12),
        plot.margin = unit(c(0.5,0.5,0,0.5), "cm"))



prow <- plot_grid(AUD_RACLAN, THIN_RACLAN, 
                  AUD_CETISL, THIN_CETISL,
                  align = 'v',
                  hjust = -1,
                  ncol=2,
                  nrow = 2)

legend_b <- get_legend(AUD_RACLAN + theme(legend.position="bottom"))
p<-plot_grid( prow, legend_b, ncol = 1, rel_heights = c(1, 0.1))
