setwd("~/Desktop/Garibaldi_Analysis-1")

library(ggplot2)
library(plyr)
library(tidyr)
library(forcats)
library(gridExtra)
library(MCMCglmm)
library(readxl)
library(vegan)
library(cowplot)

##################################################

#Calculate changes in vegetation
dataset <-read_xlsx("/Users/lizmccleary/Desktop/Garibaldi_Analysis-1/ITEX_merged_2025.xlsx")%>% 
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

##Relative cover percentage and change in cover per species 
percent_cover_summary <- coverCanopy %>%
  dplyr::group_by(DATE, SPECIES) %>%
  dplyr::summarize(TotalPercentCover = sum(cover, na.rm = TRUE))

# View the summary
print(percent_cover_summary)

#Use summary to identify which species appear in one year but not in the other
species_2022 <- c("antalp", "caraqu", "caraur", "carnig", "carspp", "casmer", "casspp", "equarv", "equvar", 
                  "galhum", "junarc", "jundru", "junmer", "kalmic", "litter", "lupspp", "moss", "other", 
                  "phygla", "pinvul", "salbar", "salgla", "triglu")

species_2024 <- c("agrhum", "bare ground", "caraur", "carlent", "carnig", "carspec", "carsty", "casmer", 
                  "chalat", "crust", "dead wood", "equarv", "equvar", "eriper", "fungi", "gauhum", "junarc", 
                  "jundru", "junmer", "kalmic", "lichen", "litter", "lupspp", "moss", "phyemp", "phygla", 
                  "pinvul", "poagla", "poaspp", "rock", "salarc", "salbar", "salcom", "salniv", "salsit", 
                  "soil", "solmul", "tomst", "triglu")

#Species in 2022 but not in 2024
species_in_2022_not_2024 <- setdiff(species_2022, species_2024)
print("Species in 2022 but not in 2024:")
print(species_in_2022_not_2024)

#Species in 2024 but not 2022
species_in_2024_not_2022 <- setdiff(species_2024, species_2022)
print("Species in 2024 but not in 2022:")
print(species_in_2024_not_2022)


ggplot(percent_cover_summary, aes(x = DATE, y = TotalPercentCover, color = SPECIES)) +
  geom_line() +
  geom_point() +
  labs(title = "Percent Cover Changes per Species Over Time",
       x = "Year", 
       y = "Total Percent Cover",
       color = "Species") +
  theme_minimal()

#####################NMDS changes in species cover
#Make wide format

SITES<- coverTopBottom %>% 
        filter(SITE %in% c("Salix", "Cassiope","Meadow"))%>% 
        filter(!SPECIES %in% "carspp")%>%
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


####CHANGE IN HEIGHT

dataset <-read_xlsx("ITEX_merged_2025.xlsx")

Height_ALL <- dataset %>%
  filter(SITE %in% c("Salix", "Meadow", "Cassiope")) %>%
  filter(HitOrder %in% c("Top")) %>%
  select(DATE, SITE, TRTMT, PLOT, SPECIES, HitOrder, CanopyHeight.mm.)

#Change PLOT to factor for randomness in MCMC
Height_ALL$PLOT <- as.factor(Height_ALL$PLOT)
Height_ALL$TRTMT <- as.factor(Height_ALL$TRTMT)

#Specify dplyr package for group_by and summarise functions 
Height_ALL <- Height_ALL %>%
  dplyr::group_by(DATE, TRTMT, PLOT, SITE) %>%
  dplyr::summarize(mean_height = mean(CanopyHeight.mm., na.rm = TRUE)) %>%
  ungroup()

#Convert to data frame before passing to MCMC 
Height_ALL_df <- as.data.frame(Height_ALL)

#prior2 <- list(R = list(V = 1, nu = 0.002), 
               #G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.v = 10000), 
                        #G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.v = 10000)))

#Using a prior without random variables 
prior2 <- list(
  R = list(V = 1, nu = 0.002)) #,
  #G = list(
    #G1 = list(V = 1, nu = 2, alpha.mu = 0, alpha.v = 10000),
            #G2 = list(V = 1, nu = 2, alpha.mu = 0, alpha.v = 10000)))


#Warning message, "Unknown or unitialized column: 'family'" --> ran with only one or none random variables
ALL_Height_OTC <- MCMCglmm(mean_height ~ I(DATE - 2022), 
                             #random = ~ DATE + PLOT, 
                             data = Height_ALL_df[Height_ALL_df$TRTMT == "warming",], 
                             family="gaussian",
                             pr = TRUE, 
                             nitt = 100000, 
                             burnin = 20000, 
                             prior = prior2)

ALL_Height_CTL <- MCMCglmm(mean_height ~ I(DATE-2022), 
                           #random = ~ DATE+PLOT, 
                           data = Height_ALL_df[Height_ALL_df$TRTMT == "control",], 
                           family ="gaussian", 
                           pr = TRUE, 
                           nitt = 100000, 
                           burnin = 20000,
                           prior = prior2)


#Test difference between treatment
canopy_m <- MCMCglmm(mean_height ~ I(DATE-2022)+TRTMT-1, 
                            #random = ~ DATE + PLOT,
                            data = Height_ALL_df, 
                            family = "gaussian", pr = TRUE, nitt = 100000, 
                            burnin = 20000, prior = prior2)



# Calculating model predictions

#CTL
nyears <- 23
niter <- length(ALL_Height_CTL$Sol[,"(Intercept)"])

ALL_Height_CTLpreds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
   ALL_Height_CTLpreds[i,j] <- ALL_Height_CTL$Sol[i,"(Intercept)"] + ALL_Height_CTL$Sol[i,"I(DATE - 2022)"]*j
  }
}


ALL_Height_CTLpreds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  ALL_Height_CTLpreds_df [i,] <- quantile(ALL_Height_CTLpreds[,i], c(0.025, 0.5, 0.975))
}

ALL_Height_CTLpreds_df <- cbind.data.frame(lower = ALL_Height_CTLpreds_df[,1], 
                                               mean = ALL_Height_CTLpreds_df[,2], 
                                               upper = ALL_Height_CTLpreds_df[,3], year = seq(1:23))
#OTC
nyears <- 23
niter <- length(ALL_Height_OTC$Sol[,"(Intercept)"])

ALL_Height_OTCpreds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    ALL_Height_OTCpreds[i,j] <- ALL_Height_OTC$Sol[i,"(Intercept)"] + ALL_Height_OTC$Sol[i,"I(DATE - 2022)"]*j
  }
}


ALL_Height_OTCpreds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  ALL_Height_OTCpreds_df [i,] <- quantile(ALL_Height_OTCpreds[,i], c(0.025, 0.5, 0.975))
}

ALL_Height_OTCpreds_df <- cbind.data.frame(lower = ALL_Height_OTCpreds_df[,1], 
                                               mean = ALL_Height_OTCpreds_df[,2], 
                                               upper = ALL_Height_OTCpreds_df[,3], year = seq(1:23))

###Plot change in height
ggplot() +
  geom_point(data= Height_ALL, 
             aes(x = DATE, y = mean_height, color=TRTMT),
             alpha = 0.8, size = 2)+
  scale_color_manual(values = c("#1B9E77", "#7570B3", "#E6AB02"), labels = c("control","warming")) +
  scale_fill_manual(values = c("#1B9E77", "#7570B3", "#E6AB02"),  labels =c("control","warming")) +
  scale_x_continuous(breaks = c(2022, 2024)) +
  geom_ribbon(data = ALL_Height_CTLpreds_df, aes(x = year + 2022, ymin = lower, ymax = upper), fill = "#1B9E77", alpha = 0.2) +
  geom_line(data = ALL_Height_CTLpreds_df, aes(x = year + 2022, y = mean), colour = "#1B9E77") +
  geom_ribbon(data = ALL_Height_OTCpreds_df, aes(x = year + 2022, ymin = lower, ymax = upper), fill = "#7570B3", alpha = 0.2) +
  geom_line(data = ALL_Height_OTCpreds_df, aes(x = year + 2022, y = mean), colour = "#7570B3",  linetype="dashed") +
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
ggplot(data=Height_ALL, aes(x=factor(SITE), y=mean_height,fill = factor(TRTMT))) + 
  geom_boxplot() +
  scale_fill_brewer("Treatment", palette="Dark2")+
  facet_wrap(~DATE)+
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

##ANOVA to determine significance of change between years and warming/control 

# Two-way ANOVA to test for interaction between Year and Treatment
anova_result <- aov(mean_height ~ DATE * TRTMT, data = Height_ALL)

# Summary of the ANOVA result
summary(anova_result)
