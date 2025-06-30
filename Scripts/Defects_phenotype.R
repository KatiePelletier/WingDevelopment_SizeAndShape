#Jun 2025
#returning to this and cleaning up the analysis. 
#The goal of this analysis is to look at the relationship between what I called wing defects (abnnormal veins, for the most part) and wing shape/size. I also want to look at the interaction with cell size (do biger cells lead to more wing defects?)

library(tidyverse)
library(lme4)
library(car)
library(emmeans)
library(effects)
library(broom)
library(geomorph)

#reading in the data. In three diffrent files. I will read in a clip then combine the data 
#this is just the phenotype of those wings. Genotype, sex, and scoring for 4 phenotypes (bianary): pcvl, acvl, extra vein, loss/shorter logintudinal vein. 
pure <- read.csv("../Data/KP_pureline_defects.csv")
pure_clean <- pure[,c(2, 4, 6:9)]

zi418ef43 <- read.csv("../Data/wingdefext_zi418ef43.csv")
zi418ef43_clean <- zi418ef43[c(4, 5, 8:10)]
zi418ef43_clean$line <- "zi418ef43"
  
zi251ef43 <- read.csv("../Data/zi251ef43defects.csv")
zi251ef43_clean <- zi251ef43[,c(4, 5, 8:10)]
zi251ef43_clean$line <- "zi251ef43"

wings <- rbind(pure_clean, zi418ef43_clean, zi251ef43_clean)

#to make variable that scores for any type of defect. Also bianary. 
wings$defect <- with(wings, 
                     ifelse((pcvl + acvl + log.vein.loss + extra.vein ) >= 1, 1, 0))

#Now I want to get the propotion of defects for each sex*line

wings_table <- with(wings, 
                  table(defect, 
                        interaction(line, sex, drop = TRUE, sep = "_") ))


wings_table <- t(wings_table) # transpose table
wings_table

wings_names <- rownames(wings_table)

wings_names <- str_split(wings_names, pattern = "_", simplify = TRUE)

wings2 <- data.frame(wings_names, wings_table[,1], wings_table[,2])
colnames(wings2) <- c("line", "sex", "normal", "defect")

wings2 <- mutate(wings2, Freq = defect/(defect+normal))
wings2

wings2$line <- factor(wings2$line, c("ef43", "ef81", "ef96", "zi192", "zi251", "zi418", "zi251ef43", "zi418ef43"))




#Now to make a plot! 

defects.plot <- ggplot(wings2, aes(x = line, y = Freq, fill = sex)) + 
  geom_col(position = "dodge") + 
  theme_bw(base_size = 18) + 
  xlab("Genotype") + 
  ylab("Frequency of Wing Defect") +
  scale_fill_grey(labels = c("Female", "Male"))

png("../Figures/freqdefects_autoFreqCutoff.png", res = 100, units = "px", width=1060, height=412)
defects.plot 
dev.off()

png("../Figures/freqdefects_FreqTo1.png", res = 100, units = "px", width=1060, height=412)
ggplot(wings2, aes(x = line, y = Freq, fill = sex)) + 
  geom_col(position = "dodge") + 
  theme_bw(base_size = 18) + 
  xlab("Genotype") + 
  ylim(0,1) + 
  ylab("Frequency of Wing Defect") +
  scale_fill_grey(labels = c("Female", "Male"))
dev.off()

  
#I don't have these images on my work computer and I want to have all the parts before I start making full figures 
# defects.wings.raw <- image_read("../Figures/wingDefectsFig.png")
# defects.wings.crop <- image_trim(defects.wings.raw)
# defects.wings <- ggdraw() + draw_image(defects.wings.crop)
# 
# png("../Figures/defectsSupFig.png", width = 3500, height = 1500, units = "px", res = 300)
# plot_grid(defects.wings, defects.plot, nrow = 1, rel_heights = c(1, 0.6))
# dev.off()

#Now to do some stats. Ovi, I need to add the shape and tricome data to this

zi418ef43_wings <- read.csv("../Data/zi418ef43_tricomesAndShape.csv")

#this is full splines. 
str(zi418ef43_wings)

#7 = tricomes, 109 = CS, 111 = size, 13:108 = shape
#I want to take the average of tricome count over the wing, but also to get a measure of variance as well. For this, I will use Levene's statistic, following the logic of the dll and MP papers. I think this still holds for looking within a wing? 

#I need this for later 
zi418ef43$flyID <-  paste(paste0(zi418ef43$plate, zi418ef43$row), 
                          zi418ef43$fly, sep = "_")


zi418ef43.catDat <- (zi418ef43_wings %>% 
                       group_by(sex, plate, Fly_ID) %>% 
                       #get the measures I want for every line 
                       summarize(tricome_mu = mean(tricomeCount, na.rm = TRUE),
                              tricome_sd = sd(tricomeCount, na.rm = TRUE), 
                              tricome_cv = tricome_sd/tricome_mu) %>% 
                       #don't need the grouping variable and makes sure there is no mistake later
                       ungroup() %>%
                       #ugly code, but working with what I have 
                       #adding back the shape and size data for the wings.
                       #because the data is identical and repeated for every observation, I can just take                         any match.
                       left_join(zi418ef43_wings, multiple = "any") %>%
                       #dropping the tricome count and wing region col because this doesn't reflect reality
                       dplyr::select(-c("tricomeCount", "wingRegion")) %>% 
                       #annoying data cleanup to match the wing defect flyID 
                       mutate(flyID = paste(plate, fly, sep = "_")) %>% 
                       #joining in the defect data 
                       left_join(zi418ef43, by = "flyID") %>% 
                       #and finally to remake this from above. This part of the code could be much cleaner.
                       #also making a log and corrected CS variable.
                       mutate(defect = 
                                if_else((pcvl + acvl + log.vein.loss + extra.vein ) >= 1, 1, 0), 
                              logCS = log(CS, base = 2), 
                              logCSc = logCS - mean(logCS), 
                              tricomeC = tricome_mu - mean(tricome_mu))
                       )
  
hist(zi418ef43.catDat$tricome_mu)
hist(zi418ef43.catDat$tricome_sd)
hist(zi418ef43.catDat$tricome_cv)
hist(zi418ef43.catDat$tricomeC)


#I forgot about the hump here, but that is real in the data, not something weird that I've done. 
hist(zi418ef43.catDat$CS)
hist(zi418ef43.catDat$logCS)
hist(zi418ef43.catDat$logCSc)

#Now for some models. 


#First, do bigger overall wings have more vein defects?
defect.size.mod <- glm(defect ~ logCSc*sex.x, 
                       data = zi418ef43.catDat, family = binomial(link = "logit"))

#males have slighly more defects and the significance is borderline 
summary(defect.size.mod)
Anova(defect.size.mod)


plot(emtrends(defect.size.mod, ~ sex.x, var = "logCSc"))


#to for the males being more likley than females to have a defect 
#the log odds ratio is pretty small, and the conf int crosses 0. 
pairs(emmeans(defect.size.mod, ~sex.x))
tidy(pairs(emmeans(defect.size.mod, ~sex.x)), conf.int = TRUE)


#What scale is this on? Looks log to me so probably the latent scale. 
#Not going to play a lot with this right now. 
plot(predictorEffect("logCSc", defect.size.mod))

#a note for later on how to get this quickplotting to work
plot(predictorEffect("logCSc", defect.size.mod), lines=list(multiline=TRUE))

#putting it back on the response scale. 
plot(predictorEffect("logCSc", defect.size.mod), axes=list(y=list(type="response", grid=TRUE)))

plot(predictorEffect("logCSc", defect.size.mod), axes=list(y=list(type="response", grid=TRUE)), 
     lines=list(multiline=TRUE))

#Now I want to ask about the cell size and defects 

#this is actually a better model than before because it accounts for the fact that larger wings could just be explained by more cells. 
defect.cellsize.mod <- glm(defect ~ (logCSc + sex.x + tricomeC)^2, 
                       data = zi418ef43.catDat, family = binomial(link = "logit"))


#now males have a bigger effect, there is also the interaction between tricome count and males, which is in the opposite direction. Wing size still doesn't seem to have an effect. 
summary(defect.cellsize.mod)

Anova(defect.cellsize.mod)

#The effect is oppisite for M and F?
plot(predictorEffect("tricomeC", defect.cellsize.mod, cov = "logCSc"))

plot(predictorEffects(defect.cellsize.mod))

#trying to make it less overwhelming 
#A little easier to look at.

#still really hard to think though. 

plot(predictorEffects(defect.cellsize.mod, ~ logCSc + tricomeC, 
                      xlevels = list(logCSc = c(-0.4, -0.09, 0.09, 0.2), 
                                     tricomeC = c(c(-30, -3, 9, 20)))),
     axes=list(grid=TRUE,
               x=list(rug=FALSE),
               y=list(type="response")),
     lines=list(multiline=TRUE)
               )


#This is a cool plot, but I gave no idea what to make of it really. 
#Having larger than overage wings and smaller than average tricome count increases the proportion of predicted defects in the wing. Not true for females (if anything its the opposite pattern)  

#I think this is a terrible way to look at this because there is no data for half of these measurements. 
#Should I do a sex based mean? That would give larger than average females and larger than average males? 
#I don't actually care that the mean of the two sexes is diffrent, just want to know if bigger than average wings tend to have more defects. 
plot(predictorEffects(defect.cellsize.mod, ~ tricomeC, 
                      xlevels = list(logCSc = c(-0.4, -0.09, 0.09, 0.2))), 
     axes=list(grid=TRUE,
               x=list(rug=FALSE),
               y=list(type="response")),
     lines=list(multiline=TRUE, z.var = "logCSc")
)

#lower tricome count = bigger cells (less cells in the little counting box). When paired with the wings also being bigger, this would mean (to me) that these wings are bigger because there are bigger cells. 

#In line with the idea that as cells get bigger, the wings get bigger. 
#Should really check that there are wings with counts this low?

#I am going to save these so I can look at multiple things at once. 
png("../Figures/zi418ef43_wingsizeC_density.png")
ggplot(zi418ef43.catDat) + 
  geom_density(aes(x = logCSc, col = sex.x))
dev.off()

png("../Figures/zi418ef43_tricomeC_density.png")
ggplot(zi418ef43.catDat) + 
  geom_density(aes(x = tricomeC, col = sex.x))
dev.off()

#Redoing this with sex specific means. 


zi418ef43.catDat2 <- (zi418ef43.catDat %>% 
                        group_by(sex.x) %>% 
                        mutate(logCSsexC = logCS - mean(logCS),
                               tricomeSexC = tricome_mu - mean(tricome_mu))
                      )

#works. 
ggplot(zi418ef43.catDat2) + 
  geom_density(aes(x = logCSsexC, col = sex.x))

ggplot(zi418ef43.catDat2) + 
  geom_density(aes(x = tricomeSexC, col = sex.x))



defect.cellsize.mod2 <- glm(defect ~ (logCSsexC + sex.x + tricomeSexC)^2, 
                           data = zi418ef43.catDat2, family = binomial(link = "logit"))


#this makes way more sense. 
summary(defect.cellsize.mod2)
Anova(defect.cellsize.mod2)

#chocse these cut points to be slightly below and above the max size. 
#maybe an interaction between biger wings, smaller tricome counts and more wing defects. 

#The CI on these are probably really big. I can make a better version of this plot later but for now this is good. 
plot(predictorEffects(defect.cellsize.mod2, ~ tricomeSexC, 
                      xlevels = list(logCSc = c(-0.15, -0.05, 0.05, 0.15))), 
     axes=list(grid=TRUE,
               x=list(rug=FALSE),
               y=list(type="response")),
     lines=list(multiline=TRUE, z.var = "logCSsexC")
)

#also want to ask if variation among regions is a good predictor.
#still a whole lot of nothing to see here. 

#keeping the two way interactions here. 
defect.cellsizeVar.mod <- glm(defect ~ (logCSsexC + sex.x + tricomeSexC + tricome_cv)^2, 
                            data = zi418ef43.catDat2, family = binomial(link = "logit"))

summary(defect.cellsizeVar.mod)
Anova(defect.cellsizeVar.mod)



#now I want to do something with shape. 
#Do I want to also centre shape? Will that matter?


#To remind myself, because I haven't worked with this data in years, I am going to take the normal data and do a PCA and just look at it for a second before I decide what to do with the analysis. 

wing.pca <- prcomp(zi418ef43.catDat2[,14:109])

wing.pca.df <- data.frame(zi418ef43.catDat2, wing.pca$x[,1:56])

#so sex is PC1, which is also size, (probably driving a lot of this) 
ggplot(wing.pca.df) + 
  geom_point(aes(x = PC1, y = PC2, col = sex.x))

ggplot(wing.pca.df) + 
  geom_point(aes(y = PC1, x = logCS, col = sex.x))


#what if I just include size in the PCA 

wing.pca2 <- prcomp(zi418ef43.catDat2[,c(14:109, 125)])

wing.pca.df2 <- data.frame(zi418ef43.catDat2, wing.pca2$x[,1:57])


#good. 
ggplot(wing.pca.df2) + 
  geom_point(aes(y = PC1, x = logCS, col = sex.x))


#this looks better in terms of the sex differences. 
ggplot(wing.pca.df2) + 
  geom_point(aes(x = PC2, y = PC3, col = sex.x))

ggplot(wing.pca.df2) + 
  geom_point(aes(x = PC4, y = PC5, col = sex.x))

ggplot(c) + 
  geom_point(aes(x = PC6, y = PC7, col = sex.x))
#seems fine to me. 

#A nicer scree plot. 
#ths ovi looks crazy 
data.frame(sdev = wing.pca2$sdev,
           PC = seq(1:length(wing.pca2$sdev))) %>% 
  ggplot(aes(x =as.factor(PC), y = sdev)) + 
    geom_col()


#this drops 1 
data.frame(sdev = wing.pca2$sdev,
           PC = seq(1:length(wing.pca2$sdev))) %>% 
  filter(PC != 1) %>% 
  ggplot(aes(x =as.factor(PC), y = sdev)) + 
  geom_col()

#this keeps to 57 (the total avail dems)
#the us unfortunately no super clear cut off here that I would use for the stats :(
data.frame(sdev = wing.pca2$sdev,
           PC = seq(1:length(wing.pca2$sdev))) %>% 
  filter(PC > 1 & PC <= 57) %>% 
  ggplot(aes(x =as.factor(PC), y = sdev)) + 
  geom_col()

#I don't really want to have to have all 56 dimentions as predictors, that is crazy. I am going to start with the first 3 and see if there is an assosiation. 


#should I also have the cell size stuff in here? 
#No interactions because I think all these effects are perpendicular by definition. 
shape.defect.mod <- glm(defect ~ logCSsexC + PC2 + PC3 + PC4, 
                        data = wing.pca.df2, family = binomial(link = "logit"))



#what do I do with missing landmarks? 
summary(shape.defect.mod)

#so there is an effect of shape but I need to remember how shape works when these things are missing. I should just drop the crossvein landmarks and recheck this result
#For the most part, these are not compleate losses and I guessed where the intersections would be.
Anova(shape.defect.mod)

#I actually thought this was way more focused on cvl cases. 
wing.pca.df2 %>% 
  group_by(sex.x) %>% 
  summarize(acvl = sum(acvl, na.rm = T), 
            pcvl = sum(pcvl, na.rm = T), 
            log.vein.loss = sum(log.vein.loss, na.rm = T),
            extra.vein = sum(extra.vein, na.rm = T))

wing.pca.df2 %>% 
  group_by(sex.x) %>% 
  summarize(acvl = sum(acvl, na.rm = T)/n(), 
            pcvl = sum(pcvl, na.rm = T)/n(), 
            log.vein.loss = sum(log.vein.loss, na.rm = T)/n(),
            extra.vein = sum(extra.vein, na.rm = T)/n())


#I thought I would plot the shape change along PC2 and PC3. 

#I need a mean shape 

#colMeans would have been a better idea here. But this is what came to my mind first. 
m <- colMeans(wing.pca.df2[,14:109])

#putting this into a matrix. 
#I made this for myself before 
#why are the first and last number 0? ths
vec2matrix <- function(x) {
  y <- unlist(x)
  m <- matrix(y, nrow = 48, ncol = 2, byrow = TRUE)
}

m2 <- vec2matrix(m)

#PC2
#tried with 2sd and then went with 3 to make the visualization easier/not have to mess with mag to see. 
#3 * sdev of PC *  rotation
pc2eff <- 3 * wing.pca2$sdev[2] * wing.pca2$rotation[,2]

pc2effM <- m2 + vec2matrix(pc2eff[1:48])


#this contains the winglinks needed for plotting 
source("KP_geomorphwingfunctions.R")

#PC2 is pretty subtle and mostly associated with the lm that would be missing in the weird ones. But there is also something weird going on with the top margin. 
plotRefToTarget(m2, pc2effM , 
                method = "points", 
                links = wing.links,
                mag = 1, 
                gridPars = wing.spec)


#PC3
pc3eff <- 3 * wing.pca2$sdev[3] * wing.pca2$rotation[,3]

pc3effM <- m2 + vec2matrix(pc3eff[1:48])


plotRefToTarget(m2, pc3effM , 
                method = "points", 
                links = wing.links,
                mag = 1, 
                gridPars = wing.spec)


#I am not convinced that the effect is just that those are the landmarks affected by the loss. I need to go back and look at the images. 



#Ideas to go forward with shape analysis:
#1 drop landmarks associated with the loss for each case (probably removes the most variable landmarks)
#Go back and annotate which landmarks are actually missing in these images. Mark these with NA. 
#I can't remove those animals from the analysis and having NAs makes doing the PCA impossible. 
#I can remove the landmark that is NA from the entire dataset. This is the most conservative option. 
