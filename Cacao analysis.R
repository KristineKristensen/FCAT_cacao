#Foreest cover

library(raster)
library(sf)
library(terra)

#creating the FCAT boudary layer
boundaryFCAT <-st_read("boundary.shp")
b <- extent(boundaryFCAT) 
boundary <- raster(b)
crs(boundary) <- "+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs +type=crs"


# forrest proportion in each plot
FCATtoCachi <- raster("data/2022084bandFCATtoCachi20230509sieved100.tif")
FCAT_veg <- crop(x=FCATtoCachi, y= boundary)

plots <- read.csv("data/plot_extent.csv", sep=";")
plotname <- plots$X
data_list <- list()
prop_list <- list()

for (p in plotname)
{kor <- as.list(subset(plots, X==p))
kant <- raster(ymx= kor$ymax, # maximum y coordinate 
               xmn=kor$xmin, # minimum x coordinate
               ymn=kor$ymin, #minimum y coordinate
               xmx=kor$xmax, #minimum y coordinate
               resolution=FCAT_veg) # defining resolution
crs(kant) <- "+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs +type=crs"
veg_plot <- crop(x=FCAT_veg,
                 y= kant)

total <- ncell(veg_plot)
forest <- sum(values(veg_plot) == 1)
proportion_forest <- forest /total
veg_data <- c(terra::values(veg_plot))
data_list[[p]] <- veg_data
prop_list[[p]] <-  proportion_forest

}

data_forrest_cover <- do.call(rbind, data_list)

forrest_proportion <- do.call(rbind, prop_list)
print(forrest_proportion)



###########################################################################
#statistical models

library(MASS)
library(DHARMa)
library(lme4)
library(car)




#load data
all <- read.csv("data/all.csv", sep= ";")

# check data
head(all, n=5)
str(all)

#################################################################################
#Healthy pods

#POISSON glm
glm_healthy <- glm(healthy ~ distance_m+ type+ forrest, data= all, family="poisson")
summary(glm_healthy)

#Quasipoisson
glm_healthy_quasi<- glm(healthy ~ distance_m+ type+ forrest, data= all, family="quasipoisson")
summary(glm_healthy_quasi)


#Negative binomial
nb_glm_healthy <- glm.nb(healthy ~ distance_m+ type+ forrest, data= all)
summary(nb_glm_healthy)


#visual investigation
# define a 2x2 plotting space
par(mfrow=c(2,2))
# Evaluate the assumptions
plot(nb_glm_healthy)

#using box-plot to see if there is homogenety of variance between plots (locations are not considered as the level in a random factor should be bigger than five to say something about variances between groups)
glm_res_healthy <- residuals(object = nb_glm_healthy, type = "pearson")
#residuals by plots
boxplot(glm_res_healthy~ plots,
        data = all)


#Testing if plots have a significant influence
glmnb_healthy_p <- MASS::glm.nb(healthy ~ distance_m+ type+forrest+ plots, data= all)
summary(glmnb_healthy_p)
anova(glmnb_healthy_p,test = "Chisq")

#mixed-effect model 
ri_healthy <- glmer(healthy ~  type+ distance_m+ forrest+ (1 |plots), family="poisson", data =all)
summary(ri_healthy)


#Random intercept model With negative binomial
ri_nb_healthy <- glmer.nb(healthy ~  type+ distance_m+forrest+ (1 |plots), data =all)
summary(ri_nb_healthy)


#check model with DHARMa
sim_healthy <- simulateResiduals(fittedModel = ri_nb_healthy, n=1000)
plot(sim_healthy, form=NULL)

#check for collinearity
vif(glm_healthy)



###################################################################################
# Pods with large dark spots


#POISSON glm
glm_manchas<- glm(manchas~ distance_m+ type+ forrest, data= all, family="poisson")
summary(glm_manchas)

#visual investigation
# define a 2x2 plotting space
par(mfrow=c(2,2))
# Evaluate the assumptions
plot(glm_manchas)



#using box-plot to see if there is homogenety of variance between plots
glm_res <- residuals(object = glm_manchas, type = "pearson") 
#residuals by plots
boxplot(glm_res~ plots, data = all)

#Testing if plots have a significant influence
glm_manchas_p<- glm(manchas~  distance_m+ type+ forrest+ plots, data= all, family="poisson")
summary(glm_manchas_p)
#error in fit

# check multicollinearity
vif(glm_manchas_p)



#Random intercept model
ri_manchas <- glmer(manchas ~ type+distance_m + forrest+(1 |plots), data =all, family = "poisson")
summary(ri_manchas)

#check model with DHARMa
all_rand_sim <- simulateResiduals(fittedModel = ri_manchas, n=1000)
plot(all_rand_sim, form=NULL)

#Random intercept model With negative binomial
ri_nb_manchas <- glmer.nb(manchas ~ type+distance_m + forrest+(1 |plots), data =all)
summary(ri_nb_manchas)

#check model with DHARMa
sim_manchas_nb <- simulateResiduals(fittedModel = ri_nb_manchas, n=1000)
plot(sim_manchas_nb, form=NULL)



#################################################################################################################

# Pods with White spots


#POISSON glm
glm_blanco<- glm(blanco~ distance_m+ type+ forrest, data= all, family="poisson")
summary(glm_blanco)

#quasipoisson
glm_blanco_q<- glm(blanco~ distance_m+ type+ forrest, data= all, family="quasipoisson")
summary(glm_blanco_q)

# check multicollinearity
vif(glm_blanco)

#Testing if plots have a significant influence
glm_blanco_place<- glm(blanco~ distance_m+ type+ forrest+plots, data= all, family="poisson")
summary(glm_blanco_place)
#error in fit, so no ANOVA is made


#using box-plot to see if there is homogenety of variance between plots
glm_res_blanco <- residuals(object = glm_blanco, type = "pearson") 
#residuals by plots
boxplot(glm_res_blanco~ plots, data = all)


#Random intercept model
ri_blanco <- glmer(blanco~ type+ distance_m + forrest+ (1 |plots), data =all, family = "poisson")
summary(ri_blanco)

#check model with DHARMa
sim_blanco <- simulateResiduals(fittedModel = glm_blanco, n=1000)
plot(sim_blanco, form=NULL)


###########################################################################################################################

# Pods with small dark spots

#POISSON glm
glm_puntitos<- glm(puntitos~ distance_m+ type + forrest, data= all, family="poisson")
summary(glm_puntitos)

#quasipoisson
glm_puntitos_q<- glm(manchas~ distance_m+ type + forrest, data= all, family="quasipoisson")
summary(glm_puntitos_q)

#visual inspection
par(mfrow=c(2,2))
plot(glm_puntitos_q)

#using box-plot to see if there is homogenety of variance between plots
glm_res_puntitos <- residuals(object = glm_puntitos, type = "pearson")
#residuals by plots
boxplot(glm_res_puntitos~ plots,data = all)

#Random intercept model
ri_puntitos <- glmer(puntitos~ distance_m+ type+ forrest+ (1|plots), data= all, family="poisson")
summary(ri_puntitos)

#check multicollinearity
vif(ri_puntitos)

#check model with DHARMa
sim_puntitos <- simulateResiduals(fittedModel = ri_puntitos, n=1000)
plot(sim_puntitos, form=NULL)

##########################################################################################################################
# trees with Witches broom

#Poisson glm
glm_bruja<- glm(bruja~ distance_m+ type+forrest, data= all, family="poisson")
summary(glm_bruja)

#visual inspection
par(mfrow=c(2,2))
plot(glm_bruja)

#quasipoisson
glm_bruja_q<- glm(bruja~ distance_m+forrest+ type, data= all, family="quasipoisson")
summary(glm_bruja_q)

#Testing if plots have a significant influence
glm_bruja_p<- glm(bruja~ distance_m+ type+ plots, data= all, family="poisson")
summary(glm_bruja_p)
#model cannot be fitted

#using box-plot to see if there is homogenety of variance between plots
glm_res_bruja <- residuals(object = glm_bruja, type = "pearson")
#residuals by plots
boxplot(glm_res_bruja~ plots, data = all)


#Random intercept
ri_witches <- glmer(bruja ~ distance_m+ type+forrest+ (1|plots), family="poisson", data=all)
summary(ri_witches)

#check multicollinearity
vif(ri_witches)

#check model with DHARMa
sim_witches <- simulateResiduals(fittedModel = ri_witches, n=1000)
plot(sim_witches, form=NULL)

#################################################################################################################################

#Binomial with 0= healthy and 1= disease

#Binomial glm
glm_disease <- glm(disease ~  distance_m + type + total+ forrest , family = "binomial", data=all)
summary(glm_disease)

#Visual inspection
par(mfrow=c(2,2))
plot(glm_disease)


#Using box-plot to see if there is homogenety of variance between plots
glm_res_disease <- residuals(object = glm_disease, type = "pearson")
#residuals by plots
boxplot(glm_res_disease~ plots, data = all)

#Random intercept
ri_disease <- glmer(disease ~ distance_m + type +forrest+ total+ (1|plots), family = "binomial", data=health )
summary(ri_disease)

#test multicollinearity
vif(ri_disease)


#check model with DHARMa
sim_disease <- simulateResiduals(fittedModel = ri_disease, n=1000)
plot(sim_disease, form=NULL)







