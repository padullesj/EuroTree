
#####
# Core Script for
# "Environmental drivers of taxonomic and functional turnover of tree assemblages in Europe"
# Journal: Oikos
#
# Author: XXXX XXXX
# Date: 14/06/2022
######

#load R libraries:
library(ecodist)
library(FD)
library(ggrepel)
library(corrplot)
library(gridExtra)
library(grid)
library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
library(plyr)
library(raster)
library(sp)
library(betapart)
library(sf) #do not load if step 7 runs into an error
library(PhyloMeasures)
library(cluster)
library(ggpubr)
library(gdm)
library(rnaturalearth) #do not load if step 7 runs into an error
library(RColorBrewer)
library(viridis)

#Set working directory if necessary:
#setwd("/Users/padulles/Documents/project/Rsession")



####
# Table of contents
#
# 1) Load and clean-up trait data
# 2) Load and clean-up tree data
# 3) Aggregate species occurrences within 100 x 100 km grid cells
# 4) Transform long-format into wide-format table (community data)
# 5) Calculate beta-diversity
# 6) Cluster analysis
# 7) Run Multiple Regression on distance Matrices (MRMs)
# 8) Individual effects of environmental variables
# 9) Mapping turnover
####



###
# 1) Load and clean-up trait data
###

#Load imputed trait data:
traits<-read.table("Processed_data/imputed_traits.csv") #I retrieved original trait data from TRY and GRooT.

#Clean-up trait data (keep selected traits and remove species with NAs):
traits<-traits[ , -which(names(traits) %in% c("leaf_P", "leaf_C", "leaf_length"))]
traits<-na.omit(traits)

#Run PCA on trait data:
res.pca <-psych::principal(traits, rotate="varimax", nfactors =4) #Select the first 4 PCA (axes of trait variation). 
spp_scores<-res.pca$scores #save scores



###
# 2) Load and clean-up tree data
###

#Load species data and prepare community matrix:
wood<-unique(fread("Raw_data/EUForestspecies.csv", select = c("X","Y", "SPECIES NAME")))
names(wood)[3]<-"species"
wood<-subset(wood, Y>1500000)
wood<-subset(wood, X<6000000)
wood<-wood[wood$species %in% unique(rownames(traits)), ] #Keep only species with values for all traits



###
# 3) Aggregate species occurrences within 100 x 100 km grid cells
###

#Define extent:
xmn=2600000
xmx=5900000
ymn=1500000
ymx=5400000

#Create raster with the specified spatial resolution (100 x 100 km):
x <- raster(ncol=(xmx-xmn)/100000, nrow=(ymx-ymn)/100000, xmn=xmn, xmx=xmx, ymn=ymn, ymx=ymx)
values(x) <- 1:ncell(x)

#Extract raster values for each initial location. This tells what initial locations are in each 100 x 100 km grid cell:
wood$xloc <-raster::extract(x, wood[,c("X", "Y")]) #X = longitude, Y = latitude

#Convert raster into data frame:
spdf <- as(x, "SpatialPixelsDataFrame") #convert into spatial object
spdf <- as.data.frame(spdf) #convert into data frame
colnames(spdf) <- c("value", "x", "y") #change column names
spdf<-spdf[spdf$value %in% unique(wood$xloc), ] #reatain only grid cells with initial locations
#write.table(spdf, "Processed_data/coordinates_cells.csv") #save to use this dataset to extract our set of predictors


###
# 4) Transform long-format into wide-format table (community data)
###

wood<-unique(wood[,c("species", "xloc")]) #select species
wood$value<-1 # add field to assign values to new columns
wood<-reshape2::dcast(wood, xloc ~ species, value.var="value") #transform into wide format
rownames(wood)<-wood$xloc #assign rownames
wood$xloc<-NULL #and remove column with names
wood <- wood  %>% mutate_all(funs(replace_na(.,0))) #replace NAs with 0
wood <- wood[-which(rowSums(wood[sapply(wood, is.numeric)]) < 2),] #discard grid cells with less than 2 species



####
# 5) Calculate beta-diversity
####

#First, determine the best agglomerative method to construct the dendrograms:
fdist<-cluster::daisy(spp_scores, metric="euclidean") #calculate Euclidean distances between species

#Test different Hierarchical Agglomerative Clustering:
h1=cophenetic(hclust(fdist,method='average'))
h2=cophenetic(hclust(fdist,method='complete'))
h3=cophenetic(hclust(fdist,method='ward.D2'))
h4=cophenetic(hclust(fdist,method='single'))

#Get correlations
out<-data.frame(upgma=cor(fdist, h1), complete=cor(fdist, h2), ward=cor(fdist, h3), single=cor(fdist, h4))

#Repeat for the different trait dimensions:
for(i in 1:4){
  
  #Calculate Euclidean distances between species:
  fdist<-cluster::daisy(as.data.frame(spp_scores[,c(i)]), metric="euclidean")
  
  #Hierarchical Agglomerative Clustering:
  h1<-cophenetic(hclust(fdist,method='average'))
  h2<-cophenetic(hclust(fdist,method='complete'))
  h3<-cophenetic(hclust(fdist,method='ward.D2'))
  h4<-cophenetic(hclust(fdist,method='single'))
  
  #Correlations:
  out2<-data.frame(upgma=cor(fdist, h1), complete=cor(fdist, h2), ward=cor(fdist, h3), single=cor(fdist, h4))
  
  #Save:
  out<-rbind(out, out2)
}
#write.table(out, "Processed_data/HAC_result.csv") #Decided to use UPGMA (average) because it showed the highest cophenetic correlation

#Taxonomic beta-diversity:
tdist<-beta.pair(wood, index.family="sorensen") #Calculate beta-diversity
#write.table(as.matrix(tdist$beta.sim), "Processed_data/taxonomic_turn.csv") #Save turnover component
#write.table(as.matrix(tdist$beta.sne), "Processed_data/taxonomic_nest.csv") #Save nestedness component

#Functional beta-diversity (all trait dimensions):
fdist<-cluster::daisy(spp_scores, metric="euclidean") #Again, get Euclidean distances between species
tree<-as.phylo(hclust(fdist, method="average")) # Use previous output to create the functional dendrogram
fdist<-phylo.beta.pair(wood, tree, index.family="sorensen") #Calculate beta-diversity
#write.table(as.matrix(fdist$phylo.beta.sim), "Processed_data/functional_turn.csv") #Save turnover component
#write.table(as.matrix(fdist$phylo.beta.sne), "Processed_data/functional_nest.csv") #Save nestedness component

#Repeat for the different trait dimensions (might take a while):
for(i in 1:4){
  print(i) #to keep track of the process
  fdist<-cluster::daisy(as.data.frame(spp_scores[,c(i)]), metric="euclidean") #get Euclidean distances between species
  tree<-as.phylo(hclust(fdist, method="average")) # Use previous output to create the functional dendrogram
  fdist<-phylo.beta.pair(wood, tree, index.family="sorensen") #Calculate beta-diversity
  #write.table(as.matrix(fdist$phylo.beta.sim), paste("Processed_data/functional_RC", i, "_turn.csv", sep="")) #Save turnover component
  #write.table(as.matrix(fdist$phylo.beta.sne), paste("Processed_data/functional_RC", i, "_turn.csv", sep="")) #Save nestedness component
} 



###
# 6) Cluster analysis
###

#First, determine the best agglomerative method to construct the dendrograms:

#Load turnover data (showing example for taxonomic turnover; should do for functional turnover too)
bd<-as.dist(read.table("Processed_data/taxonomic_turn.csv"))

# Hierarchical Agglomerative Clustering
h1<-cophenetic(hclust(bd,method='average'))
h2<-cophenetic(hclust(bd,method='complete'))
h3<-cophenetic(hclust(bd,method='ward.D2'))
h4<-cophenetic(hclust(bd,method='single'))

#Run correlations (similar results for UPGMA and Ward). After visually exploring the dendrograms, I chose Ward's method.
out<-data.frame(upgma=cor(bd, h1), complete=cor(bd, h2), ward=cor(bd, h3), single=cor(bd, h4))

#Hierarchical Agglomerative method (Ward's):

#Taxonomic turnover
bd<-as.dist(read.table("Processed_data/taxonomic_turn.csv"))
wd<-hclust(bd, method="ward.D2")
plot(wd) #plot dendrogram

#Plot silhouette to determine the "optimum" number of clusters
plot(2:10, sapply(2:10, function(i) { 
  mean(silhouette(cutree(wd, i), dmatrix=as.matrix(bd))[,"sil_width"]) }),
  xlab="Number of clusters", ylab="Average Silhouette", type="b", pch=20)

#Save group results:
groups <- as.data.frame(cutree(wd, k=2)) # cut tree into optimal number of clusters (and more):
names(groups)[1]<-"tk2" #change name of column
groups$tk3 <- cutree(wd, k=3) #3 groups
groups$tk4 <- cutree(wd, k=4) #4 groups
groups$value<-rownames(groups) #assign rownames (site IDs)

#Functional turnover
bd<-as.dist(read.table("Processed_data/functional_turn.csv"))
wd<-hclust(bd, method="ward.D2")
plot(wd) #plot dendrogram

#Plot silhouette to determine the "optimum" number of clusters
plot(2:10, sapply(2:10, function(i) { 
  mean(silhouette(cutree(wd, i), dmatrix=as.matrix(bd))[,"sil_width"]) }),
  xlab="Number of clusters", ylab="Average Silhouette", type="b", pch=20)

#Save group results:
groups$fk2 <- cutree(wd, k=2)
groups$fk3 <- cutree(wd, k=3)
groups$fk4 <- cutree(wd, k=4)

#Write.output
#write.table(groups, "Processed_data/clusters.csv")



###
# 7) Run Multiple Regression on distance Matrices (MRMs)
###

#Define routes to load data:
routes<-c("Processed_data/taxonomic_turn.csv", #taxonomic turnover
          "Processed_data/functional_turn.csv", #functional turnover
          "Processed_data/functional_RC1_turn.csv", #turnover in LES
          "Processed_data/functional_RC2_turn.csv", #turnover in RES
          "Processed_data/functional_RC3_turn.csv", #turnover in PWSM
          "Processed_data/functional_RC4_turn.csv") #turnover in TH

#Define response variables:
types<-c("TBD", "All trait dimensions", "LES", "RES", "PWSM", "TH")

#load function to normalize data:
min_max_norm <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

#Run loop for each response variable:
#WARNING! If the following code runs into an error (Error in xj[i, , drop = FALSE] : incorrect number of dimensions), 
#re-launch R and do no load packages "sf" and "rnaturalearth"
#There seems to be a problem between the functions in these packages and function MRM in package ecodist.
out<-NULL #empty object to store the result
for(i in 1:length(types)){
  print(i) #to keep track of the process
  bd<-as.matrix(read.table(routes[i])) #load response variable
  
  #load predictors:
  predictors<-read.table("Processed_data/predictors.csv") #load predictors (obtained from WorldClim and SoilGrids)
  predictors<-predictors[predictors$value %in% unique(rownames(bd)), ] #subset only grid cells for which there is data
  
  #run PCA
  res.pca <-psych::principal(predictors[,c(4:22)], rotate="varimax", nfactors =3) #perform PCA on macro-climatic variables
  predictors2<-as.data.frame(res.pca$scores) #keep PCA scores
  predictors2<-cbind(predictors2, predictors[,c(1:3,23:25)]) #merge with the other soil variables
  
  #get geographic distances:
  geo<-read.table("Processed_data/geographic_distances.csv")
  
  #run MRMs:
  m1<-MRM(as.dist(bd) ~ min_max_norm(dist(predictors2$RC1)) + 
            min_max_norm(dist(predictors2$RC2)) + 
            min_max_norm(dist(predictors2$RC3)) + 
            min_max_norm(dist(predictors2$ph)) + 
            min_max_norm(as.dist(min_max_norm(geo[,-c(1)]))), nperm=999)
  res<-as.data.frame(m1$coef) #store regression coefficients
  res$rsquared<-m1$r.squared[1] #store R-squared from MRMs
  res$type<-types[i] #identify response variable
  
  out<-rbind(out, res) #merge results
}
out$`as.dist(bd)`<-round(out$`as.dist(bd)`*10, 2) #round output and multiply by 10 the regression coefficients
out$pval<-round(out$pval, 3) #round p-value
out$rsquared<-round(out$rsquared, 2) #round R-Squared
#write.table(out, "Processed_data/mrm_results.csv") #save output



###
# 8) Individual effects of environmental variables
###

#Loop for each response variable
out<-NULL #create object to store results
out2<-NULL #create object to store results
for(i in 1:length(types)){
  print(i) #to keep track of he process 
  bd<-as.matrix(read.table(routes[i])) #load response variable
  
  #load predictors:
  predictors<-read.table("Processed_data/predictors.csv") #load predictors (obtained from WorldClim and SoilGrids)
  predictors<-predictors[predictors$value %in% unique(rownames(bd)), ]  #subset only grid cells for which there is data
  predictors$ph<-predictors$ph/10 #original soil pH units
  
  #run PCA
  res.pca <-psych::principal(predictors[,c(4:22)], rotate="varimax", nfactors =3) #perform PCA on macro-climatic variables
  predictors2<-as.data.frame(res.pca$scores) #keep PCA scores
  predictors2<-cbind(predictors2, predictors[,c(1:3,23:25)]) #merge with the other soil variables
  
  #get geographic distances:
  geo<-as.matrix(read.table("Processed_data/geographic_distances.csv"))
  
  #prepare tables for the models:
  sites <- as.matrix(predictors2$value) 
  colnames(sites) <- c("value")
  bd2 <- cbind(sites, bd)
  
  # First prepare a GDM Site-Pair table with biological and predictor data:
  exFormat3 <- formatsitepair(bd2, 3, XColumn="x", YColumn="y", predData = predictors2, siteColumn="value")
  
  # Then we add our own geographic distance matrix calculated using the format "4":
  exFormat4 <- formatsitepair(exFormat3, 4, predData=predictors2, siteColumn="value", distPreds=list(geo))
  
  #Get normalized differences:
  exFormat4$distance <- min_max_norm((exFormat4$distance))
  exFormat4$dif_RC1 <- min_max_norm(abs(exFormat4$s1.RC1- exFormat4$s2.RC1))
  exFormat4$dif_RC2 <- min_max_norm(abs(exFormat4$s1.RC2 - exFormat4$s2.RC2))
  exFormat4$dif_RC3 <- min_max_norm(abs(exFormat4$s1.RC3 - exFormat4$s2.RC3))
  exFormat4$dif_ph <- min_max_norm(abs(exFormat4$s1.ph - exFormat4$s2.ph))
  exFormat4$dif_geo <- min_max_norm(abs(exFormat4$s1.matrix_1 - exFormat4$s2.matrix_1))
  
  #Prepare output table to get slope and R-squared:
  models<-as.data.frame(matrix(nrow=4, ncol=3))
  names(models)<-c("coef", "pvalue", "r2")
  rownames(models)<-c("RC1", "RC2", "RC3", "ph")
  
  #Run models with individual variables while controlling for geographic distances:
  RC1a<-summary(lm(distance~dif_RC1 + dif_geo, data=exFormat4))
  RC2a<-summary(lm(distance~dif_RC2 + dif_geo, data=exFormat4))
  RC3a<-summary(lm(distance~dif_RC3 + dif_geo, data=exFormat4))
  pha<-summary(lm(distance~dif_ph + dif_geo, data=exFormat4))
  
  #Store coefficients, p-values and R-squares in the table:
  models[1,1]<-RC1a$coefficients[2,1]
  models[1,2]<-RC1a$coefficients[2,4]
  models[1,3]<-RC1a$r.squared
  models[2,1]<-RC2a$coefficients[2,1]
  models[2,2]<-RC2a$coefficients[2,4]
  models[2,3]<-RC2a$r.squared
  models[3,1]<-RC3a$coefficients[2,1]
  models[3,2]<-RC3a$coefficients[2,4]
  models[3,3]<-RC3a$r.squared
  models[4,1]<-pha$coefficients[2,1]
  models[4,2]<-pha$coefficients[2,4]
  models[4,3]<-pha$r.squared
  models$type<-types[i]
  
  #Generate predictions based on model fits:
  yRC1<-RC1a$coefficients[1,1] + (RC1a$coefficients[2,1]*exFormat4$dif_RC1) + (RC1a$coefficients[3,1])
  yRC2<-RC2a$coefficients[1,1] + (RC2a$coefficients[2,1]*exFormat4$dif_RC2) + (RC2a$coefficients[3,1])
  yRC3<-RC3a$coefficients[1,1] + (RC3a$coefficients[2,1]*exFormat4$dif_RC3) + (RC3a$coefficients[3,1])
  yph<-pha$coefficients[1,1] + (pha$coefficients[2,1]*exFormat4$dif_ph) + (pha$coefficients[3,1])
  
  #Merge into data frame:
  res<-data.frame(y=c(yRC1, yRC2, yRC3, yph), 
                  x=c(exFormat4$dif_RC1, exFormat4$dif_RC2, exFormat4$dif_RC3, exFormat4$dif_ph),
                  variable=c(rep("RC1", length(yRC1)), rep("RC2", length(yRC2)), rep("RC3", length(yRC3)), rep("ph", length(yph))))
  res$type<-types[i] #assign respnose variable 
  
  #Merge in each loop
  out<-rbind(res, out)
  out2<-rbind(out2, models)
}
out2$coef<-out2$coef*10 #multiply coefficients by 10
out2$pvalue<-round(out2$pvalue, 3) #round p-values
out2$r2<-round(out2$r2, 2) #round R-squared
#write.table(out, "Processed_data/regression_results.csv") #save output
#write.table(out2, "Processed_data/regression_coefficients.csv") #save output

#Example to plot results from individual regressions while controlling for geographical distances:

#This example corresponds to the effect of "Precipitation" on the different measures of turnover:
sub<-subset(out, variable=="RC1") #subset data for "Precipitation" (i.e., RC1 axis from the PCA)
bp1<-ggplot(data=sub, aes(x=x, y=y, color=type)) +
  geom_smooth(stat="smooth", method='lm', formula= y~x)+
  labs(x="\u0394 Precipitation (PC2)", y="Functional turnover") +
  scale_color_manual(values=c("brown1", "deepskyblue", "darkgreen", "gold2", "darkorchid"))+
  scale_y_continuous(limits = c(NA, 0.95), breaks = c(0.2, 0.4, 0.6, 0.8))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text = element_text(size = 15),
        legend.position="bottom",
        panel.border = element_rect(fill=NA, colour = "black", size=1),
        legend.title = element_blank())

###
# 9) Mapping turnover
###

#Example code to plot taxonomic and functional turnover across grid cells:

#Load political boundaries in Europe:
world <- ne_countries(scale = "medium", returnclass = "sf") #load world map
europe <- world[which(world$continent == "Europe"),] #select Europe
europe <- europe[which(europe$sovereignt != "Russia" & europe$sovereignt != "Ukraine" & europe$sovereignt != "Iceland"),] #remove countries with missing data
europe2 <- st_transform(europe, crs = "+init=epsg:3035") #change coordinate system

#Create a formatted ggplot theme:
theme_opts<-list(theme(panel.background = element_blank(),
                       plot.background = element_rect(fill="white"),
                       panel.border = element_blank(),
                       panel.grid.major = element_blank(),
                       axis.text.x = element_blank(),
                       axis.text.y = element_blank(),
                       axis.ticks = element_blank(),
                       axis.title.x = element_blank(),
                       axis.title.y = element_blank(),
                       legend.background=element_blank(),
                       legend.position="bottom",
                       legend.text=element_text(size=9),
                       legend.title=element_text(size=13),
                       plot.title = element_text(size=25, face="bold", hjust = 0.5)))

#Load taxonomic turnover and clean-up data:
bd<-as.matrix(read.table("Processed_data/taxonomic_turn.csv"))
diag(bd)<-NA #set diagonal to NA
mbd<-as.data.frame(rowMeans(bd, na.rm=T)) #get mean turnover across grid cells
names(mbd)[1]<-"tbd" #set name of variable

#Load functional turnover:
bd<-as.matrix(read.table("Processed_data/functional_turn.csv"))
diag(bd)<-NA #set diagonal to NA
mbd$fbd<-rowMeans(bd, na.rm=T) #get mean turnover across grid cells

#Load predictors:
predictors<-read.table("Processed_data/predictors.csv")
predictors<-predictors[predictors$value %in% unique(rownames(mbd)), ] #keep grid cells for which there is data available

#Assign coordinates:
mbd$x<-predictors$x
mbd$y<-predictors$y

#Convert raster into data frame for plotting:
df <- st_as_sf(x = mbd, coords = c("x", "y"), crs = "+init=epsg:3035")
df <- st_transform(df, 3035) #set coordinante system

#plot
tt<-ggplot() + 
  geom_sf(data = df, aes(colour = tbd), shape=15, size=2.1) +
  scale_colour_viridis(name="Mean\ntaxonomic\nturnover", direction=-1) +
  geom_sf(data = europe2, fill=NA, color="black", size=0.25) +
  xlim(c(2600000, 5900000)) +
  ylim(c(1500000, 5400000)) +
  theme_opts

#plot
ft<-ggplot() + 
  geom_sf(data = df, aes(colour = fbd), shape=15, size=2.1) +
  scale_colour_viridis(name="Mean\nfunctional\nturnover", direction=-1, breaks=c(0.2, 0.4, 0.6)) +
  geom_sf(data = europe2, fill=NA, color="black", size=0.25) +
  xlim(c(2600000, 5900000)) +
  ylim(c(1500000, 5400000)) +
  theme_opts

#Put together:
p<-ggarrange(tt, ft, labels=c("a)", "b)"), ncol = 2, nrow = 1)

#Save output:
png("Results/turnovers_map.png",
    res=600, height=4,width=8,units="in"); 
p
dev.off()

#Clean-up environment:
rm(list = ls())

###
#End of script
###