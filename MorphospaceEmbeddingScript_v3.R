# scMorphometrics morphospace embedding script | As used in Andrews et al, 2021 'Single-cell morphometrics reveals ancestral principles of notochord development'

# This script takes statistics on cell shape exported from Imaris, performs PCA to embed them in a 'morphospace', and uses Slingshot (Street et al, 2018) for trajectory inference.
# It also offers a means to export new tiff files with point centres for each cell intensity-coded for pseudotime values, which can be imported to Imaris for visualisation as a colour code

# The data input for the script are .csv files exported from Imaris, including its full suite of shape statistics. The first steps of the script use information encoded
# in the file name to allocate categorical variables (Stage, AP position etc). In this case, file names include a stage (6ss - 14ss), an AP position (0% - 100%) and DV layer
# (D, M, V). eg. '20191016 12ss Embryo 2 100% M_Detailed.csv'. All files are then read in together, and compiled into a single dataframe.

#

library(dplyr)
library(data.table)
library(ggplot2)
library(plotly)
library(tidyverse)
library(factoextra)
library(viridis)
library(reshape2)


########## 1. Import all cell shape data in .csv format ##########

NotochordlistM <- list() #Import notochord cells classified as 'middle layer'
listcsvM <- dir(pattern = "* M_")
for (k in 1:length(listcsvM)){
  
  NotochordlistM[[k]] <- read.csv(listcsvM[k], skip=3)
  
  NotochordlistM[[k]]$Identity <- listcsvM[k]
  
}

NotochordlistD <- list() #Import notochord cells classified as 'dorsal (Müller) layer'
listcsvD <- dir(pattern = "* D_")
for (k in 1:length(listcsvD)){
  
  NotochordlistD[[k]] <- read.csv(listcsvD[k], skip=3)
  
  NotochordlistD[[k]]$Identity <- listcsvD[k]
  
}

NotochordlistV <- list() #Import notochord cells classified as 'ventral (Müller) layer'
listcsvV <- dir(pattern = "* V_")
for (k in 1:length(listcsvV)){
  
  NotochordlistV[[k]] <- read.csv(listcsvV[k], skip=3)
  
  NotochordlistV[[k]]$Identity <- listcsvV[k]
  
}

########## 2. Classify DV layers ##########

NotochordDFD <- rbindlist(NotochordlistD, fill = TRUE) #Group data into dorsal, middle and ventral dataframes and classify layers as '1', '2', '3' respectively
NotochordDFD$Layer <- '1'

NotochordDFM <- rbindlist(NotochordlistM, fill = TRUE)
NotochordDFM$Layer <- '2'

NotochordDFV <- rbindlist(NotochordlistV, fill = TRUE)
NotochordDFV$Layer <- '3'


Data <- rbind(NotochordDFM, NotochordDFD, NotochordDFV, fill = TRUE) #Merge dorsal, middle and ventral dataframes
  

Data$Unit <- NULL #Delete unused columns
Data$Category <- NULL
Data$Collection <- NULL
Data$Image <- NULL
Data$Time <- NULL
Data$X <- NULL
Data$ReferenceFrame <- NULL

########## 3. Classify AP positions ##########

Notochord0 <- Data[grep(" 0%", Data$Identity), ]
Notochord0$Region <- '0'

Notochord25 <- Data[grep("25%", Data$Identity), ]
Notochord25$Region <- '25'

Notochord50 <- Data[grep("50%", Data$Identity), ]
Notochord50$Region <- '50'

Notochord75 <- Data[grep("75%", Data$Identity), ]
Notochord75$Region <- '75'

Notochord100 <- Data[grep("100%", Data$Identity), ]
Notochord100$Region <- '100'

NotochordALL <- Data[grep("ALL", Data$Identity), ]
NotochordALL$Region <- 'All'

Data <- rbind(Notochord0,
              Notochord25,
              Notochord50,
             Notochord100,
             NotochordALL, fill = TRUE)


########## 4. Classify somite stages ##########


Notochord6 <- Data[grep("6ss", Data$Identity), ]
Notochord6$Stage <- '6'

Notochord8 <- Data[grep("8ss", Data$Identity), ]
Notochord8$Stage <- '8'

Notochord10 <- Data[grep("10ss", Data$Identity), ]
Notochord10$Stage <- '10'

Notochord12 <- Data[grep("12ss", Data$Identity), ]
Notochord12$Stage <- '12'

Notochord14 <- Data[grep("14ss", Data$Identity), ]
Notochord14$Stage <- '14'

Data <- rbind(Notochord6,
              Notochord8,
              Notochord10,
              Notochord12,
              Notochord14)


########## 5. Isolate useful shape variables ##########

Area <-subset(Data, Variable=="Area")
BB_AA_X <-subset(Data, Variable=="BoundingBoxAA Length X")
BB_AA_Y <-subset(Data, Variable=="BoundingBoxAA Length Y")
BB_AA_Z <-subset(Data, Variable=="BoundingBoxAA Length Z")
BB_OO_A <-subset(Data, Variable=="BoundingBoxOO Length A")
BB_OO_B <-subset(Data, Variable=="BoundingBoxOO Length B")
BB_OO_C <-subset(Data, Variable=="BoundingBoxOO Length C")
CentreHomogenousMassX <-subset(Data, Variable=="Center of Homogeneous Mass X")
CentreHomogenousMassY <-subset(Data, Variable=="Center of Homogeneous Mass Y")
CentreHomogenousMassZ <-subset(Data, Variable=="Center of Homogeneous Mass Z")
EllipsoidA_X <-subset(Data, Variable=="Ellipsoid Axis A X")
EllipsoidA_Y <-subset(Data, Variable=="Ellipsoid Axis A Y")
EllipsoidA_Z <-subset(Data, Variable=="Ellipsoid Axis A Z")
EllipsoidB_X <-subset(Data, Variable=="Ellipsoid Axis B X")
EllipsoidB_Y <-subset(Data, Variable=="Ellipsoid Axis B Y")
EllipsoidB_Z <-subset(Data, Variable=="Ellipsoid Axis B Z")
EllipsoidC_X <-subset(Data, Variable=="Ellipsoid Axis C X")
EllipsoidC_Y <-subset(Data, Variable=="Ellipsoid Axis C Y")
EllipsoidC_Z <-subset(Data, Variable=="Ellipsoid Axis C Z")
EllipsoidLength_A <-subset(Data, Variable=="Ellipsoid Axis Length A")
EllipsoidLength_B <-subset(Data, Variable=="Ellipsoid Axis Length B")
EllipsoidLength_C <-subset(Data, Variable=="Ellipsoid Axis Length C")
Ellipticity_Oblate <-subset(Data, Variable=="Ellipticity (oblate)")
Ellipticity_Prolate <-subset(Data, Variable=="Ellipticity (prolate)")
Triangles <-subset(Data, Variable=="Number of Triangles")
X <-subset(Data, Variable=="Position X")
Y <-subset(Data, Variable=="Position Y")
Z <-subset(Data, Variable=="Position Z")
Sphericity <-subset(Data, Variable=="Sphericity")
Volume <-subset(Data, Variable=="Volume")

DAPI <- subset(Data, Channel=='1')
DAPI_Centre_X <- subset(DAPI, Variable=='Center of Image Mass X')
DAPI_Centre_Y <- subset(DAPI, Variable=='Center of Image Mass Y')
DAPI_Centre_Z <- subset(DAPI, Variable=='Center of Image Mass Z')

Sum <- rbind(Area,
             BB_AA_X,
             BB_AA_Y,
             BB_AA_Z,
             BB_OO_A,
             BB_OO_B,
             BB_OO_C,
             CentreHomogenousMassX,
             CentreHomogenousMassY,
             CentreHomogenousMassZ,
             EllipsoidA_X,
             EllipsoidA_Y,
             EllipsoidA_Z,
             EllipsoidB_X,
             EllipsoidB_Y,
             EllipsoidB_Z,
             EllipsoidC_X,
             EllipsoidC_Y,
             EllipsoidC_Z,
             EllipsoidLength_A,
             EllipsoidLength_B,
             EllipsoidLength_C,
             Ellipticity_Oblate,
             Ellipticity_Prolate,
             Triangles,
             X,
             Y,
             Z,
             Sphericity,
             Volume, 
             DAPI_Centre_X,
             DAPI_Centre_Y,
             DAPI_Centre_Z)

Sum$Channel <- NULL #Remove unused columns


########## 6. Convert dataframe to wide format, with each shape variable as a distinct column ##########

DataW <- dcast(Sum, ID + Identity + Layer + Region + Stage ~ Variable, value.var='Value')

colnames(DataW) <- c("ID",
                     "Identity",
                     "Layer",
                     "Region",
                     "Stage",
                     "Area",
                     "BoundingBoxAA_X",
                     "BoundingBoxAA_Y",
                     "BoundingBoxAA_Z",
                     "BoundingBoxOO_A",
                     "BoundingBoxOO_B",
                     "BoundingBoxOO_C",
                     "CenterHomogeneousMass_X",
                     "CenterHomogeneousMass_Y",
                     "CenterHomogeneousMass_Z",
                     "CenterImageMassX",
                     "CenterImageMassY",
                     "CenterImageMassZ",
                     "EllipsoidAxisA_X",
                     "EllipsoidAxisA_Y",
                     "EllipsoidAxisA_Z",
                     "EllipsoidAxisB_X",
                     "EllipsoidAxisB_Y",
                     "EllipsoidAxisB_Z",
                     "Major_ellipsoid_axis_AP",
                     "Major_ellipsoid_axis_DV",
                     "Major_ellipsoid_axis_ML",
                     "EllipsoidAxisLength_A",
                     "EllipsoidAxisLength_B",
                     "EllipsoidAxisLength_C",
                     "Ellipticity_oblate",
                     "Ellipticity_prolate",
                     "Number_of_Triangles",
                     "PositionX",
                     "PositionY",
                     "PositionZ",
                     "Sphericity",
                     "Volume") #Label columns appropriately



########## 7. Calculate further shape metrics ##########

DataW$EllipsoidAxisA_X <- abs(DataW$EllipsoidAxisA_X) #Use absolute values for ellipsoid axes
DataW$EllipsoidAxisA_Y <- abs(DataW$EllipsoidAxisA_Y)
DataW$EllipsoidAxisA_Z <- abs(DataW$EllipsoidAxisA_Z)
DataW$EllipsoidAxisB_X <- abs(DataW$EllipsoidAxisB_X)
DataW$EllipsoidAxisB_Y <- abs(DataW$EllipsoidAxisB_Y)
DataW$EllipsoidAxisB_Z <- abs(DataW$EllipsoidAxisB_Z)
DataW$Major_ellipsoid_axis_AP <- abs(DataW$Major_ellipsoid_axis_AP)
DataW$Major_ellipsoid_axis_DV <- abs(DataW$Major_ellipsoid_axis_DV)
DataW$Major_ellipsoid_axis_ML <- abs(DataW$Major_ellipsoid_axis_ML)

DataW$APvsDV <- DataW$BoundingBoxAA_X / DataW$BoundingBoxAA_Y #Calculate AP:DV, AP:ML, DV:ML ratios (uncouples cell shape from cell size)
DataW$APvsML <- DataW$BoundingBoxAA_X / DataW$BoundingBoxAA_Z
DataW$DVvsML <- DataW$BoundingBoxAA_Y / DataW$BoundingBoxAA_Z

DataW$MeanLat <- ((DataW$BoundingBoxAA_Z + DataW$BoundingBoxAA_Y) / 2) #Calculate mean cell diameter from ML and DV bounding box lengths

DataW$Elongation <- DataW$BoundingBoxAA_X / DataW$MeanLat #Calculate cell elongation/AP anisotropy using AP length and mean diameter

DataW$NuclearDispX <- abs(DataW$CenterImageMassX - DataW$PositionX) #Calculate true distance of nucleus from centre of homogenous mass on X, Y, Z planes
DataW$NuclearDispY <- abs(DataW$CenterImageMassY - DataW$PositionY)
DataW$NuclearDispZ <- abs(DataW$CenterImageMassZ - DataW$PositionZ)

DataW$SurfaceConvolution <- DataW$Area / DataW$Volume #Calculate surface area:volume ratio

DataW$Cuboidness <- (DataW$Volume / (DataW$BoundingBoxOO_A * DataW$BoundingBoxOO_B * DataW$BoundingBoxOO_C)) #Calculate cuboidness

DataW$Flatness <- DataW$BoundingBoxOO_B/DataW$BoundingBoxOO_A #Calculate flatness

DataW$Sphericity <- (DataW$Volume)^(2/3) * 6 #Calculate sphericity (3 steps)
DataW$Sphericity <- DataW$Sphericity * pi^(1/3)
DataW$Sphericity <- DataW$Sphericity / DataW$Area

DataW <- as.data.frame(DataW)


########## 7. Screen for redundancy ##########

library(Hmisc) #Quantify correlation with rcorr()

Variables <- (DataW[,c(6, 10:12, 25:27, 37, 38, 39, 40, 43, 44, 45, 46)]) #Include column numbers of variables to screen

rcorr(as.matrix(Variables)) #Gives correlation scores

library(corrplot) #Visualise correlation with corrplot()

Variables <- cor(DataW[,c(6, 10:12, 25:27, 37, 38, 39, 40, 43, 44, 45, 46)]) #Include column numbers of variables to screen
round(Variables, 2)

corrplot(Variables, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)

########## 8. Perform PCA ##########
  
library(factoextra)
library(FactoMineR)
library(viridis)

res.pca <- prcomp(DataW[,c(25:27, 33, 37, 38, 39, 40, 41, 43, 44, 45, 46, 48, 49)], scale = TRUE) #Feed shape metrics to PCA (no ID, position, type)

fviz_eig(res.pca, addlabels = TRUE, fill = 'black') #Show scree plot (PC weightings)

fviz_contrib(res.pca, choice = "var", axes = 5, top = 10) #Show correlates of specific PCs ('axes' defines PC of interest)

quali.sup <- as.numeric(DataW[, 1]) #Set DV position as colour code for plotting PCA results

fviz_pca_ind(res.pca, #Plot PCA with colour code for DV position
             col.ind = quali.sup,
             addEllipses = FALSE, 
             alpha = 1,
             label = FALSE) +
  scale_colour_viridis() +
  theme_minimal() +
  coord_fixed() 

fviz_pca_biplot(res.pca, axes = c(2, 3), col.ind = quali.sup, #Plot compass plot with PCA in background
                addEllipses = FALSE, 
                alpha = 0.2, 
                label = "var", 
                repel = TRUE, 
                colour = 'grey') +
  scale_colour_viridis() +
  theme_minimal() +
  coord_fixed()


########## 9. Extract PC coordinates and merge with raw values ##########

PCs <- as.data.frame(res.pca$x)
AllData <- cbind(PCs, DataW)

AllData$Layer <- as.numeric(AllData$Layer)
AllData$Stage <- as.numeric(AllData$Stage)

########## 10. Plot PCA in ggplot ##########

ggplot(AllData, aes(x = PC1, y = PC2, colour = PC3)) +
  geom_point() +
  scale_colour_viridis(option = 'magma') +
  theme_minimal() +
  coord_fixed() 

########## 11. Filter data by stage ##########

ggplot(AllData, aes(x = PC1, y = PC2, colour = Layer)) +
  geom_point() +
  scale_colour_viridis() +
  theme_minimal() +
  coord_fixed() +
  facet_wrap(~Stage)

########## 12. Follow region-specific population over time ##########


Background6 <- AllData
Background6$Stage <- 6
Background8 <- AllData
Background8$Stage <- 8
Background10 <- AllData
Background10$Stage <- 10
Background12 <- AllData
Background12$Stage <- 12
Background14 <- AllData
Background14$Stage <- 14
Background <- rbind(Background6, 
                    Background8,
                    Background10,
                    Background12,
                    Background14)

sample <- subset(AllData, Region==50 & Layer==2)

sample$Stage <- as.factor(sample$Stage)
sample$Stage <- factor(sample$Stage, levels = c("6", "8", "10", "12", "14"))

ggplot(Background, aes(x = PC1, y = PC2)) + #Plot cells from a specific region/cell layer over time
  geom_point(colour = 'grey90', alpha = 0.2) +
  geom_point(data = sample, aes(x = PC1, y = PC2, colour = Stage)) +
  scale_colour_brewer(palette = "Reds") +
  theme_minimal() +
  coord_fixed()



########## 13. Infer pseudotemporal axis ##########

library('slingshot')

sample <- subset(AllData, Region==50 & Layer==2) #Adjust subsetting for preference of position/cell layer

sample$Stage <- as.numeric(sample$Stage)

sample <- as.matrix(sapply(sample, as.numeric))  
rownames(sample) <- paste0('C',1:344)

PC <- sampleDF[,1:3]
PCR <- as.matrix(sapply(PC, as.numeric))  
rownames(PCR) <- paste0('C',1:344)

cl <- sampleDF$Stage
colData(sim)$Stage <- cl

sim <- slingshot(sim, clusterLabels = 'Stage', PCR)

colors <- colorRampPalette(brewer.pal(5,'Reds')[-2])(100)
plotcol <- colors[cut(cl, breaks=100)]

lin1 <- getLineages(PCR, cl, start.clus = '3')

plot(PCR[,1:2], col = cl, asp = 1, pch = 16) +
  lines(lin1, lwd = 3, col = 'black')

crv1 <- getCurves(lin1, thresh = 1)

plot(PCR[,1:2], col = plotcol, asp = 1, pch = 16) +
lines(crv1, lwd = 3, col = 'black')

p <- plot(AllData[,1:2], col = 'gray95', pch=16, asp = 1, cex = 0.8, xlim=range(-2:5), ylim=range(-2.5, 8))
+ points(PCR[,1:2], col = plotcol, pch=16, asp = 1, cex = 0.9, xlim=range(-2:5), ylim=range(-2.5, 8))
+ lines(crv1, lwd=2, col='black', dims = 1:3)
+ lines(lin1, lwd = 2, cex = 1.2, col = 'blue3', dims = 1:2)

########## 14. Plot shape transformations across pseudotemporal axis ##########

TrajVec <- as.vector(slingCurves(crv1)[[1]]$ord)

Traj <- sampleDF[sampleDF$IDs %in% TrajVec,] 

sampleDFT <- sampleDF[match(TrajVec, sampleDF$IDs),]

sampleDFT$traj <- seq.int(nrow(sampleDFT))

sampleDFT$trajN <- sampleDFT$traj / max(sampleDFT$traj)

ggplot(sampleDFT, aes(x = trajN, y = Volume)) +
  geom_smooth(span = 1, colour = 'black') +
  theme_minimal() 

########## 15. Export pseudotime colour code for Imaris heatmaps ##########

sub12 <- subset(AllData, Stage==12) #Isolate embryos of a given stage (6, 8, 10, 12, 14)

unique(sub12$Identity) #Find names of individuals of the given stage - copy one into line 432.

one12 <-  sub12[grep("20191016 12ss Embryo 2 ALL", sub12$Identity), ] #Isolate cells for one embryo

one12$ReducedX <- one12$PositionX - min(one12$PositionX) #Normalise X axis (2 steps)
one12$NormalisedX <- one12$ReducedX / max(one12$ReducedX)

one12$Slice <- cut(one12$PositionZ, breaks = seq(0, 85.78, #Max Z position
                                                 by = 0.739), #Z step size
                   labels = 1:116) #Number of slices

var_list <- split(one12, one12$Slice) #Makes list of slices to generate images for

plot_list = list() #Makes dot plots for each slice, with dot intensity scaled to normalised pseudotime position 
for (i in 1:116) { #Change according to Z slice number
  p = ggplot(var_list[[i]], aes(x = PositionX, y = PositionY, alpha = trajN)) +
    geom_point(size = 5, shape = 16) + #Optimise for cell size - points must lie inside each cell
    scale_alpha(range = c(0.0029, 0.753)) + #Defines the range for pseudotime values. Use (0, 1) for a single embryo, or limits when normalised across population for comparison between embryos
    theme_bw() +
    xlim(0, 345.94) + #Image dimensions - must correspond to raw data!
    ylim(151.41, 0) +
    scale_x_continuous(limits = c(0, 345.94), expand = c(0, 0)) +
    scale_y_reverse(limits = c(151.41, 0), expand = c(0, 0)) +
    coord_fixed() +
    theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),legend.position="none",
          panel.background=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank())
  
  plot_list[[i]] = p
}


# Save plots to tiff. Makes a separate file for each plot.
for (i in 1:116) { #Adjust according to slice number - 1 graph per slice
  file_name = paste("12ss_pseudotime", i, ".tiff", sep="")
  tiff(file_name, units="px", width=(738), height=(323), res=50) #Change according to image pixel dimensions - refer again to raw data
  print(plot_list[[i]])
  dev.off()
  
}


