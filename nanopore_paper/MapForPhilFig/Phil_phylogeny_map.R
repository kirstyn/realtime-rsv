################################################################################
# Brunker et al 2020, Wellcome Open Research
# Map Philippines by region                        
################################################################################

rm(list=ls())
# Note: set working directory to wellcome2020-paper/WGS_tree_Easia_SEA4/MapForPhilFig
# Load libraries
library(dplyr)
library(rgdal)
library(ggplot2)
library(broom)
library(rgeos)
library(prettymapr)
library(RColorBrewer)
library(wesanderson)
# Load data
phl_province <- readOGR("PHL", "PHL_province")
phl_province_data <- fortify(phl_province) # Process shapefiles into dataframes - needed for using ggplot2 on spatial data!

# 83 provinces are split into 17 regions - but the regions need adding!
admin <- read.csv("PHL/phl_adminboundaries.csv"); head(admin)

test1 = match(admin$admin2Name_en, phl_province@data$NAME_1)
admin$admin2Name_en[which(is.na(test1))]

test2 = match(phl_province@data$NAME_1, admin$admin2Name_en)
phl_province@data$NAME_1[which(is.na(test2))]
phl_province@data$NAME_1 <- as.character(phl_province@data$NAME_1 )
phl_province@data$NAME_1[grep("Metropolitan Manila", phl_province@data$NAME_1)] <- "NCR, City of Manila, First District"
phl_province@data$NAME_1[grep("North Cotabato", phl_province@data$NAME_1)] <- "Cotabato City"
phl_province@data$NAME_1[grep("Shariff Kabunsuan", phl_province@data$NAME_1)] <- "Maguindanao"

phl_province@data$Province <- admin$admin2Name_en[match(phl_province@data$NAME_1, admin$admin2Name_en)]
phl_province@data$Province <- as.character(phl_province@data$Province)
phl_province@data$Province[grep("NCR, City of Manila", phl_province@data$NAME_1)] <- "NCR, City of Manila"
phl_province@data$Province[grep("Cotabato", phl_province@data$NAME_1)] <- "Region XII"
phl_province@data$Province[grep("Maguindanao", phl_province@data$NAME_1)] <- "Autonomous Region in Muslim Mindanao "

phl_province@data$REGION <- admin$admin1Name_en[match(phl_province@data$NAME_1, admin$admin2Name_en)]
phl_province@data$reg <- "None"
phl_province@data$reg[which(phl_province@data$REGION == "Region I ")] <- "1"
phl_province@data$reg[which(phl_province@data$REGION == "Region III ")] <- "3"
phl_province@data$reg[which(phl_province@data$REGION == "Region IV-A ")] <- "4A"
phl_province@data$reg[which(phl_province@data$REGION == "Region IV-B ")] <- "4B"
phl_province@data$reg[which(phl_province@data$REGION == "Region V ")] <- "5"
phl_province@data$reg[which(phl_province@data$REGION == "National Capital Region ")] <- "NCR"
regions <- sort(unique(phl_province@data$reg))
phl_province$reg = factor(phl_province@data$reg, levels=regions)

#-------------------------------------------------------------------------------
# Create Philippines map
gCentroid(phl_province)
bbox(phl_province)

# Colour scheme
reg_cases = c("1", "3", "4A", "4B", "5", "NCR", "None")
pal=wes_palette(7, name = "FantasticFox1", type = "continuous")
reg_colours = c(pal[c(1:4,6:7)], "transparent")

# pdf("figs/PhilSeqMap.pdf", height=8, width=6)
plot(phl_province, border="darkgrey", col = "grey", lwd=.5, xlim = c(120,125), ylim = c(6,19))
plot(phl_province, col=reg_colours[match(phl_province@data$reg, reg_cases)],  border="darkgrey", lwd=.5, add=TRUE)
#legend("topright", legend=reg_cases[1:6], fill=reg_colours[1:6], bty = "n", title="Region", cex=1.3)
addnortharrow(pos="bottomleft", scale=0.8)
# dev.off()



