################################################################################
#                           Map Philippines by region                         #
################################################################################
rm(list=ls())

# If ggplot2 hangs - reinstall!:
# deps <- tools::package_dependencies("ggplot2", recursive = TRUE)$ggplot2
# for (dep in deps)
#   try(install.packages(dep))

# Load libraries
library(dplyr)
library(rgdal)
library(ggplot2)
library(broom)
library(rgeos)
library(prettymapr)
library(RColorBrewer)
# Load data
phl_province <- readOGR("PHL", "PHL_province")
phl_province_data <- fortify(phl_province) # Process shapefiles into dataframes - needed for using ggplot2 on spatial data!

# 83 provinces are split into 17 regions - but the regions need adding!
admin <- read.csv("PHL/phl_adminboundaries.csv"); head(admin)

test1 = match(admin$admin2Name_en, phl_province@data$NAME_1)
admin$admin2Name_en[which(is.na(test1))]
# "Davao" = "Region XI" # Davao Occidental
# "Isabela" = "Region II" # City of Isabela, Cagayan Valley
# "Cotabato" = "Region XII" # Soccsksargen (North Cotabato) # Cotabato; Cotabato City
# "NCR" = "National Capital Region" # NCR, City of Manila, First District; NCR, Fourth District; NCR, Second District; NCR, Third District
# "Shariff Kabunsuan" = "Autonomous Region in Muslim Mindanao "

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
#reg_colours = c(brewer.pal(6,"Spectral"), "transparent")
#reg_colours  <- c(wes_palette(6, name = "FantasticFox1", type = "continuous"), "transparent")
pal=wes_palette(7, name = "FantasticFox1", type = "continuous")
reg_colours = c(pal[c(1:4,6:7)], "transparent")

pdf("PhilSeqMap3.pdf", height=8, width=6)
plot(phl_province, border="darkgrey", col = "grey", lwd=.5, xlim = c(120,125), ylim = c(6,19))
plot(phl_province, col=reg_colours[match(phl_province@data$reg, reg_cases)],  border="darkgrey", lwd=.5, add=TRUE)
legend("topright", legend=reg_cases[1:6], fill=reg_colours[1:6], bty = "n", title="Region", cex=1.3)
addnortharrow(pos="bottomleft", scale=0.8)
dev.off()


# Alternative in gglpot
gg <- ggplot(data = phl_province) # Make gg object
values <- data.frame("group" = factor(1:length(phl_province$reg)), "color" = color) # Make auxilliary data.frame
gg$data <- merge(gg$data, values, by = c("group")) # Add color data to gg object by merging
# Plot gg object

gg + 
  geom_polygon(data=phl_province, aes(x=long, y=lat), fill="grey", col="white") +
  geom_polygon(aes(x = long, y = lat, group = group),  fill = gg$data$color, colour = 1) + 
  theme_void() +
  theme(legend.title = element_blank()) +
  coord_equal()
ggsave("test2.pdf", height=8, width=8)


deps <- tools::package_dependencies("ggplot2", recursive = TRUE)$ggplot2
for (dep in deps)
  try(install.packages(dep))


ggplot() +
  geom_polygon(data=phl_province, aes(x=long, y=lat), fill="grey", col="white") +
#  geom_polygon(data=phl_province, aes(x=long, y=lat, group=reg, fill = reg_colours), col="white", size=0.1) +
  theme_void() +
  theme(legend.title = element_blank()) +
  coord_equal()
ggsave("test2.pdf", height=8, width=8)


# Subset data
prov_seq_data <- phl_province[which(phl_province$reg!="None"),]

# Set factor levels
prov_seq <- sort(unique(prov_seq_data$reg))
prov_seq <- drop.levels(prov_seq, reorder=TRUE)
prov_seq_data$regions <- factor(prov_seq_data$reg, levels=prov_seq)

ggplot() +
  geom_polygon(data=phl_province, aes(x=long, y=lat, group=group, fill=factor(group)), col="white") +
  # geom_polygon(data=prov_seq_data, aes(x=long, y=lat, group=group, fill=regions), col="white", size=0.1) +
  #scale_fill_manual(values=viridis(n=length(prov_seq)), na.value="grey", breaks=prov_seq, labels=prov_seq) +
  theme_void() +
  theme(legend.title = element_blank()) +
  coord_equal()
ggsave("test2.pdf", height=8, width=8)





