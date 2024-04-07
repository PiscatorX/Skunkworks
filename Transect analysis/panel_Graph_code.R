library(readxl)
library(ggplot2)
library(dplyr)
#install this package
#https://teunbrand.github.io/ggh4x/index.html
library(ggh4x)


#Graphs:
#Import Datasets:
CHNS_Dataset_SM_INT <- readxl::read_excel("CHNS_Dataset_SM&INT.xlsx")

MasterDataFrame2.0 <- CHNS_Dataset_SM_INT %>%
  select(c("Habitat", "Veg", "Site", "Transect",
           "Distance", "Wet_mass", "Dry_mass", "Volume", 
           "Depth_interval","C_percent", 
           "N_percent")) %>%
           mutate(DBD = Dry_mass / Volume) %>%
           #This orders Veg in a way we want them to appear on the graph 
           dplyr::mutate(Veg = factor(Veg, levels = c("Saltmarsh","IntertidalSeagrass")),
                        Site = factor(Site, levels = c("Upper","Middle","Lower"))) 
           
            
    
#Code from Andrew for distance sampled:
#Set the distance for saltmarsh as negative
MasterDataFrame2.0$Distance[MasterDataFrame2.0$Veg == "Saltmarsh"] <- -MasterDataFrame2.0$Distance[MasterDataFrame2.0$Veg == "Saltmarsh"]

unique(MasterDataFrame2.0$Site)


#setup colour for vegetation types
Veg_cols <- c("Orange4", "Springgreen1")
names(Veg_cols) <- c("Saltmarsh", "IntertidalSeagrass")
#labels for legend
veg_labels = c("Salt marsh", "Seagrass")
names(veg_labels) <- c("Saltmarsh", "IntertidalSeagrass")

veg_labels_INT = c("Intertidal seagrass")
names(veg_labels_INT) <- c("IntertidalSeagrass")

legend_labels = c('Upper','Middle','Lower',veg_labels)
names(legend_labels)[1:3] <- c('Upper','Middle','Lower')


#Transect:
#SOC:

Transect_Graph <- ggplot(data = MasterDataFrame2.0, aes(x = Distance, y = C_percent, shape = Transect, colour  = Veg, group = Transect),) + 
  geom_line(aes(color=Transect), alpha = 0.75, linewidth = 1) +
  geom_point(size = 3,  alpha = 0.5) +
  labs(y = "SOC (%)", x  = "Distance from seagrass edge (m)") +
  scale_colour_manual(values = Veg_cols, labels = veg_labels, , name = "Vegetation type") +
  theme(panel.grid.minor = element_blank(),
        legend.key.size = unit(1,'cm'),
        panel.background = element_blank(),
        strip.background = element_rect(colour = "black", fill=NA, size=1),
        strip.text = element_text(size = 12, colour = "black"),
        axis.text = element_text(size = 12, colour = "#434343"),
        axis.title = element_text(size = 12,colour = "black"),
        panel.spacing.x = unit(0.75, "lines"),
        legend.text = element_text(colour = "black"),
        legend.title = element_text(colour = "black"),
        legend.background  = element_rect(colour = "black", fill=NA, size = 0.75),
        panel.border = element_rect(colour = "black", fill=NA, linewidth = 0.75)) +
  facet_wrap(~Site + Veg, scales = "free", nrow = 3, ncol = 2, labeller = as_labeller(legend_labels)) + 
facetted_pos_scales(
    y = list(
      Veg == "IntertidalSeagrass" ~ scale_y_continuous(limits = c(0, 2)),
      Veg == "Saltmarsh" ~ scale_y_continuous(limits = c(0, 18))

)
)
Transect_Graph

ggsave("Transect_Graph.pdf")

