#Analysis for: Bacterial, but not fungal, communities show spatial heterogeneity in European beech (Fagus sylvatica L.) deadwood
#Jason Bosch

#####SET UP#####

###Reproducibility###

library(groundhog)

groundhog_day <- "2022-09-01"
meta.groundhog(groundhog_day)

###Libraries###

#CRAN libraries
groundhog.library(c("vegan",
                    "stringr",
                    "ggplot2",
                    "ggsignif",
                    "ggrepel",
                    "RColorBrewer",
                    "gridExtra",
                    "rnaturalearth",
                    "rnaturalearthdata",
                    "ggspatial",
                    "sf",
                    "cluster",
                    "NbClust",
                    "ggdendro",
                    "ggvenn",
                    "shipunov"),
                  groundhog_day)

#GIT Libraries
groundhog.library(c("github::JLSteenwyk/ggpubfigs",
                    "github::NicolasH2/ggdendroplot"),
                  groundhog_day)

# #Problem installing with groundhog
library(NetCoMi)

#BioConductor Libraries
library(phyloseq)

###Directory###

setwd("~/PostDoc/02_Pojects/01_Microbial_co-occurrence/04_Analysis_Results/R_FINAL/")

###Functions###

#Import QIIME2 taxonomy function
import_qiime2_taxonomy <- function(qiime2_taxonomy) {
  taxonomy_cleaning <- as.data.frame(qiime2_taxonomy[,1])
  taxonomy_cleaning[,1] <- gsub("[a-z]__","",taxonomy_cleaning[,1])
  taxonomy_split <- as.data.frame(matrix(nrow = 1,ncol = 7))
  for (n in 1:nrow(taxonomy_cleaning)) {
    taxonomy_split[n,] <- str_split_fixed(taxonomy_cleaning[n,1],"; |;", n = 7)
  }
  colnames(taxonomy_split) <- c("Domain","Phylum","Class","Order","Family","Genus","Species")
  rownames(taxonomy_split) <- rownames(qiime2_taxonomy)
  tax_table(as.matrix(taxonomy_split))
}

#####LOAD DATA#####

#Bacteria

#Import the data
bact_asv_table_raw <- read.table("~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Bacteria_POOLED/export/OTU.tsv", 
                                 header = TRUE, row.names = 1, sep = "\t", comment.char = "#")
bact_tree_final <-  read_tree(treefile = "~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Bacteria_POOLED/export/tree.nwk", 
                              errorIfNULL = FALSE)
bact_QIIME2_taxonomy_final <- import_qiime2_taxonomy(read.table("~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Bacteria_POOLED/export/taxonomy.tsv", header = TRUE, row.names = 1, sep = "\t"))
bact_metadata_raw <- read.table("~/PostDoc/02_Pojects/01_Microbial_co-occurrence/06_Miscellaneous/Metadata.tsv", 
                                header = TRUE, row.names = 1, sep = "\t")
bact_metadata_final <- sample_data(bact_metadata_raw)

#Filtering
#Remove any ASV with fewer 0.1% of the mean number of reads per sample
bact_asv_final <- otu_table(bact_asv_table_raw[(rowSums(bact_asv_table_raw)>round(mean(colSums(bact_asv_table_raw))*0.001)),],taxa_are_rows = TRUE)
#Removed 4097 taxa as too rare, kept 3521

#Create phyloseq object
#No tree to enable merging with fungal dataset later
physeq_bact <- phyloseq(bact_asv_final, bact_metadata_final, bact_QIIME2_taxonomy_final)

#Filter to remove Eukrayotic sequences as well as chloroplasts and mitochondria and unassigned
#Without tree
physeq_bact <- subset_taxa(physeq_bact, tax_table(physeq_bact)[,"Domain"]!="Eukaryota")
physeq_bact <- subset_taxa(physeq_bact, tax_table(physeq_bact)[,"Domain"]!="Unassigned")
table(as.vector(physeq_bact@tax_table[,"Genus"]=="Mitochondria"))
physeq_bact <- subset_taxa(physeq_bact, tax_table(physeq_bact)[,"Genus"] != "Mitochondria")
table(as.vector(physeq_bact@tax_table[,"Genus"]=="Chloroplast"))
physeq_bact <- subset_taxa(physeq_bact, tax_table(physeq_bact)[,"Genus"] != "Chloroplast")
table(physeq_bact@tax_table[,"Domain"])
#Removed 224 taxa; kept 1 Archaea, 3296 Bacteria

#Fungi

#Import the data
fung_asv_table_raw <- read.table("~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Fungi_POOLED//export/OTU.tsv", 
                                 header = TRUE, row.names = 1, sep = "\t", comment.char = "#")
fung_QIIME2_taxonomy_final <- import_qiime2_taxonomy(read.table("~/PostDoc/02_Pojects/01_Microbial_co-occurrence/02_Cleaned_Data/QIIME_Fungi_POOLED//export/taxonomy.tsv", header = TRUE, row.names = 1, sep = "\t"))
fung_metadata_raw <- read.table("~/PostDoc/02_Pojects/01_Microbial_co-occurrence/06_Miscellaneous/Metadata.tsv", 
                                header = TRUE, row.names = 1, sep = "\t")
fung_metadata_final <- sample_data(fung_metadata_raw)

#Filtering
#Remove any ASV with fewer 0.1% of the mean number of reads per sample
fung_asv_final <- otu_table(fung_asv_table_raw[(rowSums(fung_asv_table_raw)>round(mean(colSums(fung_asv_table_raw))*0.001)),],taxa_are_rows = TRUE)
#Removed 2327 taxa as too rare, kept 956

#Create phyloseq object
physeq_fung <- phyloseq(fung_asv_final, fung_metadata_final, fung_QIIME2_taxonomy_final)

#All detected taxa are fungi, so no need for further filtering steps
table(physeq_fung@tax_table[,"Domain"])

#####ANALYSIS#####

###Constant colours for figures###
scale_colours <- friendly_pal("contrast_three",2)
names(scale_colours) <- unique(physeq_bact@sam_data$Scale)
tree_colours <- brewer.pal(10, "Paired")
names(tree_colours) <- unique(physeq_bact@sam_data$Tree)

###Study Site Map###

#Load data
world <- ne_countries(scale = "medium", returnclass = "sf")
roads <- ne_download(scale = "large",category = "cultural",type = "roads",returnclass = "sf")

#Remove unnecessary roads
roads <- st_crop(roads, xmin = 11.50, xmax = 19.25, ymin = 48.25, ymax = 51.25)

#Extract locations and labels
world_points<- st_centroid(st_make_valid(world)) #Make valid fixes some weird broken geometry but the fix also messes up the plotting
world_points <- cbind(st_make_valid(world), st_coordinates(st_centroid(st_make_valid(world)$geometry)))
world_points_used <- world_points[world_points$name%in%c("Czech Rep.","Austria","Germany","Poland","Slovakia"),c("name","Y","X")]
world_points_used$geometry <- rep("black",5)
colnames(world_points_used) <- c("Name","Lat","Lon","Colour")

#To aid in labelling, all points should be in one object.
map_objects <- as.data.frame(matrix(data = c("Žofín Forest",48.66369499053803,14.706172165867045,"red",
                                             "Prague",50.083333,14.416667,"black",
                                             "Brno",49.1925,16.608333,"black",
                                             "Ostrava",49.835556,18.2925,"black",
                                             "Plzeň",49.7475,13.3775,"black",
                                             "Liberec",50.766667,15.066667,"black",
                                             "Olomouc",49.593889,17.250833,"black"),
                                    ncol = 4,dimnames = list(NULL,c("Name","Lat","Lon","Colour")),byrow = TRUE))
map_objects_names <- rbind(map_objects,world_points_used)
map_objects_names$Name <- gsub("Czech Rep.","Czechia",map_objects_names$Name)
map_objects_names$Size <- c(rep(1,nrow(map_objects)),rep(2,nrow(world_points_used)))

#Generate annotated world map
sf::sf_use_s2(FALSE) #Fixes some geometry problem.
#May need to redo map to get labels in best spot
map <- 
  ggplot(data = world) +
  geom_sf(fill= "gray90") +
  geom_sf(data = roads,colour = "light grey") +
  labs( x = "Longitude", y = "Latitude") +
  coord_sf(xlim = c(11.50, 19.25), ylim = c(48.25, 51.25), expand = FALSE) +
  annotation_scale(location = "br", width_hint = 0.5) +
  annotation_north_arrow(location = "bl", which_north = "true", style = north_arrow_fancy_orienteering) +
  theme(panel.grid.major = element_line(color = "gray60", linetype = "dashed", size = 0.25), panel.background = element_rect(fill = "aliceblue")) +
  geom_point(data = map_objects, mapping = aes(x=as.numeric(Lon), y= as.numeric(Lat), color=Colour)) +
  scale_color_manual(values = c("black" = "black", "red" = "red")) +
  theme(legend.position = "None") +
  geom_text_repel(data = map_objects_names, mapping = aes(x=as.numeric(Lon), y= as.numeric(Lat), label = Name, size=Size)) +
  scale_size(range = c(4, 5.5))

#Figure SX7
ggsave("Figure_S1.svg", map, width = 8, height = 5)
ggsave("Figure_S1_Supplementary_Data.tiff", map, width = 8, height = 5, dpi = 600)
ggsave("Figure_S1_Supplementary_Data.jpg", map, width = 8, height = 5, dpi = 400)


###Alpha diversity###

#Rarefy communities
physeq_bact_rarefied <- rarefy_even_depth(physeq_bact,sample.size = min(colSums(physeq_bact@otu_table)),rngseed = 123456, replace = FALSE)
#410 ASVs were removed during rarefaction
physeq_fung_rarefied <- rarefy_even_depth(physeq_fung,sample.size = min(colSums(physeq_fung@otu_table)),rngseed = 123456, replace = FALSE)
#No ASVs were removed

#Calculate the alpha diversity values
alpha_diversity_bact <- estimate_richness(physeq_bact_rarefied, measures = c("Observed","Simpson","Shannon"))
alpha_diversity_fung <- estimate_richness(physeq_fung_rarefied, measures = c("Observed","Simpson","Shannon"))
#Add in metadata to allow groupings
table_alpha_diversity_bact <- merge.data.frame(alpha_diversity_bact, bact_metadata_final, by = "row.names")
table_alpha_diversity_fung <- merge.data.frame(alpha_diversity_fung, fung_metadata_final, by = "row.names")

#Draw richness plots

Alpha_observed_bacteria <- 
  ggplot(table_alpha_diversity_bact, aes(x=Scale, y=Observed)) +
  geom_boxplot(aes(colour = Scale)) +
  geom_point (aes(colour = Scale), size = 0.5, position = position_jitterdodge(jitter.width = 0.1))  +
  scale_y_continuous(limits = c(0,950),breaks = c(0,100,200,300,400,500,600,700,800,900)) +
  theme_classic() +
  scale_color_manual(values = scale_colours) +
  labs(title = "Bacteria: Richness") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("Composite","Fine")), map_signif_level=TRUE, colour="black", y_position = 890)

Alpha_observed_fungi <- 
  ggplot(table_alpha_diversity_fung, aes(x=Scale, y=Observed)) +
  geom_boxplot(aes(colour = Scale)) +
  geom_point (aes(colour = Scale), size = 0.5, position = position_jitterdodge(jitter.width = 0.1))  +
  scale_y_continuous(limits = c(0,200),breaks = c(0,25,50,75,100,125,150,175)) +
  theme_classic() +
  scale_color_manual(values = scale_colours) +
  labs(title = "Fungi: Richness") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("Composite","Fine")), map_signif_level=TRUE, colour="black", y_position = 180)

#Draw Simpsons plots

Alpha_simpson_bacteria <- 
  ggplot(table_alpha_diversity_bact, aes(x=Scale, y=Simpson)) +
  geom_boxplot(aes(colour = Scale)) +
  geom_point (aes(colour = Scale), size = 0.5, position = position_jitterdodge(jitter.width = 0.1))  +
  scale_y_continuous(limits = c(0,1.1),breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)) +
  theme_classic() +
  scale_color_manual(values = scale_colours) +
  labs(title = "Bacteria: Simpson (1-D)") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("Composite","Fine")), map_signif_level=TRUE, colour="black", y_position = 1)

Alpha_simpson_fungi <- 
  ggplot(table_alpha_diversity_fung, aes(x=Scale, y=Simpson)) +
  geom_boxplot(aes(colour = Scale)) +
  geom_point (aes(colour = Scale), size = 0.5, position = position_jitterdodge(jitter.width = 0.1))  +
  scale_y_continuous(limits = c(0,1.1),breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)) +
  theme_classic() +
  scale_color_manual(values = scale_colours) +
  labs(title = "Fungi: Simpson (1-D)") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("Composite","Fine")), map_signif_level=TRUE, colour="black", y_position = 1)

#Draw Shannon plots

Alpha_shannon_bacteria <- 
  ggplot(table_alpha_diversity_bact, aes(x=Scale, y=Shannon)) +
  geom_boxplot(aes(colour = Scale)) +
  geom_point (aes(colour = Scale), size = 0.5, position = position_jitterdodge(jitter.width = 0.1))  +
  scale_y_continuous(limits = c(0,6.7),breaks = c(0,1,2,3,4,5,6)) +
  theme_classic() +
  scale_color_manual(values = scale_colours) +
  labs(title = "Bacteria: Shannon") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("Composite","Fine")), map_signif_level=TRUE, colour="black", y_position = 6.2)

Alpha_shannon_fungi <- 
  ggplot(table_alpha_diversity_fung, aes(x=Scale, y=Shannon)) +
  geom_boxplot(aes(colour = Scale)) +
  geom_point (aes(colour = Scale), size = 0.5, position = position_jitterdodge(jitter.width = 0.1))  +
  scale_y_continuous(limits = c(0,4.2),breaks = c(0,0.5,1,1.5,2,2.5,3,3.5,4)) +
  theme_classic() +
  scale_color_manual(values = scale_colours) +
  labs(title = "Fungi: Shannon") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("Composite","Fine")), map_signif_level=TRUE, colour="black", y_position = 3.9)

#Figure X1
Fig_X1 <- grid.arrange(Alpha_observed_bacteria + labs(tag = "A"), Alpha_simpson_bacteria + labs(tag = "B"), Alpha_shannon_bacteria + labs(tag = "C"), 
                          Alpha_observed_fungi + labs(tag = "D"), Alpha_simpson_fungi + labs(tag="E"), Alpha_shannon_fungi + labs(tag="F"), ncol=3)
ggsave(filename = paste("Figure_1.svg",sep=""), Fig_X1, width = 9, height = 6)
ggsave(filename = paste("Figure_1.tiff",sep=""), Fig_X1, width = 9, height = 6, dpi = 600)
ggsave(filename = paste("Figure_1.jpg",sep=""), Fig_X1, width = 9, height = 6, dpi = 400)


###Relative abundance###

##Bacteria

#Normalisation
physeq_bact_proportion <- transform_sample_counts(physeq_bact, function (x) x/sum(x))
physeq_fung_proportion <- transform_sample_counts(physeq_fung, function (x) x/sum(x))

#Phylum level (Anything less than 1% is merged together)
Phylum_plots <- list()
Phylum_RA_tables <- list()
for (BF in c("bact","fung")) {
  physeq_phylum <- tax_glom(get(paste("physeq_",BF,"_proportion",sep = "")), taxrank = "Phylum",NArm = FALSE)
  physeq_phylum_dominant <- merge_taxa(physeq_phylum,taxa_names(filter_taxa(physeq_phylum, function(x) mean(x) < 0.01, TRUE)))
  for (ROW in rownames(physeq_phylum_dominant@tax_table[!is.na(physeq_phylum_dominant@tax_table[,"Phylum"]),])) {
    if (physeq_phylum_dominant@tax_table[ROW,"Phylum"] == "") {
      physeq_phylum_dominant@tax_table[ROW,"Phylum"] <- paste("Unidentified_",physeq_phylum_dominant@tax_table[ROW,"Domain"],sep = "")
    }
  }
  phylum_colours <- colorRampPalette(brewer.pal(12, "Paired"))(length(unique(tax_table(physeq_phylum_dominant)[,"Phylum"])))
  names(phylum_colours) <- sort(unique(tax_table(physeq_phylum_dominant)[,"Phylum"]))
  plot <- 
    plot_bar(physeq_phylum_dominant, fill = "Phylum") + 
    geom_bar(stat = "identity", position = "fill", colour = "black") + 
    labs(x="Tree", y = "Abundance") + 
    facet_wrap(~Tree, scales = "free_x", nrow = 1) + 
    scale_fill_manual(values = phylum_colours, labels = c(sort(unique(physeq_phylum_dominant@tax_table[,"Phylum"]),decreasing = F),"Rare_Taxa")) + 
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  Phylum_plots[BF] <- list(plot)
  
  RA_table <- as.data.frame(rowMeans(otu_table(physeq_phylum_dominant)))
  RA_table["Composite"] <- rowMeans(otu_table(physeq_phylum_dominant)[,rownames(physeq_phylum_dominant@sam_data[physeq_phylum_dominant@sam_data[,"Scale"]=="Composite",])])
  RA_table["Fine"] <- rowMeans(otu_table(physeq_phylum_dominant)[,rownames(physeq_phylum_dominant@sam_data[physeq_phylum_dominant@sam_data[,"Scale"]=="Fine",])])
  RA_table["Phylum"] <- tax_table(physeq_phylum_dominant)[,"Phylum"]
  colnames(RA_table)[1] <- "Overall"
  RA_table <- RA_table[,c("Phylum","Overall","Composite","Fine")]
  RA_table <- RA_table[order(RA_table$Phylum),]
  Phylum_RA_tables[BF] <- list(RA_table)
  write.table(RA_table, file = paste("RA_table_",BF,"_phylum.csv",sep = ""),quote = T,row.names = F, col.names = T, sep = ",")
}

#Order level (Anything less than 1% is merged together)
Order_plots <- list()
Order_RA_tables <- list()
for (BF in c("bact","fung")) {
  physeq_order <- tax_glom(get(paste("physeq_",BF,"_proportion",sep = "")), taxrank = "Order",NArm = FALSE)
  physeq_order_dominant <- merge_taxa(physeq_order,taxa_names(filter_taxa(physeq_order, function(x) mean(x) < 0.01, TRUE)))
  for (ROW in rownames(physeq_order_dominant@tax_table[!is.na(physeq_order_dominant@tax_table[,"Order"]),])) {
    if (physeq_order_dominant@tax_table[ROW,"Order"] == "") {
      if (physeq_order_dominant@tax_table[ROW,"Phylum"] == "") {
        physeq_order_dominant@tax_table[ROW,"Order"] <- paste("Unidentified_",physeq_order_dominant@tax_table[ROW,"Domain"],sep = "")
      }
      else
        physeq_order_dominant@tax_table[ROW,"Order"] <- paste("Unidentified_",physeq_order_dominant@tax_table[ROW,"Phylum"],sep = "")
    }
  }
  for (ROW in rownames(physeq_order_dominant@tax_table[!is.na(physeq_order_dominant@tax_table[,"Order"]),])) {
    if (physeq_order_dominant@tax_table[ROW,"Order"] == "uncultured") {
      physeq_order_dominant@tax_table[ROW,"Order"] <- paste("Uncultured_",physeq_order_dominant@tax_table[ROW,"Phylum"],sep = "")
    }
  }
  for (ROW in rownames(physeq_order_dominant@tax_table[!is.na(physeq_order_dominant@tax_table[,"Order"]),])) {
    if (physeq_order_dominant@tax_table[ROW,"Order"] == "WD260") {
      physeq_order_dominant@tax_table[ROW,"Order"] <- "Gammaproteobacteria WD260"
    }
  }
  order_colours <- colorRampPalette(brewer.pal(12, "Paired"))(length(unique(tax_table(physeq_order_dominant)[,"Order"])))
  names(order_colours) <- sort(unique(tax_table(physeq_order_dominant)[,"Order"]))
  plot <- 
    plot_bar(physeq_order_dominant, fill = "Order") + 
    geom_bar(stat = "identity", position = "fill", colour = "black") + 
    labs(x="Tree", y = "Abundance") + 
    facet_wrap(~Tree, scales = "free_x", nrow = 1) + 
    scale_fill_manual(values = order_colours, labels = c(sort(unique(physeq_order_dominant@tax_table[,"Order"]),decreasing = F),"Rare_Taxa")) + 
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  Order_plots[BF] <- list(plot)
  
  RA_table <- as.data.frame(rowMeans(otu_table(physeq_order_dominant)))
  RA_table["Composite"] <- rowMeans(otu_table(physeq_order_dominant)[,rownames(physeq_order_dominant@sam_data[physeq_order_dominant@sam_data[,"Scale"]=="Composite",])])
  RA_table["Fine"] <- rowMeans(otu_table(physeq_order_dominant)[,rownames(physeq_order_dominant@sam_data[physeq_order_dominant@sam_data[,"Scale"]=="Fine",])])
  RA_table["Order"] <- tax_table(physeq_order_dominant)[,"Order"]
  colnames(RA_table)[1] <- "Overall"
  RA_table <- RA_table[,c("Order","Overall","Composite","Fine")]
  RA_table <- RA_table[order(RA_table$Order),]
  Phylum_RA_tables[BF] <- list(RA_table)
  write.table(RA_table, file = paste("RA_table_",BF,"_order.csv",sep = ""),quote = T,row.names = F, col.names = T, sep = ",")
}

#Genus level (Anything less than 1% is merged together)
Genus_plots <- list()
Genus_RA_tables <- list()
for (BF in c("bact","fung")) {
  physeq_genus <- tax_glom(get(paste("physeq_",BF,"_proportion",sep = "")), taxrank = "Genus",NArm = FALSE)
  physeq_genus_dominant <- merge_taxa(physeq_genus,taxa_names(filter_taxa(physeq_genus, function(x) mean(x) < 0.01, TRUE)))
  for (ROW in rownames(physeq_genus_dominant@tax_table[!is.na(physeq_genus_dominant@tax_table[,"Genus"]),])) {
    if (physeq_genus_dominant@tax_table[ROW,"Genus"] == "") {
      if (physeq_genus_dominant@tax_table[ROW,"Phylum"] == "") {
        physeq_genus_dominant@tax_table[ROW,"Genus"] <- paste("Unidentified_",physeq_genus_dominant@tax_table[ROW,"Domain"],sep = "")
      }
      else if (physeq_genus_dominant@tax_table[ROW,"Order"] == "") {
        physeq_genus_dominant@tax_table[ROW,"Genus"] <- paste("Unidentified_",physeq_genus_dominant@tax_table[ROW,"Phylum"],sep = "")
      }
      else
        physeq_genus_dominant@tax_table[ROW,"Genus"] <- paste("Unidentified_",physeq_genus_dominant@tax_table[ROW,"Order"],sep = "")
    }
  }
  for (ROW in rownames(physeq_genus_dominant@tax_table[!is.na(physeq_genus_dominant@tax_table[,"Genus"]),])) {
    if (physeq_genus_dominant@tax_table[ROW,"Genus"] == "uncultured") {
      if (physeq_genus_dominant@tax_table[ROW,"Order"] == "uncultured") {
        physeq_genus_dominant@tax_table[ROW,"Genus"] <- paste("Uncultured_",physeq_genus_dominant@tax_table[ROW,"Phylum"],sep = "")
      }
      else
        physeq_genus_dominant@tax_table[ROW,"Genus"] <- paste("Uncultured_",physeq_genus_dominant@tax_table[ROW,"Order"],sep = "")
    }
    for (ROW in rownames(physeq_genus_dominant@tax_table[!is.na(physeq_genus_dominant@tax_table[,"Genus"]),])) {
      if (physeq_genus_dominant@tax_table[ROW,"Genus"] == "67-14") {
        physeq_genus_dominant@tax_table[ROW,"Genus"] <- "Solirubrobacterales 64-14"
      }
      if (physeq_genus_dominant@tax_table[ROW,"Genus"] == "WD260") {
        physeq_genus_dominant@tax_table[ROW,"Genus"] <- "Gammaproteobacteria WD260"
      }
    }
    
  }
  genus_colours <- colorRampPalette(brewer.pal(12, "Paired"))(length(unique(tax_table(physeq_genus_dominant)[,"Genus"])))
  names(genus_colours) <- sort(unique(tax_table(physeq_genus_dominant)[,"Genus"]))
  plot <- 
    plot_bar(physeq_genus_dominant, fill = "Genus") + 
    geom_bar(stat = "identity", position = "fill", colour = "black") + 
    labs(x="Tree", y = "Abundance") + 
    facet_wrap(~Tree, scales = "free_x", nrow = 1) + 
    scale_fill_manual(values = genus_colours, labels = c(sort(unique(physeq_genus_dominant@tax_table[,"Genus"]),decreasing = F),"Rare_Taxa")) + 
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  Genus_plots[BF] <- list(plot)
  
  RA_table <- as.data.frame(rowMeans(otu_table(physeq_genus_dominant)))
  RA_table["Composite"] <- rowMeans(otu_table(physeq_genus_dominant)[,rownames(physeq_genus_dominant@sam_data[physeq_genus_dominant@sam_data[,"Scale"]=="Composite",])])
  RA_table["Fine"] <- rowMeans(otu_table(physeq_genus_dominant)[,rownames(physeq_genus_dominant@sam_data[physeq_genus_dominant@sam_data[,"Scale"]=="Fine",])])
  RA_table["Genus"] <- tax_table(physeq_genus_dominant)[,"Genus"]
  colnames(RA_table)[1] <- "Overall"
  RA_table <- RA_table[,c("Genus","Overall","Composite","Fine")]
  RA_table <- RA_table[order(RA_table$Genus),]
  Phylum_RA_tables[BF] <- list(RA_table)
  write.table(RA_table, file = paste("RA_table_",BF,"_genus.csv",sep = ""),quote = T,row.names = F, col.names = T, sep = ",")
}

#Figure X3
Fig_X3<- grid.arrange(Genus_plots$bact + labs(tag = "A"), Genus_plots$fung + labs(tag = "B"), ncol=1)
ggsave(filename = paste("Figure_3.svg",sep=""), Fig_X3, width = 16, height = 10)
ggsave(filename = paste("Figure_3.tiff",sep=""), Fig_X3, width = 16, height = 10, dpi = 600)
ggsave(filename = paste("Figure_3.jpg",sep=""), Fig_X3, width = 16, height = 10, dpi = 400)


#Figure SX5
Fig_SX5 <- grid.arrange(Phylum_plots$bact + labs(tag = "A"), Phylum_plots$fung + labs(tag = "B"), Order_plots$bact + labs(tag = "C"), 
                        Order_plots$fung + labs(tag = "D"), ncol=2)
ggsave(filename = paste("Figure_S4_Supplementary_Data.svg",sep=""), Fig_SX5, width = 32, height = 10)
ggsave(filename = paste("Figure_S4_Supplementary_Data.tiff",sep=""), Fig_SX5, width = 32, height = 10, dpi = 600)
ggsave(filename = paste("Figure_S4_Supplementary_Data.jpg",sep=""), Fig_SX5, width = 32, height = 10, dpi = 400)


###Beta diversity/Ordination###

#Normalisation
physeq_bact_proportion <- transform_sample_counts(physeq_bact, function (x) x/sum(x))
physeq_fung_proportion <- transform_sample_counts(physeq_fung, function (x) x/sum(x))

#Distances
JaccardDist_bact <- distance(physeq_bact_proportion,method = "jaccard")
HellingerDist_bact <- distance(otu_table(decostand(physeq_bact_proportion@otu_table,"hellinger",MARGIN = 2),taxa_are_rows = TRUE),method = "euclidean")
JaccardDist_fung <- distance(physeq_fung_proportion,method = "jaccard")
HellingerDist_fung <- distance(otu_table(decostand(physeq_fung_proportion@otu_table,"hellinger",MARGIN = 2),taxa_are_rows = TRUE),method = "euclidean")

#Calculate the ordinations
set.seed(123456)
ordination_Jaccard_bact <- ordinate(physeq_bact_proportion, method="NMDS", distance=JaccardDist_bact)
set.seed(123456)
ordination_Hellinger_bact <- ordinate(physeq_bact, method="NMDS", distance=HellingerDist_bact)
set.seed(123456)
ordination_Jaccard_fung <- ordinate(physeq_fung_proportion, method="NMDS", distance=JaccardDist_fung)
set.seed(123456)
ordination_Hellinger_fung <- ordinate(physeq_fung, method="NMDS", distance=HellingerDist_fung)

#Plotting
NMDS_phyloseq_Jaccard_bact <- 
  plot_ordination(physeq_bact_proportion, ordination_Jaccard_bact, color = "Tree", shape = "Scale") +
  geom_point(size=3) + 
  scale_colour_manual(values = c(tree_colours)) + 
  theme_classic() +
  labs(title = "Bacteria: NMDS (Jaccard)") +
  annotate(geom="label", x = max(ordination_Jaccard_bact$points[,"MDS1"]), y = max(ordination_Jaccard_bact$points[,"MDS2"])+0.1, label = paste("Stress: ", round(ordination_Jaccard_bact$stress, digits = 3)), vjust = "inward", hjust = "inward") +
  theme(plot.title = element_text(hjust = 0.5))

NMDS_phyloseq_Hellinger_bact <- 
  plot_ordination(physeq_bact, ordination_Hellinger_bact, color = "Tree", shape = "Scale") + 
  geom_point(size=3) + 
  theme(aspect.ratio=1) + 
  scale_colour_manual(values = tree_colours) + 
  theme_classic() +
  labs(title = "Bacteria: NMDS (Hellinger)") +
  annotate(geom="label", x = max(ordination_Hellinger_bact$points[,"MDS1"]), y = max(ordination_Hellinger_bact$points[,"MDS2"])+0.1, label = paste("Stress: ", round(ordination_Hellinger_bact$stress, digits = 3)), vjust = "inward", hjust = "inward") +
  theme(plot.title = element_text(hjust = 0.7))

NMDS_phyloseq_Jaccard_fung <- 
  plot_ordination(physeq_fung_proportion, ordination_Jaccard_fung, color = "Tree", shape = "Scale") +
  geom_point(size=3) + 
  scale_colour_manual(values = c(tree_colours)) + 
  theme_classic() +
  labs(title = "Fungi: NMDS (Jaccard)") +
  annotate(geom="label", x = max(ordination_Jaccard_fung$points[,"MDS1"]), y = max(ordination_Jaccard_fung$points[,"MDS2"])+0.1, label = paste("Stress: ", round(ordination_Jaccard_fung$stress, digits = 3)), vjust = "inward", hjust = "inward") +
  theme(plot.title = element_text(hjust = 0.5))

NMDS_phyloseq_Hellinger_fung <- 
  plot_ordination(physeq_fung, ordination_Hellinger_fung, color = "Tree", shape = "Scale") + 
  geom_point(size=3) + 
  theme(aspect.ratio=1) + 
  scale_colour_manual(values = tree_colours) + 
  theme_classic() +
  labs(title = "Fungi: NMDS (Hellinger)") +
  annotate(geom="label", x = max(ordination_Hellinger_fung$points[,"MDS1"]), y = max(ordination_Hellinger_fung$points[,"MDS2"])+0.1, label = paste("Stress: ", round(ordination_Hellinger_fung$stress, digits = 3)), vjust = "inward", hjust = "inward") +
  theme(plot.title = element_text(hjust = 0.5))

#Figure X2
Figure_X2 <- grid.arrange(NMDS_phyloseq_Jaccard_bact + labs(tag = "A"), NMDS_phyloseq_Jaccard_fung + labs(tag = "B"), ncol=1)
ggsave(filename = paste("Figure_2.svg",sep=""), Figure_X2, width = 8, height = 8)
ggsave(filename = paste("Figure_2.tiff",sep=""), Figure_X2, width = 8, height = 8, dpi = 600)
ggsave(filename = paste("Figure_2.jpg",sep=""), Figure_X2, width = 8, height = 8, dpi = 400)


#Figure SX2
Figure_SX2 <- grid.arrange(NMDS_phyloseq_Hellinger_bact + labs(tag = "A"), NMDS_phyloseq_Hellinger_fung + labs(tag = "B"), ncol=1)
ggsave(filename = paste("Figure_S3.svg",sep=""), Figure_SX2, width = 8, height = 8)
ggsave(filename = paste("Figure_S3.tiff",sep=""), Figure_SX2, width = 8, height = 8, dpi = 600)
ggsave(filename = paste("Figure_S3.jpg",sep=""), Figure_SX2, width = 8, height = 8, dpi = 400)


###Cluster the samples###

#Normalisation
physeq_bact_proportion <- transform_sample_counts(physeq_bact, function (x) x/sum(x))
physeq_fung_proportion <- transform_sample_counts(physeq_fung, function (x) x/sum(x))

#Distance
JaccardDist_bact <- distance(physeq_bact_proportion,method = "jaccard")
JaccardDist_fung <- distance(physeq_fung_proportion,method = "jaccard")

#Find the best Linkage Method to use for clustering 
#define linkage methods
m <- c( "average", "single", "complete", "ward")
names(m) <- c( "average", "single", "complete", "ward")
#function to compute agglomerative coefficient
ac_bact <- function(x) {
  agnes(JaccardDist_bact, diss = TRUE, method = x)$ac
}
ac_fung <- function(x) {
  agnes(JaccardDist_fung, diss = TRUE, method = x)$ac
}
#calculate agglomerative coefficient for each clustering linkage method
sapply(m, ac_bact)
sapply(m, ac_fung)
#Ward has highest ac in both cases, so best method to use

#Cluster samples
final_clust_bact <- hclust(JaccardDist_bact, method = "ward.D2")
final_clust_fung <- hclust(JaccardDist_fung, method = "ward.D2")

#Determine the Optimal Number of Clusters by choosing the one returned most often by a list of 20 indices
inds <- c("kl","ch","hartigan","cindex","db","silhouette","duda","pseudot2","ratkowsky","ball","ptbiserial","gap","frey","mcclain","gamma","gplus","tau","dunn","sdindex","sdbw")
#In previous agnes function, the "ward" option is equivalent to "ward.D2" in the hclust function. So will use the "ward.D2" method in NbClust as well for equivalence.
#https://link.springer.com/article/10.1007/s00357-014-9161-z
restable_bact <- data.frame()
for (i in (1:length(inds))) {
  restable1 <- data.frame()
  restable1[1,1] <- inds[i]
  restable1[1,2] <- NbClust(as.data.frame(t(otu_table(physeq_bact_proportion@otu_table,taxa_are_rows = TRUE))),diss=JaccardDist_bact,distance = NULL,method="ward.D2",index=inds[i])$Best.nc[1]
  restable_bact <- rbind(restable_bact,restable1)
}
restable_fung <- data.frame()
for (i in (1:length(inds))) {
  restable1 <- data.frame()
  restable1[1,1] <- inds[i]
  restable1[1,2] <- NbClust(as.data.frame(t(otu_table(physeq_fung_proportion@otu_table,taxa_are_rows = TRUE))),diss=JaccardDist_fung,distance = NULL,method="ward.D2",index=inds[i])$Best.nc[1]
  restable_fung <- rbind(restable_fung,restable1)
}

#plot histogram
Clustering_inds_bact <- ggplot(restable_bact, aes(x=V2)) +
  geom_histogram(bins = max(restable_bact$V2),col="black",fill="dark grey",binwidth = 1) +
  scale_x_continuous(n.breaks = max(restable_bact$V2)) +
  scale_y_continuous(breaks = c(0,2,4,6,8,10)) +
  theme_classic() +
  labs(title = "Optimal number of clusters", y="Count", x = "Clusters") +
  theme(plot.title = element_text(hjust = 0.5))
Clustering_inds_fung <- ggplot(restable_fung, aes(x=V2)) +
  geom_histogram(bins = max(restable_fung$V2),col="black",fill="dark grey",binwidth = 1) +
  scale_x_continuous(n.breaks = max(restable_fung$V2)) +
  scale_y_continuous(breaks = c(0,2,4,6,8,10)) +
  theme_classic() +
  labs(title = "Optimal number of clusters", y="Count", x = "Clusters") +
  theme(plot.title = element_text(hjust = 0.5))

#Colour clusters by vector is not officially supported but seems to work for now.
bact_clust_colours <- tree_colours[str_extract(final_clust_bact$labels,".{10}")][final_clust_bact$order]
fung_clust_colours <- tree_colours[str_extract(final_clust_fung$labels,".{10}")][final_clust_fung$order]

#Perform bootstrapping of clusters
set.seed(123456)
Bootstrap_clust_bact <- Bclust(t(as.data.frame(otu_table(physeq_bact_proportion@otu_table, taxa_are_rows = TRUE))), FUN = function(.x) hclust(distance(otu_table(t(.x),taxa_are_rows = TRUE), method = "jaccard"), method = "ward.D2"), iter=1000,mc.cores = 7)
set.seed(123456)
Bootstrap_clust_fung <- Bclust(t(as.data.frame(otu_table(physeq_fung_proportion@otu_table, taxa_are_rows = TRUE))), FUN = function(.x) hclust(distance(otu_table(t(.x),taxa_are_rows = TRUE), method = "jaccard"), method = "ward.D2"), iter=1000,mc.cores = 7)

#Get labels and co-ordinates 
#Not how this is meant to be used but it seems to work
plot(Bootstrap_clust_bact)
Bootstrap_labels_bact <- Bclabels(hcl = final_clust_bact,values = Bootstrap_clust_bact$values)
plot(Bootstrap_clust_fung)
Bootstrap_labels_fung <- Bclabels(hcl = final_clust_fung,values = Bootstrap_clust_fung$values)

#Final set up
dendro_clust_bact <- dendro_data(final_clust_bact)
dendro_clust_fung <- dendro_data(final_clust_fung)
dendro_clust_bact_labels <- as.data.frame(cbind(Bootstrap_labels_bact$coords[,"x"],cbind(Bootstrap_labels_bact$coords[,"y"],Bootstrap_labels_bact$labels)))
colnames(dendro_clust_bact_labels) <- c("x","y","label")
dendro_clust_fung_labels <- as.data.frame(cbind(Bootstrap_labels_fung$coords[,"x"],cbind(Bootstrap_labels_fung$coords[,"y"],Bootstrap_labels_fung$labels)))
colnames(dendro_clust_fung_labels) <- c("x","y","label")

#Plot clustering
#Dendrocut needs to be tried manually
Cluster_graph_bact <- 
  ggdendrogram(dendro_clust_bact) + 
  geom_dendro(final_clust_bact, dendrocut=26, groupCols=c(colorRampPalette(brewer.pal(8, "Dark2"))(14),"black")) +
  geom_text(data = dendro_clust_bact_labels, aes(x = x, y = y-0.05, label = paste(round(label*100),"%",sep="")), size = 2, vjust = 0) +
  labs(y = "", x = "") +
  scale_y_continuous(breaks = NULL) +
  coord_cartesian(xlim = c(0, 41)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, colour = bact_clust_colours))
Cluster_graph_fung <- 
  ggdendrogram(dendro_clust_fung) + 
  geom_dendro(final_clust_fung, dendrocut=24, groupCols=c(colorRampPalette(brewer.pal(8, "Dark2"))(14),"black")) +
  geom_text(data = dendro_clust_fung_labels, aes(x = x, y = y-0.05, label = paste(round(label*100),"%",sep="")), size = 2, vjust = 0) +
  labs(y = "", x = "") +
  scale_y_continuous(breaks = NULL) +
  coord_cartesian(xlim = c(0, 41)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, colour = fung_clust_colours))

#Figure SX3
Figure_SX3 <- grid.arrange(Cluster_graph_bact + labs(tag = "A"), Clustering_inds_bact + labs(tag = "B"), Cluster_graph_fung + labs(tag = "C"), Clustering_inds_fung + labs(tag = "D"), ncol=2, widths = c(2,1))
ggsave(filename = paste("Figure_S2.svg",sep=""), Figure_SX3, width = 10, height = 8)
ggsave(filename = paste("Figure_S2.tiff",sep=""), Figure_SX3, width = 10, height = 8, dpi = 600)
ggsave(filename = paste("Figure_S2.jpg",sep=""), Figure_SX3, width = 10, height = 8, dpi = 400)


###Compare Bulk and Fine association networks
#Doing comparison at genus level through CCREPE
#Do bacteria and fungi separately for maximum power

#Merge to an appropriate level
physeq_bact_Genus <- tax_glom(physeq_bact, taxrank = "Genus",NArm = FALSE)
physeq_fung_Genus <- tax_glom(physeq_fung, taxrank = "Genus",NArm = FALSE)

#Remove taxa with fewer than 60 reads
physeq_bact_Genus_filtered <- subset_taxa(physeq_bact_Genus, rowSums(otu_table(physeq_bact_Genus)) >= 60)
physeq_fung_Genus_filtered <- subset_taxa(physeq_fung_Genus, rowSums(otu_table(physeq_fung_Genus)) >= 60)

# #Make sure that both groups contain the same taxa
physeq_bact_Genus_filtered_common <- subset_taxa(physeq_bact_Genus_filtered, rowSums(otu_table(physeq_bact_Genus_filtered)[,sample_data(physeq_bact_Genus_filtered)[,"Scale"] == "Fine"]) > 1 & rowSums(otu_table(physeq_bact_Genus_filtered)[,sample_data(physeq_bact_Genus_filtered)[,"Scale"] == "Composite"]) > 1)
physeq_fung_Genus_filtered_common <- subset_taxa(physeq_fung_Genus_filtered, rowSums(otu_table(physeq_fung_Genus_filtered)[,sample_data(physeq_fung_Genus_filtered)[,"Scale"] == "Fine"]) > 1 & rowSums(otu_table(physeq_fung_Genus_filtered)[,sample_data(physeq_fung_Genus_filtered)[,"Scale"] == "Composite"]) > 1)

#Generate the networks
network_bacteria <- netConstruct(physeq_bact_Genus_filtered_common, group = sample_data(physeq_bact_Genus_filtered_common)$Scale%in%"Fine", measurePar = list(min.subj = 10), filtTax = "none", measure = "ccrepe", normMethod = "fractions", sparsMethod = "softThreshold", seed = 123456)
network_fungi <- netConstruct(physeq_fung_Genus_filtered_common, group = sample_data(physeq_fung_Genus_filtered_common)$Scale%in%"Fine",  measurePar = list(min.subj = 10), filtTax = "none", measure = "ccrepe", normMethod = "fractions", sparsMethod = "softThreshold", seed = 123456)

#Analyse partial network and plot
netprops_bacteria <- netAnalyze(network_bacteria, clustMethod = "cluster_fast_greedy")
netprops_fungi <- netAnalyze(network_fungi, clustMethod = "cluster_fast_greedy")

#Make sure group nodes are colourblind safe. Edges are assigned in the functions.
netcolours_Bacteria <- brewer.pal(max(length(unique(netprops_bacteria$clustering$clust1)),length(unique(netprops_bacteria$clustering$clust2))), "Set3")
netcolours_Fungi <- brewer.pal(max(length(unique(netprops_fungi$clustering$clust1)),length(unique(netprops_fungi$clustering$clust2))), "Set3")

#Create labels

labels_bacteria <- as.vector(tax_table(physeq_bact_Genus_filtered_common)[,"Genus"])
names(labels_bacteria) <- (rownames(tax_table(physeq_bact_Genus_filtered_common)))

labels_fungi <- as.vector(tax_table(physeq_fung_Genus_filtered_common)[,"Genus"])
names(labels_fungi) <- (rownames(tax_table(physeq_fung_Genus_filtered_common)))

#Plot networks
jpeg(filename = "Network_Bacteria.jpg", width = 2880, height = 1610)
plot(netprops_bacteria,
     shortenLabels = "none",
     layoutGroup = "union",
     sameLayout = TRUE,
     groupNames = c("",""),
     labels = labels_bacteria,
     labelScale = FALSE,
     rmSingles = "inboth",
     nodeTransp = 40,
     borderCol = "gray60",
     hubTransp = 20,
     hubBorderWidth = 2,
     hubBorderCol = "gray90",
     edgeTranspLow = 80,
     edgeTranspHigh = 30,
     cexLabels = 0,
     cexTitle = 5,
     colorVec = netcolours_Bacteria,
     posCol = "#2c7bb6", 
     negCol = "#d7191c")
dev.off()

jpeg(filename = "Network_Fungi.jpg", width = 2880, height = 1610)
plot(netprops_fungi,
     shortenLabels = "none",
     layoutGroup = "union",
     sameLayout = TRUE,
     groupNames = c("",""),
     labels = labels_fungi,
     labelScale = FALSE,
     rmSingles = "inboth",
     nodeTransp = 40,
     borderCol = "gray60",
     hubTransp = 20,
     hubBorderWidth = 2,
     hubBorderCol = "gray90",
     edgeTranspLow = 80,
     edgeTranspHigh = 30,
     cexLabels = 0,
     cexTitle = 5,
     colorVec = netcolours_Fungi,
     posCol = "#2c7bb6", 
     negCol = "#d7191c")
dev.off()

#Compare the networks
#Can't do permutation tests because sparsity means that too many ASVs are dropped and the dimensions do not match.
Network_comparison_bacteria <- netCompare(netprops_bacteria,
                                          cores = 7L,
                                          seed = 12345)
Network_comparison_fungi <- netCompare(netprops_fungi,
                                       cores = 7L,
                                       seed = 12345)

#Extract associations (Only keep those with high edge value > 0.5)
Net_Associations_Bact_Fine <- data.frame()
for (x in rownames(network_bacteria$assoMat1)) {
  edge_weight <- as.vector(network_bacteria$assoMat1[x,])
  int_partner <- colnames(network_bacteria$assoMat1)
  int_res <- data.frame("ASV_1" = rep(x,nrow(network_bacteria$assoMat1)),"ASV_2" = int_partner,"Weight" = edge_weight)
  int_res <- int_res[int_res$Weight > 0.5,]
  Net_Associations_Bact_Fine <- rbind(Net_Associations_Bact_Fine,int_res)
}
Net_Associations_Bact_Bulk <- data.frame()
for (x in rownames(network_bacteria$assoMat2)) {
  edge_weight <- as.vector(network_bacteria$assoMat2[x,])
  int_partner <- colnames(network_bacteria$assoMat2)
  int_res <- data.frame("ASV_1" = rep(x,nrow(network_bacteria$assoMat2)),"ASV_2" = int_partner,"Weight" = edge_weight)
  int_res <- int_res[int_res$Weight > 0.5,]
  Net_Associations_Bact_Bulk <- rbind(Net_Associations_Bact_Bulk,int_res)
}
Net_Associations_Fung_Fine <- data.frame()
for (x in rownames(network_fungi$assoMat1)) {
  edge_weight <- as.vector(network_fungi$assoMat1[x,])
  int_partner <- colnames(network_fungi$assoMat1)
  int_res <- data.frame("ASV_1" = rep(x,nrow(network_fungi$assoMat1)),"ASV_2" = int_partner,"Weight" = edge_weight)
  int_res <- int_res[int_res$Weight > 0.5,]
  Net_Associations_Fung_Fine <- rbind(Net_Associations_Fung_Fine,int_res)
}
Net_Associations_Fung_Bulk <- data.frame()
for (x in rownames(network_fungi$assoMat2)) {
  edge_weight <- as.vector(network_fungi$assoMat2[x,])
  int_partner <- colnames(network_fungi$assoMat2)
  int_res <- data.frame("ASV_1" = rep(x,nrow(network_fungi$assoMat2)),"ASV_2" = int_partner,"Weight" = edge_weight)
  int_res <- int_res[int_res$Weight > 0.5,]
  Net_Associations_Fung_Bulk <- rbind(Net_Associations_Fung_Bulk,int_res)
}

#Venn of interactions new or shared
Net_Association_Bact_Overlap <- ggvenn(list(Fine = paste(Net_Associations_Bact_Fine$ASV_1,Net_Associations_Bact_Fine$ASV_2,sep = "_"), 
                                            Bulk = paste(Net_Associations_Bact_Bulk$ASV_1,Net_Associations_Bact_Bulk$ASV_2,sep = "_")),
                                       fill_alpha = 0, text_size = 3, stroke_size = 0.5, set_name_size = 5)
Net_Association_Fung_Overlap <- ggvenn(list(Fine = paste(Net_Associations_Fung_Fine$ASV_1,Net_Associations_Fung_Fine$ASV_2,sep = "_"), 
                                            Bulk = paste(Net_Associations_Fung_Bulk$ASV_1,Net_Associations_Fung_Bulk$ASV_2,sep = "_")),
                                       fill_alpha = 0, text_size = 3, stroke_size = 0.5, set_name_size = 5)

#Make a figure
Net_Association_Overlap_Fig <- grid.arrange(Net_Association_Bact_Overlap + labs(tag = "A"), Net_Association_Fung_Overlap + labs(tag = "B"), nrow=1)
ggsave(filename = "Net_Association_Overlap_Fig.svg", Net_Association_Overlap_Fig, width = 8, height = 4)

#Manually construct network comparison table

##Information for text

#Number of fungal ASVs
summary(colSums(physeq_fung_Genus@otu_table>1))

#Genus variation
physeq_bact_Genus_RA <- transform_sample_counts(physeq_bact_Genus, function (x) sqrt(x/sum(x)))
physeq_bact_Genus_RA@tax_table[physeq_bact_Genus_RA@tax_table[,"Genus"] == "Pseudomonas",]
physeq_bact_Genus_RA@otu_table["08f55b866c61fd6c3104d94c26230721",c("sample_049A","sample_049B","sample_049C")]
physeq_fung_Genus_RA <- transform_sample_counts(physeq_fung_Genus, function (x) sqrt(x/sum(x)))
physeq_fung_Genus_RA@tax_table[physeq_fung_Genus_RA@tax_table[,"Genus"]=="Kretzschmaria",]
physeq_fung_Genus_RA@otu_table["d56ca28652a1b8b0c1495a60fa09ef5d",c("sample_007A","sample_007B","sample_007C")]