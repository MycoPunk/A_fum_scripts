#A fum pan genome analysis, first pass data investigation of PIRATE results
#4.Oct.2020

#set packages 
library(data.table)
library(tidyverse)
library(hrbrthemes)
library(forcats)
library(ggforce)
library(ggimage)
library(ggtree)
library(ape)
library(phytools)



setwd("~/Desktop/Project_Afum_pangenome/")
PIRATE.gene_families<-as.data.frame(fread("PIRATE.gene_families.tsv"))

##basic stats:
#how many gene fmailies are there?
n_gene_fams<- nrow(PIRATE.gene_families)
#19016

#how many strains?
#names(PIRATE.gene_families)
#names of cols to exclude
cols_to_exclude<- c("allele_name",
                    "gene_family",
                    "consensus_gene_name",
                    "consensus_product",
                    "threshold",
                    "alleles_at_maximum_threshold",
                    "number_genomes",
                    "average_dose",
                    "min_dose",
                    "max_dose",
                    "genomes_containing_fissions",
                    "genomes_containing_duplications",
                    "number_fission_loci",
                    "number_duplicated_loci",
                    "no_loci",
                    "products",
                    "gene_names",
                    "min_length(bp)",
                    "max_length(bp)",
                    "average_length(bp)")

strains_only<- PIRATE.gene_families[,!names(PIRATE.gene_families) %in% cols_to_exclude]
strain_names<- names(strains_only)

length(strain_names)
length(unique(strain_names))
#print to cross ref for tree building
#write.table(strain_names, "strain_names.txt", sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

#how many of the gene families are in every genome (of 308)
n_gene_fams_core<- sum(PIRATE.gene_families$number_genomes == 308)
#1492
#that's what percent out of the total?
(n_gene_fams_core*100)/n_gene_fams
#7.846024 %

#how many of the gene families are singletons (accessory)
n_gene_fams_singletons<- sum(PIRATE.gene_families$number_genomes == 1)
#6165 - core genome is smaller than the accessory genome 
#that's what percent out of the total?
(n_gene_fams_singletons*100)/n_gene_fams
#32.42007 %

#graph the distribution of gene presence in a gene family (distribution of core to accessory genes)
#plot
gene_fam_totals<-as.data.frame(PIRATE.gene_families$number_genomes)
colnames(gene_fam_totals) <- 'count'

#set groups
gene_fam_totals$group = 0                        
for (i in 1:nrow(gene_fam_totals)){
  if (gene_fam_totals$count[i] == 1) {
   gene_fam_totals$group[i] = "Singleton"
   } else if (gene_fam_totals$count[i] == 308) {
   gene_fam_totals$group[i] = "Core"
   } else {
   gene_fam_totals$group[i] = "Accessory"
   }
}

#plot
#pallet: Zissou1 = c("#3B9AB2", "#78B7C5", "#EBCC2A", "#E1AF00", "#F21A00")
#throughout: "#78B7C5" = accessory, #F21A00 = singelton, #E1AF00 = core
p <- gene_fam_totals %>%
ggplot(aes(x = count)) +
  geom_bar(aes(color = group), fill = "white",
                 position = "identity") +
  ggtitle("Distribution of gene families across pangenome") +
  ylab("n gene families") + xlab("n genomes") +
  theme_ipsum() +
  theme(
    plot.title = element_text(size=15), legend.title = element_blank()) +
  scale_color_manual(values = c("#78B7C5","#E1AF00", "#F21A00")) 
p


###plot as donut chart 
#make df of totals
fam_dist_df <- data.frame(
  category=c("Singleton", "Accessory", "Core"),
  count=c(n_gene_fams_singletons, n_gene_fams -(n_gene_fams_singletons + n_gene_fams_core), n_gene_fams_core)
)
# Compute percentages
fam_dist_df$fraction <- fam_dist_df$count / sum(fam_dist_df$count)
# Compute the cumulative percentages (top of each rectangle)
fam_dist_df$ymax <- cumsum(fam_dist_df$fraction)
# Compute the bottom of each rectangle
fam_dist_df$ymin <- c(0, head(fam_dist_df$ymax, n=-1))
# Compute label position
fam_dist_df$labelPosition <- (fam_dist_df$ymax + fam_dist_df$ymin) / 2
# Compute a good label
fam_dist_df$label <- paste0(fam_dist_df$category, "\n", fam_dist_df$count)
#plot
ggplot(fam_dist_df, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
  geom_rect() +
  geom_text( x=1.5, aes(y=labelPosition, label=label, color=category), size=4) + # x here controls label position (inner / outer)
  scale_fill_manual(values=c("#78B7C5","#E1AF00", "#F21A00"))+
  scale_color_manual(values=c("#78B7C5","#E1AF00", "#F21A00"))+
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_void() +
  theme(legend.position = "none")



##which strain has the highest/lowest number of singletons and accessory gene fams? 
gene_fam_by_strain<-as.data.frame(PIRATE.gene_families[,21:ncol(PIRATE.gene_families)])
##make binary (if gene = 1, if not = 0)
#fill in zeros
gene_fam_by_strain_zeros<- sapply(gene_fam_by_strain, gsub, pattern = "^\\s*$" , replacement = 0 )
#fill in ones
gene_fam_by_strain_ones<- as.data.frame(replace(gene_fam_by_strain_zeros, gene_fam_by_strain_zeros!="0", 1))
#change to numeric
gene_fam_by_strain_ones_num <- mutate_all(gene_fam_by_strain_ones, function(x) as.numeric(as.character(x)))


#subset to remove core genes form accessory and singletons 
all_accessory<- gene_fam_by_strain_ones_num[rowSums(gene_fam_by_strain_ones_num) != 308,]
#subset to get only singletons 
singletons_only<- gene_fam_by_strain_ones_num[rowSums(gene_fam_by_strain_ones_num) == 1,]

#Which strain has the largest accessory genome (genes not in the core?) 
accessory_by_strain<- as.data.frame(colSums(all_accessory))
colnames(accessory_by_strain) <- "totals"
accessory_by_strain$strain<- row.names(accessory_by_strain)
#get max
accessory_by_strain[which.max(accessory_by_strain$totals),]
#get min
accessory_by_strain[which.min(accessory_by_strain$totals),]
#sort 
accessory_by_strain<- accessory_by_strain[order(accessory_by_strain$totals),]


#graph accessory genome size by strain
p<- accessory_by_strain %>%
  mutate(name = fct_reorder(strain, totals)) %>%
  ggplot( aes(x=name, y=totals)) +
  geom_bar(stat="identity", fill="#78B7C5", alpha=2, width=1) +
  xlab("") + ylab("n acessory gene families") +
  ggtitle("accessory genome size by strain") +
  theme(text=element_text(size=9), 
        axis.text.x = element_text(size = 2, angle=90, hjust=1), legend.position = "none")+
  facet_zoom(ylim = c(6000, max(accessory_by_strain$totals)), zoom.data = ifelse(a <= 6000,  FALSE))

p


#for singletons 
#Which strain has the largest accessory genome (genes not in the core?) 
singletons_only_by_strain<- as.data.frame(colSums(singletons_only))
colnames(singletons_only_by_strain) <- "totals"
singletons_only_by_strain$strain<- row.names(singletons_only_by_strain)
#get max
singletons_only_by_strain[which.max(singletons_only_by_strain$totals),]
#get min
singletons_only_by_strain[which.min(singletons_only_by_strain$totals),]
#how many are zero?
no_singeltons<-data.frame(singletons_only_by_strain[singletons_only_by_strain$totals == 0,])
nrow(no_singeltons)

#sort 
singletons_only_by_strain<- singletons_only_by_strain[order(singletons_only_by_strain$totals),]


#graph accessory genome size by strain
p<- singletons_only_by_strain %>%
  mutate(name = fct_reorder(strain, totals)) %>%
  ggplot( aes(x=name, y=totals)) +
  geom_bar(stat="identity", fill="#F21A00", alpha=2, width=1) +
  xlab("") + ylab("n singleton gene families") +
  ggtitle("singleton genome size by strain") +
  theme(text=element_text(size=9), 
        axis.text.x = element_text(size = 2, angle=90, hjust=1), legend.position = "none")
p


###graph presence / absence matrix on to big tree
tree <- read.tree("pan_genome.SNP.fasttree.tre")

#circ <- ggtree(tree, layout = "circular")
#circ
#tree$tip.label
#length(strain_names)

##trim tree to only the strain used in the pan genome analysis 
#in this case - we need to remove the 20 AF100 strains that are redundant 

redundant_strains<- c("Afum_DMC_AF1001_15",
                      "Afum_DMC_AF100_12_42",
                      "Afum_DMC_AF100_12_9",
                      "Afum_DMC_AF100_1_11",
                      "Afum_DMC_AF100_1_14",
                      "Afum_DMC_AF100_1_18",
                      "Afum_DMC_AF100_1_20_C",
                      "Afum_DMC_AF100_1_20_OE",
                      "Afum_DMC_AF100_1_24",
                      "Afum_DMC_AF100_1_3",
                      "Afum_DMC_AF100_1_4",
                      "Afum_DMC_AF100_1_8",
                      "Afum_DMC_AF100_2B",
                      "Afum_DMC_AF100_3B",
                      "Afum_DMC_AF100_4B",
                      "Afum_DMC_AF100_5B",
                      "Afum_DMC_AF100_6B",
                      "Afum_DMC_AF100_7B",
                      "Afum_DMC_AF100_8B",
                      "Afum_DMC_AF100_9B")

strain_names_wo_redundant<-strain_names[!strain_names %in% redundant_strains]
#subset input data to remove these strains
accessory_by_strain_matched<- accessory_by_strain[rownames(accessory_by_strain) %in% strain_names_wo_redundant,]
singletons_only_by_strain_matched<- singletons_only_by_strain[rownames(singletons_only_by_strain) %in% strain_names_wo_redundant,]

##use mapping file to rename the strains to match the way they appear in the tree

#mapping file to map old names to new names + clade annotations
name_map<-read.delim("namemappingfile_12Oct2020.tab", header = TRUE, sep = "\t", fill = TRUE, strip.white = TRUE)

#accessory_by_strain_matched_w_names<- accessory_by_strain_matched
row.names(accessory_by_strain_matched) <- name_map$name_pop_genome[match(row.names(accessory_by_strain_matched), name_map$name_Pan_genome)]
row.names(singletons_only_by_strain_matched) <- name_map$name_pop_genome[match(row.names(singletons_only_by_strain_matched), name_map$name_Pan_genome)]

#attach clade annotations - note- update this later when you designate clades
singletons_only_by_strain_matched$clade<- name_map$clade[match(row.names(singletons_only_by_strain_matched), name_map$name_pop_genome)]
accessory_by_strain_matched$clade<- name_map$clade[match(row.names(accessory_by_strain_matched), name_map$name_pop_genome)]


##plot tree
#match tree (for the couple sp. that are  miss-named)
tree_me$tip.label[tree_me$tip.label=="F18149"] <- "F18149-JCVI"
tree_me$tip.label[tree_me$tip.label=="Afu_343-P/11"] <- "Afu_343-P-11"
tree_me$tip.label[tree_me$tip.label=="AFIS1435CDC_6"] <- "AFIS1435_CDC_6"

#split by clade
grA_me<- split(row.names(accessory_by_strain_matched), accessory_by_strain_matched$clade)

#set colors by group
tree_grA_me <- ggtree::groupOTU(tree_me, grA_me)
str(tree_grA_me)
levels(attributes(tree_grA_me)$group) 
levels(attributes(tree_grA_me)$group)[1] <- "clade1"
# Reorder factor levels if needed (this controls the order in the legend)
attributes(tree_grA_me)$group <- factor(x = attributes(tree_grA_me)$group, 
                                        levels = c("clade1", "clade2", "clade3", "clade4", "clade5", "clade6", "clade7", "clade8"))

#set colots
my_cols_me <- c(clade1 = "#382147",
                clade2 ="#7D7A70",
                clade3 = "#ABA778",
                clade4 = "#F2CC35",
                clade5 = "#F2CC85",
                #clade7 = "#FFCA98",
                clade6 = "#F7A583",
                clade7 = "#ED7F6F",
                clade8 = "#D4494E")
 


#"#A6669A",
#"#6B3F8A" 

names(my_cols_me) <- levels(attributes(tree_grA_me)$group)
scales::show_col(my_cols_me); my_cols_me

#optional ultrametric tree 
#tree_grA_me_ultra<- force.ultrametric(tree_grA_me)

#simple plot
tree_plot_me <- 
  ggtree(tr = tree_grA_me, 
         # color by group attribute, check str(tree_grA_me)
         mapping = aes(color = group), 
         layout  = 'circular') + 
  geom_treescale(x=3, y=NULL, color = "white") +
  # set line thickness
#  size = 1)
  # adjust coloring of main groups
  scale_color_manual(name = 'Clade', values = my_cols_me) + 
  theme(legend.title=element_text(size=9), # The title of legend 
      legend.text=element_text(size=7)) +
  guides(color = guide_legend(override.aes = list(size = 4))) 


# plot and ddd the tip labels
tree_plot_me + geom_tiplab(size = 1, align = TRUE, linesize = .25, linetype = 0)

#shrink the the min val of each data set
accessory_by_strain_matched$totals<- accessory_by_strain_matched$totals - min(accessory_by_strain_matched$totals)
singletons_only_by_strain_matched$totals<- singletons_only_by_strain_matched$totals - min(singletons_only_by_strain_matched$totals)

##to visualize, re-scale both data sets on a scale of 0-1
#function
scale <- function(x){(x-min(x))/(max(x)-min(x))}
accessory_by_strain_matched$totals<- scale(accessory_by_strain_matched$totals)
singletons_only_by_strain_matched$totals<- scale(singletons_only_by_strain_matched$totals)

#join annotation data (accessory)
accessory_by_strain_matched_input<- data.table(cbind((tip_lbs = row.names(accessory_by_strain_matched)), clade = (accessory_by_strain_matched$clade), val = (accessory_by_strain_matched$totals)))
tree_plot_me <- tree_plot_me %<+% data.table(accessory_by_strain_matched_input) 

#join annotation data (singleton)
singletons_only_by_strain_matched_input<- data.table(cbind((tip_lbs = row.names(singletons_only_by_strain_matched)), clade = (singletons_only_by_strain_matched$clade), val2 = (singletons_only_by_strain_matched$totals)))
tree_plot_me <- tree_plot_me %<+% data.table(singletons_only_by_strain_matched_input) 

#process updated tree view data
tree_dt_me <- data.table(tree_plot_me$data)
head(tree_dt_me)

# select only the tip labels and order by coord y
tree_dt_me <- tree_dt_me[isTip == TRUE][order(y)]

##add circular barplot
# Define variable to control the x coordinates of bars (segments)
#for accessory
my_factor <- 1
x_base_me <- max(tree_dt_me$x) + abs(min(as.numeric(tree_dt_me$val), na.rm = TRUE))*my_factor + 1.5
# Define variable to control the x coordinates of segments & labels
my_x_me <- x_base_me + max(as.numeric(tree_dt_me$val), na.rm = TRUE)*my_factor + 1.5

#for singeltons
x_base_me2 <- max(tree_dt_me$x) + abs(min(as.numeric(tree_dt_me$val2), na.rm = TRUE))*my_factor + 2.6
# Define variable to control the x coordinates of segments & labels
my_x_me2 <- x_base_me + max(as.numeric(tree_dt_me$val2), na.rm = TRUE)*my_factor + 2.6

# Need to add a value (usually a small amount) to `ymax` 
# to force a complete circle (otherwise a stripe of white remains).
# This value needs to be just right:
# -  not too big because affects the visual which results in strange angles on labels;
# -  not too small because otherwise the remaining strip of white is not completely filled.
# Is a matter of try and error until you get it right. 
# Also, better to be biased towards smaller values since big value can lead to 
# big displacements when drawing the tree.
fill_in_value <- 0.0005


# Plot the tree with circular barplot
tree_bars_me <- 
  #tree_labeled_me +
  tree_plot_me +
  #Add a disc to plot bars on top of it
  geom_rect(data = tree_dt_me,
            aes(xmin = x_base_me + min(as.numeric(val), na.rm = TRUE)*my_factor,
                ymin = 0,
                xmax = x_base_me + max(as.numeric(val)*my_factor, na.rm = TRUE)*my_factor,
                # Add also fill_in_value to `ymax` to force a complete circle
                ymax = max(y) + fill_in_value), 
            color = NA, # set NA so to avoid coloring borders
            fill = "#F0F1F2",
            alpha = 0.5) +
  # Add bars for each tip label
  geom_rect(data = tree_dt_me,
            aes(xmin = x_base_me,
                ymin = y - 0.1,
                xmax = x_base_me + as.numeric(val)*my_factor,
                ymax = y + 0.1),
            fill = "#78B7C5",
            show.legend = FALSE,
            # no borders? color for now 
            color = "#78B7C5") +
  #add second set of rings for singletons
  #ring
  geom_rect(data = tree_dt_me,
            aes(xmin = x_base_me2 + min(as.numeric(val2), na.rm = TRUE)*my_factor,
                ymin = 0,
                xmax = x_base_me2 + max(as.numeric(val2)*my_factor, na.rm = TRUE)*my_factor,
                # Add also fill_in_value to `ymax` to force a complete circle
                ymax = max(y) + fill_in_value), 
            color = NA, # set NA so to avoid coloring borders
            fill = "#F0F1F2",
            alpha = 0.5) +
  #bars
  geom_rect(data = tree_dt_me,
            aes(xmin = x_base_me2,
                ymin = y - 0.1,
                xmax = x_base_me2 + as.numeric(val2)*my_factor,
                ymax = y + 0.1),
            #width = 1,
            fill = "#F21A00",
            show.legend = FALSE,
            # no borders 
            color = "#F21A00") +
theme(
  # Set font size & family - affects legend only 
  # "sans" = "Arial" and is the default on Windows OS; check windowsFonts()
  text = element_text(size = 10, family = "sans"),
  # Grab bottom-right (x=1, y=0) legend corner 
  legend.justification = c(1,0),
  # and position it in the bottom-right plot area.
  legend.position = c(1.1, 0.29),
  legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
  # Set height of legend items (keys). This affects spacing between them as well.
  legend.key.height = unit(4, "mm"),
  # Set margin around entire plot.
  plot.margin = unit(c(t = .8, r = 2.6, b = -0.6, l = 0), "cm") 
)

tree_bars_me + geom_tiplab2(size = 1, offset = .03, align = TRUE, linesize = .04, linetype = 1)

#ggsave(file="pan_genome_tree.png",device="png")



###are there unique gene families (or unique losses) by clade?
View(accessory_by_strain_matched)
#count data per family


#need to match annotation data to the all_accessory data 
###modify this:
##which strain has the highest/lowest number of singletons and accessory gene fams? 
gene_fam_by_strain_w_anno<-as.data.frame(PIRATE.gene_families[,c(1:4,21:ncol(PIRATE.gene_families))])
##make binary (if gene = 1, if not = 0)
#fill in zeros
gene_fam_by_strain_w_anno_zeros<- sapply(gene_fam_by_strain_w_anno, gsub, pattern = "^\\s*$" , replacement = 0 )
#fill in ones
gene_fam_by_strain_w_anno_ones<- as.data.frame(replace(gene_fam_by_strain_w_anno_zeros[,5:ncol(gene_fam_by_strain_w_anno_zeros)], gene_fam_by_strain_zeros!="0", 1))
#change to numeric
gene_fam_by_strain_w_anno_ones_num <- mutate_all(gene_fam_by_strain_w_anno_ones, function(x) as.numeric(as.character(x)))
#bind annotations
gene_fam_by_strain_w_anno_ones_num2<- cbind(gene_fam_by_strain_w_anno_zeros[,1:4], gene_fam_by_strain_w_anno_ones_num)


#subset to remove core genes and singletons 
#remove core
all_accessory_anno<- gene_fam_by_strain_w_anno_ones_num2[rowSums(gene_fam_by_strain_w_anno_ones_num2[,5:ncol(gene_fam_by_strain_w_anno_ones_num2)]) != 308,]
#remove singeltons
all_accessory_anno_no_sing<- all_accessory_anno[rowSums(all_accessory_anno[,5:ncol(all_accessory_anno)]) != 1,]
dim(all_accessory_anno_no_sing)

#split list of strain names by clade (as factors)
accessory_by_strain_matched_factor_list<- as.list(accessory_by_strain_matched$strain, accessory_by_strain_matched$clade)
accessory_by_strain_matched_list <- setNames(accessory_by_strain_matched$clade, accessory_by_strain_matched$strain)
#accessory_by_strain_matched_list_factor<- as.factor(accessory_by_strain_matched_list)

#split the count data frame by the factor list
#accessory_by_strain_matched_list_factor[accessory_by_strain_matched_list_factor == "clade1"]
clade1_names<- accessory_by_strain_matched_list[accessory_by_strain_matched_list == "clade1"]
clade1<- cbind(all_accessory_anno_no_sing[,1:4], all_accessory_anno_no_sing[,colnames(all_accessory_anno_no_sing) %in% names(clade1_names)])

clade2_names<- accessory_by_strain_matched_list[accessory_by_strain_matched_list == "clade2"]
clade2<- cbind(all_accessory_anno_no_sing[,1:4], all_accessory_anno_no_sing[,colnames(all_accessory_anno_no_sing) %in% names(clade2_names)])

clade3_names<- accessory_by_strain_matched_list[accessory_by_strain_matched_list == "clade3"]
clade3<- cbind(all_accessory_anno_no_sing[,1:4], all_accessory_anno_no_sing[,colnames(all_accessory_anno_no_sing) %in% names(clade3_names)])

clade4_names<- accessory_by_strain_matched_list[accessory_by_strain_matched_list == "clade4"]
clade4<- cbind(all_accessory_anno_no_sing[,1:4], all_accessory_anno_no_sing[,colnames(all_accessory_anno_no_sing) %in% names(clade4_names)])

clade5_names<- accessory_by_strain_matched_list[accessory_by_strain_matched_list == "clade5"]
clade5<- cbind(all_accessory_anno_no_sing[,1:4], all_accessory_anno_no_sing[,colnames(all_accessory_anno_no_sing) %in% names(clade5_names)])

clade6_names<- accessory_by_strain_matched_list[accessory_by_strain_matched_list == "clade6"]
clade6<- cbind(all_accessory_anno_no_sing[,1:4], all_accessory_anno_no_sing[,colnames(all_accessory_anno_no_sing) %in% names(clade6_names)])

clade7_names<- accessory_by_strain_matched_list[accessory_by_strain_matched_list == "clade7"]
clade7<- cbind(all_accessory_anno_no_sing[,1:4], all_accessory_anno_no_sing[,colnames(all_accessory_anno_no_sing) %in% names(clade7_names)])

clade8_names<- accessory_by_strain_matched_list[accessory_by_strain_matched_list == "clade8"]
clade8<- cbind(all_accessory_anno_no_sing[,1:4], all_accessory_anno_no_sing[,colnames(all_accessory_anno_no_sing) %in% names(clade8_names)])

#are there gene families that are exclusive to clade 1?
#clade1_rowsums<- clade1[rowSums(clade1[,5:ncol(clade1)]) >0,]

#calculate rowsums

#QC make sure they match 
#length(names(clade7)) -5 
#length(clade7_names)  #yep

clade1$totals<-  rowSums(clade1[,5:ncol(clade1)])
clade2$totals<-  rowSums(clade2[,5:ncol(clade2)])
clade3$totals<-  rowSums(clade3[,5:ncol(clade3)])
clade4$totals<-  rowSums(clade4[,5:ncol(clade4)])
clade5$totals<-  rowSums(clade5[,5:ncol(clade5)])
clade6$totals<-  rowSums(clade6[,5:ncol(clade6)])
clade7$totals<-  rowSums(clade7[,5:ncol(clade7)])
#clade8$totals<-  sum(clade8[,5:ncol(clade8)]) #because this is really only one... 

#get fams exclusive to clade1
exclusive_to_clade1<- clade1[(clade1$totals > 0) & 
                               (clade2$totals == 0) &
                               (clade3$totals == 0) &
                               (clade4$totals == 0) &
                               (clade5$totals == 0) &
                               (clade6$totals == 0) &
                               (clade7$totals == 0),]

dim(exclusive_to_clade1) 
#there are 261 accessory gene fams exclusive to clade 1
length(clade1_names) #there are 96 strains in clade 1
sort(exclusive_to_clade1$totals) # top hit is distributed in 46 of the 96 isolates in clade1
length(clade1_names)*.33 #none in over half the isolates in clade 1
#get all present in more than 1/3 of the isolates 
exclusive_to_clade1_of_note<- exclusive_to_clade1[exclusive_to_clade1$totals > length(clade1_names)*.33,]
dim(exclusive_to_clade1_of_note)
#there are 5 of note- all "hypothetical protein 

#spot check
#test<-exclusive_to_clade1[exclusive_to_clade1$totals == 46,]
#test$gene_family #g010930
#test2<- clade2[clade2$gene_family == "g010930",]




#exclusive to clade2
exclusive_to_clade2<- clade2[(clade2$totals > 0) & 
                               (clade1$totals == 0) &
                               (clade3$totals == 0) &
                               (clade4$totals == 0) &
                               (clade5$totals == 0) &
                               (clade6$totals == 0) &
                               (clade7$totals == 0),]

dim(exclusive_to_clade2) 
length(clade2_names) #there are 40 isolatres in clade 2
sort(exclusive_to_clade2$totals)
#there are 25 accessory gene fams exclusive to clade2, but not ver wide spread
exclusive_to_clade2_of_note<- exclusive_to_clade2[exclusive_to_clade2$totals > length(clade2_names)*.33,]
dim(exclusive_to_clade2_of_note)
#none of note


#exclusive to clade3
exclusive_to_clade3<- clade3[(clade3$totals > 0) & 
                               (clade1$totals == 0) &
                               (clade2$totals == 0) &
                               (clade4$totals == 0) &
                               (clade5$totals == 0) &
                               (clade6$totals == 0) &
                               (clade7$totals == 0),]

dim(exclusive_to_clade3) 
length(clade3_names)
sort(exclusive_to_clade3$totals)
#there are 19 accessory gene fams exclusive to clade3, but not ver wide spread
exclusive_to_clade3_of_note<- exclusive_to_clade3[exclusive_to_clade3$totals > length(clade3_names)*.33,]
dim(exclusive_to_clade3_of_note)
#none of note



#exclusive to clade4
exclusive_to_clade4<- clade4[(clade4$totals > 0) & 
                               (clade1$totals == 0) &
                               (clade2$totals == 0) &
                               (clade3$totals == 0) &
                               (clade5$totals == 0) &
                               (clade6$totals == 0) &
                               (clade7$totals == 0),]

dim(exclusive_to_clade4) 
length(clade4_names)
sort(exclusive_to_clade4$totals)
#there are 7 accessory gene fams exclusive to clade4, but not ver wide spread
exclusive_to_clade4_of_note<- exclusive_to_clade4[exclusive_to_clade4$totals > length(clade4_names)*.33,]
dim(exclusive_to_clade4_of_note)
#none of note


#exclusive to clade5
exclusive_to_clade5<- clade5[(clade5$totals > 0) & 
                               (clade1$totals == 0) &
                               (clade2$totals == 0) &
                               (clade3$totals == 0) &
                               (clade4$totals == 0) &
                               (clade6$totals == 0) &
                               (clade7$totals == 0),]

dim(exclusive_to_clade5) 
length(clade5_names)
sort(exclusive_to_clade5$totals)
#there are 165 accessory gene fams exclusive to clade2, but not ver wide spread
exclusive_to_clade5_of_note<- exclusive_to_clade5[exclusive_to_clade5$totals > length(clade5_names)*.33,]
dim(exclusive_to_clade5_of_note)
#none of note


#exclusive to clade6
exclusive_to_clade6<- clade6[(clade6$totals > 0) & 
                               (clade1$totals == 0) &
                               (clade2$totals == 0) &
                               (clade3$totals == 0) &
                               (clade4$totals == 0) &
                               (clade5$totals == 0) &
                               (clade7$totals == 0),]

dim(exclusive_to_clade6) 
length(clade6_names)
sort(exclusive_to_clade6$totals)
#there are 69 accessory gene fams exclusive to clade2, but not ver wide spread
length(clade6_names)*.5 #none in more than half
exclusive_to_clade6_of_note<- exclusive_to_clade6[exclusive_to_clade6$totals > length(clade6_names)*.33,]
dim(exclusive_to_clade6_of_note)
#one of note - a hypothetical protein 

#exclusive to clade7
exclusive_to_clade7<- clade7[(clade7$totals > 0) & 
                               (clade1$totals == 0) &
                               (clade2$totals == 0) &
                               (clade3$totals == 0) &
                               (clade4$totals == 0) &
                               (clade5$totals == 0) &
                               (clade6$totals == 0),]

dim(exclusive_to_clade7) 
length(clade7_names)
sort(exclusive_to_clade7$totals)
#there are 3 accessory gene fams exclusive to clade2, but not very wide spread
exclusive_to_clade7_of_note<- exclusive_to_clade7[exclusive_to_clade7$totals > length(clade7_names)*.33,]
dim(exclusive_to_clade7_of_note)
#2 of note - both "hypothetical protein"


#graph the fams of interest onto the tree 
#clade 1: (5) exclusive_to_clade1_of_note
#clade 6: (1) exclusive_to_clade6_of_note
#clade 7: (2) exclusive_to_clade7_of_note

#tree of clade1 w/ presence absence for the fams of interest


clade1_names



##subset the clade-trees from the larger tree
#fixnames calde1
clade1_names_fixed <- name_map$name_pop_genome[match(names(clade1_names), name_map$name_Pan_genome)]
#subset tree to just clade1
clade1_tree<- keep.tip(tree_grA_me, clade1_names_fixed)
#fixnames calde6
clade6_names_fixed <- name_map$name_pop_genome[match(names(clade6_names), name_map$name_Pan_genome)]
#subset tree to just clade1
clade6_tree<- keep.tip(tree_grA_me, clade6_names_fixed)
#fixnames calde7
clade7_names_fixed <- name_map$name_pop_genome[match(names(clade7_names), name_map$name_Pan_genome)]
#subset tree to just clade1
clade7_tree<- keep.tip(tree_grA_me, clade7_names_fixed)


#clade 1 plot
clade1_tree_plot <- 
  ggtree(tr = clade1_tree, 
         # color by group attribute, check str(tree_grA_me)
         #mapping = aes(color = group), 
         layout  = 'circular', 
         branch.length = "none") + 
  geom_treescale(x=3, y=NULL, color = NA) +
  scale_color_manual(name = 'Clade', values = my_cols_me) + 
  theme(legend.title=element_text(size=9), # The title of legend 
        legend.text=element_text(size=7)) +
  guides(color = guide_legend(override.aes = list(size = 4))) 


clade1_tree_plot + geom_tiplab(size = 1, align = TRUE, linesize = .25, linetype = 0)

##plot presence / absence of gene fams 

#mutate cols into rows
row.names(exclusive_to_clade1_of_note)<- exclusive_to_clade1_of_note$gene_family
clade_1_binary<- data.frame(t(exclusive_to_clade1_of_note[,5:(ncol(exclusive_to_clade1_of_note) -1)]))
#turn into presence absence 
clade_1_binary_presence<- data.frame(sapply(clade_1_binary, gsub, pattern = "1", replacement = "present"))
clade_1_binary_presence_absence<- sapply(clade_1_binary_presence, gsub, pattern = "0", replacement = "absent")
View(clade_1_binary_presence_absence)
row.names(clade_1_binary_presence_absence) <- row.names(clade_1_binary)
#fix rownames
clade_1_binary_presence_absence_fixed <- name_map$name_pop_genome[match(rownames(clade_1_binary_presence_absence), name_map$name_Pan_genome)]
row.names(clade_1_binary_presence_absence) <- clade_1_binary_presence_absence_fixed










#what is "traits" here?
clade1_tree_plot2 <-  gheatmap(clade1_tree_plot, clade_1_binary_presence_absence, offset=0.2, width=0.2, low="white", high="black", colnames_position = "top", font.size=2, color="black") +
  scale_fill_manual(values=c("white", "black"))

clade1_tree_plot2
clade1_tree_plot2 + geom_tiplab(size = 1, align = TRUE, linesize = .25, offset = 15, linetype = 0)

#add plot title Clade 1 specific families 




###where do the genes fall on each chr?
#graph where the core genes fall on the chromosome 

#graph where the accessory genes fall on the chromosome

#how many genes fall at chromosome ends? (in the last X BPs). 

