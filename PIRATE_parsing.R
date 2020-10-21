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

#mapping file to map old names to new names 
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


ggsave(file="pan_genome_tree.png",device="png")

help(ggsave)


###where do the genes fall on each chr?
#graph where the core genes fall on the chromosome 

#graph where the accessory genes fall on the chromosome

#how many genes fall at chromosome ends? (in the last X BPs). 

