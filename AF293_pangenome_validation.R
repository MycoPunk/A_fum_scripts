#A fum pan genome analysis: 
#validate parameters for PIRATE pangenome 
#here, AF293 should have approximately the same number of gene families as there are gene models for AF293
#there are currently 9781 gene models for AF293 (https://mycocosm.jgi.doe.gov/Aspfu1/Aspfu1.info.html)


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
library(ggplot2)
library(stringr)



setwd("~/Desktop/Project_Afum_pangenome/")
PIRATE.gene_families<-as.data.frame(fread("PIRATE.gene_families.ordered_ver13.tsv")) 




##basic stats:
#how many gene fmailies are there?
n_gene_fams<- nrow(PIRATE.gene_families)
n_gene_fams

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
                    "average_length(bp)", 
                    "cluster", 
                    "cluster_order")

strains_only<- PIRATE.gene_families[,!names(PIRATE.gene_families) %in% cols_to_exclude]
strain_names<- names(strains_only)

length(strain_names)
length(unique(strain_names)) 

#266 strains

#how many of the gene families are in every genome
n_gene_fams_core<- sum(PIRATE.gene_families$number_genomes == 266)


#that's what percent out of the total?
(n_gene_fams_core*100)/n_gene_fams

#w/wobble = 1
n_gene_fams_core_w1<- sum(PIRATE.gene_families$number_genomes >= 265)
(n_gene_fams_core_w1*100)/n_gene_fams


#w/wobble = 2
n_gene_fams_core_w2<- sum(PIRATE.gene_families$number_genomes >= 264)
(n_gene_fams_core_w2*100)/n_gene_fams


#by percentage: 
#present in 99% of genomes 
.99*266
n_gene_fams_core_w99per<- sum(PIRATE.gene_families$number_genomes >= .99*266)
(n_gene_fams_core_w99per*100)/n_gene_fams


#present in 95% of genomes
.95*266
n_gene_fams_core_w95per<- sum(PIRATE.gene_families$number_genomes >= .95*266)
(n_gene_fams_core_w95per*100)/n_gene_fams

n_gene_fams - n_gene_fams_core_w95per

#how many of the gene families are singletons (accessory)
n_gene_fams_singletons<- sum(PIRATE.gene_families$number_genomes == 1)

#that's what percent out of the total?
(n_gene_fams_singletons*100)/n_gene_fams

##which strain has the highest/lowest number of singletons and accessory gene fams? 
gene_fam_by_strain<-as.data.frame(PIRATE.gene_families[,23:ncol(PIRATE.gene_families)])
ncol(gene_fam_by_strain)


##make binary (if gene = 1, if not = 0)
#fill in zeros
gene_fam_by_strain_zeros<- sapply(gene_fam_by_strain, gsub, pattern = "^\\s*$" , replacement = 0 )
#fill in ones
gene_fam_by_strain_ones<- as.data.frame(replace(gene_fam_by_strain_zeros, gene_fam_by_strain_zeros!="0", 1))
#change to numeric
gene_fam_by_strain_ones_num <- mutate_all(gene_fam_by_strain_ones, function(x) as.numeric(as.character(x)))

#there were a total of how many gene fams found in AF293?
sum(gene_fam_by_strain_ones_num$Afum_AF293)
#out of how many gene fams?
nrow(gene_fam_by_strain_ones_num)


AF293<- as.data.frame(PIRATE.gene_families$Afum_AF293)
dim(AF293)

multiple_genes_in_fam<- AF293 %>% 
  filter(str_detect(AF293$`PIRATE.gene_families$Afum_AF293`, pattern = ":"))


# Count the number of 'a's in each element of string
multiple_genes_in_fam$number.of.multiples <- str_count(multiple_genes_in_fam$`PIRATE.gene_families$Afum_AF293`, ":")
View(multiple_genes_in_fam)
sum(multiple_genes_in_fam$number.of.multiples)
nrow(multiple_genes_in_fam)


#over all
AF293<- as.data.frame(PIRATE.gene_families$Afum_AF293)
dim(AF293)

#replace ; w/ :
AF293$`PIRATE.gene_families$Afum_AF293`<- gsub("\\;", ":", AF293$`PIRATE.gene_families$Afum_AF293`)

#remove null values (fams not in AF293)
AF293_subset<-  data.frame(AF293[!(is.na(AF293$`PIRATE.gene_families$Afum_AF293`) | AF293$`PIRATE.gene_families$Afum_AF293`==""), ])
dim(AF293)
#15623
dim(AF293_subset)
#9583

#get totals for how many genes are in each fam
AF293_subset$number_of_genes_in_fam <- str_count(AF293_subset$AF293...is.na.AF293..PIRATE.gene_families.Afum_AF293.....AF293..PIRATE.gene_families.Afum_AF293....., ":")
table(AF293_subset$number_of_genes_in_fam)

#which means there are how many gene fams without a representative in AF293? 
nrow(AF293) - nrow(AF293_subset)
#6040

#did any genes in AF293 get assigned to more than one family?

#remove parens for matching
AF293_subset$AF293...is.na.AF293..PIRATE.gene_families.Afum_AF293.....AF293..PIRATE.gene_families.Afum_AF293.....<- gsub("\\(", "", AF293_subset$AF293...is.na.AF293..PIRATE.gene_families.Afum_AF293.....AF293..PIRATE.gene_families.Afum_AF293.....)
AF293_subset$AF293...is.na.AF293..PIRATE.gene_families.Afum_AF293.....AF293..PIRATE.gene_families.Afum_AF293.....<- gsub("\\)", "", AF293_subset$AF293...is.na.AF293..PIRATE.gene_families.Afum_AF293.....AF293..PIRATE.gene_families.Afum_AF293.....)

#remove doubbles and tripples
AF293_singles_only<-  AF293_subset %>% 
  filter(str_detect(AF293_subset$number_of_genes_in_fam, pattern = "0"))

AF293_doubbles_only<-  AF293_subset %>% 
  filter(str_detect(AF293_subset$number_of_genes_in_fam, pattern = "1"))

dim(AF293_doubbles_only)

AF293_tripples_only<-  AF293_subset %>% 
  filter(str_detect(AF293_subset$number_of_genes_in_fam, pattern = "2"))

AF293_quadruple_only<-  AF293_subset %>% 
  filter(str_detect(AF293_subset$number_of_genes_in_fam, pattern = "3"))
View(AF293_quadruple_only)

#split the strings for doubbles and tripples strings at ":" 
doubbles_to_append<- AF293_doubbles_only %>% separate(AF293...is.na.AF293..PIRATE.gene_families.Afum_AF293.....AF293..PIRATE.gene_families.Afum_AF293....., c("A","B"), sep = ":")
tripples_to_append<- AF293_tripples_only %>% separate(AF293...is.na.AF293..PIRATE.gene_families.Afum_AF293.....AF293..PIRATE.gene_families.Afum_AF293....., c("D","E", "F"), sep = ":")
quadruples_to_append<- AF293_quadruple_only %>% separate(AF293...is.na.AF293..PIRATE.gene_families.Afum_AF293.....AF293..PIRATE.gene_families.Afum_AF293....., c("G","H", "I", "J"), sep = ":")

#appned 
All_AF293_genes<- c(AF293_singles_only$AF293...is.na.AF293..PIRATE.gene_families.Afum_AF293.....AF293..PIRATE.gene_families.Afum_AF293....., 
                    doubbles_to_append$A, 
                    doubbles_to_append$B, 
                    tripples_to_append$D,
                    tripples_to_append$E,
                    tripples_to_append$F,
                    quadruples_to_append$G,
                    quadruples_to_append$H,
                    quadruples_to_append$I,
                    quadruples_to_append$J)

#check that no genes are being assigned to more than one family
length(All_AF293_genes) #this is correct. 
length(unique(All_AF293_genes)) #no overlaps - no genes are being assigned to more than one fam. 

