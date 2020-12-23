#this script identifies the number of clusters present in the A fum dataset. 
#using DAPC in the the adegnet package

#set wd
setwd("~/Desktop/Project_Afum_pangenome/")

library("poppr")
library("ape") # To visualize the tree using the "nj" function
library("magrittr")
library("adegenet")
library("reshape")
library('vcfR')
library('parallel')
library('parallel')

#set seed for reproducibility
set.seed(666)

#read in vcf

#vcf_me<- read.vcfR("pop_for_pan.SNP.combined_selected.NO_TEs.vcf.gz") #old version

vcf_me<- read.vcfR("pop_for_pan3_.SNP.combined_selected.NO_TEs.vcf.gz") #new version W/o TE's and w 267 strains


#convert vcf to genlight object
gl_Afum <- vcfR2genlight(vcf_me)
rm(vcf_me)

#run find clusters -note this takes about 90 minutes at ~61,000 SNPS
grp1<-find.clusters(gl_Afum, max.n.clust=20) #note n-clusters needs to be less than the n of individuals
#Choose the number PCs to retain (>=1): 
#  200 #for the selection of k, there is no benefit to reducing the number of PCs
#Choose the number of clusters (>=2): 
#  3 #3 is where the 'elbow' in the curve is


###run DAPC

#first, determine the number of PC's to retain in the analysis 
#run dapc with default settings
dapc1 <- dapc(gl_Afum, grp1$grp, n.da=100, n.pca=20)
temp <- optim.a.score(dapc1)
temp #a-score recommends 2 PCs #5 with 3 groups with 20 possible


#re-run dapc with 
#dapc2 <- dapc(gl_Afum, grp1$grp, n.da=100, n.pca=6)
dapc2 <- dapc(gl_Afum, grp1$grp, n.da=100, n.pca=3)

#plot the k=4 population groups 
#colors for this project: "#382147","#7D7A70","#ABA778","#F2CC35", "#F2CC85", "#FFCA98", "#F7A583","#ED7F6F","#D4494E"
myCol <- c(clade1 = "#56326E",
              clade2 ="#ED7F6F",
                clade3 = "#ABA778")
                #clade4 = "#F7C165")

#scales::show_col(myCol); myCol

scatter(dapc2, posi.da="bottomright", bg="white",
        pch=20, 
        #cell=0, 
        cstar=0, 
        col=myCol, 
        scree.pca=TRUE,
        posi.pca="bottomleft",
        solid=.4,
        cex = 1.5,
        clab=0,
        leg=TRUE,
        txt.leg=paste("Cluster",1:3))


#in a single dimension 
scatter(dapc2,1,1, col=myCol, bg="white",
        scree.da=FALSE, legend=TRUE, solid=.4)

##extract group membership
#summary(dapc2)
#how many strains / groups?
dim(dapc2$posterior)
##  267   3
#assignplot(dapc2, subset=1:150,cex.lab=.10,pch=1)
#round(head(dapc2$posterior),3)

#get groups
groups<- data.frame(clade = dapc2$grp)

#STRUCTURE-like plot (complot)
#look for the most admixed individuals
admixed <- which(apply(dapc2$posterior,1, function(e) all(e<0.9))) #those having no more than 0.9 probability of membership to any group
length(admixed)
#there are no individuals likely to be admixed between the three populations


#graph the four admixed strains STRUCTURE style
par(mar=c(8,4,5,1), xpd=TRUE)
#compoplot(dapc2, subset=admixed, cleg=.6, posi=list(x=0,y=1.2), lab=names(dapc2$grp), show.lab = TRUE)

#plot all
compoplot(dapc2, cleg=.6, posi=list(x=0,y=1.2), lab=names(dapc2$grp), col.pal = myCol)


##print clade assignments
#write.table(groups, "DAPC_groupings_no_TEs_20Dec2020.txt", sep="\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
name_map<-read.delim("namemappingfile_12Oct2020.tab", header = TRUE, sep = "\t", fill = TRUE, strip.white = TRUE)

#clades<-read.delim("DAPCA_groupings.txt", header = TRUE, sep = "\t", fill = TRUE, strip.white = TRUE)
clades<- data.frame(cbind(name_pop_genome= rownames(groups), groups))
#match pan_genome_names w/ name_Pan_genome

#change oddly formatted strains
clades$name_pop_genome[clades$name_pop_genome=="AFIS1435CDC_6"] <- "AFIS1435_CDC_6"
clades$name_pop_genome[clades$name_pop_genome=="Afu_343-P/11"] <- "Afu_343-P-11"

#match names
clades$name_Pan_genome <- name_map$name_Pan_genome[match(clades$name_pop_genome, name_map$name_pop_genome)]

#looks good - print it
write.table(clades, "clade_map_20Dec2020_3groups.tab", sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
View(clades)
rm(dapc1)
#####try to determine number of PC's with xval estimation

#xval<- xvalDapc(tab(gl_Afum, NA.method = "mean"), grp1$grp,
#         n.pca = 1:8, n.rep = 80,
#         parallel = "multicore", ncpus = 6L)


#gl_Afum <- vcfR2genlight(vcf_me)
#rm(vcf_me)
#grp1<-find.clusters(gl_Afum, max.n.clust=20) #choose to retain 200 PCs
#xval<- xvalDapc(tab(gl_Afum, NA.method = "mean"), grp1$grp,
#                 n.pca.max = 200, training.set = 0.9,
#                 result = "groupMean", center = TRUE, scale = FALSE,
#                 n.pca = NULL, n.rep = 30, xval.plot = TRUE,
#                 parallel = "multicore", ncpus = 6L)

#xval
#$`Number of PCs Achieving Lowest MSE`
#[1] "180"

#note you can make a graph directly from here, rather than re-running dapc
#xval$DAPC

#number of PCs
#numPCs = as.numeric(xval$`Number of PCs Achieving Lowest MSE`)

#now run DACP
#note n.da is the number of axes retained in the discriminant analysis- leave off to identify interactively
#dapc1 = dapc(gl_Afum, grp1$grp, n.pca = numPCs)



#re-run dapc with 
#dapc2 <- dapc(gl_Afum, grp1$grp, n.da=100, n.pca=6)
#dapc2 <- dapc(gl_Afum, grp1$grp, n.da=100, n.pca=<>) #change n.pca to match the output of xval, lowest MSE
