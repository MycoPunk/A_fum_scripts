
#load packages
library(pophelper)
library(grid)
library(gridExtra)

setwd("~/Desktop/Project_Afum_pangenome/")

#first get the group assignments
inds <- read.table("strain_names_in_order.csv",header=FALSE,stringsAsFactors=F)
#inds$V1

ffiles <- list.files(path=".",pattern="*.meanQ",full.names=T)
flist <- readQ(files=ffiles)

rownames(flist[[1]]) <- (inds$V1)
if(length(unique(sapply(flist,nrow)))==1) flist <- lapply(flist,"rownames<-",inds$V1)
# show row names of all runs and all samples
#lapply(flist, rownames)

# view head of first converted file
head(flist[[1]])
#view file atributes like this
attributes(flist[[1]])

#tabulateQ
tr1 <- tabulateQ(qlist=flist)
tabulateQ(flist)
tabulateQ(flist, writetable=TRUE, sorttable=TRUE)

#summarize Q
sr1 <- summariseQ(tr1)
summariseQ(tr1, writetable=TRUE)

#plot
p1 <- plotQ(alignK(sortQ(flist)[1:5]),imgoutput="join",returnplot=T,exportplot=F,quiet=T,basesize=11,
            clustercol=c("#56326E",
                         "#ED7F6F",
                         "#ABA778",
                         "#F2CC85", 
                         "#D4494E",
                         "#7D7A70"),
            showlegend=T, legendlab=c("Clade 1","Clade 2","Clade 3","Clade 4","Clade 5","Clade 6"),
            legendkeysize=10,legendtextsize=10,legendmargin=c(2,2,2,0),
            showindlab=F,useindlab=T,showyaxis=T, indlabsize = 4,
            sortind="all", sharedindlab = F, panelspacer=0.4,
            barbordersize=1)
grid.arrange(p1$plot[[1]])



#my_cols_me<-c("#56326E","#ED7F6F","#ABA778","#F2CC85", "#D4494E","#7D7A70")
#scales::show_col(my_cols_me); my_cols_me


#get likely admixed individuals 

admixed<- flist2[apply(flist2[1:3], 1, function(x) all(x < .75)), ] 
#here there are 6 of them 

#Afum_266.3.Cluster1 Afum_266.3.Cluster2 Afum_266.3.Cluster3
#AF293             0.605493            0.372659            0.021848
#CF098             0.605721            0.372531            0.021748
#CM7632            0.621613            0.345604            0.032783
#F18304            0.599702            0.378497            0.021801
#F7763             0.699989            0.146177            0.153834
#RSF2S8            0.363823            0.636174            0.000004


#plot admixed individuals only
admixed_names<- rownames(admixed)
admixed_names
#subset input 
only_admixed<- sapply(flist,function(x) x[match(c(admixed_names),rownames(x)),])

#plot
p1 <- plotQ(alignK(sortQ(only_admixed)[2]),returnplot=T,exportplot=F,quiet=T,basesize=20,barsize = .5,
            clustercol=c("#56326E",
                         "#ED7F6F",
                         "#ABA778",
                         "#F2CC85", 
                         "#D4494E",
                         "#7D7A70"),
            showlegend=T, legendlab=c("Clade 1","Clade 2","Clade 3"),
            legendkeysize=10,legendtextsize=10,legendmargin=c(2,2,2,0),
            showindlab=T,useindlab=T,showyaxis=T, indlabsize = 10,
            sortind="all", sharedindlab = F, panelspacer=0.4,
            barbordersize=0)
grid.arrange(p1$plot[[1]])
