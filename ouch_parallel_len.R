library(ouch)
library(geiger)
library(phytools)

key <- 'totalNoGap'
args <- commandArgs(trailingOnly=TRUE)
direct <- '/data1/home/dbarth/LinkerProject/Newicks'

#some string processing to unpack info we need
a <- strsplit(args[1],'/')[[1]][2]
a <- strsplit(a,'_') #cmd line arg (gt_id, linker, .props.txt)

if(identical(a[[1]][2],'Linker')){
    parent <- a[[1]][1]
    gt <- parent
    nolengths <- strsplit(a[[1]][4],'.lengths.txt')[[1]][1]
    linker <- paste(a[[1]][2],a[[1]][3],nolengths,sep='_')
} else {
    parent <- a[[1]][1]
    gt <- paste(a[[1]][1],a[[1]][2],sep='_')
    nolengths <- strsplit(a[[1]][5],'.lengths.txt')[[1]][1]
    linker <- paste(a[[1]][3],a[[1]][4],nolengths,sep='_')
    }

print(a[[1]]) 
#load in data
tree <- read.tree(paste(direct,paste(gt,'newick',sep='.'),sep='/')) # load tree
data <- read.table(args[1],header=TRUE,row.names=1)
data <- data[c(key)]
treedata <- list(tree, data)
sapply(treedata,class)
nc <- with(treedata,name.check(tree,data))
tree <- with(treedata,drop.tip(tree,nc$tree_not_data))

#if the tree has 10+ leaf nodes
if((length(tree$tip.label)>10)){

     ### make an ouchtree out of the phy-format tree
     ot <- ape2ouch(tree)
     
     ### merge data with tree info
     otd <- as(ot,"data.frame")
     ### in these data, it so happens that the rownames correspond to node names
     ### we will exploit this correspondence in the 'merge' operation:
     data$labels <- rownames(data)
     otd <- merge(otd,data,by="labels",all=TRUE)
     rownames(otd) <- otd$nodes
     otd$totalNoGap<-as.numeric(otd$totalNoGap)
     ### this data-frame now contains the data as well as the tree geometry

     ### now remake the ouch tree
     ot <- with(otd,ouchtree(nodes=nodes,ancestors=ancestors,times=times,labels=labels))

     setwd("/data1/home/dbarth/LinkerProject/len_residues_kco/") #changed for testing purposes, orig: /data1/home/dbarth/LinkerProject/ModelLengths/
     sink(file=paste('len_model_fit',gt,paste(linker,'.txt',sep=''),sep='_'),append=FALSE)
     b1 <- brown(tree=ot,data=otd['totalNoGap'])
     print("brown")
     print(summary(b1))
     ### evaluate an OU model with a single, global selective regime
     otd$regimes <- as.factor("global")
     h1 <- hansen(
                  tree=ot,
                  data=otd['totalNoGap'],
                  regimes=otd["regimes"],
                  sqrt.alpha=c(1),
                  sigma=c(1),
                  maxit=10000
                  )
     print("hansen")
     print(summary(h1))

     wn <- fitContinuous(phy=tree,dat=data['totalNoGap'],model='white')
     print("white_noise")
     print(wn$opt)

     sink()
     ## End(Not run)
}

