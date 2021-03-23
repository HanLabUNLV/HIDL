library(ouch)
library(geiger)
library(phytools)

residue <- 'E' #special residues P, G, A
args <- commandArgs(trailingOnly=TRUE)
direct <- '/data1/home/dbarth/LinkerProject/Newicks'

#some string processing to unpack info we need
a <- strsplit(args[1],'/')[[1]][2]
a <- strsplit(a,'_') #cmd line arg (gt_id, linker, .props.txt)

if(identical(a[[1]][2],'Linker')){
    parent <- a[[1]][1]
    gt <- parent
    noprops <- strsplit(a[[1]][4],'.props.txt')[[1]][1]
    linker <- paste(a[[1]][2],a[[1]][3],noprops,sep='_')
} else { 
    parent <- a[[1]][1]
    gt <- paste(a[[1]][1],a[[1]][2],sep='_')
    noprops <- strsplit(a[[1]][5],'.props.txt')[[1]][1]
    linker <- paste(a[[1]][3],a[[1]][4],noprops,sep='_')
    }
 
#load in data
tree <- read.tree(paste(direct,paste(gt,'newick',sep='.'),sep='/')) # load tree
data <- read.table(args[1],header=TRUE)

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
     #print(otd)
     ### this data-frame now contains the data as well as the tree geometry

     ### now remake the ouch tree
     ot <- with(otd,ouchtree(nodes=nodes,ancestors=ancestors,times=times,labels=labels))

     setwd(paste("/data1/home/dbarth/LinkerProject/",paste(residue,'residues',sep='_'),sep='/'))
     sink(file=paste('model_fit',gt,paste(linker,'.txt',sep=''),sep='_'),append=TRUE)
     b1 <- brown(tree=ot,data=otd[c(residue)])
     print("brown")
     print(summary(b1))
     
     ### evaluate an OU model with a single, global selective regime
     otd$regimes <- as.factor("global")
     h1 <- hansen(
                  tree=ot,
                  data=otd[c(residue)],
                  regimes=otd["regimes"],
                  sqrt.alpha=c(1),
                  sigma=c(1),
                  maxit=10000
                  )
     print("hansen")
     print(summary(h1))
     sink()
     ## End(Not run)
}

