library(ouch)
library(geiger)
library(phytools)

args <- commandArgs(trailingOnly=TRUE)
residue <- args[2]
#args <- 'LinkerFastasFeatures/ENSGT00550000075016_0_Linker_4_5.AAC.tsv' 
direct <- '/data2/han_lab/dbarth/LinkerProject'
#some string processing to unpack info we need
a <- strsplit(args[1],'/')[[1]][2]
if(identical(a[[1]][2],'Linker')){
    parent <- a[[1]][1]
    gt <- parent
    nolengths <- strsplit(a[[1]][4],split='.',fixed=TRUE)[[1]][1]
    linker <- paste(a[[1]][2],a[[1]][3],nolengths,sep='_')
} else {
    split <- strsplit(a,split='_',fixed=TRUE)
    if(length(split[[1]])>3){
      #subtree
      gt <- paste(split[[1]][1], split[[1]][2],sep='_')
      linker <- paste(split[[1]][3], strsplit(split[[1]][4], split='.AAC')[[1]][1],sep='_')  
    } else {
      gt <- split[[1]][1]
      linker <- paste(split[[1]][2], strsplit(split[[1]][3], split='.AAC')[[1]][1],sep='_')  
    }
}
#load in tree
tree <- read.tree(paste(direct,'Newicks',paste(gt,'newick',sep='.'),sep='/')) # load tree


#find short branch lengths, delete them
tip.branches = match(1:length(tree$tip), tree$edge[,2])
bls = cbind.data.frame(taxa=tree$tip, tree$edge[tip.branches,], bl=tree$edge.length[tip.branches])
drop2 = bls[bls$bl<0.0001,]
drop2$taxa = as.character(drop2$taxa)

#change domain/linker here
#get rownames from file
rows <- read.table(paste('KCO_Domain_pids/',gt,'_',linker,'_pids.txt',sep=''), header=FALSE)
#need 10 tips to do analysis
if((dim(rows)[1]<10)){
    stop(paste("Error: not enough tips for tree ",gt,sep=''))
}

#read in data, filter it
data <- read.table(args[1],header=FALSE)
data <- data[,-1]
colnames(data) <- c('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y')
data=data[,colSums(data)!=0]

rownames(data) <- t(rows)
treedata <- list(tree, data)
sapply(treedata,class)
nc <- with(treedata,name.check(tree,data))

if (!is.atomic(nc)) { #if there's no tips to drop, then pass
    tree <- with(treedata,drop.tip(tree,nc$tree_not_data))
}

tree <- with(treedata,drop.tip(tree, drop2$taxa)) #delete short bl
nc <- with(treedata,name.check(tree,data))

if (!is.atomic(nc)) { #if there's no data to drop, then pass
    data <- as.data.frame(data[-which(rownames(data) %in% nc$data_not_tree),,drop=FALSE])
}
#this may seem like the dumbest line of code ever, but I'm pretty sure without it, everything will fail
tree <- read.tree(text=write.tree(tree))
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
     varnames <- colnames(otd)
     varnames <- varnames[-which(varnames %in% c("labels", "nodes", "ancestors", "times", "regimes"))]      
     lowtri = sum(1:length(varnames))
     otd[varnames] <- sapply(otd[varnames], as.numeric, simplify=FALSE)
     ### this data-frame now contains the data as well as the tree geometry
     ### now remake the ouch tree
     ot <- with(otd,ouchtree(nodes=nodes,ancestors=ancestors,times=times,labels=labels))

     #set a file to print output
     setwd("/local/han_lab/AAC_uv_results_domain") #changed for testing purposes, orig: /data1/home/dbarth/LinkerProject/ModelLengths/
     sink(file=paste(residue,'model_fit',gt,paste(linker,'.txt',sep=''),sep='_'),append=FALSE)

     #BM model
     b1 <- brown(tree=ot,data=otd[residue])
     print("brown")
     print(summary(b1))

     ### evaluate an OU model with a single, global selective regime
     otd$regimes <- as.factor("global")
     h1 <- hansen(
                  tree=ot,
                  data=otd[residue],
                  regimes=otd["regimes"],
                  sqrt.alpha=c(1), 
                  method="Nelder-Mead",
                  sigma=c(1), #t(rep.int(1,lowtri)),
                  maxit=100000
                  )
     print("hansen")
     print(summary(h1))

     #WN model
     wndata = as.matrix(data[residue])
     wn <- fitContinuous(phy=tree,dat=wndata,model='white')
     print("white_noise")
     print(wn)

     sink()
     ## End(Not run)
}






