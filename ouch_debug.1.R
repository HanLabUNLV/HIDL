library(ouch)
library(geiger)
library(phytools)

key <- 'totalNoGap'
#args <- commandArgs(trailingOnly=TRUE)
args <- 'LinkerFastasFeatures/ENSGT00550000075016_0_Linker_4_5.GAAC.tsv' 
direct <- '/data2/han_lab/dbarth/LinkerProject'

#some string processing to unpack info we need
a <- strsplit(args[1],'/')[[1]][2]
a <- strsplit(a,'_') #cmd line arg (gt_id, linker, .props.txt)

if(identical(a[[1]][2],'Linker')){
    parent <- a[[1]][1]
    gt <- parent
    nolengths <- strsplit(a[[1]][4],split='.',fixed=TRUE)[[1]][1]
    linker <- paste(a[[1]][2],a[[1]][3],nolengths,sep='_')
} else {
    parent <- a[[1]][1]
    gt <- paste(a[[1]][1],a[[1]][2],sep='_')
    nolengths <- strsplit(a[[1]][5],split='.',fixed=TRUE)[[1]][1]
    linker <- paste(a[[1]][3],a[[1]][4],nolengths,sep='_')
    }

#load in tree
tree <- read.tree(paste(direct,paste(gt,'newick',sep='.'),sep='/')) # load tree
#tree <-  di2multi(tree, 0.0001)

#find short branch lengths, delete them
tip.branches = match(1:length(tree$tip), tree$edge[,2])
bls = cbind.data.frame(taxa=tree$tip, tree$edge[tip.branches,], bl=tree$edge.length[tip.branches])
drop2 = bls[bls$bl<0.0001,]
drop2$taxa = as.character(drop2$taxa)

#get rownames from file
rows <- read.table(paste('KCO_Linker_pids/',gt,'_',linker,'_pids.txt',sep=''), header=FALSE)

#read in data, filter it
data <- read.table(args[1],header=FALSE)
data=data[,colSums(data)!=0]
data <- data[,-1]
rownames(data) <- t(rows)


#limit to 15 rows so that it's easy to read in stack overflow
treedata <- list(tree, data[1:15,])
sapply(treedata,class)
nc <- with(treedata,name.check(tree,data))
tree <- with(treedata,drop.tip(tree,nc$tree_not_data))
tree <- with(treedata,drop.tip(tree, drop2$taxa)) #delete short bl
nc <- with(treedata,name.check(tree,data))
data <- data[-which(rownames(data) %in% nc$data_not_tree),]

#write tree to a file
write.tree(tree,file='test_tree.newick')
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
     otd[varnames] <- sapply(otd[varnames], as.numeric)
     ### this data-frame now contains the data as well as the tree geometry
     ### now remake the ouch tree
     ot <- with(otd,ouchtree(nodes=nodes,ancestors=ancestors,times=times,labels=labels))

     #This will still give an error, but copying and pasting the tree written to the file in produces no error..
     b1 <- brown(tree=ot,data=otd[varnames])
     print("brown")
     print(summary(b1))

     ### evaluate an OU model with a single, global selective regime
     otd$regimes <- as.factor("global")
     h1 <- hansen(
                  tree=ot,
                  data=otd[varnames],
                  regimes=otd["regimes"],
                  sqrt.alpha=t(rep.int(1,lowtri)),#c(1, 1, 1, 1, 1, 3, 3, 3, 3, 6, 6, 6, 6, 10, 10, 10, 15, 15, 21), 
                  method="Nelder-Mead",
                  sigma=t(rep.int(1,lowtri)),
                  maxit=50000
                  )
     print("hansen")
     print(summary(h1))

     #wn <- fitContinuous(phy=tree,dat=data['totalNoGap'],model='white')
     #print("white_noise")
     #print(wn$opt)

     #sink()
     ## End(Not run)
}






