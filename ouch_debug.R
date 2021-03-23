library(ouch)
library(geiger)
library(phytools)

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
#tree = read.tree(text='((ENSCAFP00000019418:0.0735803,((((((ENSMAUP00000014127:0.0328701,ENSPEMP00000016993:0.02382)Cricetidae:0.0023201,ENSMOCP00000002175:0.04084)Cricetidae:0.0110901,((MGP_SPRETEiJ_P0059033:0.0066001,MGP_CAROLIEiJ_P0056799:0.00968)Mus:0.0228902,ENSRNOP00000073536:0.02747)Murinae:0.0229401)Muroidea:0.0542403,ENSDORP00000005057:0.08843)Rodentia:0.0077101,ENSCLAP00000005438:0.0839002)Rodentia:0.0383304,((ENSPCOP00000010808:0.01859,ENSMICP00000026445:0.02551)Lemuriformes:0.0343502,((((((ENSPPAP00000002002:0.0022401,ENSP00000266079:0.00405)Homininae:0.0013401,ENSGGOP00000016478:0.00271)Homininae:0.0080702,ENSNLEP00000032061:0.00933)Hominoidea:0.0066201,(((ENSRBIP00000040386:0.00057,ENSRROP00000034917:0.02862)Rhinopithecus:0.0047301,ENSCANP00000029315:0.00782)Colobinae:0.0049701,(((ENSCATP00000028349:0.00281,ENSPANP00000014391:0.00945)Cercopithecinae:0.0004201,ENSMLEP00000016334:0.04489)Cercopithecinae:0.0025902,((ENSMMUP00000038297:0.02311,ENSMNEP00000046191:0.10245)Macaca:0.0254401,ENSMFAP00000002296:0.00114)Macaca:0.0034501)Cercopithecinae:0.0059301)Cercopithecidae:0.0122701)Catarrhini:0.0125401,((ENSCCAP00000013808:0.0413301,ENSCJAP00000019620:0.02311)Cebidae:0.0010901,ENSANAP00000001861:0.01478)Platyrrhini:0.0200901)Simiiformes:0.0315901,ENSTSYP00000007234:0.04373)Haplorrhini:0.0094701)Primates:0.0060601)Euarchontoglires:0.0119501)Boreoeutheria:0.0888803,ENSGALP00000061667:0.2738806)Amniota;')
tree = read.tree(file='test_tree.newick')

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
treedata <- list(tree, data[1:15,])
sapply(treedata,class)
nc <- with(treedata,name.check(tree,data))
tree <- with(treedata,drop.tip(tree,nc$tree_not_data))
tree <- with(treedata,drop.tip(tree, drop2$taxa)) #delete short bl
nc <- with(treedata,name.check(tree,data))
data <- data[-which(rownames(data) %in% nc$data_not_tree),]

print(data)
#write.tree(tree,file='test_tree.newick')
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

     #setwd("/data1/home/dbarth/LinkerProject/len_residues_kco/") #changed for testing purposes, orig: /data1/home/dbarth/LinkerProject/ModelLengths/
     #sink(file=paste('len_model_fit',gt,paste(linker,'.txt',sep=''),sep='_'),append=FALSE)
     #print(otd[c('V1','V2','V3','V4','V5','V6')])
     #print(typeof(otd['V1']))
     b1 <- brown(tree=ot,data=otd[varnames])
     #print("brown")
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






