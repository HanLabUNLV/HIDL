library(phytools)
library(geiger)


args = commandArgs(trailingOnly=TRUE)
linkerdirname = args[1]
domaindirname = args[2]
linkerD2ndirname = args[3]
domainD2ndirname = args[4]

getDist <- function(chardir, D2ndir, seqname) {

  maxbl = NULL
  maxD2n = NULL
  if (grepl("Linker",seqname)) {
    treename = unlist(strsplit(seqname, "_Linker_"))[1]
  } else {
    treename = unlist(strsplit(seqname, "_Domain_"))[1]
  }

  print(seqname)
  tree1 <- read.tree(paste0('Newicks/', treename, '.newick')) # load tree
  data1 <- read.table(paste0(chardir, '/', seqname, '.lengths.txt'),header=TRUE)

  dat <- data.frame(length=data1$totalNoGap, row.names=data1$seqId)

  treedata <- list(tree1, dat)
  sapply(treedata,class)
  nc <- with(treedata,name.check(tree1,dat))
  if ( nc == "OK") { 
    tree <- tree1
  }
  else {
    tree <- with(treedata,drop.tip(tree1,nc$tree_not_data))
  }
  if (is.null(tree)) {
    return (NULL)
  }
  if (length(tree$tip.label) < 2) {
    return (NULL)
  }


  D <- cophenetic(tree)
  D[!lower.tri(D)] <- 0
  
  maxbl <- max(D) # maximum distance between nodes = 2* tree height
  sumbl <- sum(D) # maximum distance between nodes = 2* tree height

  jD2name = paste0(D2ndir,'/', seqname, ".matrix")
  D2n = NULL
  maxD2n = NULL 
  sumD2n = NULL
  if (file.exists(jD2name)) {
      colnum <- as.numeric(readLines(jD2name, n=1))
      print(colnum)
      if (is.na(colnum)) {
        return (NULL)
      }
      if (colnum < length(tree$tip.label)) {
        return (NULL)
      }
      jD2 <- read.table(jD2name, skip = 1, row.names=1, fill=TRUE, col.names = paste0("V",seq_len(colnum+1))) 
      colnames(jD2) <- rownames(jD2)
      jD2 <- jD2[rownames(jD2) %in% tree$tip.label,]
      jD2 <- jD2[,colnames(jD2) %in% tree$tip.label]
      limit <- abs(jD2-11.5129254649702)<1e-10
      rmidx <- apply(limit, 1, any, na.rm = FALSE) & apply(limit, 2, any, na.rm = FALSE)
      jD2 <- jD2[is.na(rmidx),]
      jD2 <- jD2[,is.na(rmidx)]
      D2n <- as(jD2, "matrix")
      maxD2n <-max(D2n, na.rm=TRUE)
      sumD2n <-sum(D2n, na.rm=TRUE)
  }

  D <- D[rownames(D) %in% rownames(D2n),]
  D <- D[,colnames(D) %in% colnames(D2n)]
  return( c(maxbl, sumbl, maxD2n, sumD2n, sumD2n/sum(D)))

}


files <- list.files(path=linkerdirname, pattern="*.txt", full.names=TRUE, recursive=FALSE)
linkerlist <- sub(".lengths.txt", "", basename(files))
result.linker = as.data.frame(matrix(ncol=5, nrow=length(linkerlist)))
names(result.linker) = c("maxbl", "totalbl", "maxD2n", "totalD2n", "ratio")
rownames(result.linker) = linkerlist
for (i in 1:length(linkerlist)) {
  ret <-  getDist(chardir= linkerdirname, D2ndir=linkerD2ndirname, seqname=linkerlist[i])  
  if (! is.null(ret)) {
    result.linker[i,1:length(ret)] =  ret
  }
}
write.table(result.linker, "distance.linker.bl.D2n.txt", sep="\t", quote=FALSE)

files <- list.files(path=domaindirname, pattern="*.txt", full.names=TRUE, recursive=FALSE)
domainlist <- sub(".lengths.txt", "", basename(files))
result.domain = as.data.frame(matrix(ncol=5, nrow=length(domainlist)))
names(result.domain) = c("maxbl", "totalbl", "maxD2n", "totalD2n", "ratio")
rownames(result.domain) = domainlist
for (i in 1:length(domainlist)) {
  ret <-  getDist(chardir= domaindirname, D2ndir=domainD2ndirname, seqname=domainlist[i])  
  if (! is.null(ret)) {
    result.domain[i,] =  ret
  }
}
write.table(result.domain, "distance.domain.bl.D2n.txt", sep="\t", quote=FALSE)



linker <- result.linker
linker <- linker[!is.na(linker$maxbl),]
linker <- linker[linker$maxbl > 0.00001,]
linker <- linker[!is.na(linker$linkermaxD2n) & linker$linkermaxD2n < 11.51292,]
linker <- linker[order(linker$linkermaxD2n),]
domain <- result.domain
domain <- domain[!is.na(domain$maxbl),]
domain <- domain[domain$maxbl > 0.00001,]
domain <- domain[!is.na(domain$domainmaxD2n) & domain$domainmaxD2n < 11.51292,]

pdf("distance.pdf")
plot(linker[,c(1,3)])
plot(linker[,c(2,4)])
hist(linker[,5])
boxplot(linker[,c(1,3,5)])
boxplot(linker[,c(2,4)])
plot(domain[c(1,3)]
plot(domain[,c(1,3)]
plot(domain[,c(1,3)])
plot(domain[,c(2,4)])
hist(domain[,5])
boxplot(domain[,c(1,3,5)])
boxplot(domain[,c(2,4)])
dev.off()



ou <- read.table("AAC_uv_results_linker.txt", header=TRUE)
ou_A <- ou[ou[,1]=="A",]
rownames(ou_A) <- ou_A$genetree_linker
ou_A_linker <- ou_A[rownames(linker),]
ou_A_linker <- cbind.data.frame(linker, ou_A_linker)
ou_A_linker <- ou_A_linker[!is.na(ou_A_linker$residue),]
