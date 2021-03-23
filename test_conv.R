library(geiger)
library(phytools)
library(mvMORPH)
library(parallel)
nb_cores = 4L

args = commandArgs(trailingOnly=TRUE)
dirname = args[1]
seqname = args[2]
dirname = 'LinkerLengthsKCO/'
seqname = 'ENSGT00390000001446_0_Linker_2_3'


bootstrap_mv_BM_OU <- function(fit_mvBM, fit_mvOU, vcv="fixedRoot", nsim=100 ) {

  obs_diff = -2 * (logLik(fit_mvBM) - logLik(fit_mvOU) )
  
  simul_mvBM <- simulate(fit_mvBM, tree=tree, nsim=nsim) 
  simul_mvOU <- simulate(fit_mvOU, tree=tree, nsim=nsim) 

  #bootstrapBMBM <- lapply(simul_mvBM, function(x) mvBM(tree, x, echo=F, diagnostic=F, param=list(constraint="equal"))) 
  #bootstrapBMOU <- lapply(simul_mvBM, function(x) mvOU(tree, x, echo=F, diagnostic=F, param=list(constraint="equal", vcv=vcv))) 
  #bootstrapOUBM <- lapply(simul_mvOU, function(x) mvBM(tree, x, echo=F, diagnostic=F, param=list(constraint="equal"))) 
  #bootstrapOUOU <- lapply(simul_mvOU, function(x) mvOU(tree, x, echo=F, diagnostic=F, param=list(constraint="equal", vcv=vcv))) 

  bootstrapBMBM <- mclapply(1:nsim, function(x){ mvBM(tree, simul_mvBM[[x]], echo=F, diagnostic=F, param=list(constraint="equal")) },   mc.cores = getOption("mc.cores", nb_cores))
  bootstrapBMOU <- mclapply(1:nsim, function(x){ mvOU(tree, simul_mvBM[[x]], echo=F, diagnostic=F, param=list(constraint="equal", vcv=vcv), optimization="subplex" ) },   mc.cores = getOption("mc.cores", nb_cores))
  bootstrapOUBM <- mclapply(1:nsim, function(x){ mvBM(tree, simul_mvOU[[x]], echo=F, diagnostic=F, param=list(constraint="equal")) },   mc.cores = getOption("mc.cores", nb_cores))
  bootstrapOUOU <- mclapply(1:nsim, function(x){ mvOU(tree, simul_mvOU[[x]], echo=F, diagnostic=F, param=list(constraint="equal", vcv=vcv), optimization="subplex" ) },   mc.cores = getOption("mc.cores", nb_cores))
  null_dist = -2 * (sapply(bootstrapBMBM, logLik) - sapply(bootstrapBMOU, logLik))
  test_dist = -2 * (sapply(bootstrapOUBM, logLik) - sapply(bootstrapOUOU, logLik))
 
  lrs <- data.frame(null=null_dist, test=test_dist)
  p1 <- ggplot(lrs)+
  geom_density(mapping=aes(x=null, fill='null'))+
  geom_density(mapping=aes(x=test, fill='test'))+
  geom_vline(xintercept=obs_diff, col="red") +
  xlab('Delta')+
  ylab('Density')
  ggsave('bootstrap.pdf',p1)
 
}

if (grepl("Linker",seqname)) {
  treename = unlist(strsplit(seqname, "_Linker_"))[1]
} else {
  treename = unlist(strsplit(seqname, "_Domain_"))[1]
}

tree <- read.tree(paste0('Newicks/',treename, '.newick')) # load tree
cntdata <- read.table(paste0(dirname, '/',seqname, '.lengths.txt'),header=TRUE, row.names = 1)
treedata <- list(tree, cntdata)
sapply(treedata,class)
nc <- with(treedata,name.check(tree,cntdata))
tree <- with(treedata,drop.tip(tree,nc$tree_not_data))
cntdata <- cntdata[tree$tip.label,]
AAcodes <- colnames(cntdata)[4:23]
dat <- cntdata[,AAcodes]/cntdata$totalNoGap

cnt <- dat[nrow(dat)+1:nrow(dat)*2,]
rownames(cnt) <- rownames(dat)
cnt[is.na(cnt)] <- 0
total <- data.frame(total=rep(0, nrow(dat)), row.names = rownames(dat))
files <- list.files(path='LinkerLengthsKCO/', pattern=paste0( treename, "_*"), full.names=TRUE, recursive=FALSE)
files <- c(files, list.files(path='DomainLengthsKCO/', pattern=paste0( treename, "_*"), full.names=TRUE, recursive=FALSE))
for (f in files) {
  data <- read.table(f,header=TRUE, row.names = 1)
  data <- data[rownames(dat),]
  data[is.na(data)] <- 0
  cnt <- cnt + data[,AAcodes]
  total <- total + data$totalNoGap
}
allprops <- cnt/total$total

asinTransform <- function(p) { asin(sqrt(p)) }
dat <- asinTransform(dat)
allprops <- asinTransform(allprops)


for (aa in AAcodes) {

  print(aa)
  aa_dat <- dat[aa]
  if (sum(aa_dat) == 0) {
    next;
  }
  mv_dat <- cbind.data.frame(dat[aa], allprops[aa]) 
  colnames(mv_dat) <- c("seq", "total")
  print(mv_dat)

#
#  ### make an ouchtree out of the phy-format tree
#  ot <- ape2ouch(tree)
#  ### merge data with tree info
#  otd <- as(ot,"data.frame")
#  ### in these data, it so happens that the rownames correspond to node names
#  ### we will exploit this correspondence in the 'merge' operation:
#  aa_dat$labels <- rownames(aa_dat)
#  m <- merge(otd,aa_dat,by="labels",all=TRUE)
#  otd <- m[order(as.numeric(as.character(m$nodes))+0),]
#  rownames(otd) <- otd$nodes
#  print(otd)
#  ### this data-frame now contains the data as well as the tree geometry
#
#  ### now remake the ouch tree
#  ot <- with(otd,ouchtree(nodes=nodes,ancestors=ancestors,times=times,labels=labels))
#  otd$regimes <- as.factor("global")
#
#  bm <- brown(tree=ot, data=otd[c(aa)])
#  c(summary(bm)$loglik, summary(bm)$aic)
#  ou_NM <- hansen(tree=ot, data=otd[c(aa)], regimes = otd["regimes"], sqrt.alpha = 1, sigma = c(1), method = "Nelder-Mead", maxit = 200000)
#  c(summary(ou_NM)$loglik, summary(ou_NM)$aic, summary(ou_NM)$conv.code
#  ou_BFGS <- hansen(tree=ot, data=otd[c(aa)], regimes = otd["regimes"], sqrt.alpha = 1, sigma = c(1), method = "BFGS", maxit = 200000)
#  c(summary(ou_BFGS)$loglik, summary(ou_BFGS)$aic, summary(ou_BFGS)$conv.code)
#  ou_L_BFGS_B <- hansen(tree=ot, data=otd[c(aa)], regimes = otd["regimes"], sqrt.alpha = 1, sigma = c(1), method = "L-BFGS-B", maxit = 200000)
#  c(summary(ou_L_BFGS_B)$loglik, summary(ou_L_BFGS_B)$aic, summary(ou_L_BFGS_B)$conv.code)
#  ou_NM <- hansen(tree=ot, data=otd[c(aa)], regimes = otd["regimes"], sqrt.alpha = 1, sigma = coef(bm)$sigma.sq.matrix, method = "Nelder-Mead", maxit = 200000)
#  c(summary(ou_NM)$loglik, summary(ou_NM)$aic, summary(ou_NM)$conv.code)
#  ou_BFGS <- hansen(tree=ot, data=otd[c(aa)], regimes = otd["regimes"], sqrt.alpha = 1, sigma = coef(bm)$sigma.sq.matrix, method = "BFGS", maxit = 200000)
#  c(summary(ou_BFGS)$loglik, summary(ou_BFGS)$aic, summary(ou_BFGS)$conv.code)
#  ou_L_BFGS_B <- hansen(tree=ot, data=otd[c(aa)], regimes = otd["regimes"], sqrt.alpha = 1, sigma = coef(bm)$sigma.sq.matrix, method = "L-BFGS-B", maxit = 200000, parscale = c(1000, 1), factr=0)
#  c(summary(ou_L_BFGS_B)$loglik, summary(ou_L_BFGS_B)$aic, summary(ou_L_BFGS_B)$conv.code)
#

  # geiger
  #fitBM<-fitContinuous(tree,aa_dat[aa])
  #fitBM
  #fitOU<-fitContinuous(tree,aa_dat[aa], model="OU", control=list(niter=1000))
  #fitOU

  # mvMORPH
  fit_BM <- mvBM(tree, aa_dat[aa], model="BM1", scale.height=TRUE)
  fit_OUfixed <- mvOU(tree, aa_dat[aa], model="OU1", scale.height=TRUE, param=list(vcv="fixedRoot"))
  fit_OUrandom <- mvOU(tree, aa_dat[aa], model="OU1", scale.height=TRUE, param=list(vcv="randomRoot"))

  fit_mvBM <- mvBM(tree, mv_dat, model="BM1", scale.height=TRUE, method="pic", param=list(constraint="equal"))
  fit_mvOUfixed <- vector("list", 4)
  fit_mvOUfixed[[1]] <- mvOU(tree, mv_dat, model="OU1", scale.height=TRUE, param=list(decomp="diagonal", decompSigma="equal", vcv="fixedRoot"), optimization="BFGS", control = list(maxit = 100000, parscale=c(1000,10,1,1)) , echo=F)
  fit_mvOUfixed[[2]] <- mvOU(tree, mv_dat, model="OU1", scale.height=TRUE, param=list(decomp="diagonal", decompSigma="equal", vcv="fixedRoot"), optimization="SANN", control = list(maxit = 100000, parscale=c(1000,10,1,1)) , echo=F)
  fit_mvOUfixed[[3]] <- mvOU(tree, mv_dat, model="OU1", scale.height=TRUE, param=list(decomp="diagonal", decompSigma="equal", sigma=fit_mvBM$sigma[1,], vcv="fixedRoot"), optimization="BFGS", control = list(maxit = 100000, parscale=c(1000,10,1,1)) , echo=F)
  fit_mvOUfixed[[4]] <- mvOU(tree, mv_dat, model="OU1", scale.height=TRUE, param=list(decomp="diagonal", decompSigma="equal", sigma=fit_mvBM$sigma[1,], vcv="fixedRoot"), optimization="SANN", control = list(maxit = 100000, parscale=c(1000,10,1,1)) , echo=F)
  conv <- unlist(sapply(fit_mvOUfixed, "[[", "convergence"))
  hess <- unlist(sapply(fit_mvOUfixed, "[[", "hess.values"))
  goodlist <- fit_mvOUfixed[!conv & !hess]
  if (length(goodlist)) {
    goodAICc <- unlist(sapply(goodlist, "[[", "AICc")) 
    bestidx <- which (goodAICc == min(goodAICc))
    fit <- goodlist[[bestidx]]
    print(fit)
    print(fit$param$optimization)
  }

  fit_mvOUrandom <- vector("list", 4)
  fit_mvOUrandom[[1]] <- mvOU(tree, mv_dat, model="OU1", scale.height=TRUE, param=list(decomp="diagonal", decompSigma="equal", vcv="randomRoot"), optimization="BFGS", control = list(maxit = 100000, parscale=c(1000,10,1,1)) , echo=F)
  fit_mvOUrandom[[2]] <- mvOU(tree, mv_dat, model="OU1", scale.height=TRUE, param=list(decomp="diagonal", decompSigma="equal", vcv="randomRoot"), optimization="SANN", control = list(maxit = 100000, parscale=c(1000,10,1,1)) , echo=F)
  fit_mvOUrandom[[3]] <- mvOU(tree, mv_dat, model="OU1", scale.height=TRUE, param=list(decomp="diagonal", decompSigma="equal", sigma=fit_mvBM$sigma[1,], vcv="randomRoot"), optimization="BFGS", control = list(maxit = 100000, parscale=c(1000,10,1,1)) , echo=F)
  fit_mvOUrandom[[4]] <- mvOU(tree, mv_dat, model="OU1", scale.height=TRUE, param=list(decomp="diagonal", decompSigma="equal", sigma=fit_mvBM$sigma[1,], vcv="randomRoot"), optimization="SANN", control = list(maxit = 100000, parscale=c(1000,10,1,1)) , echo=F)
  conv <- unlist(sapply(fit_mvOUrandom, "[[", "convergence"))
  hess <- unlist(sapply(fit_mvOUrandom, "[[", "hess.values"))
  goodlist <- fit_mvOUrandom[!conv & !hess]
  if (length(goodlist)) {
    goodAICc <- unlist(sapply(goodlist, "[[", "AICc")) 
    bestidx <- which (goodAICc == min(goodAICc))
    fit <- goodlist[[bestidx]]
    print(fit)
    print(fit$param$optimization)
  }

  #bootstrap_mv_BM_OU(fit_mvBM, fit_mvOUfixed)
  
}


