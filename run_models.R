library(geiger)
library(phytools)
library(mvMORPH)
library(parallel)
library(ggplot2)
nb_cores = 4L

args = commandArgs(trailingOnly=TRUE)
dirname = args[1]
seqname = args[2]
dirname = "LinkerLengthsKCO/"
seqname = "ENSGT00390000000018_0_Linker_5_6"
outdir = "model_results/"


bootstrap_mv_BM_OU <- function(aa, seqname, fit_mvBM, fit_mvOU, vcv="fixedRoot", opt_method="BFGS", nsim=100 ) {

  obs_lr = -2 * (logLik(fit_mvBM) - logLik(fit_mvOU) )
  bootstrapBMBM = vector("list", nsim*2)
  bootstrapBMOU = vector("list", nsim*2)
  bootstrapOUBM = vector("list", nsim*2)
  bootstrapOUOU = vector("list", nsim*2)

  rep = 0 
  total = 1
  repeat{
    simul_mvBM <- simulate(fit_mvBM, tree=tree, nsim=nsim) 
    tmpBMBM <- lapply(simul_mvBM, function(x) mvBM(tree, x, echo=F, diagnostic=F, model="BM1", scale.height=TRUE, param=list(constraint="equaldiagonal"))) 
    tmpBMOU <- lapply(simul_mvBM, function(x) mvOU(tree, x, echo=F, diagnostic=F, model="OU1", scale.height=TRUE, param=list(decomp="diagonal", decompSigma="equaldiagonal", vcv=vcv), optimization=opt_method, control = list(maxit = 100000, parscale=c(1000,10,1) ))) 
    goodsamples = !sapply(tmpBMBM, inherits, "try-error") & !sapply(tmpBMOU, inherits, "try-error")
    goodsamples = goodsamples & !sapply(tmpBMBM, "[[", "convergence") & !sapply(tmpBMBM, "[[", "hess.values") & !sapply(tmpBMOU, "[[", "convergence") & !sapply(tmpBMOU, "[[", "hess.values")
    if (sum(goodsamples)) {
      bootstrapBMBM[total:(total+sum(goodsamples)-1)] = tmpBMBM[goodsamples]
      bootstrapBMOU[total:(total+sum(goodsamples)-1)] = tmpBMOU[goodsamples]
    }
    total = total+sum(goodsamples)
    #bootstrapBMBM <- mclapply(1:nsim, function(x){ mvBM(tree, simul_mvBM[[x]], echo=F, diagnostic=F, model="BM1", scale.height=TRUE, param=list(constraint="equaldiagonal")) },   mc.cores = getOption("mc.cores", nb_cores))
    #bootstrapBMOU <- mclapply(1:nsim, function(x){ mvOU(tree, simul_mvBM[[x]], echo=F, diagnostic=F, model="OU1", scale.height=TRUE, param=list(decomp="diagonal", decompSigma="equaldiagonal", vcv=vcv), optimization=opt_method, control = list(maxit = 100000, parscale=c(1000,10,1) )) },   mc.cores = getOption("mc.cores", nb_cores))
    if( total >= nsim | rep > 10){
        break
    }
    rep = rep+1
  } 
  rep = 0 
  total = 1
  repeat{  
    simul_mvOU <- simulate(fit_mvOU, tree=tree, nsim=nsim) 
    tmpOUBM <- mclapply(1:nsim, function(x){ mvBM(tree, simul_mvOU[[x]], echo=F, diagnostic=F, model="BM1", scale.height=TRUE, param=list(constraint="equaldiagonal")) },   mc.cores = getOption("mc.cores", nb_cores))
    tmpOUOU <- mclapply(1:nsim, function(x){ mvOU(tree, simul_mvOU[[x]], echo=F, diagnostic=F, model="OU1", scale.height=TRUE, param=list(decomp="diagonal", decompSigma="equaldiagonal", vcv=vcv), optimization=opt_method, control = list(maxit = 100000, parscale=c(1000,10,1) )) },   mc.cores = getOption("mc.cores", nb_cores))
    goodsamples = !sapply(tmpOUBM, inherits, "try-error") & !sapply(tmpOUOU, inherits, "try-error")
    goodsamples = goodsamples & !sapply(tmpOUBM, "[[", "convergence") & !sapply(tmpOUBM, "[[", "hess.values") & !sapply(tmpOUOU, "[[", "convergence") & !sapply(tmpOUOU, "[[", "hess.values")
    if (sum(goodsamples)) {
      bootstrapOUBM[total:(total+sum(goodsamples)-1)] = tmpOUBM[goodsamples]
      bootstrapOUOU[total:(total+sum(goodsamples)-1)] = tmpOUOU[goodsamples]
    }
    total = total+sum(goodsamples)
    #bootstrapOUBM <- lapply(simul_mvOU, function(x) mvBM(tree, x, echo=F, diagnostic=F, model="BM1", scale.height=TRUE, param=list(constraint="equaldiagonal"))) 
    #bootstrapOUOU <- lapply(simul_mvOU, function(x) mvOU(tree, x, echo=F, diagnostic=F, model="OU1", scale.height=TRUE, param=list(decomp="diagonal", decompSigma="equaldiagonal", vcv=vcv), optimization=opt_method, control = list(maxit = 100000, parscale=c(1000,10,1) ))) 
    if( total >= nsim | rep > 10){
      break
    }
    rep = rep+1
  }

  nsim <- min (sum(!sapply(bootstrapBMBM, is.null)), sum(!sapply(bootstrapOUBM, is.null))) 
  if (nsim < 100) {
    print("bootstrap failed")
    return (NULL)
  }
  bootstrapBMBM <- bootstrapBMBM[1:nsim]
  bootstrapBMOU <- bootstrapBMOU[1:nsim]
  bootstrapOUBM <- bootstrapOUBM[1:nsim]
  bootstrapOUOU <- bootstrapOUOU[1:nsim]
  
  null_lr_dist = -2 * (sapply(bootstrapBMBM, logLik) - sapply(bootstrapBMOU, logLik))
  test_lr_dist = -2 * (sapply(bootstrapOUBM, logLik) - sapply(bootstrapOUOU, logLik))

  lrs <- data.frame(null=null_lr_dist, test=test_lr_dist)
  p1 <- ggplot(lrs)+
  geom_density(mapping=aes(x=null, fill='null'), alpha=0.7)+
  geom_density(mapping=aes(x=test, fill='test'), alpha=0.7)+
  geom_vline(xintercept=obs_lr, col="red") +
  xlab('Delta')+
  ylab('Density')
  ggsave(paste0(outdir, aa, '_', seqname, '.lr.bootstrap.pdf'),p1)

  sigma_hat = fit_mvOU$sigma[1,1]
  sigmas = sapply(bootstrapOUOU, "[[", "sigma")[1,]
  deltas = sigmas-sigma_hat
  d = quantile(deltas, c(0.05,0.95))
  ci_sigma = sigma_hat - c(d[2], d[1]) 
  sigma_hat = fit_mvOU$sigma[1,1]
  alphas = sapply(bootstrapOUOU, "[[", "alpha")[1,]
  deltas = alphas-alpha_hat
  d = quantile(deltas, c(0.05,0.95))
  ci_alpha = alpha_hat - c(d[2], d[1]) 
  alpha_total_hat = fit_mvOU$alpha[2,2]
  alpha_totals = sapply(bootstrapOUOU, "[[", "alpha")[4,]
  deltas = alpha_totals-alpha_total_hat
  d = quantile(deltas, c(0.05,0.95))
  ci_alpha_total = alpha_total_hat - c(d[2], d[1]) 
  pvalue = (1 + sum(null_lr_dist >= obs_lr)) / (nsim + 1)
  ret = c(ci_sigma, ci_alpha, ci_alpha_total, pvalue)
  names(ret) = c("sigma5%", "sigma95%","alpha5%", "alpha95%","alphatotal5%", "alphatotal95%","pvalue")
  return (ret)

}






if (grepl("Linker",seqname)) {
  treename = unlist(strsplit(seqname, "_Linker_"))[1]
} else {
  treename = unlist(strsplit(seqname, "_Domain_"))[1]
}

tree <- read.tree(paste0('Newicks/', treename, '.newick')) # load tree
cntdata <- read.table(paste0(dirname, '/',seqname, '.lengths.txt'),header=TRUE, row.names = 1)
treedata <- list(tree, cntdata)
sapply(treedata,class)
nc <- with(treedata,name.check(tree,cntdata))
tree <- with(treedata,drop.tip(tree,nc$tree_not_data))
cntdata <- cntdata[tree$tip.label,]
AAcodes <- colnames(cntdata)[4:23]
dat <- cntdata[,AAcodes]/cntdata$totalNoGap
if (nrow(dat) < 10) {
  print("not enough tips in the phylogeny")
  quit()
}

cnt <- dat[nrow(dat)+1:nrow(dat)*2,]
rownames(cnt) <- rownames(dat)
cnt[is.na(cnt)] <- 0
total <- data.frame(total=rep(0, nrow(dat)), row.names = rownames(dat))
gtlenname = paste0("GTLengthsKCO/",treename,".lengths.txt")
if (file.exists(gtlenname)) {
  gt <- read.table(gtlenname, header=TRUE)
  total <- gt$total
  cnt <- gt[,-1]
} else {
  # get all linkers and domains from same gene tree
  files <- list.files(path='LinkerLengthsKCO/', pattern=paste0( treename, "_*"), full.names=TRUE, recursive=FALSE)
  files <- c(files, list.files(path='DomainLengthsKCO/', pattern=paste0( treename, "_*"), full.names=TRUE, recursive=FALSE))
  for (f in files) {
    data <- read.table(f,header=TRUE, row.names = 1)
    data <- data[rownames(dat),]
    data[is.na(data)] <- 0
    cnt <- cnt + data[,AAcodes]
    total <- total + data$totalNoGap
  }
  gt <- cbind.data.frame(total, cnt)
  write.table(gt, gtlenname, sep="\t", quote=FALSE)
}
allprops <- cnt/total[,1]


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

  results = data.frame(matrix(nrow=1, ncol=158))
  r_i = 1
  # mvMORPH
  # uv
  fit_BM <- mvBM(tree, aa_dat[aa], model="BM1", scale.height=TRUE)
  fit_OUfixed <- mvOU(tree, aa_dat[aa], model="OU1", scale.height=TRUE, param=list(vcv="fixedRoot"))
  fit_OUrandom <- mvOU(tree, aa_dat[aa], model="OU1", scale.height=TRUE, param=list(vcv="randomRoot"))
  names(results)[r_i:(r_i+3)] = paste0("BM_", c("sigma", "theta", "lnL", "AICc"))
  results[r_i:(r_i+3)] = c(fit_BM$sigma, fit_BM$theta, fit_BM$LogLik, fit_BM$AICc)
  r_i = r_i+4
  names(results)[r_i:(r_i+6)] = paste0("OUr_", c("sigma", "alpha", "theta", "lnL", "AICc", "conv", "hessian"))
  results[r_i:(r_i+6)] = c(fit_OUrandom$sigma, fit_OUrandom$alpha, fit_OUrandom$theta, fit_OUrandom$LogLik, fit_OUrandom$AICc, fit_OUrandom$convergence, fit_OUrandom$hess.values)
  r_i = r_i+7
  names(results)[r_i:(r_i+7)] = paste0("OUf_", c("sigma", "alpha", "theta0", "theta1", "lnL", "AICc", "conv", "hessian"))
  results[r_i:(r_i+7)] = c(fit_OUfixed$sigma, fit_OUfixed$alpha, fit_OUfixed$theta, fit_OUfixed$LogLik, fit_OUfixed$AICc, fit_OUfixed$convergence, fit_OUfixed$hess.values)
  r_i = r_i+8

  # mv
  fit_mvBM <- mvBM(tree, mv_dat, model="BM1", scale.height=TRUE, param=list(constraint="equaldiagonal"))
  names(results)[r_i:(r_i+5)] = paste0("mvBM_", c("sigma0", "sigma1", "theta0", "theta1", "lnL", "AICc"))
  results[r_i:(r_i+5)] = c(fit_mvBM$sigma[1,], fit_mvBM$theta, fit_mvBM$LogLik, fit_mvBM$AICc)
  r_i = r_i+6
  fit_mvOU_list <- vector("list", 8)
  try({fit_mvOU_list[[1]] <- mvOU(tree, mv_dat, model="OU1", scale.height=TRUE, param=list(decomp="diagonal", decompSigma="equaldiagonal", vcv="fixedRoot"), optimization="BFGS", control = list(maxit = 100000, parscale=c(1000,10,1)) , echo=F)})
  try({fit_mvOU_list[[2]] <- mvOU(tree, mv_dat, model="OU1", scale.height=TRUE, param=list(decomp="diagonal", decompSigma="equaldiagonal", vcv="fixedRoot"), optimization="SANN", control = list(maxit = 100000, parscale=c(1000,10,1)) , echo=F)})
  try({fit_mvOU_list[[3]] <- mvOU(tree, mv_dat, model="OU1", scale.height=TRUE, param=list(decomp="diagonal", decompSigma="equaldiagonal", sigma=fit_mvBM$sigma[1,1], vcv="fixedRoot"), optimization="BFGS", control = list(maxit = 100000, parscale=c(1000,10,1)) , echo=F)})
  try({fit_mvOU_list[[4]] <- mvOU(tree, mv_dat, model="OU1", scale.height=TRUE, param=list(decomp="diagonal", decompSigma="equaldiagonal", sigma=fit_mvBM$sigma[1,1], vcv="fixedRoot"), optimization="SANN", control = list(maxit = 100000, parscale=c(1000,10,1)) , echo=F)})
  try({fit_mvOU_list[[5]] <- mvOU(tree, mv_dat, model="OU1", scale.height=TRUE, param=list(decomp="diagonal", decompSigma="equaldiagonal", vcv="randomRoot"), optimization="BFGS", control = list(maxit = 100000, parscale=c(1000,10,1)) , echo=F)})
  try({fit_mvOU_list[[6]] <- mvOU(tree, mv_dat, model="OU1", scale.height=TRUE, param=list(decomp="diagonal", decompSigma="equaldiagonal", vcv="randomRoot"), optimization="SANN", control = list(maxit = 100000, parscale=c(1000,10,1)) , echo=F)})
  try({fit_mvOU_list[[7]] <- mvOU(tree, mv_dat, model="OU1", scale.height=TRUE, param=list(decomp="diagonal", decompSigma="equaldiagonal", sigma=fit_mvBM$sigma[1,1], vcv="randomRoot"), optimization="BFGS", control = list(maxit = 100000, parscale=c(1000,10,1)) , echo=F)})
  try({fit_mvOU_list[[8]] <- mvOU(tree, mv_dat, model="OU1", scale.height=TRUE, param=list(decomp="diagonal", decompSigma="equaldiagonal", sigma=fit_mvBM$sigma[1,1], vcv="randomRoot"), optimization="SANN", control = list(maxit = 100000, parscale=c(1000,10,1)) , echo=F)})
 
  sigma <- t(sapply(fit_mvOU_list, "[[", "sigma")[1:2,])
  colnames(sigma) <- c("sigma", "sigmacov")
  alpha <- t(sapply(fit_mvOU_list, "[[", "alpha")[c(1,4),])
  colnames(alpha) <- c("alpha", "alphatotal")
  theta <- sapply(fit_mvOU_list, "[[", "theta")
  theta_f =  t(matrix(unlist(theta[1:4]), nrow=4))
  theta_r =  t(matrix(unlist(theta[5:8]), nrow=2))
  theta_r <- cbind(matrix(rep(0,4), nrow=4), theta_r[,1], matrix(rep(0,4),nrow=4), theta_r[,2] )
  theta <- rbind(theta_f, theta_r)
  colnames(theta) <- c("theta0", "theta1", "thetatotal0", "thetatotal1")
  lnL <- unlist(sapply(fit_mvOU_list, "[[", "LogLik"))
  AICc <- unlist(sapply(fit_mvOU_list, "[[", "AICc"))
  conv <- unlist(sapply(fit_mvOU_list, "[[", "convergence"))
  hess <- unlist(sapply(fit_mvOU_list, "[[", "hess.values"))
  params <- sapply(fit_mvOU_list, "[[", "param")
  opt <- unlist(params["optimization",])
  root <- unlist(params["root",])
  mvOUresults <- cbind.data.frame(root, opt, sigma, alpha, theta, lnL, AICc, conv, hess)
  mvOUvector <- unlist(mvOUresults)
  names(results)[r_i:(r_i+111)] = names(mvOUvector)
  results[r_i:(r_i+111)] = mvOUvector
  r_i = r_i+112

  goodfitidx <- !conv & !hess
  if (any(goodfitidx)) {
    minAICc <- min(mvOUresults[goodfitidx,"AICc"])
    if (minAICc < fit_mvBM$AICc) {
      bestfitidx <- which(goodfitidx & mvOUresults[,"AICc"]==minAICc)
      fit_mvOU <- fit_mvOU_list[[bestfitidx]]
      bestresults = mvOUresults[bestfitidx,]
      names(bestresults) = paste0("best_", names(bestresults))
      names(results)[r_i:(r_i+13)] = names(bestresults)
      results[r_i:(r_i+13)] = bestresults
      r_i = r_i+14
  
      if (bestresults$best_root) {
        vcv = "fixedRoot"
      } else {
        vcv = "randomRoot"
      }
      ret = bootstrap_mv_BM_OU(aa, seqname, fit_mvBM, fit_mvOU, vcv)
      if (ret) {
        names(results)[r_i:(r_i+length(ret)-1)] = names(ret)
        results[r_i:(r_i+length(ret)-1)] = ret
      }
    }
  }
  write.table(results, paste0(outdir, aa, '_', seqname, ".txt"), sep="\t", quote=FALSE, row.names=FALSE)
  
}


