library(tidyverse)
library(effsize)

data <- read_delim('total_AAC_data.tsv', delim='\t')
data <- mutate(data, 'small' = length<6, 'medium' = (length>=6)&(length<=14), 'large' = (length>14)&(length<=30))

sizes <- c('small','medium', 'large')
varnames <- c('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y')


for (size in sizes) {
    cat('\n ------------------------- \n')
    print(size)
    test <- data.frame(matrix(ncol=20, nrow=2)) 
    colnames(test) <- varnames
    rownames(test) <- c('p-value', 'weighed.f.effect')
  
    for (var in varnames){
      l_filtered <- as.numeric(unlist(filter(data, get(size)==TRUE, linker==TRUE)[,var]))
      d_filtered <- as.numeric(unlist(filter(data, get(size)==TRUE, linker==FALSE)[,var]))
      wtest <- wilcox.test(l_filtered, d_filtered)
      test['p-value', var] <- wtest$p.value

      #calculate cohen's d, first transform to z
      U <- wtest$statistic
      n1 <- as.double(length(l_filtered))
      n2 <- as.double(length(d_filtered))
      #mu <- n1*n2/2
      #sig <- sqrt(n1*n2*(n1+n2+1)/12)
      #z <- (U-mu)/sig
      #test['cohens-d', var] <- z/sqrt(n1*n2)

      #caluclate f-effect
      test['weighed.f.effect',var] <- (U/(n1*n2))*sign(mean(l_filtered)-mean(d_filtered))
      
      #create plot
      #plt <- ggplot() +
      #      geom_violin(data=data[data[,size]==TRUE,], mapping=aes(x=linker, y=get(var)))+
      #      geom_line(mapping=aes(x=c(1,2), y=c(1.05,1.05)))+
      #      geom_text(mapping=aes(x=1.5, y=1.1, label=paste('f=',round(test['f-effect',var],4))))+
      #      scale_x_discrete(labels=c('Domain','Linker'))+
      #      xlab('')+
      #      ylab(paste('proportion of ', var))
      #      ggtitle(paste('Distributions of ', size, ' size'))
      #ggsave(paste(size,'_',var,'_residue_violins.png'), plt, path='linkerstatsfigs')
     
    }
  test <- as_tibble(t(test))
  test['residue'] <- varnames
  print(arrange(test, weighed.f.effect))
  #print(arrange(test, weighed.f.effect))
  cat('\n ------------------------- ')
}

