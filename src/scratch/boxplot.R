result_dir = '~/Documents/workspace/phospho_network/script_files/analysis_results'
one_file   = 'matrix_1up.csv'
aup_file   = 'matrix_aup.csv'


aup <- read.csv(paste(result_dir,aup_file,sep = '/'),header = F,as.is = T)
one <- read.csv(paste(result_dir,one_file,sep = '/'),header = F,as.is = T)

par(mfrow=c(1,2))

one2 <- one[-1,-1]
colnames(one2) <- c(seq(0,20,1),1000,Inf)
one3 <- apply(one2,2,as.numeric)
rownames(one3) <- one[-1,1]

aup2 <- aup[-1,-1]
colnames(aup2) <- c(seq(0,20,1),1000,Inf)
aup3 <- apply(aup2,2,as.numeric)
rownames(aup3) <- aup[-1,1]

# boxplot(one3,xlab='penalty factor',ylab='overall q2 for 4 best prediction',main = 'one level up',las=2)
# boxplot(aup3,xlab='penalty factor',ylab='overall q2 for 4 best prediction',main = 'all possible upstreams',las=2)

boxplot(one3[apply(one3,1,mean)>0,],xlab='penalty factor',ylab='overall q2 for 4 best prediction',main = 'one level up',las=2,ylim=c(-0.3,0.8))
boxplot(aup3[apply(aup3,1,mean)>0,],xlab='penalty factor',ylab='overall q2 for 4 best prediction',main = 'all possible upstreams',las=2,ylim=c(-0.3,0.8))
