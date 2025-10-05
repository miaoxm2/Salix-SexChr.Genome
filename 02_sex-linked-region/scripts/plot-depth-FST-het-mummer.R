library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpattern)
library("ggpubr")
library("ggrepel")

## The multiplot function is currently defined as
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

pai <-read.csv("/0Salix_Project/1Project_genome/5-call_variants/0final/3-ind/h1_swi_7f7m.chr15.f-m.windowed.pi",sep = " ",header = FALSE,col.names = c("coord","f","m"))
pai <-read.csv("h1_swi_7f7m.chr15.f-m.1k.windowed.pi",sep = " ",header = FALSE,col.names = c("coord","f","m"))
fst<- read.csv("/0Salix_Project/1Project_genome/5-call_variants/0final/3-ind/h1_swi_7f7m.chr15.f-m.windowed.weir.fst",sep = "\t")
#unique and filtered coord
coord<-read.csv("../../../10-assembly/1hic/mummer/final.chr15.h1-h2.uniq.coords.plot",sep = "\t")
#raw coord.not in use
#coord<-read.csv("../../../4-mummer_result/02-version.final/h1-h2.chr15.mum.coords.iden",,sep=" ",header = F)
#colnames(coord)<-c("start.h1","end.h1","start.h2","end.h2","iden")
head(coord)
faf<- read.csv("allelefq7f.chr15.frq",sep = " ",col.names = c("chr","coord","a1","f1","a2","f2"))
maf<- read.csv("allelefq7m.chr15.frq",sep = " ",col.names = c("chr","coord","a1","f1","a2","f2"))
fmaf<- read.csv("allelefq7f7m.chr15.frq",sep = " ",col.names = c("chr","coord","fa","ff","ma","mf"))

#2024-01-18 new round###
#1-fst####

#read files
setwd("/0Salix_Project/1Project_genome/5-call_variants/0final/3-ind/fst/")

#fst for h1. whole genome or chr15
#window size=10kb, sliding window=10kb
fst.h1 <- read.csv("h1mask_She1-20ind.f-m.10k.windowed.weir.fst",sep = "\t")
#window size=10kb, sliding window=1kb
fst.h1 <- read.csv("h1mask_She1-20ind.f-m.10k2.windowed.weir.fst",sep = "\t")

#fst for h2. whole genome or chr15
#window size=10kb, sliding window=10kb
fst.h2 <- read.csv("h2mask_She1-20ind.f-m.10k.windowed.weir.fst",sep = "\t")
#window size=10kb, sliding window=1kb
fst.h2 <- read.csv("h2mask_She1-20ind.f-m.10k2.windowed.weir.fst",sep = "\t")
head(fst.h2)

fst.h1<-fst.h1%>%mutate(FST=case_when(
  WEIGHTED_FST<=0 ~ 0,
  WEIGHTED_FST>0 ~ WEIGHTED_FST
))%>%filter(grepl("chr",CHROM))%>%
  mutate(CHR=as.numeric(sub("[a-z]1_chr","",CHROM)))

fst.h2<-fst.h2%>%mutate(FST=case_when(
  WEIGHTED_FST<=0 ~ 0,
  WEIGHTED_FST>0 ~ WEIGHTED_FST
))%>%filter(grepl("chr",CHROM))%>%
  mutate(CHR=as.numeric(sub("[a-z]2_chr","",CHROM)))

#1-1 plot fst for whole genome####
length(unique(fst.h1$CHROM))

#h1
line=quantile(fst.h1$FST,0.99)

g_list = list()
i=1
for (id in sort(unique(fst.h1$CHROM))) {
 # line<-quantile(fst.h1[which(fst.h1$CHROM==id),]$FST,0.95)
  plot<-ggplot(fst.h1[which(fst.h1$CHROM==id),])+geom_hline(yintercept = line,colour="red",linetype = 'dashed')+
    geom_line(aes(x=BIN_START,y=FST))+
    ggtitle(id)+
    xlab("")+labs(y = (~F[ST]~  FvsM))+scale_x_continuous(breaks =seq(0,max(fst.h1[which(fst.h1$CHROM==id),2]),2e6),labels =seq(0,max(fst.h1[which(fst.h1$CHROM==id),2])/1e6,2))+
    scale_y_continuous(limits = c(-0.1,1))+theme_bw()+
    theme(axis.title = element_text(size=30),axis.text = element_text(size = 28),plot.title= element_text(size=30))
  print(plot)
  g_list[[i]] = plot
  i=i+1
}
pdf('h1-She1.FM-Fst.line2.pdf', height =30, width = 60)
multiplot(g_list[[1]],g_list[[2]], g_list[[3]],g_list[[4]],g_list[[5]],g_list[[6]],g_list[[7]],g_list[[8]],g_list[[9]],g_list[[10]],g_list[[11]],g_list[[12]],g_list[[13]],g_list[[14]],g_list[[15]],g_list[[16]],g_list[[17]],g_list[[18]],g_list[[19]], cols=5)
dev.off()

install.packages("qqman")
install.packages("Cairo")
library(qqman)

fst.h1$CHR<-as.numeric(sub("[a-z]1_chr","",fst.h1$CHROM))
fst.h1$Fst<-fst.h1$WEIGHTED_FST
fst.h1$POS<-fst.h1$BIN_START
fst.h1$SNP<-fst.h1$BIN_START
fst.h1%>%select("SNP","CHR","POS","Fst")



colorset<-c("#003f5c","#58508d","#bc5090","#ff6361","#ffa600")
pdf("h1.fst.manhattan.pdf",width = 16.8,height = 8)
manhattan(fst.h1%>%filter(grepl("chr",CHROM))%>%select("SNP","CHR","POS","Fst"),chr="CHR",bp="POS",p="Fst",snp="SNP",col=colorset,logp=F,suggestiveline = F,genomewideline = F,ylab="Fst",ylim=c(0,1),font.lab=10,main="S.herbacea Female vs Male",cex.lab=2, cex.axis=2, cex.main=4, cex.sub=2)
dev.off() 

#h2
fst.h2$CHR<-as.numeric(sub("[a-z]2_chr","",fst.h2$CHROM))
fst.h2$Fst<-fst.h2$WEIGHTED_FST
fst.h2$POS<-fst.h2$BIN_START
fst.h2$SNP<-fst.h2$BIN_START
fst.h2%>%select("SNP","CHR","POS","Fst")

colorset<-c("#003f5c","#58508d","#bc5090","#ff6361","#ffa600")
pdf("h2.fst.manhattan.pdf",width = 16.8,height = 8)
manhattan(fst.h2%>%filter(grepl("chr",CHROM))%>%select("SNP","CHR","POS","Fst"),
          chr="CHR",bp="POS",p="Fst",snp="SNP",
          col=colorset,logp=F,suggestiveline = F,genomewideline = F,ylab="Fst",ylim=c(0,1),
          font.lab=10,main="S.herbacea Female vs Male",
          cex.lab=2, cex.axis=2, cex.main=4, cex.sub=2)
dev.off() 




#h2
line=quantile(fst.h2$FST,0.99)
g_list = list()
i=1
for (id in sort(unique(fst.h2$CHROM))) {
  plot<-ggplot(fst.h2[which(fst.h2$CHROM==id),])+geom_hline(yintercept = line,colour="red",linetype = 'dashed')+
    geom_line(aes(x=BIN_START,y=FST))+
    ggtitle(id)+
    xlab("")+labs(y = (~F[ST]~  FvsM))+scale_x_continuous(breaks =seq(0,max(fst.h2[which(fst.h2$CHROM==id),2]),2e6),labels =seq(0,max(fst.h2[which(fst.h2$CHROM==id),2])/1e6,2))+
    scale_y_continuous(limits = c(-0.1,1))+theme_bw()+
    theme(axis.title = element_text(size=30),axis.text = element_text(size = 28),plot.title= element_text(size=30))
  print(plot)
  g_list[[i]] = plot
  i=i+1
}
pdf('h2-She1.FM-Fst.pdf', height =30, width = 60)
multiplot(g_list[[1]],g_list[[2]], g_list[[3]],g_list[[4]],g_list[[5]],g_list[[6]],g_list[[7]],g_list[[8]],g_list[[9]],g_list[[10]],g_list[[11]],g_list[[12]],g_list[[13]],g_list[[14]],g_list[[15]],g_list[[16]],g_list[[17]],g_list[[18]],g_list[[19]], cols=5)
dev.off()

#plot for only one chromosome
pdf('h1-She1.chr15.FM-Fst.pdf', height =6, width = 15)
ggplot(fst.h1[which(fst.h1$CHROM=="h1_chr15"),])+geom_line(aes(x=BIN_START,y=WEIGHTED_FST))+
  ggtitle("h1-she1 Fm Fst")+
  xlab("h1_chr15")+ylab("Fst")+scale_x_continuous(breaks =seq(0,max(fst.h1[which(fst.h1$CHROM=="h1_chr15"),2]),1e6),labels =seq(0,max(fst.h1[which(fst.h1$CHROM=="h1_chr15"),2])/1e6,1))+
  scale_y_continuous(limits = c(-0.1,1))+theme_bw()+
  theme(axis.title = element_text(size=30),axis.text = element_text(size = 28),plot.title= element_text(size=30))
dev.off()


pdf('h2-She1.chr15.FM-Fst.pdf', height =6, width = 15)
ggplot(fst.h2[which(fst.h2$CHROM=="h2_chr15"),])+geom_line(aes(x=BIN_START,y=WEIGHTED_FST))+
  ggtitle("h2-she1 Fm Fst")+
  xlab("h2_chr15")+ylab("Fst")+scale_x_continuous(breaks =seq(0,17e6,1e6),labels =seq(0,17,1))+
  scale_y_continuous(limits = c(-0.1,1))+theme_bw()+
  theme(axis.title = element_text(size=30),axis.text = element_text(size = 28),plot.title= element_text(size=30))
dev.off()


ggplot(fst.she)+geom_line(aes(x=BIN_START,y=WEIGHTED_FST))+
  ggtitle("She-32ind F-M Fst")+
  xlab("h1_chr15")+ylab("Fst")+scale_x_continuous(breaks =seq(0,15300000,1e6),labels =seq(0,15.3,1))+
  theme(axis.title = element_text(size=30),axis.text = element_text(size = 28),plot.title= element_text(size=30))
fst.she1 <- read.csv("fst/h1mask_She1-20ind.f-m.10k.chr15.windowed.weir.fst",sep = "\t")
ggplot(fst.she1)+geom_line(aes(x=BIN_START,y=WEIGHTED_FST))+
  ggtitle("She1-20ind F-M Fst")+
  xlab("h1_chr15")+ylab("Fst")+scale_x_continuous(breaks =seq(0,15300000,1e6),labels =seq(0,15.3,1))+
  theme(axis.title = element_text(size=30),axis.text = element_text(size = 28),plot.title= element_text(size=30))
fst.she2 <- read.csv("fst/h1mask_She2-12ind.f-m.10k.chr15.windowed.weir.fst",sep = "\t")
ggplot(fst.she2)+geom_line(aes(x=BIN_START,y=WEIGHTED_FST))+
  ggtitle("She2-12ind F-M Fst")+
  xlab("h1_chr15")+ylab("Fst")+scale_x_continuous(breaks =seq(0,15300000,1e6),labels =seq(0,15.3,1))+
  theme(axis.title = element_text(size=30),axis.text = element_text(size = 28),plot.title= element_text(size=30))

fst.she12 <-read.csv("fst/h1mask_She-32ind.1-2.10k.chr15.windowed.weir.fst",sep = "\t")
ggplot(fst.she12)+geom_line(aes(x=BIN_START,y=WEIGHTED_FST))+
  ggtitle("She 1-2 Fst")+
  xlab("h1_chr15")+ylab("Fst")+scale_x_continuous(breaks =seq(0,15300000,1e6),labels =seq(0,15.3,1))+
  theme(axis.title = element_text(size=30),axis.text = element_text(size = 28),plot.title= element_text(size=30))

fst.h1.chr15<-fst.h1%>%filter(CHR==15 & WEIGHTED_FST>0)
fst.h1.chr15<-fst.h1%>%filter(CHR==15)%>%mutate(WEIGHTED_FST=case_when(
  WEIGHTED_FST<=0 ~ 0,
  WEIGHTED_FST>0 ~ WEIGHTED_FST
))

quantile(fst.h1.chr15$WEIGHTED_FST,0.5)
quantile(fst.h1.chr15$WEIGHTED_FST,0.95)
quantile(fst.h1$FST,0.95, na.rm = TRUE)

plot(FST_values, type = "l", main = "FST across the chromosome")
abline(h = median_fst, col = "blue", lty = 2)  # 中位数cutoff线
abline(h = high_cutoff_fst, col = "red", lty = 2)  # 95百分位数cutoff线
legend("topright", legend=c("Median FST", "95th Percentile FST"), col=c("blue", "red"), lty=2)


#1-fst-pool####
setwd("/0Salix_Project/1Project_genome/5-call_variants/0final/2-pool/")
#h1
fst<-read.csv("h1_pool_f_m.fst",sep="\t",header = F)
colnames(fst)<-c("chr","coord","snp","cov","mincov","fst")
fstp1<-fst%>%mutate(fst=as.numeric(gsub("1:2=","",fst)))%>%mutate(coord=coord-5000)%>%filter(grepl("chr",chr))
line=quantile(fstp1$fst,0.99)
g_list = list()
i=1
for (id in sort(unique(fstp1$chr))) {
  # line<-quantile(fstp1[which(fstp1$chr==id),]$fst,0.95)
  plot<-ggplot(fstp1[which(fstp1$chr==id),])+geom_hline(yintercept = line,colour="red",linetype = 'dashed')+
    geom_line(aes(x=coord,y=fst))+
    ggtitle(id)+
    xlab("")+labs(y = (~F[ST]~  FvsM))+scale_x_continuous(breaks =seq(0,max(fstp1[which(fstp1$chr==id),2]),2e6),labels =seq(0,max(fstp1[which(fstp1$chr==id),2])/1e6,2))+
    scale_y_continuous(limits = c(0,0.4))+theme_bw()+
    theme(axis.title = element_text(size=30),axis.text = element_text(size = 28),plot.title= element_text(size=30))
  print(plot)
  g_list[[i]] = plot
  i=i+1
}
pdf('h1-She1.pool.FM-Fst.line.pdf', height =30, width = 60)
multiplot(g_list[[1]],g_list[[2]], g_list[[3]],g_list[[4]],g_list[[5]],g_list[[6]],g_list[[7]],g_list[[8]],g_list[[9]],g_list[[10]],g_list[[11]],g_list[[12]],g_list[[13]],g_list[[14]],g_list[[15]],g_list[[16]],g_list[[17]],g_list[[18]],g_list[[19]], cols=5)
dev.off()

#h2
fst<-read.csv("h2_pool_f_m.fst",sep="\t",header = F)
colnames(fst)<-c("chr","coord","snp","cov","mincov","fst")
fstp2<-fst%>%mutate(fst=as.numeric(gsub("1:2=","",fst)))%>%mutate(coord=coord-5000)%>%filter(grepl("chr",chr))
line=quantile(fstp2$fst,0.99)
head(fstp2)
g_list = list()
i=1
for (id in sort(unique(fstp2$chr))) {
  # line<-quantile(fstp2[which(fstp2$chr==id),]$fst,0.95)
  plot<-ggplot(fstp2[which(fstp2$chr==id),])+geom_hline(yintercept = line,colour="red",linetype = 'dashed')+
    geom_line(aes(x=coord,y=fst))+
    ggtitle(id)+
    xlab("")+labs(y = (~F[ST]~  FvsM))+scale_x_continuous(breaks =seq(0,max(fstp2[which(fstp2$chr==id),2]),2e6),labels =seq(0,max(fstp2[which(fstp2$chr==id),2])/1e6,2))+
    scale_y_continuous(limits = c(0,0.4))+theme_bw()+
    theme(axis.title = element_text(size=30),axis.text = element_text(size = 28),plot.title= element_text(size=30))
  print(plot)
  g_list[[i]] = plot
  i=i+1
}
pdf('h2-She1.pool.FM-Fst.line.pdf', height =30, width = 60)
multiplot(g_list[[1]],g_list[[2]], g_list[[3]],g_list[[4]],g_list[[5]],g_list[[6]],g_list[[7]],g_list[[8]],g_list[[9]],g_list[[10]],g_list[[11]],g_list[[12]],g_list[[13]],g_list[[14]],g_list[[15]],g_list[[16]],g_list[[17]],g_list[[18]],g_list[[19]], cols=5)
dev.off()


#2-pca####
pca<- read.csv("pca/h1mask_She-32ind.plink.pca.eigenvec",header = F,sep = " ")
ggplot(pca,aes(x=V3,y=V4))+geom_point()+geom_text_repel(aes(label=V1))


#3-dp####
#read files
setwd("/0Salix_Project/1Project_genome/5-call_variants/0final/3-ind/cov")

# individuals modified.
## h1
dpall <- read.csv("h1mask_She1.20ind.ave",sep = "\t")

dpall1<-dpall%>%replace(is.na(.), 0)%>%rowwise()%>%
  mutate(f=mean(c(h1mask_She1.F05,h1mask_She1.F06,h1mask_She1.F08,h1mask_She1.F09,h1mask_She1.F11,h1mask_She1.F16,h1mask_She1.F22,h1mask_She1.F32,h1mask_She1.F36,h1mask_She1.F40,h1mask_She1.M05)),
         m=mean(c(h1mask_She1.M06,h1mask_She1.M08,h1mask_She1.M09,h1mask_She1.M11,h1mask_She1.M16,h1mask_She1.M22,h1mask_She1.M32,h1mask_She1.M36,h1mask_She1.M40)))%>%
  select(chr,start,f,m)%>%rename(coord=start)%>%mutate(ratio=f/(f+m))%>%filter(grepl("chr",chr))


g_list = list()
i=1
for (id in sort(unique(dpall1$chr))) {
  plot<-ggplot(dpall1[which(dpall1$chr==id),])+geom_point(aes(x=coord,y=ratio))+
    ggtitle(id)+
    xlab("")+ylab("Depth F/(F+M)")+scale_x_continuous(breaks =seq(0,max(dpall1[which(dpall1$chr==id),2]),2e6),labels =seq(0,max(dpall1[which(dpall1$chr==id),2])/1e6,2))+
    scale_y_continuous(limits = c(-0.1,1))+theme_bw()+
    theme(axis.title = element_text(size=30),axis.text = element_text(size = 28),plot.title= element_text(size=30))
  print(plot)
  g_list[[i]] = plot
  i=i+1
}
pdf('h1-ind.FM-dp.pdf', height =30, width = 60)
multiplot(g_list[[1]],g_list[[2]], g_list[[3]],g_list[[4]],g_list[[5]],g_list[[6]],g_list[[7]],g_list[[8]],g_list[[9]],g_list[[10]],g_list[[11]],g_list[[12]],g_list[[13]],g_list[[14]],g_list[[15]],g_list[[16]],g_list[[17]],g_list[[18]],g_list[[19]], cols=5)
dev.off()

## h2
dpall.h2<-read.csv("h2mask_She1.20ind.ave",sep = "\t")
dpall1.h2<-dpall.h2%>%replace(is.na(.), 0)%>%rowwise()%>%
  mutate(f=mean(c(h2mask_She1.F05,h2mask_She1.F06,h2mask_She1.F08,h2mask_She1.F09,h2mask_She1.F11,h2mask_She1.F16,h2mask_She1.F22,h2mask_She1.F32,h2mask_She1.F36,h2mask_She1.F40,h2mask_She1.M05)),
         m=mean(c(h2mask_She1.M06,h2mask_She1.M08,h2mask_She1.M09,h2mask_She1.M11,h2mask_She1.M16,h2mask_She1.M22,h2mask_She1.M32,h2mask_She1.M36,h2mask_She1.M40)))%>%
  select(chr,start,f,m)%>%rename(coord=start)%>%mutate(ratio=f/(f+m))%>%filter(grepl("chr",chr))

g_list = list()
i=1
for (id in sort(unique(dpall1.h2$chr))) {
  plot<-ggplot(dpall1.h2[which(dpall1.h2$chr==id),])+geom_point(aes(x=coord,y=ratio))+
    ggtitle(id)+
    xlab("")+ylab("Depth F/(F+M)")+scale_x_continuous(breaks =seq(0,max(dpall1.h2[which(dpall1.h2$chr==id),2]),2e6),labels =seq(0,max(dpall1.h2[which(dpall1.h2$chr==id),2])/1e6,2))+
    scale_y_continuous(limits = c(-0.1,1))+theme_bw()+
    theme(axis.title = element_text(size=30),axis.text = element_text(size = 28),plot.title= element_text(size=30))
  print(plot)
  g_list[[i]] = plot
  i=i+1
}
pdf('h2-ind.FM-dp.pdf', height =30, width = 60)
multiplot(g_list[[1]],g_list[[2]], g_list[[3]],g_list[[4]],g_list[[5]],g_list[[6]],g_list[[7]],g_list[[8]],g_list[[9]],g_list[[10]],g_list[[11]],g_list[[12]],g_list[[13]],g_list[[14]],g_list[[15]],g_list[[16]],g_list[[17]],g_list[[18]],g_list[[19]], cols=5)
dev.off()

#test
test<-dpall.h2%>%replace(is.na(.), 0)%>%rowwise()%>%
  mutate(f=mean(c(h2mask_She1.F05,h2mask_She1.F06,h2mask_She1.F08,h2mask_She1.F09,h2mask_She1.F11,h2mask_She1.F16,h2mask_She1.F22,h2mask_She1.F32,h2mask_She1.F36,h2mask_She1.F40,h2mask_She1.M05)),
          m=mean(c(h2mask_She1.M06,h2mask_She1.M08,h2mask_She1.M09,h2mask_She1.M11,h2mask_She1.M16,h2mask_She1.M22,h2mask_She1.M32,h2mask_She1.M36,h2mask_She1.M40)))%>%mutate(ratio=f/(f+m))%>%
  filter(chr=="h2_chr15")%>%pivot_longer(cols = h2mask_She1.F05:h2mask_She1.M40,names_to = "id",values_to = "dp")%>%mutate(sex=case_when(
    grepl("F",id)~"F", grepl("M",id)~"M"))
  
# pools
dpp <- read.table("../../2-pool/h1_pool_fm.normalised.cov")
dpp.h2 <- read.table("../../2-pool/h2_pool_fm.cov")

dp.15 <- dp%>%rename("chr"=h1_chr19,"coord"=X0)%>%filter(chr=="h1_chr15")%>%rowwise(chr,coord)%>%summarise(f=mean(c(h1_swi_f05:h1_swi_f40)),m=mean(h1_swi_m05:h1_swi_m40))

test<-dp%>%rename("chr"=h1_chr19,"coord"=X0)%>%filter(chr=="h1_chr15")%>%
  select(h1_swi_F06:h1_swi_F11,h1_swi_F22,h1_swi_F32,h1_swi_F36,h1_swi_M06:h1_swi_M11,h1_swi_M22,h1_swi_M32,h1_swi_M36)
library(RColorBrewer)
pheatmap(as.matrix(test), cluster_row = FALSE,cluster_cols = FALSE, scale = "row",labels_row = seq(0,15.8,0.5),main = "Read depth on h1_chr15",fontsize_col = 20,color = colorRampPalette(rev(brewer.pal(n = 10, name ="Spectral")))(100))

dpp15 <-dpp%>%rename("pool_f"=V4,"pool_m"=V5,"coord"=V2)%>%
  filter(V1=="h1_chr15")%>%select(coord,pool_f,pool_m)%>%mutate(fratio=pool_f/(pool_m+pool_f))
ggplot(dpp15)+geom_line(aes(x=coord,y=pool_f))+scale_y_continuous(limits = c(0,5))+xlab("")+ylab("female depth")+ggtitle("")+
  theme(axis.title = element_text(size=30),axis.text = element_text(size = 28),plot.title= element_text(size=30))+scale_x_continuous(breaks =seq(0,15300000,1e6),labels =seq(0,15.3,1) )

#plot dp for whole genome####

dpp2<-dpp%>%rename("CHROM"=V1,"start"=V2,"end"=V3,"pool_f"=V4,"pool_m"=V5)%>%
  mutate(CHR=as.numeric(sub("[a-z]1_chr","",CHROM)))%>%
  mutate(fratio=pool_f/(pool_m+pool_f))%>%filter(grepl("chr",CHROM))%>%filter(!is.na(fratio))

dpp2.h2<-dpp.h2%>%rename("CHROM"=V1,"start"=V2,"end"=V3,"pool_f"=V4,"pool_m"=V5)%>%
  mutate(CHR=as.numeric(sub("[a-z]2_chr","",CHROM)))%>%
  mutate(fratio=pool_f/(pool_m+pool_f))%>%filter(grepl("chr",CHROM))%>%
  #filter(pool_f<=1 & pool_m<=1)%>%
  filter(!is.na(fratio))

colorset<-c("#003f5c","#58508d","#bc5090","#ff6361","#ffa600")
manhattan(dpp2%>%filter(grepl("chr",CHROM)),chr="CHR",bp="start",p="fratio",snp="fratio",
          col=colorset,logp=F,suggestiveline = F,genomewideline = F,ylab="Depth F/(F+M)",ylim=c(0,1),font.lab=10,main="S.herbacea Depth ratio",
          cex.lab=2, cex.axis=2, cex.main=4, cex.sub=2)

length(unique(dpp2$CHROM))
#h1
g_list = list()
i=1
for (id in sort(unique(dpp2$CHROM))) {
  plot<-ggplot(dpp2[which(dpp2$CHROM==id),])+geom_point(aes(x=start,y=fratio))+
    ggtitle(id)+
    xlab("")+ylab("Depth F/(F+M)")+scale_x_continuous(breaks =seq(0,max(dpp2[which(dpp2$CHROM==id),2]),2e6),labels =seq(0,max(dpp2[which(dpp2$CHROM==id),2])/1e6,2))+
    scale_y_continuous(limits = c(-0.1,1))+theme_bw()+
    theme(axis.title = element_text(size=30),axis.text = element_text(size = 28),plot.title= element_text(size=30))
  print(plot)
  g_list[[i]] = plot
  i=i+1
}
pdf('h1-pool.FM-dp.pdf', height =30, width = 60)
multiplot(g_list[[1]],g_list[[2]], g_list[[3]],g_list[[4]],g_list[[5]],g_list[[6]],g_list[[7]],g_list[[8]],g_list[[9]],g_list[[10]],g_list[[11]],g_list[[12]],g_list[[13]],g_list[[14]],g_list[[15]],g_list[[16]],g_list[[17]],g_list[[18]],g_list[[19]], cols=5)
dev.off()

length(unique(dpp2.h2$CHROM))
#h1
g_list = list()
i=1
for (id in sort(unique(dpp2.h2$CHROM))) {
  plot<-ggplot(dpp2.h2[which(dpp2.h2$CHROM==id),])+geom_point(aes(x=start,y=fratio))+
    ggtitle(id)+
    xlab("")+ylab("Depth F/(F+M)")+scale_x_continuous(breaks =seq(0,max(dpp2.h2[which(dpp2.h2$CHROM==id),2]),2e6),labels =seq(0,max(dpp2.h2[which(dpp2.h2$CHROM==id),2])/1e6,2))+
    scale_y_continuous(limits = c(-0.1,1))+theme_bw()+
    theme(axis.title = element_text(size=30),axis.text = element_text(size = 28),plot.title= element_text(size=30))
  print(plot)
  g_list[[i]] = plot
  i=i+1
}
pdf('h2-pool.FM-dp.pdf', height =30, width = 60)
multiplot(g_list[[1]],g_list[[2]], g_list[[3]],g_list[[4]],g_list[[5]],g_list[[6]],g_list[[7]],g_list[[8]],g_list[[9]],g_list[[10]],g_list[[11]],g_list[[12]],g_list[[13]],g_list[[14]],g_list[[15]],g_list[[16]],g_list[[17]],g_list[[18]],g_list[[19]], cols=5)
dev.off()



#4-gene####

#allele frequency###
p1<- fmaf%>%filter(!is.na(mf))%>%filter(ff>0.45&mf<0.1)%>%ggplot()+geom_point(aes(x=coord,y="f"),color="red")
p2<-fmaf%>%filter(!is.na(ff))%>%filter(mf>0.45&ff<0.1)%>%ggplot()+geom_point(aes(x=coord,y="male"),color="blue")

fmaf1<-fmaf%>%filter(!is.na(mf))%>%mutate(af=ff-mf)%>%mutate(color=case_when(af>0.4 ~"F-bias",af< -0.4 ~"M-bias"))
r2<-ggplot(fmaf1)+geom_point(aes(x=coord,y=af,color=color))+scale_x_continuous(breaks =seq(0,15300000,1e6),labels =seq(0,15.3,1) )+xlab("h1_chr15")+ylab("F-M")+ggtitle("Allele frequency difference")+scale_y_continuous(breaks = seq(-0.5,0.5,0.1))+
  scale_color_discrete(guide="none")+
  theme(axis.title = element_text(size=30),axis.text = element_text(size = 28),plot.title= element_text(size=30))

png(filename = "/0Salix_Project/1Project_genome/5-call_variants/0final/3-ind/h1chr15.mummer.afd.png",width = 2800,height = 2000,res = 80)
ggarrange(p1,r2,ncol = 1,nrow = 2,heights= c(1.5,1),align = "v")
dev.off()


#5-heterozygosity rate of average f and m ####
##5-1 read files
setwd("/0Salix_Project/1Project_genome/5-call_variants/0final/3-ind/heter")
### h1 chr15
het3 <- read.csv("h1mask_She1-20ind.snps.filter.chr15.bi.hetden10kb.ave",sep = " ")
colnames(het3)<-c("coord","female","male")
head(het3)

### h1 all chr
hetall.h1<-read.csv("h1_mask.She1-20ind.snp.vcf.hetden10kb.ave2",sep = " ")
colnames(hetall.h1)<-c("chr","coord","female","male")
head(hetall.h1)

#### h2 all chr
hetall.h2<-read.csv("h2_mask.She1-20ind.snp.vcf.hetden10kb.ave2",sep = " ")
colnames(hetall.h2)<-c("chr","coord","female","male")
head(hetall.h2)

##5-2 normalize by mean of each sex
### all chr
hetall2.h1<-hetall.h1%>%mutate(f_mod=female/mean(hetall.h1$female),m_mod=male/mean(hetall.h1$male))
hetall2.h2<-hetall.h2%>%mutate(f_mod=female/mean(hetall.h2$female),m_mod=male/mean(hetall.h2$male))


##5-3 calculate the ratio of f/m and mark bias
### all chr
hetall3.h1<-hetall2.h1%>%mutate(ratio=case_when(
  f_mod==0 & m_mod>0 ~ 0, 
  f_mod>0 & m_mod==0 ~ 100, 
  f_mod>0 & m_mod>0 ~ f_mod/m_mod))%>%mutate(color=case_when(
    ratio>2 ~ "f-bias",
    ratio<0.5 ~ "m-bias", ratio==100 ~ "infinite"
  ))

hetall4.h1<-hetall.h1%>%mutate(ratio=case_when(
  female==0 & male>0 ~ 0, 
  female>0 & male==0 ~ 100, 
  female>0 & male>0 ~ female/male))%>%
  filter(ratio>0 & ratio <100)%>%mutate(color=case_when(
    ratio>2 ~ "f-bias",
    ratio<0.5 ~ "m-bias"))


hetall3.h2<-hetall2.h2%>%mutate(ratio=case_when(
  f_mod==0 & m_mod>0 ~ 0, 
  f_mod>0 & m_mod==0 ~ 100, 
  f_mod>0 & m_mod>0 ~ f_mod/m_mod))%>%mutate(color=case_when(
    ratio>2 ~ "f-bias",
    ratio<0.5 ~ "m-bias", ratio==100 ~ "infinite"
  ))
hetall4.h2<-hetall.h2%>%mutate(ratio=case_when(
  female==0 & male>0 ~ 0, 
  female>0 & male==0 ~ 100, 
  female>0 & male>0 ~ female/male))%>%
  filter(ratio>0 & ratio <100)%>%mutate(color=case_when(
    ratio>2 ~ "f-bias",
    ratio<0.5 ~ "m-bias"))

##5-4 plot

### plot loop for each chr
g_list = list()
i=1
for (id in sort(unique(hetall4.h1$chr))) {
  plot<-ggplot(hetall4.h1[which(hetall4.h1$chr==id),])+
    geom_point(aes(x=coord,y=log10(ratio),color=color))+
    ggtitle(id)+
    xlab("")+labs(y=("Heter F/M" ~log[10]))+
    scale_x_continuous(breaks =seq(0,max(hetall4.h1[which(hetall4.h1$chr==id),2]),2e6),labels =seq(0,max(hetall4.h1[which(hetall4.h1$chr==id),2])/1e6,2))+
    scale_y_continuous(limits = c(-2,2))+theme_bw()+
    theme(axis.title = element_text(size=30),axis.text = element_text(size = 28),plot.title= element_text(size=30),legend.position = 'none')
  print(plot)
  g_list[[i]] = plot
  i=i+1
}
pdf('h1-She1.FM-het4.pdf', height =30, width = 60)
multiplot(g_list[[1]],g_list[[2]], g_list[[3]],g_list[[4]],g_list[[5]],g_list[[6]],g_list[[7]],g_list[[8]],g_list[[9]],g_list[[10]],g_list[[11]],g_list[[12]],g_list[[13]],g_list[[14]],g_list[[15]],g_list[[16]],g_list[[17]],g_list[[18]],g_list[[19]], cols=5)
dev.off()
g_list = list()
i=1
for (id in sort(unique(hetall4.h2$chr))) {
  plot<-ggplot(hetall4.h2[which(hetall4.h2$chr==id),])+
    geom_point(aes(x=coord,y=log10(ratio),color=color))+
    ggtitle(id)+
    xlab("")+labs(y=("Heter F/M" ~log[10]))+
    scale_x_continuous(breaks =seq(0,max(hetall4.h2[which(hetall4.h2$chr==id),2]),2e6),labels =seq(0,max(hetall4.h2[which(hetall4.h2$chr==id),2])/1e6,2))+
    scale_y_continuous(limits = c(-2,2))+theme_bw()+
    theme(axis.title = element_text(size=30),axis.text = element_text(size = 28),plot.title= element_text(size=30),legend.position = 'none')
  print(plot)
  g_list[[i]] = plot
  i=i+1
}
pdf('h2-She1.FM-het4.pdf', height =30, width = 60)
multiplot(g_list[[1]],g_list[[2]], g_list[[3]],g_list[[4]],g_list[[5]],g_list[[6]],g_list[[7]],g_list[[8]],g_list[[9]],g_list[[10]],g_list[[11]],g_list[[12]],g_list[[13]],g_list[[14]],g_list[[15]],g_list[[16]],g_list[[17]],g_list[[18]],g_list[[19]], cols=5)
dev.off()



#alternative ploting
library(patchwork)
wrap_plots(g_list[[1]],g_list[[2]], g_list[[3]],g_list[[4]],g_list[[5]],g_list[[6]],g_list[[7]],g_list[[8]],g_list[[9]],g_list[[10]],g_list[[11]],g_list[[12]],g_list[[13]],g_list[[14]],g_list[[15]],g_list[[16]],g_list[[17]],g_list[[18]],g_list[[19]],ncol = 5,nrow = 4,byrow = F)+
  plot_layout(guides = "collect") &theme(legend.position = "right")



#plot the ratio of female and male
png(filename = "heterzygosity_count_fmratio_10kb.chr15.png",width = 3000,height = 1500,res = 250)
ggplot(het5)+ 
  annotate("rect",xmin=3e6,xmax=5e6,ymin=-Inf,ymax=Inf,fill="#FFC36D", colour = "grey",alpha=0.5)+
  annotate("rect",xmin=6e6,xmax=11.37e6,ymin=-Inf,ymax=Inf,fill="grey70", colour = "grey",alpha=0.5)+
  annotate("rect",xmin=10.8e6,xmax=11.4e6,ymin=-Inf,ymax=Inf,fill="#FFC36D", colour = "grey",alpha=0.5)+
  geom_point(aes(x=coord,y=log10(ratio),color=color))+
  scale_x_continuous(limits = c(0,16e6),breaks = seq(0,16e6,2e6))+
  scale_y_continuous(limits = c(-1,2),breaks = seq(-2,2,0.5))+
  xlab("h1_chr15")+ylab("F/M ratio heterozygous(10kb)")+
  geom_line(aes(x=coord,y=0))+
  theme_bw()+theme(text = element_text(size=20))
dev.off()
                   
ggplot(het5)
  geom_line(aes(x=coord,y=0.5))+
  geom_point(aes(x=coord,y=ratio,color=color))+
  scale_x_continuous(limits = c(0,16e6),breaks = seq(0,16e6,2e6))+
  scale_y_continuous(limits=c(0,1),breaks = seq(0,1,0.2))+
  xlab("h1_chr15")+ylab("F/F+M heterozygosity number per 10kb")+
  theme_bw()+theme(text = element_text(size=20))

#directly plot for each sex
geom_line(aes(x=coord,y=female),color="#FFC36D", linewidth=0.3)+geom_area(aes(x=coord,y=female),fill="#FFC36D", alpha=0.2)+
  geom_line(aes(x=coord,y=male),color="#69b3a2", linewidth=0.3)+geom_area(aes(x=coord,y=male),fill="#69b3a2", alpha=0.2)+
  scale_y_continuous(limits = c(0,200))+scale_x_continuous(limits = c(0,16e6),breaks = seq(0,16e6,2e6))+
  xlab("h1_chr15")+ylab("Average heterozygosity number within sex per 10kb")

## female-specific regions coordination ####
setwd("/0Salix_Project/1Project_genome/5-call_variants/0final/3-ind/allele_diff/")
hemi<-read.table("h1mask_She1-dp.M0F1.bed")
colnames(hemi)<-c("chr","start","end")
head(hemi)

alle<-read.table("h1mask_She1-20ind.all_variant.allele_fheter_mhomo.strict.geno.male11.pos")
colnames(alle)<-c("chr","pos")
head(alle)
pdf("h1.f-specific.coord.pdf",width = 30,height = 4)
ggplot()+
  annotate("rect",xmin=3.14e6,xmax=6e6,ymin=-Inf,ymax=Inf,fill="#fdbe64", colour = "grey",alpha=0.5)+
  annotate("rect",xmin=6e6,xmax=11.53e6,ymin=-Inf,ymax=Inf,fill="#ffe7b8",alpha=0.5)+
  annotate("rect",xmin=10.8e6,xmax=11.53e6,ymin=-Inf,ymax=Inf,fill="#FF8C00", colour = "grey",alpha=0.5)+
  geom_segment(data=hemi,aes(x=start,xend=end,y=0,yend=1),color="#619d3c")+
  # geom_segment(data=alle,aes(x=pos,xend=pos,y=0,yend=1),color="#008080")+
  scale_x_continuous(limits=c(0,16e6),breaks =seq(0,16000000,1e6),labels =seq(0,16,1))+scale_y_continuous(breaks =seq(0,1,1),labels =seq(0,1,1))+theme_bw()+xlab("h1 chr15")+
  theme(axis.title = element_text(size=30),axis.text = element_text(size = 28),plot.title= element_text(size=30),axis.title.y=element_blank(),axis.text.y = element_blank())
dev.off()



#6-merge####
#alignment
a1<-
  ggplot(coord)+
  annotate("rect",xmin=3.14e6,xmax=6e6,ymin=-Inf,ymax=Inf,fill="#fdbe64", colour = "grey",alpha=0.5)+
  annotate("rect",xmin=6e6,xmax=11.53e6,ymin=-Inf,ymax=Inf,fill="#ffe7b8",alpha=0.5)+
  annotate("rect",xmin=10.8e6,xmax=11.53e6,ymin=-Inf,ymax=Inf,fill="#FF8C00", colour = "grey",alpha=0.5)+
  geom_segment(aes(x=h1_chr15_start,y=h2_chr15_start,xend=h1_chr15_end,yend=h2_chr15_end),size=3)+
   xlab("")+ylab("h2 chr15 (Mb)")+ #ggtitle("Alignment")+
  scale_x_continuous(breaks =seq(0,16000000,1e6),labels =seq(0,16,1),minor_breaks = NULL)+
  scale_y_continuous(breaks =seq(0,17000000,1e6),labels =seq(0,17,1),minor_breaks = NULL)+theme_bw()+
  theme(axis.title = element_text(size=30),axis.text = element_text(size = 28),
        plot.title= element_text(size=30),axis.title.x=element_blank(), axis.text.x=element_blank())



#fst
line<-quantile(fst.h1[which(fst.h1$CHROM=="h1_chr15"),]$FST,0.95)
#a2<-ggplot(fst.h1[which(fst.h1$CHROM=="h1_chr15"),])+
  geom_hline(yintercept = line,colour="red",linetype = 'dashed')+
  annotate("rect",xmin=3.14e6,xmax=6e6,ymin=-Inf,ymax=Inf,fill="#fdbe64", colour = "grey",alpha=0.5)+
  annotate("rect",xmin=6e6,xmax=11.53e6,ymin=-Inf,ymax=Inf,fill="#ffe7b8",alpha=0.5)+
  annotate("rect",xmin=10.8e6,xmax=11.53e6,ymin=-Inf,ymax=Inf,fill="#FF8C00", colour = "grey",alpha=0.5)+
  geom_line(aes(x=BIN_START,y=WEIGHTED_FST))+
  xlab("")+ #labs(y=(~F[ST]), title = (~F[ST]~  " F vs M"))+
  scale_x_continuous(breaks =seq(0,16e6,1e6),labels =seq(0,16,1))+
  scale_y_continuous(limits = c(-0.1,1))+theme_bw()+
  theme(axis.title = element_text(size=30),axis.text = element_text(size = 28),plot.title= element_text(size=30),axis.title.x=element_blank(), axis.text.x=element_blank())

line<-quantile(fstp1$fst,0.95)
a2<-ggplot(fstp1[which(fstp1$chr=="h1_chr15"),])+
  geom_hline(yintercept = line,colour="red",linetype = 'dashed')+
  annotate("rect",xmin=3.14e6,xmax=6e6,ymin=-Inf,ymax=Inf,fill="#fdbe64", colour = "grey",alpha=0.5)+
  annotate("rect",xmin=6e6,xmax=11.53e6,ymin=-Inf,ymax=Inf,fill="#ffe7b8",alpha=0.5)+
  annotate("rect",xmin=10.8e6,xmax=11.53e6,ymin=-Inf,ymax=Inf,fill="#FF8C00", colour = "grey",alpha=0.5)+
  geom_line(aes(x=coord,y=fst))+
  xlab("")+ #labs(y=(~F[ST]), title = (~F[ST]~  " F vs M"))+
  scale_x_continuous(breaks =seq(0,16e6,1e6),labels =seq(0,16,1),minor_breaks = NULL)+
  scale_y_continuous(limits = c(0,0.4),minor_breaks = NULL)+theme_bw()+
  theme(axis.title = element_text(size=30),axis.text = element_text(size = 28),plot.title= element_text(size=30),axis.title.x=element_blank(), axis.text.x=element_blank())


#dp
a3<-ggplot(dpp15)+
  annotate("rect",xmin=3.14e6,xmax=6e6,ymin=-Inf,ymax=Inf,fill="#fdbe64", colour = "grey",alpha=0.5)+
  annotate("rect",xmin=6e6,xmax=11.53e6,ymin=-Inf,ymax=Inf,fill="#ffe7b8",alpha=0.5)+
  annotate("rect",xmin=10.8e6,xmax=11.53e6,ymin=-Inf,ymax=Inf,fill="#FF8C00", colour = "grey",alpha=0.5)+
  geom_point(aes(x=coord,y=pool_f/(pool_m+pool_f)))+
  ylab("dp")+#ggtitle("Depth F/(F+M)")+
  theme_bw()+
  theme(axis.title = element_text(size=30),axis.text = element_text(size = 28),plot.title= element_text(size=30),
        axis.title.x=element_blank(), axis.text.x=element_blank())+
  scale_x_continuous(breaks =seq(0,16e6,1e6),labels =seq(0,16,1),minor_breaks = NULL )+
  scale_y_continuous(limits = c(0,1),minor_breaks = NULL)

#hemy
a5<-ggplot()+
  annotate("rect",xmin=3.14e6,xmax=6e6,ymin=-Inf,ymax=Inf,fill="#fdbe64", colour = "grey",alpha=0.5)+
  annotate("rect",xmin=6e6,xmax=11.53e6,ymin=-Inf,ymax=Inf,fill="#ffe7b8",alpha=0.5)+
  annotate("rect",xmin=10.8e6,xmax=11.53e6,ymin=-Inf,ymax=Inf,fill="#FF8C00", colour = "grey",alpha=0.5)+
  geom_segment(data=hemi,aes(x=start,xend=end,y=0,yend=1),color="#619d3c")+
  xlab("")+ylab("Hemi")+#ggtitle("Hemizygous sites")+
  # geom_segment(data=alle,aes(x=pos,xend=pos,y=0,yend=1),color="#008080")+
  scale_x_continuous(limits=c(0,16e6),breaks =seq(0,16000000,1e6),labels =seq(0,16,1),minor_breaks = NULL)+
    scale_y_continuous(breaks = NULL)+theme_bw()+
  theme(axis.title = element_text(size=30),axis.text = element_text(size = 28),plot.title= element_text(size=30),
          axis.text.x=element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank())


#heter
a4<-ggplot(het5)+ 
  annotate("rect",xmin=3.14e6,xmax=6e6,ymin=-Inf,ymax=Inf,fill="#fdbe64", colour = "grey",alpha=0.5)+
  annotate("rect",xmin=6e6,xmax=11.53e6,ymin=-Inf,ymax=Inf,fill="#ffe7b8",alpha=0.5)+
  annotate("rect",xmin=10.8e6,xmax=11.53e6,ymin=-Inf,ymax=Inf,fill="#FF8C00", colour = "grey",alpha=0.5)+
  geom_point(aes(x=coord,y=log10(ratio),color=color))+
  xlab("h1 chr15 (Mb)")+labs(y ="Heter")+#title=("Number of heterozygous sites F/M" ~log[10]),
  scale_x_continuous(breaks =seq(0,16e6,1e6),labels =seq(0,16,1),minor_breaks = NULL)+guides(color=FALSE)+
  scale_y_continuous(limits = c(-1,2),breaks = seq(-2,2,1),minor_breaks = NULL)+theme_bw()+
  theme(axis.title = element_text(size=30),axis.text = element_text(size = 28),plot.title= element_text(size=30))

setwd("/0Salix_Project/1Project_genome/5-call_variants/0final/")
pdf("h1_chr15.align-fst-dp-heter.pdf",width = 14,height = 20)
ggarrange(a1,a3,a2,a5,a4,ncol = 1,nrow = 5,heights= c(2,0.8,0.8,0.5,0.8),align = "v")
dev.off()
