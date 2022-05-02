# -- Polygenic counts with standardized length ------------

setwd('/workdir/Therkildsen/Angsd/AllsamplesMinusMixUps/Top0.99DiffStatGeno/Beagle/')
GenoPhenoAll=read.table('D_PolygenicDosage_beagle_Top0.99DiffStat_2.txt',header=T)
head(GenoPhenoAll)

# Compute standardized length
table(GenoPhenoAll$Population)
lines=c('R1Gen5','R2Gen5','D1Gen5','D2Gen5','U1Gen5','U2Gen5','RGen0')
gen0mean = mean(GenoPhenoAll[which(GenoPhenoAll$Population == 'RGen0'),'STD_Length'])
gen5mean = mean(GenoPhenoAll[which(GenoPhenoAll$Population %in% c('R1Gen5', 'R2Gen5')),'STD_Length'])
#GenoPhenoAll$Population[(which(GenoPhenoAll$Population %in% c('R1Gen5', 'R2Gen5')))]
gen10mean = mean(GenoPhenoAll[which(GenoPhenoAll$Population %in% c('R1Gen10', 'R2Gen10')),'STD_Length'])


GenoPhenoAll$NormLength=rep(NA,nrow(GenoPhenoAll))
GenoPhenoAll[grep('Gen0', GenoPhenoAll$Population), 'NormLength'] = GenoPhenoAll[grep('Gen0', GenoPhenoAll$Population), 'STD_Length'] - genmean
GenoPhenoAll[grep('Gen5', GenoPhenoAll$Population), 'NormLength'] = GenoPhenoAll[grep('Gen5', GenoPhenoAll$Population), 'STD_Length'] - gen5mean
GenoPhenoAll[grep('Gen10', GenoPhenoAll$Population), 'NormLength'] = GenoPhenoAll[grep('Gen10', GenoPhenoAll$Population), 'STD_Length'] - gen10mean

             
 
colrs=c(brewer.pal(9,'Greens')[4],brewer.pal(10,'Paired')[4],'darkgoldenrod1',brewer.pal(9,'Oranges')[c(6)],'steelblue1','dodgerblue4')
colors=c(colrs,'gray36',colrs)
#lines=c('R1Gen5','R2Gen5','D1Gen5','D2Gen5','U1Gen5','U2Gen5','RGen0')

Reldata=read.table('/workdir/Therkildsen/MAPGD/AllExperimental_MinusMixups_AllContigs_RGen0_subsamp20K_stitchNames.rel',header=T)
Reldata$X.SAMPLE_X=as.character(Reldata$X.SAMPLE_X)
Reldata$SAMPLE_Y=as.character(Reldata$SAMPLE_Y)

rel=diag(1,nrow=length(which(Reldata$X.SAMPLE_X==Reldata$X.SAMPLE_X[1]))+1)
rel[lower.tri(rel)]=Reldata$θ_XY 
rel=t(rel)
rel[lower.tri(rel)]=Reldata$θ_XY
rownames(rel)=c(Reldata$X.SAMPLE_X[1],Reldata$SAMPLE_Y[1:length(which(Reldata$X.SAMPLE_X==Reldata$X.SAMPLE_X[1]))])
colnames(rel)=c(Reldata$X.SAMPLE_X[1],Reldata$SAMPLE_Y[1:length(which(Reldata$X.SAMPLE_X==Reldata$X.SAMPLE_X[1]))])
covmat=rel[which(row.names(rel)%in%GenoPhenoAll$SampleID),which(colnames(rel)%in%GenoPhenoAll$SampleID)]
covmat=covmat[order(match(row.names(covmat), GenoPhenoAll$SampleID)),order(match(colnames(covmat), GenoPhenoAll$SampleID))]
covmat=as.matrix(covmat)
#covmat[covmat<0]=0 #I don't think it's necessary to set negative values to 0
covmat=make.positive.definite(covmat) #use this if lmekin complains "In random term 1: A variance matrix is not non-negative definite"

rel=diag(1,nrow=length(which(Reldata$X.SAMPLE_X==Reldata$X.SAMPLE_X[1]))+1)
rel[lower.tri(rel)]=Reldata$γ_XY
rel=t(rel)
rel[lower.tri(rel)]=Reldata$γ_XY
rownames(rel)=c(Reldata$X.SAMPLE_X[1],Reldata$SAMPLE_Y[1:length(which(Reldata$X.SAMPLE_X==Reldata$X.SAMPLE_X[1]))])
colnames(rel)=c(Reldata$X.SAMPLE_X[1],Reldata$SAMPLE_Y[1:length(which(Reldata$X.SAMPLE_X==Reldata$X.SAMPLE_X[1]))])
covmat2=rel[which(row.names(rel)%in%GenoPhenoAll$SampleID),which(colnames(rel)%in%GenoPhenoAll$SampleID)]
covmat2=covmat2[order(match(row.names(covmat2),GenoPhenoAll$SampleID)),order(match(colnames(covmat2),GenoPhenoAll$SampleID))]
covmat2=as.matrix(covmat2)
#covmat[covmat<0]=0 #I don't think it's necessary to set negative values to 0
covmat2=make.positive.definite(covmat2) #use this if lmekin complains "In random term 1: A variance matrix is not non-negative definite"

rel=diag(1,nrow=length(which(Reldata$X.SAMPLE_X==Reldata$X.SAMPLE_X[1]))+1)
rel[lower.tri(rel)]=Reldata$δ
rel=t(rel)
rel[lower.tri(rel)]=Reldata$δ
rownames(rel)=c(Reldata$X.SAMPLE_X[1],Reldata$SAMPLE_Y[1:length(which(Reldata$X.SAMPLE_X==Reldata$X.SAMPLE_X[1]))])
colnames(rel)=c(Reldata$X.SAMPLE_X[1],Reldata$SAMPLE_Y[1:length(which(Reldata$X.SAMPLE_X==Reldata$X.SAMPLE_X[1]))])
covmat3=rel[which(row.names(rel)%in%GenoPhenoAll$SampleID),which(colnames(rel)%in%GenoPhenoAll$SampleID)]
covmat3=covmat3[order(match(row.names(covmat3),GenoPhenoAll$SampleID)),order(match(colnames(covmat3),GenoPhenoAll$SampleID))]
covmat3=as.matrix(covmat3)
#covmat[covmat<0]=0 #I don't think it's necessary to set negative values to 0
covmat3=make.positive.definite(covmat3) #use this if lmekin complains "In random term 1: A variance matrix is not non-negative definite"

rel=diag(1,nrow=length(which(Reldata$X.SAMPLE_X==Reldata$X.SAMPLE_X[1]))+1)
rel[lower.tri(rel)]=Reldata$δ
rel=t(rel)
rel[lower.tri(rel)]=Reldata$δ
rownames(rel)=c(Reldata$X.SAMPLE_X[1],Reldata$SAMPLE_Y[1:length(which(Reldata$X.SAMPLE_X==Reldata$X.SAMPLE_X[1]))])
colnames(rel)=c(Reldata$X.SAMPLE_X[1],Reldata$SAMPLE_Y[1:length(which(Reldata$X.SAMPLE_X==Reldata$X.SAMPLE_X[1]))])
covmat4=rel[which(row.names(rel)%in%GenoPhenoAll$SampleID),which(colnames(rel)%in%GenoPhenoAll$SampleID)]
covmat4=covmat4[order(match(row.names(covmat4),GenoPhenoAll$SampleID)),order(match(colnames(covmat4),GenoPhenoAll$SampleID))]
covmat4=as.matrix(covmat4)
#covmat[covmat<0]=0 #I don't think it's necessary to set negative values to 0
covmat4=make.positive.definite(covmat4) #use this if lmekin complains "In random term 1: A variance matrix is not non-negative definite"

######Gen0-5 & Individual Pops and combined plot



# With regression lines only going through the relevant area for each pop

Break1=21300
Break2=26000
ABwidth=5

setEPS()
postscript(file='Dcopies_Top0.99DiffStat_LocalizedRegLines.eps',height=6,width=6)

par(mai=c(1,1,0.5,0.5))
lines=c('R1Gen5','R2Gen5','D1Gen5','D2Gen5','U1Gen5','U2Gen5','RGen0')
for (focal in c('R1Gen5','R2Gen5','D1Gen5','D2Gen5','U1Gen5','U2Gen5','RGen0')){
  GenoPhenoSub= GenoPhenoAll[which(GenoPhenoAll$Population==focal),]
  smpeind=which(colnames(covmat) %in% GenoPhenoSub$SampleID)
  lmefitPop=lmekin(GenoPhenoSub$NormLength~GenoPhenoSub$IncrSums+(1|GenoPhenoSub$SampleID),varlist= coxmeMlist(list(covmat[smpeind, smpeind],covmat2[smpeind, smpeind],covmat3[smpeind, smpeind])),method='ML')
  print(focal)
  print(lmefitPop)
  line=focal
  #Uncomment these lines to output individual plots too
  # png(file=paste0(focal,'_Dcopies_Top0.99DiffStat.png'),height=8,width=8,units="in",res=300)
  # plot(GenoPhenoSub $IncrSums[which(GenoPhenoSub $Population==line)], GenoPhenoSub $NormLength[which(GenoPhenoSub $Population==line)],col=colors[which(lines==focal)],pch=16,xlab=paste(c('Polygenic dosage of D allele'),collapse=''),ylab='NormLength',main=paste(c(focal,' Top 1% Diffstat SNPs'),collapse=''),cex.axis=1.5,cex.lab=2)
  # abline(lmefitPop$coef$fixed,lwd=2,col=colors[which(lines==focal)])
  # #abline(lm(GenoPhenoAll $NormLength[which(GenoPhenoAll $Population==line)]~ GenoPhenoAll $IncrSums[which(GenoPhenoAll $Population==line)]),col=colors[which(lines==focal)],lwd=1.5)#
  # dev.off()
  
  if (focal=='R1Gen5'){
    plot(GenoPhenoAll $IncrSums[which(GenoPhenoAll $Population==line)], GenoPhenoAll $NormLength[which(GenoPhenoAll $Population==line)],col=colors[which(lines==focal)],pch=16,xlab=paste(c('Count of putative slow-growing alleles'),collapse=''),ylab='Individual length (mm)',main=paste(c(''),collapse=''),xlim=c(min(GenoPhenoAll$IncrSums)-500,max(GenoPhenoAll$IncrSums)+500),ylim=c(min(GenoPhenoAll$NormLength)-10,max(GenoPhenoAll$NormLength)+10),cex.axis=1.2,cex.lab=1.3,cex=1.5)
    #clip(21000,26000,20,90)
    ablineclip(lmefitPop$coef$fixed,lwd=ABwidth,col=colors[which(lines==focal)], x1=Break1, x2=Break2)
  }
  else {
    if (focal=='RGen0'){
      points(GenoPhenoAll $IncrSums[which(GenoPhenoAll $Population==line)], GenoPhenoAll $NormLength[which(GenoPhenoAll $Population==line)],col=colors[which(lines==focal)],pch=1,cex=1.5)
      ablineclip(lmefitPop$coef$fixed,lwd=ABwidth,col=colors[which(lines==focal)], x1=Break1, x2=Break2)
    }
    if (focal=='R2Gen5'){
      points(GenoPhenoAll $IncrSums[which(GenoPhenoAll $Population==line)], GenoPhenoAll $NormLength[which(GenoPhenoAll $Population==line)],col=colors[which(lines==focal)],pch=16,cex=1.5)
      ablineclip(lmefitPop$coef$fixed,lwd=ABwidth,col=colors[which(lines==focal)], x1=Break1, x2=Break2)
    }
    if (focal=='U1Gen5' | focal=='U2Gen5'){
      points(GenoPhenoAll $IncrSums[which(GenoPhenoAll $Population==line)], GenoPhenoAll $NormLength[which(GenoPhenoAll $Population==line)],col=colors[which(lines==focal)],pch=16,cex=1.5)
      ablineclip(lmefitPop$coef$fixed,lwd=ABwidth,col=colors[which(lines==focal)], x1=15800, x2=Break1)
    }
    if (focal=='D1Gen5' | focal=='D2Gen5'){
      points(GenoPhenoAll $IncrSums[which(GenoPhenoAll $Population==line)], GenoPhenoAll $NormLength[which(GenoPhenoAll $Population==line)],col=colors[which(lines==focal)],pch=16,cex=1.5)
      ablineclip(lmefitPop$coef$fixed,lwd=ABwidth,col=colors[which(lines==focal)], x1=Break2, x2=31500)
    }}
}
legend('bottomleft',legend=c('Gen0','Up1','Up2','Control1','Control2','Down1','Down2'),col=colors[c(7,5,6,1:4)],pch=c(1,rep(16,6)),bg='white',bty='n',cex=1.1,pt.cex=1.6, inset=0.02)
#lwd=4
dev.off()


