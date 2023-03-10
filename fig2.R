#Author: Cao Wei
#update: 2023/03/09
#content: figure2

#2A
temp<-subset(CRCKUL,cells=Cells(CRCKUL)[which(CRCKUL$celltype %in% c('Cancer','Epithelial'))])
temp <- NormalizeData(temp, verbose = FALSE)
temp <- FindVariableFeatures(temp)
temp <- ScaleData(temp, verbose = FALSE)
temp <- RunPCA(temp, verbose = FALSE)
temp <- FindNeighbors(temp, reduction = "pca", dims = 1:20,k.param = 20)
temp <- FindClusters(temp, resolution = 1)
temp <- RunUMAP(temp, dims = 1:20,return.model = T)
temp$patient_CNV<-temp$patient
temp$patient_CNV[which(temp$celltype=='Epithelial')]<-'All patients'
temp$patient_CNV<-factor(temp$patient_CNV,levels=c('All patients',sort(unique(CRCKUL$patient))))
DimPlot(temp,group.by='patient_CNV',label=F,raster=T)+ggtitle('')+scale_color_manual(name='Epithelial',labels = levels(temp$patient_CNV),values=c('red',scales::hue_pal()(12)))+theme(text = element_text(size=20))
#2B
library(RCircos)
data("UCSC.HG38.Human.CytoBandIdeogram")
RCircos.Set.Core.Components(UCSC.HG38.Human.CytoBandIdeogram, chr.exclude=c('chrY'),tracks.inside=12, tracks.outside=0 )
RCircos.Set.Plot.Area()     
RCircos.Chromosome.Ideogram.Plot()
p=1
temp.poly<-unique(temp[which(temp$celltype=='Cancer'&temp$stat!='noCNV'),c('patient','region','stat')])
for(i in names(cancer_list)){
  temp.data<-cancer_list[[i]]@assays$inferCNV@meta.features
  temp.data$plot<-0
  temp.data$plot[which(temp.data$segment %in% which(CRCKUL_CNV_ploidy[[x]]=='amplified'))]<-1
  temp.data$plot[which(temp.data$segment %in% which(CRCKUL_CNV_ploidy[[x]]=='deleted'))]<-(-1)
  RCircos.Heatmap.Plot(temp.data, data.col=5, track.num=p, side='in')
  p=p+1
}
#2C
gene_UMI_count<-rowSums(CRCKUL@assays$RNA@counts[,which(CRCKUL$celltype=='Cancer'|(CRCKUL$celltype=='Epithelial'&CRCKUL$tissue=='N'))])
CRCKUL_patient_gene_foldchange<-as.matrix(as.data.frame(lapply(names(cancer_list),function(x){
  temp<-CRCKUL_Epi_Cancer_list[[x]]
  cancer_exp<-rowMeans(temp@assays$RNA@data[,which(temp$celltype=='Cancer')])
  normal_exp<-rowMeans(temp@assays$RNA@data[,which(temp$celltype=='Epithelial')])
  return(cancer_exp/normal_exp)
})))
colnames(CRCKUL_patient_gene_foldchange)<-names(cancer_list)
CRCKUL_patient_gene_foldchange<-replace(CRCKUL_patient_gene_foldchange,is.nan(CRCKUL_patient_gene_foldchange),NA)
CRCKUL_patient_gene_foldchange<-replace(CRCKUL_patient_gene_foldchange,CRCKUL_patient_gene_foldchange>5,5)
CRCKUL_patient_gene_foldchange<-replace(CRCKUL_patient_gene_foldchange,CRCKUL_patient_gene_foldchange<0.2,0.2)
temp<-as.data.frame(t(as.data.frame(lapply(names(cancer_list),function(x){
  gene_in_region<-cancer_list[[x]]@assays$RNA@meta.features$CNV_region
  ratio_data<-data.frame(ratio=CRCKUL_patient_gene_foldchange[,x])
  ratio_data$stat<-NA
  ratio_data$stat[which(gene_in_region %in% which(CRCKUL_CNV_ploidy[[x]]=='amplified'))]<-'amplified'
  ratio_data$stat[which(gene_in_region %in% which(CRCKUL_CNV_ploidy[[x]]=='deleted'))]<-'deleted'
  ratio_data$stat[which(gene_in_region %in% which(CRCKUL_CNV_ploidy[[x]]=='noCNV'))]<-'noCNV'
  ratio_data$patient<-x
  ratio_data$gene_UMI_count<-gene_UMI_count
  return(t(ratio_data))
}))))
temp$gene_UMI_count<-as.numeric(temp$gene_UMI_count)
temp$ratio<-as.numeric(temp$ratio)
temp$stat<-factor(temp$stat,levels = c('amplified','noCNV','deleted'))
temp<-temp[which(!(is.na(temp$stat)|is.na(temp$ratio)|temp$gene_UMI_count<100)),]
temp<-reshape2::melt(reshape2::dcast(temp[,c('patient','stat','ratio')],patient~stat,fun.aggregate = mean,value.var = 'ratio',fill = NA_real_),id.vars=c('patient'),value.name = 'ratio',variable.name = "stat")
ggplot(temp,aes(x=stat,y=ratio))+geom_boxplot(aes(fill=stat))+geom_line(aes(group=patient),linetype='dashed',color='grey50')+scale_fill_manual(values=c(noCNV='white',amplified='red',deleted='blue'))+geom_hline(aes(yintercept=1),color='red',linetype='dotted')+theme_bw()+xlab('CNV region status')+ylab('Cancer/Epithelial expression ratio')+theme(text = element_text(size = 20),axis.text.x = element_text(size=16),legend.position = 'none')+ggprism::add_pvalue(rstatix::wilcox_test(temp, ratio ~ stat,paired = T),y.position = c(3.7,3.5,3.7),bracket.shorten = c(0.05, 0, 0.05),label.size = 5,label = "p = {p}")
#2D
temp.patient<-c('CRC1','KUL19','KUL24','KUL27','KUL29','KUL30')
temp.used.gene<-c('GTF3A','PSMA7','MRPS26','COMMD6','RPL7','SMIM26')
cowplot::plot_grid(plotlist = lapply(1:6,function(i){
  datas<-cancer_list[[temp.patient[i]]]@assays$RNA
  gene<-temp.used.gene[i]
  gene_in_region<-datas@meta.features[gene,'CNV_region']
  plot_obj<-data.frame(expression=expm1(datas@data[gene,]),CNV=cancer_list[[temp.patient[i]]]@assays$inferCNV@data[gene_in_region,])
  r_value<-signif(cor(plot_obj$expression,plot_obj$CNV),2)
  p_value<-signif(summary(lm(expression~CNV, plot_obj))$coefficients[2,4],2)
  fig<-ggplot(plot_obj, aes(CNV, expression)) +
    geom_point(color=scales::hue_pal()(12)[match(temp.patient[i],names(cancer_list))]) +
    geom_smooth(method = "lm",inherit.aes = F,aes(CNV, expression),formula = y~x) +
    theme_bw()+xlab('')+ylab('')+theme(axis.text.x = element_text(size=16),axis.text.y = element_text(size=16))+
    ggtitle(paste0(gene,'  R=',r_value,'  p.val=',p_value))
  return(fig)
}),ncol = 3)
#2E
ITH_list<-data.frame(patient=names(cancer_list))
rownames(ITH_list)<-ITH_list$patient
ITH_list$ITHgex<-unlist(lapply(cancer_list,function(x){
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, nfeatures=1000)
  pearson_value<-cor(as.matrix(x@assays$RNA@data[x@assays$RNA@var.features,]), method = "pearson", use = "pairwise.complete.obs")
  ITHgex <- quantile(pearson_value[lower.tri(pearson_value)],probs = 0.75)-quantile(pearson_value[lower.tri(pearson_value)],probs = 0.25)
  return(ITHgex)
}))
ITH_list$ITHcnv<-unlist(lapply(cancer_list[ITH_list$patient],function(x){
  pearson_value<-cor(as.matrix(x@assays$inferCNV@counts), method = "pearson", use = "pairwise.complete.obs")
  ITHcnv <- quantile(pearson_value[lower.tri(pearson_value)],probs = 0.75)-quantile(pearson_value[lower.tri(pearson_value)],probs = 0.25)
  return(ITHcnv)
}))
r.val<-cor(ITH_list$ITHcnv,ITH_list$ITHgex)
p.val<-summary(lm(ITHcnv~ITHgex, ITH_list))$coefficients[2,4]
fig<-ggplot(ITH_list,aes(x=ITHcnv,y=ITHgex))+geom_point(aes(size=5,color=patient))+stat_smooth(method="lm")+theme_bw()+geom_text(x=0.15, y=0.23, label=paste0('R=',signif(r.val,2),"\np.val=",signif(p.val,2)),size=8,color='black')+theme(text = element_text(size = 20),legend.position = 'none')
LabelPoints(plot=fig, points=ITH_list$patient,repel = TRUE,xnudge = 0, ynudge = 0)
#2F
Highvar_CNV_gene_corlist<-lapply(names(cancer_list),function(set){
  dd<-dealMerged(cancer_list[[set]])
  high_var_gene_list<-dd@assays$RNA@var.features
  gene_in_region<-cancer_list[[set]]@assays$RNA@meta.features[high_var_gene_list,'CNV_region']
  names(gene_in_region)<-high_var_gene_list
  gene_in_region<-gene_in_region[which(!is.na(gene_in_region))]
  d1<-cancer_list[[set]]@assays$RNA@data
  d2<-cancer_list[[set]]@assays$inferCNV@data
  cor_val<-unlist(lapply(names(gene_in_region),function(x){
    cor(expm1(d1[x,]),d2[gene_in_region[[x]],])
  }))
  df<-data.frame(gene=names(gene_in_region),Rval=cor_val,region=gene_in_region)
  df$pct<-rowSums(cancer_list[[set]]@assays$RNA@data[df$gene,]>0)/length(Cells(cancer_list[[set]]))
  return(df)
})
names(Highvar_CNV_gene_corlist)<-names(cancer_list)
temp<-as.data.frame(t(as.data.frame(lapply(names(cancer_list),function(x){
  temp<-Highvar_CNV_gene_corlist[[x]]
  temp$patient<-x
  temp$stat<-CRCKUL_CNV_variability[[x]][temp$region]
  temp<-temp[which(temp$stat %in% c('lowvar','highvar')),c('gene','stat','Rval','patient')]
  return(t(temp))
}))))
temp$Rval<-as.numeric(temp$Rval)
ggplot(temp,aes(x=stat,y=Rval))+geom_boxplot(aes(fill=stat))+scale_fill_manual(values=c(lowvar='orchid1',highvar='purple'))+facet_wrap('patient')+xlab('CNV region grouped by variation')+ylab('R value of gene with located CNV region')+theme_bw()+theme(text = element_text(size = 25),axis.text.x = element_text(size=25))+ggprism::add_pvalue(temp %>% rstatix::group_by(patient) %>% rstatix::wilcox_test(Rval ~ stat) %>% rstatix::add_xy_position(),y.position = 0.9,label.size = 5)
#2G
CRCKUL_CNV_variability<-as.data.frame(lapply(names(cancer_list),function(x){
  var.limit=0.001
  temp<-cancer_list[[x]]@assays$inferCNV@data
  CNV_stat<-unlist(apply(temp,1,var))
  region_stat<-rep('noCNV',dim(temp)[1])
  CNV_region<-which(CRCKUL_CNV_ploidy[[x]] %in% c('amplified','deleted'))
  region_stat[CNV_region]<-'lowvar'
  region_stat[intersect(which(CNV_stat>var.limit),CNV_region)]<-'highvar'
  return(region_stat)
}))
colnames(CRCKUL_CNV_variability)<-names(cancer_list)
temp0<-lapply(names(cancer_list),function(x){
  pat<-subset(CRCKUL,cells=Cells(CRCKUL)[which(CRCKUL$patient==x&CRCKUL$celltype=='Cancer')])
  pat <- NormalizeData(pat, verbose = FALSE)
  stat_ratio<-apply(pat@assays$RNA@data,1,function(y){
    mean_val<-mean(y)
    if(mean_val==0){return(NA)}
    return(sd(y)/mean(y))
  })
  names(stat_ratio)<-rownames(pat@assays$RNA@data)
  stat_ratio<-stat_ratio[!is.na(stat_ratio)]
  return(stat_ratio)
})
names(temp0)<-names(cancer_list)
temp<-as.data.frame(lapply(names(cancer_list),function(x){
  stat_ratio<-temp0[[x]]
  hvg_site<-order(stat_ratio,decreasing = T)[1:2000]
  all_VAR<-CRCKUL_CNV_variability[[x]][cancer_list[[x]]@assays$RNA@meta.features[names(stat_ratio),'CNV_region']]
  stat_ratio<-as.numeric(table(all_VAR[hvg_site])[c('highvar','lowvar','noCNV')]/table(all_VAR)[c('highvar','lowvar','noCNV')])
  return(stat_ratio)
}))
colnames(temp)<-names(cancer_list)
temp$stat<-c('highvar','lowvar','noCNV')
temp<-reshape2::melt(temp,id.vars=c('stat'),value.name = 'ratio',variable.name = "patient")
ggplot(temp,aes(x=stat,y=ratio))+geom_boxplot(aes(fill=stat))+scale_fill_manual(values=c(noCNV='white',lowvar='orchid1',highvar='purple'))+geom_line(aes(group=patient),linetype='dashed',color='grey50')+ scale_y_continuous(labels = scales::percent)+xlab('CNV region grouped by variation')+ylab('High CV gene ratio')+theme_bw()+theme(text = element_text(size = 18),axis.text.x = element_text(size = 15))+ggprism::add_pvalue(rstatix::wilcox_test(temp, ratio ~ stat,paired = T) %>% rstatix::add_xy_position(x = "stat", dodge = 0.8),label.size = 5)
#S5A,S5B
temp<-subset(CRCKUL,cell=Cells(CRCKUL)[which((CRCKUL$celltype=='Epithelial'&CRCKUL$tissue=='N')|(CRCKUL$celltype=='Cancer'&CRCKUL$tissue=='T'))])
temp<-temp@meta.data[which(!is.na(temp$entropy)),]
df_p_val1 <- temp %>%
  rstatix::group_by(patient) %>%
  rstatix::wilcox_test(entropy ~ celltype) %>%
  rstatix::adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  rstatix::add_significance(p.col = "p.adj") %>% 
  rstatix::add_xy_position(x = "dose", dodge = 0.8) 
ggplot(temp)+geom_boxplot(aes(x=patient,y=entropy,fill=celltype),color='black')+scale_fill_manual(values=c('red','blue'))+scale_y_log10()+theme_bw()+theme(axis.text.x = element_text(angle = -60, hjust = 0,vjust=0.7),text = element_text(size = 20))+ ggprism::add_pvalue(df_p_val1,x=1:12,label.size = 5)
stem_list<-c("DNMT3B","PFAS","XRCC5","HAUS6","TET1","IGF2BP1","PLAA","TEX10","MSH6","DLGAP5","SKIV2L2","SOHLH2","RRAS2","PAICS","CPSF3","LIN28B","IPO5","BMPR1A","ZNF788","ASCC3","FANCB","HMGA2","TRIM24","ORC1","HDAC2","HESX1","INHBE","MIS18A","DCUN1D5","MRPL3","CENPH","MYCN","HAUS1","GDF3","TBCE","RIOK2","BCKDHB","RAD1","NREP","ADH5","PLRG1","ROR1","RAB3B","DIAPH3","GNL2","FGF2","NMNAT2","KIF20A","CENPI","DDX1","XXYLT1","GPR176","BBS9","C14orf166","BOD1","CDC123","SNRPD3","FAM118B","DPH3","EIF2B3","RPF2","APLP1","DACT1","PDHB","C14orf119","DTD1","SAMM50","CCL26","MED20","UTP6","RARS2","ARMCX2","RARS","MTHFD2","DHX15","HTR7","MTHFD1L","ARMC9","XPOT","IARS","HDX","ACTRT3","ERCC2","TBC1D16","GARS","KIF7","UBE2K","SLC25A3","ICMT","UGGT2","ATP11C","SLC24A1","EIF2AK4","GPX8","ALX1","OSTC","TRPC4","HAS2","FZD2","TRNT1","MMADHC","SNX8","CDH6","HAT1","SEC11A","DIMT1","TM2D2","FST","GBE1")
temp$stemness_score<-colMeans(CRCKUL$RNA@data[intersect(stem_list,rownames(CRCKUL$RNA@data)),rownames(temp)])
df_p_val1 <- temp %>%
  rstatix::group_by(patient) %>%
  rstatix::wilcox_test(stemness_score ~ celltype) %>%
  rstatix::adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  rstatix::add_significance(p.col = "p.adj") %>% 
  rstatix::add_xy_position(x = "dose", dodge = 0.8) 
ggplot(temp)+geom_boxplot(aes(x=patient,y=stemness_score,fill=celltype),color='black')+scale_fill_manual(values=c('red','blue'))+scale_y_log10()+theme_bw()+theme(axis.text.x = element_text(angle = -60, hjust = 0,vjust=0.7),text = element_text(size = 20))+ ggprism::add_pvalue(df_p_val1,x=1:12,label.size = 5)+ylab('Stemness score')
#S5C
temp<-as.data.frame(t(as.data.frame(lapply(names(cancer_list),function(x){
  temp<-Highvar_CNV_gene_corlist[[x]]
  temp$patient<-x
  temp1<-cancer_list[[x]]@assays$inferCNV@data
  temp$CNV_var<-unlist(apply(temp1,1,var))[temp$region]
  temp$stat<-CRCKUL_CNV_variability[[x]][temp$region]
  temp$cor_class<-'lowcor'
  temp$cor_class[which(temp$Rval>=0.3)]<-'highcor'
  return(t(temp))
}))))
temp$Rval<-as.numeric(temp$Rval)
temp$CNV_var<-as.numeric(temp$CNV_var)
temp<-temp[which(temp$Rval>=0),]
ggplot(temp)+geom_density(aes(x=Rval,fill=stat),alpha=0.6)+geom_vline(xintercept = 0.3,colour='red',linetype='dashed',size=1.5)+xlab('Correlation coefficient of gene and located region')+ylab('Gene number density')+scale_fill_manual('CNV region',values=c(noCNV='white',lowvar='orchid1',highvar='purple'))+theme_bw()+theme(text = element_text(size = 20))+geom_text(x=0.2, y=7, label='low',size=8,color='black')+geom_text(x=0.4, y=7, label='high',size=8,color='black')
