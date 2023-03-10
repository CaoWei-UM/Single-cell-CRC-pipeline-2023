#Author: Cao Wei
#update: 2023/03/09
#content: figure1

#1B
DimPlot(CRCKUL,group.by='tissue',label=F)+ggtitle('')+scale_color_manual(labels = c("Normal", "Tumor"),values=scales::hue_pal()(2))+theme(legend.text = element_text(size=20),axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank(),legend.position = 'none')+xlab('')+ylab('')
#1C
DimPlot(CRCKUL,group.by='patient')+ggtitle('')+theme(legend.text = element_text(size=20),axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank(),legend.position = 'none')+xlab('')+ylab('')
#1D
DimPlot(CRCKUL,group.by='celltype',label=T,label.size=6,repel=T)+theme(legend.position = "none",axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank())+ggtitle('')+xlab('')+ylab('')
#S1B
DimPlot(CRCKUL,group.by='sample',label=F)+theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank())+ggtitle('')
#S4A
temp<-subset(CRCKUL,cell=Cells(CRCKUL)[which(CRCKUL$celltype=='Epithelial'&CRCKUL$tissue=='N')])
DimPlot(temp,group.by='epithelial_subtype',label=F)+ggtitle('')
#S1A
ggplot(CRCKUL@meta.data,aes(x=sample,y=nCount_RNA,fill=sample))+geom_violin()+facet_wrap('patient')+ scale_y_log10()+ylab('UMI count')
ggplot(CRCKUL@meta.data,aes(x=sample,y=nFeature_RNA,fill=sample))+geom_violin()+facet_wrap('patient')+ scale_y_log10()+ylab('Gene count')
#S1C
ggplot(CRCKUL@meta.data,aes(x=patient,fill=tissue))+stat_count()+coord_flip()
#1E,S1D
DotPlot(CRCKUL,group.by = 'celltype',features = c('PTPRC','CD3D','CD3E','CD68','CSF1R','MS4A2','TPSAB1','CD79A','CD79B','IGHA1','IGHG1','PECAM1','VWF','COL1A1','DCN','EPCAM','MMP7','CEACAM6'),cols = c("grey", "red"),dot.scale = 10)+  theme(text = element_text(size = 20),axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5,size = 16),axis.text.y = element_text(size = 20))+ylab('Cell type')+xlab('Marker gene')
FeaturePlot(CRCKUL,c('CD3D','CD68','TPSAB1','CD79A','IGHA1','PECAM1','DCN','EPCAM','CEACAM6'),cols = c("grey", "red"))
#S4B
temp<-subset(CRCKUL,cell=Cells(CRCKUL)[which(CRCKUL$celltype=='Epithelial'&CRCKUL$tissue=='N')])
DotPlot(temp,group.by = 'epithelial_subtype',features =c("TMSB4X", "CCL5", "CREM", "EEF1A1", "CA4", "LYPD8", "MT1H", "OTOP2", "BEST4", "CA7", "AQP8", "SLC26A3", "GUCA2A", "CEACAM7", "CEACAM1", "CA1", "SELENBP1", "S100A6", "PHGR1", "KRT19", "MT1G", "MT2A", "TUBA1A", "PTMS", "PCSK1N", "SCGN", "CRYBA2", "UGT2B17", "ADH1C", "MGST1", "SPINK4", "ITLN1", "MUC2", "WFDC2", "TFF3", "SNHG5", "HES1", "TMSB10", "MLEC", "RARRES2", "NUPR1",'SH2D6','PLCG2','CCL23','PTPRC'),cols = c("grey", "red"),dot.scale = 10)+  theme(text = element_text(size = 20),axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5,size = 16),axis.text.y = element_text(size = 20))+ylab('Cell type')+xlab('Marker gene')
#1F,S1E
temp<-subset(CRCKUL,cell=Cells(CRCKUL)[which((CRCKUL$celltype=='Epithelial'&CRCKUL$tissue=='N')|(CRCKUL$celltype=='Cancer'&CRCKUL$tissue=='T'))])
library(ggprism)
df<-reshape2::melt(reshape2::dcast(temp@meta.data[,c('celltype','patient','nCount_RNA')],celltype~patient,fun.aggregate=mean))
colnames(df)<-c('celltype','patient','value')
p.values<-wilcox.test(df$value[which(df$celltype=='Epithelial')],df$value[which(df$celltype=='Cancer')],paired = T)$p.value
p.data<-data.frame(group1='Epithelial',group2='Cancer',p.adj=paste0('p.val=',signif(p.values,3)),y.position=max(df$nCount_RNA)*1.05)
ggplot(df,aes(x=celltype,y=value))+geom_boxplot(aes(fill=celltype))+geom_point()+geom_line(aes(group=patient))+theme_bw()+add_pvalue(p.data,label.size = 4)+ylab('Number of detected genes')
df<-reshape2::melt(reshape2::dcast(temp@meta.data[,c('celltype','patient','nFeature_RNA')],celltype~patient,fun.aggregate=mean))
colnames(df)<-c('celltype','patient','value')
p.values<-wilcox.test(df$value[which(df$celltype=='Epithelial')],df$value[which(df$celltype=='Cancer')],paired = T)$p.value
p.data<-data.frame(group1='Epithelial',group2='Cancer',p.adj=paste0('p.val=',signif(p.values,3)),y.position=max(df$nFeature_RNA)*1.05)
ggplot(df,aes(x=celltype,y=value))+geom_boxplot(aes(fill=celltype))+geom_point()+geom_line(aes(group=patient))+theme_bw()+add_pvalue(p.data,label.size = 4)+ylab('Number of detected UMIs')
#1G,S1F
temp<-subset(CRCKUL,cell=Cells(CRCKUL)[which(CRCKUL$celltype %in% c('Epithelial','Cancer'))],invert=T)
PlotSampleGroupPairBox<-function(measure){
  library(ggprism)
  raw_data<-temp@meta.data[,c('celltype',measure,'tissue','patient')]
  pairs<-unique(as.character(raw_data$tissue))
  df<-data.frame(value=tapply(raw_data[[measure]],factor(paste(raw_data$celltype,raw_data$patient,raw_data$tissue)),mean))
  df$group<-gsub('\\s.*','',rownames(df))
  df$pair<-gsub('.*\\s','',rownames(df))
  df$id<-sapply(strsplit(rownames(df),'\\s'),function(x){paste(x[1],x[2],sep = '.')})
  p.values<-sapply(unique(df$group),function(x){
    temp.df<-reshape2::dcast(df[which(df$group==x),c('value','pair','id')],id~pair,value.var = 'value')
    return(wilcox.test(temp.df[[pairs[1]]],temp.df[[pairs[2]]],paired = T)$p.value)
  })
  p.values<-p.adjust(p.values)
  p.data<-data.frame(group1=pairs[1],group2=pairs[2],p.adj=paste0('p.val=',signif(p.values,3)),y.position=max(df$value)*1.05,group=names(p.values))
  fig<-ggplot(df,aes(x=pair,y=value))+geom_boxplot(aes(fill=pair))+geom_point()+geom_line(aes(group=id))+facet_wrap('group',nrow = 1)+theme_bw()+add_pvalue(p.data,label.size = 4)
  return(fig)
}
PlotSampleGroupPairBox('nCount_RNA')+ylab('Number of detected UMIs')+xlab('')+scale_x_discrete(labels=c('N'='Normal tissue', 'T'='Tumor tissue'))+theme(text = element_text(size = 20),axis.title.x = element_blank(),axis.text.x = element_blank(),legend.position="none")
PlotSampleGroupPairBox('nFeature_RNA')+ylab('Number of detected genes')+xlab('')+scale_x_discrete(labels=c('N'='Normal tissue', 'T'='Tumor tissue'))+theme(text = element_text(size = 20),axis.title.x = element_blank(),axis.text.x = element_blank(),legend.position="none")
#1H,1I
CRCKUL_EpiCan<-subset(CRCKUL,cell=Cells(CRCKUL)[which((CRCKUL$celltype=='Epithelial'&CRCKUL$tissue=='N')|(CRCKUL$celltype=='Cancer'&CRCKUL$tissue=='T'))])
types<-c('Epithelial',"Undiff","Colonocytes","CT_colonocytes","Goblet","BEST4_OTOP2","EECs",'Cancer')
CRCKUL_EpiCan<-subset(CRCKUL_EpiCan,cell=Cells(CRCKUL_EpiCan)[which(CRCKUL_EpiCan$epithelial_subtype %in% types)])
CRCKUL_EpiCan<-NormalizeData(CRCKUL_EpiCan)
types<-c('Epithelial',"Undiff","Colonocytes","CT_colonocytes","Goblet","BEST4_OTOP2","EECs",'Cancer')
temp_list<-lapply(types,function(type){
  if(type=='Epithelial'){
    temp<-CRCKUL_EpiCan
    site_1<-which(temp$celltype==type)
    site_2<-which(temp$celltype=='Cancer')
  }else{
    temp<-subset(CRCKUL_EpiCan,cells=Cells(CRCKUL_EpiCan)[which(CRCKUL_EpiCan$epithelial_subtype %in% c(type,'Cancer'))])
    temp<-NormalizeData(temp)
    temp$epithelial_subtype<-factor(as.character(temp$epithelial_subtype),levels=c('Cancer',type))
    site_1<-which(temp$epithelial_subtype==type)
    site_2<-which(temp$epithelial_subtype=='Cancer')
  }
  ratio1<-rowSums(temp@assays$RNA@data[,site_1]>0)/length(site_1)
  ratio2<-rowSums(temp@assays$RNA@data[,site_2]>0)/length(site_2)
  kept.gene<-rownames(temp@assays$RNA@data)[which(ratio1>=0.5&ratio2>=0.5)]
  CV1<-apply(temp@assays$RNA@data[kept.gene,site_1],1,function(y){sd(y)/mean(y)})
  CV2<-apply(temp@assays$RNA@data[kept.gene,site_2],1,function(y){sd(y)/mean(y)})
  p.val<-wilcox.test(CV1,CV2,paired=T)$p.value
  return(list(p=p.val,CV=CV1,CCV=CV2))
})
names(temp_list)<-types
p.dat<-unlist(lapply(temp_list,function(x){x$p}))
dat<-lapply(names(temp_list),function(x){t(data.frame(celltype=x,CV=temp_list[[x]]$CV,CCV=temp_list[[x]]$CCV))})
dat<-as.data.frame(t(as.data.frame(dat)))
dat$CV<-as.numeric(dat$CV)
dat$CCV<-as.numeric(dat$CCV)
dat$celltype<-factor(dat$celltype,levels = c('Epithelial',setdiff(levels(CRCKUL_EpiCan$epithelial_subtype),'Cancer')))
dat$percent<-(1:dim(dat)[1])/dim(dat)[1]
dat2<-reshape2::melt(dat,id='celltype')
dat2$variable<-factor(dat2$variable,levels=c('CV','CCV'))
library(gghalves)
ggplot(dat,aes(x=celltype))+geom_half_violin(aes(y=CV,fill=celltype),colour='black',side = "l")+geom_half_violin(aes(y=CCV),colour='black',fill="blue",side = "r")+geom_boxplot(data=dat2,aes(x=celltype,y=value,fill=variable),width=0.1,position = position_dodge(width = 0.3))+geom_vline(xintercept = 1.5,color='grey40',linetype='dashed')+theme_bw()+xlab('Cell type')+ylab('Coefficient of variation')+theme(axis.text.x = element_blank(),axis.text.y = element_text(size=16),text = element_text(size=20),legend.position = 'none')+stat_compare_means(data=dat2,aes(x=celltype,y=value,fill=variable),method='wilcox.test',label = 'p.signif',paired = T,size=5)
ggplot(dat,aes(x=CV,y=percent))+geom_line(aes(color=group))+scale_y_continuous(labels = scales::percent)+ylab('Gene percent')+xlab('Coefficient of variation')+theme_bw()
#S4C,S4D
for(i in unqiue(CRCKUL$patient)){
  temp<-subset(CRCKUL_EpiCan,patient=i)
  temp<-NormalizeData(temp)
  types<-names(which(table(temp0$epithelial_subtype)>=100))#remove epithelial subtypes with low cell count
  site_1<-which(temp$celltype=='Epithelial')
  site_2<-which(temp$celltype=='Cancer')
  ratio1<-rowSums(temp@assays$RNA@data[,site_1]>0)/length(site_1)
  ratio2<-rowSums(temp@assays$RNA@data[,site_2]>0)/length(site_2)
  kept.gene<-rownames(temp@assays$RNA@data)[which(ratio1>=0.5&ratio2>=0.5)]
  CV1<-apply(temp@assays$RNA@data[kept.gene,site_1],1,function(y){sd(y)/mean(y)})
  CV2<-apply(temp@assays$RNA@data[kept.gene,site_2],1,function(y){sd(y)/mean(y)})
  p.val<-wilcox.test(CV1,CV2,paired=T)$p.value
  percent<-(1:length(CV1))/length(CV1)
  dat<-data.frame(CV=CV1,CCV=CV2,p.val=p.val,percent=percent)
  dat2<-reshape2::melt(dat,id='celltype')
  dat2$variable<-factor(dat2$variable,levels=c('CV','CCV'))
  library(gghalves)
  ggplot(dat2,aes(x=celltype,y=CV,fill=celltype))+geom_violin()+geom_boxplot(data=dat2,aes(x=celltype,y=value,fill=variable),width=0.1,position = position_dodge(width = 0.3))+theme_bw()+xlab('Cell type')+ylab('Coefficient of variation')+theme(axis.text.x = element_blank(),axis.text.y = element_text(size=16),text = element_text(size=20),legend.position = 'none')+stat_compare_means(method='wilcox.test',label = 'p.signif',paired = T,size=5)
  ggplot(dat2,aes(x=CV,y=percent))+geom_line(aes(color=celltype))+scale_y_continuous(labels = scales::percent)+ylab('Gene percent')+xlab('Coefficient of variation')+theme_bw()
}
#S2A,S3A,S3C
igGenes = c("IGHA1", "IGHA2", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHD", "IGHE", 
            "IGHM", "IGLC1", "IGLC2", "IGLC3", "IGLC4", "IGLC5", "IGLC6", "IGLC7", "IGKC","JCHAIN")
hgGenes = c("HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", 
            "HBZ")
PlotDiVolcano(object,group.by,ident.1,ident.2){
  DEG<-FindMarkers(object, group.by=group.by,ident.1 = ident.1,ident.2 = ident.2)
  DEG <- DEG[DEG$p_val_adj<p.value.threshold,]
  DEG<-DEG[which(!(rownames(DEG) %in% c(igGenes,hgGenes))),]
  DEG$pct.diff<-DEG$pct.1-DEG$pct.2
  DEG$gene_stat<-'non-significant'
  DEG$gene_stat[which(DEG$avg_log2FC>logFC&DEG$pct.diff>pct.diff)]<-ident.1
  DEG$gene_stat[which(DEG$avg_log2FC<(-logFC)&DEG$pct.diff<(-pct.diff))]<-ident.2
  point_colors<-c("#F8766D","grey","#00BFC4")
  names(point_colors)<-c(ident.1,ident.2, "non-significant")
  fig<-ggplot(DEG, aes(x=pct.diff, y=avg_log2FC)) +
    geom_point(aes(colour = gene_stat), size=1) +
    scale_colour_manual(values = point_colors) +
    scale_x_continuous(labels = scales::percent) +
    xlab('Percentage of cells difference') +
    ylab('Average log2 fold change') +
    theme_bw()+ theme(legend.title = element_blank())
  topn.label<-as.integer(topn.label)
  top_genes1<-DEG[which(DEG$gene_stat==ident1_name),]
  top_genes2<-DEG[which(DEG$gene_stat==ident2_name),]
  top_genes1<-head(rownames(top_genes1)[order(abs(top_genes1$pct.diff)*abs(top_genes1$avg_log2FC), decreasing=T)],topn.label)
  top_genes2<-head(rownames(top_genes2)[order(abs(top_genes2$pct.diff)*abs(top_genes2$avg_log2FC), decreasing=T)],topn.label)
  fig<-LabelPoints(plot=fig, points=c(top_genes1,top_genes2),repel = TRUE,xnudge = 0, ynudge = 0)
  return(fig)
}
temp<-subset(CRCKUL,cell=Cells(CRCKUL)[which((CRCKUL$celltype=='Epithelial'&CRCKUL$tissue=='N')|(CRCKUL$celltype=='Cancer'&CRCKUL$tissue=='T'))])
PlotDiVolcano(temp,'celltype','Epithelial','Cancer')
temp<-subset(CRCKUL,cell=Cells(CRCKUL)[which(CRCKUL$celltype=='Endothelial')])
PlotDiVolcano(temp,'tissue','N','T')
temp<-subset(CRCKUL,cell=Cells(CRCKUL)[which(CRCKUL$celltype=='Stromal')])
PlotDiVolcano(temp,'tissue','N','T')
