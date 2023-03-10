#Author: Cao Wei
#update: 2023/03/09
#content: figure4

#4A
CRCKUL_CNVcluster_DEG<-as.data.frame(t(as.data.frame(lapply(names(cancer_list),function(x){
  obj<-NormalizeData(cancer_list[[x]]) #use normlized data for each patient's cancer cell
  Idents(obj)<-obj$CNV_cluster
  temp<-FindAllMarkers(obj)
  temp<-temp[which(temp$p_val_adj<0.05),]
  temp<-temp[which(!(temp$gene %in% c(hgGenes,igGenes))),]
  dat<-temp[,c('gene','avg_log2FC','p_val_adj')]
  dat$patient<-x
  return(t(dat))
}))))
CRCKUL_CNVcluster_DEG$avg_log2FC<-as.numeric(CRCKUL_CNVcluster_DEG$avg_log2FC)
CRCKUL_CNVcluster_DEG$p_val_adj<-as.numeric(CRCKUL_CNVcluster_DEG$p_val_adj)
CRCKUL_CNVcluster_DEG<-CRCKUL_CNVcluster_DEG[which(!is.na(CRCKUL_CNVcluster_DEG$gene)),]
CRCKUL_CNVcluster_DEG$range<-as.numeric(CRCKUL_CNVcluster_DEG$range)
CRCKUL_CNVcluster_DEG<-CRCKUL_CNVcluster_DEG[which(CRCKUL_CNVcluster_DEG$range>log2(1.5)),]
CRCKUL_CNVcluster_DEG<-CRCKUL_CNVcluster_DEG[which(CRCKUL_CNVcluster_DEG$gene %in% unique(unlist(lapply(CRCKUL_CNVcluster_DEG_list,function(x){x$gene})))),]
temp<-ggplot(CRCKUL_CNVcluster_DEG, aes(x=avg_log2FC, y=-log10(p_val_adj))) +
  ggrastr::geom_point_rast(aes(colour = patient), size=1,raster.dpi = getOption("ggrastr.default.dpi", 300)) +
  ylab('-log10 adjust P-value')+geom_vline(aes(xintercept=log2(1.8)),color='red',linetype='dotted')+theme_bw()+theme(text = element_text(size = 20))+xlab('Average log2 fold change between subclones')
temp1<-rownames(CRCKUL_CNVcluster_DEG)[which(abs(CRCKUL_CNVcluster_DEG$avg_log2FC)>1.5)]
LabelPoints(plot=temp, points=temp1, label=CRCKUL_CNVcluster_DEG[temp1,'gene'],repel = TRUE,xnudge = 0, ynudge = 0,size=5)
#4B
CRCKUL_CNVcluster_DEG_GO<-as.data.frame(t(as.data.frame(lapply(names(CRCKUL_CNVcluster_DEG_list),function(x){
  temp<-CRCKUL_CNVcluster_DEG_list[[x]]
  temp$stat<-CRCKUL_CNV_variability[[x]][as.numeric(temp$region)]
  temp$patient<-x
  return(t(temp[,c('gene','region_cor','stat')]))
}))))
colnames(CRCKUL_CNVcluster_DEG_GO)<-c('gene','avg_log2FC','cor','stat','patient')
CRCKUL_CNVcluster_DEG_GO<-CRCKUL_CNVcluster_DEG_GO[which(!is.na(CRCKUL_CNVcluster_DEG_GO$cor)),]
CRCKUL_CNVcluster_DEG_GO$cor<-as.numeric(CRCKUL_CNVcluster_DEG_GO$cor)
CRCKUL_CNVcluster_DEG_GO$avg_log2FC<-as.numeric(CRCKUL_CNVcluster_DEG_GO$avg_log2FC)
CRCKUL_CNVcluster_DEG_GO$cor_level<-'Low'
CRCKUL_CNVcluster_DEG_GO$cor_level[which(CRCKUL_CNVcluster_DEG_GO$cor>0.3)]<-'High'
CRCKUL_CNVcluster_DEG_GO$cor_level<-factor(CRCKUL_CNVcluster_DEG_GO$cor_level,levels=c('Low','High'))
ggplot(CRCKUL_CNVcluster_DEG_GO,aes(x=cor,fill=stat))+geom_histogram(color='black')+geom_vline(aes(xintercept=0.3),colour='red',linetype='dashed',size=1.5)+scale_fill_manual(name='CNV region',values = c(noCNV='white',lowvar='orchid1',highvar='purple'))+xlab('R value of gene and region')+ylab('CNV subclones\' DEG count')+theme_bw()+theme(text = element_text(size = 20))+geom_text(x=0.2, y=20, label='low',size=8,color='black')+geom_text(x=0.4, y=20, label='high',size=8,color='black')
#4C
library(clusterProfiler)
library(org.Hs.eg.db)
high_cor_gene<-unique(CRCKUL_CNVcluster_DEG_GO$gene[which(CRCKUL_CNVcluster_DEG_GO$cor_level=='High')])
low_cor_gene<-unique(CRCKUL_CNVcluster_DEG_GO$gene[which(CRCKUL_CNVcluster_DEG_GO$cor_level=='Low'&CRCKUL_CNVcluster_DEG_GO$avg_log2FC>log2(1.8))])
temp1<-list(high_cor=bitr(high_cor_gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID,
            low_cor=bitr(low_cor_gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID)
temp.gomf<-compareCluster(temp1,fun='enrichGO',qvalueCutoff=0.05,OrgDb='org.Hs.eg.db',ont = "MF")
dotplot(temp.gomf,showCategory=30,font.size = 12)
#4D
ggplot(CRCKUL_CNVcluster_DEG_GO[which(!(CRCKUL_CNVcluster_DEG_GO$cor_level=='Low'&CRCKUL_CNVcluster_DEG_GO$avg_log2FC<=log2(1.8))),],aes(x=cor_level,fill=stat))+stat_count(color='black')+scale_fill_manual(name='CNV region',values = c(noCNV='white',lowvar='orchid1',highvar='purple'))+xlab('Genes grouped by correlation coefficient value')+ylab('CNV subclones\' DEG count')+theme_bw()+theme(text = element_text(size = 20))
#4E,4F
temp<-CRCKUL_CNVcluster_DEG_GO_involved_patients[which(CRCKUL_CNVcluster_DEG_GO_involved_patients$pathway %in% temp.gomf@compareClusterResult$Description[which(temp.gomf@compareClusterResult$Cluster=='low_cor')]),]
temp$pathway<-factor(temp$pathway,levels=names(sort(table(temp$pathway))))
ggplot(temp,aes(x=pathway))+stat_count(aes(fill=Patient),color='black')+theme_bw()+scale_fill_manual(values=scales::hue_pal()(12)[match(unique(temp$Patient),names(cancer_list))])+theme(text = element_text(size = 20),axis.text.x = element_text(size=16))+xlab('Involved pathway')+ylab('Patient count')+coord_flip()

