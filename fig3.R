#Author: Cao Wei
#update: 2023/03/09
#content: figure3

#3A
temp<-subset(CRCKUL,cells=Cells(CRCKUL)[which(CRCKUL$celltype %in% c('Epithelial','Cancer'))])
temp$type<-as.character(temp$celltype)
temp$type[which(temp$type=='Cancer')]<-paste0('Cancer.',temp$patient[which(temp$type=='Cancer')])
library(infercnv)
options(stringsAsFactors=T)
temp_dir<-file.path(tempdir(), stringi::stri_rand_strings(1,8))
dir.create(temp_dir)
anno_file<-cbind(rownames(temp@meta.data),temp$type)
data.table::fwrite(anno_file,paste0(temp_dir,'/infercnv_anno.loc'), quote=F, row.names=F, col.names=F,sep="\t")
count_file<-as.data.frame(temp@assays[[DefaultAssay(temp)]]@counts)
count_file<-count_file[order(rownames(count_file)),]
data.table::fwrite(count_file,paste0(temp_dir,'/infercnv.matrix'), quote=F,row.names=T,sep = "\t")
infercnv_obj <- CreateInfercnvObject(raw_counts_matrix=paste0(temp_dir,'/infercnv.matrix'),
                                     annotations_file=paste0(temp_dir,'/infercnv_anno.loc'),
                                     delim="\t",
                                     gene_order_file='/yourpath/gene_order_hg38.txt',
                                     ref_group_names='Epithelial')
infercnv_obj <- infercnv::run(infercnv_obj, out_dir=temp_dir,
                              cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                              cluster_by_groups=T, cluster_references =T,
                              denoise=T, HMM=F, HMM_type='i6', 
                              analysis_mode="samples", no_prelim_plot=T, 
                              num_threads=16, plot_steps = F,output_format='pdf')
#3B
mapCNR2Seurat <- function(CNR, Seurat_obj, CNV.ratio=NULL){
  Seurat_meta<-Seurat_obj@assays$inferCNV@meta.features
  chr<-Seurat_meta$chromosome
  start<-Seurat_meta$start_pos
  end<-Seurat_meta$end_pos
  InferCNV<-as.matrix(Seurat_obj@assays$inferCNV@counts)
  meta<-CNR$metadata
  match_pos<-lapply(1:length(chr),function(x){
    which(meta$chromosome==chr[x]&meta$end>start[x]&meta$start<end[x])
  })
  DNA_CNV<-t(as.data.frame(lapply(match_pos,function(x){colMeans(CNR$log2[x,])})))
  if(!is.numeric(CNV.ratio)){
    CNV.ratio<-Matrix::nnzero(InferCNV)/prod(dim(InferCNV))
    print(paste0("Estimated CNV percent is ",100*CNV.ratio,'%.'))
  }
  max_thes<-quantile(DNA_CNV,1-CNV.ratio/2)
  min_thes<-quantile(DNA_CNV,CNV.ratio/2)
  DNA_CNV <- replace(DNA_CNV,DNA_CNV>min_thes&DNA_CNV<max_thes,0)
  DNA_value<-as.numeric(DNA_CNV)
  RNA_value<-as.numeric(InferCNV)
  range_diff<-mean(abs(RNA_value[which(RNA_value!=0)]))/mean(abs(DNA_value[which(DNA_value!=0)]))
  DNA_CNV<-DNA_CNV*range_diff
  DNA_CNV <- replace(DNA_CNV,DNA_CNV>max(InferCNV),max(InferCNV))
  DNA_CNV <- replace(DNA_CNV,DNA_CNV<min(InferCNV),min(InferCNV))
  rownames(DNA_CNV)<-rownames(InferCNV)
  CNR$mapped<-as.data.frame(DNA_CNV)
  return(CNR)
}
RunCNVPCA<-function(object,dims=20){
  CNV_mat <- object@assays$inferCNV@data
  CNR_info<-object@assays$inferCNV@meta.features
  if(!all(c('segment','chromosome','start_pos','end_pos') %in% colnames(CNR_info))){stop('segment, chromosome, start_pos, end_pos must in meta.features in inferCNV assay')}
  CNV_list<-split(CNR_info,cut(CNR_info$segment,0:max(CNR_info$segment)))
  CNV_info <- data.frame(t(data.frame(lapply(CNV_list,function(x){
    if(length(unique(x$segment))>1||length(unique(x$chromosome))>1){
      stop('segment or chromosome not unique')
    }
    return(c(min(x$start_pos),max(x$end_pos),length(x$segment)))
  }))))
  colnames(CNV_info)<-c('start_pos','end_pos','n_probe')
  segment_weight <- CNV_info$end_pos-CNV_info$start_pos
  pca.results<-prcomp(t(CNV_mat*segment_weight), rank. = dims, center = F)
  sdev <- pca.results$sdev
  feature.loadings <- pca.results$rotation
  cell.embeddings <- pca.results$x
  rownames(x = feature.loadings) <- rownames(x = CNV_mat)
  colnames(x = feature.loadings) <- paste0('PC_', 1:dims)
  rownames(x = cell.embeddings) <- colnames(x = CNV_mat)
  colnames(x = cell.embeddings) <- colnames(x = feature.loadings)
  reduction.data <- CreateDimReducObject(
    embeddings = cell.embeddings,
    loadings = feature.loadings,
    assay = 'inferCNV',
    stdev = sdev,
    key = 'PC_'
  )
  object[['CNVpca']]<-reduction.data
  return(object)
}
CRC1_cancer_CNR<-list()
CNR_cells <- colnames(CRC1_DNA_CNR$log2)
kept <- match(names(grep('^subclone',CRC1_DNA_CNR$CNV_cluster,value = T)),CNR_cells)
CRC1_cancer_CNR$log2<-CRC1_DNA_CNR$log2[,kept]
CRC1_cancer_CNR$depth<-CRC1_DNA_CNR$depth[,kept]
CRC1_cancer_CNR$weight<-CRC1_DNA_CNR$weight[,kept]
CRC1_cancer_CNR$segment<-CRC1_DNA_CNR$segment[,kept]
CRC1_cancer_CNR$CNV_cluster<-CRC1_DNA_CNR$CNV_cluster[kept]
CRC1_cancer_CNR$ploidy_segment<-CRC1_DNA_CNR$ploidy_segment[kept,]
CRC1_cancer_CNR<-mapCNR2Seurat(CRC1_cancer_CNR,cancer_list$CRC1)$mapped
CRC1_cancer_DRNA <- CreateSeuratObject(counts = cbind(cancer_list$CRC1@assays$inferCNV@counts,CRC1_cancer_CNR),assay = 'inferCNV', project = 'CRC1_DNA')
CRC1_cancer_DRNA@assays$inferCNV@meta.features<-cancer_list$CRC1@assays$inferCNV@meta.features
raw_CNV <- as.matrix(GetAssayData(CRC1_cancer_DRNA, assay='inferCNV', slot='counts'))
raw_CNV_feature <- CRC1_cancer_DRNA@assays$inferCNV@meta.features
loc <- round((raw_CNV_feature$start_pos+raw_CNV_feature$end_pos)/2)
chr <- raw_CNV_feature$chromosome
raw_CNV <- cbind.data.frame(loc,raw_CNV)
chr_index <- unique(chr)
raw_CNV <- cbind.data.frame(match(chr,chr_index),raw_CNV)
CNR_smooth <- winsorize(data=raw_CNV, tau=2.5, k=15, assembly='hg38',verbose=FALSE)
CNR_segment <- multipcf(data=CNR_smooth, Y=raw_CNV, gamma=40, assembly='hg38',verbose=FALSE)
seg <- CNR_segment[,-c(1:5)]
info <- CNR_segment[,c(1:5)]
info$chrom <- chr_index[info$chrom]
seg_list <- list(meta=info,segments=seg)
CRC1_cancer_DRNA@assays$inferCNV@data <- as.matrix(seg_list[['segments']])
seg <- seg_list[['meta']]$n.probes
segment_region <- unlist(lapply(1:length(seg), function(i){rep(i,seg[i])}))
CRC1_cancer_DRNA@assays$inferCNV@meta.features$segment <- segment_region
CRC1_cancer_DRNA<-RunCNVPCA(CRC1_cancer_DRNA)
CRC1_cancer_DRNA<-FusionCNV::RunCNVUMAP(CRC1_cancer_DRNA,1:20)
CRC1_cancer_DRNA$CNV_subclone<-c(cancer_list$CRC1$CNV_cluster,CRC1_DNA_CNR$CNV_cluster)[Cells(CRC1_cancer_DRNA)]
CRC1_cancer_DRNA$dataset<-'RNA'
CRC1_cancer_DRNA$dataset[which(Cells(CRC1_cancer_DRNA) %in% names(CRC1_DNA_CNR$CNV_cluster))]<-'DNA'
CRC1_cancer_DRNA@reductions$CNVumap@cell.embeddings
temp<-cbind.data.frame(CRC1_cancer_DRNA@reductions$CNVumap@cell.embeddings,CRC1_cancer_DRNA@meta.data)
ggplot(temp,aes(x=CNVumap_1,y=CNVumap_2))+ggrastr::geom_point_rast(aes(color=CNV_subclone,size=dataset,stroke=dataset),raster.dpi = getOption("ggrastr.default.dpi", 300))+theme_bw()+scale_color_manual(values=c(ggsci::pal_startrek("uniform")(5)))+scale_size_manual(values=c(DNA=4,RNA=1))+scale_stroke_manual(values=c(DNA=1,RNA=1))
#S6A
for(i in cancer_list){DimPlot(i,group.by='CNV_cluster')}
#3C
library(slingshot)
clusterLabels <- as.character(cancer_list$CRC1@meta.data$CNV_cluster)
data <- Embeddings(cancer_list$CRC1$CNVumap)
slingshot_obj <- slingshot(data, clusterLabels, start.clus = 'subclone_1')
slingshot_obj@slingParams$group <- clusterLabels
dat<-cbind.data.frame(slingshot_obj@reducedDim[,1:2],slingshot_obj@slingParams$group)
colnames(dat)[3]<-'group'
dat$group <- factor(as.character(dat$group), levels = sort(unique(dat$group)))
curve_dat<-as.data.frame(t(as.data.frame(lapply(names(slingshot_obj@curves),function(x){
  c<-slingshot_obj@curves[[x]]
  return(t(cbind(c$s[c$ord,1:2],x)))
}))))
colnames(curve_dat)[3]<-'curves'
curve_dat[[1]]<-as.numeric(curve_dat[[1]])
curve_dat[[2]]<-as.numeric(curve_dat[[2]])
centers<-as.data.frame(apply(dat[,1:2],2,function(x){tapply(x,group,mean)}))
directions<-as.data.frame(unique(unlist(lapply(slingshot_obj@lineages,function(x){
  lapply(1:(length(x)-1),function(y){c(x[y],x[y+1])})
}),recursive = F)))
pos_s<-centers[as.character(directions[1,]),]
pos_e<-centers[as.character(directions[2,]),]
centers$name=levels(group)
names(color.use)<-c(levels(group),unique(curve_dat$curves))
ggplot()+ggrastr::geom_point_rast(data=dat,aes_string(x=colnames(dat)[1],y=colnames(dat)[2],color=group.by),size=0.5,raster.dpi = getOption("ggrastr.default.dpi", 300))+geom_path(data=curve_dat,aes_string(x=colnames(dat)[1],y=colnames(dat)[2],color='curves'),size=2)+geom_point(data=centers,aes_string(x=colnames(dat)[1],y=colnames(dat)[2]),color='grey70',size=10)+geom_segment(aes(x=pos_s[,1], y=pos_s[,2], xend=pos_e[,1], yend=pos_e[,2]), arrow = arrow(length=unit(.5, 'cm')),color='black', lwd=1)+ggrepel::geom_text_repel(data=centers,aes_string(x=colnames(dat)[1],y=colnames(dat)[2],label='name'),size=5)
#3D
stem_list<-c("DNMT3B","PFAS","XRCC5","HAUS6","TET1","IGF2BP1","PLAA","TEX10","MSH6","DLGAP5","SKIV2L2","SOHLH2","RRAS2","PAICS","CPSF3","LIN28B","IPO5","BMPR1A","ZNF788","ASCC3","FANCB","HMGA2","TRIM24","ORC1","HDAC2","HESX1","INHBE","MIS18A","DCUN1D5","MRPL3","CENPH","MYCN","HAUS1","GDF3","TBCE","RIOK2","BCKDHB","RAD1","NREP","ADH5","PLRG1","ROR1","RAB3B","DIAPH3","GNL2","FGF2","NMNAT2","KIF20A","CENPI","DDX1","XXYLT1","GPR176","BBS9","C14orf166","BOD1","CDC123","SNRPD3","FAM118B","DPH3","EIF2B3","RPF2","APLP1","DACT1","PDHB","C14orf119","DTD1","SAMM50","CCL26","MED20","UTP6","RARS2","ARMCX2","RARS","MTHFD2","DHX15","HTR7","MTHFD1L","ARMC9","XPOT","IARS","HDX","ACTRT3","ERCC2","TBC1D16","GARS","KIF7","UBE2K","SLC25A3","ICMT","UGGT2","ATP11C","SLC24A1","EIF2AK4","GPX8","ALX1","OSTC","TRPC4","HAS2","FZD2","TRNT1","MMADHC","SNX8","CDH6","HAT1","SEC11A","DIMT1","TM2D2","FST","GBE1")
cancer_list$CRC1$stem<-colMeans(cancer_list$CRC1@assays$RNA[stem_list,])
library(ggnewscale)
data<-as.matrix(cancer_list$CRC1@meta.data[,c('curve1','curve2'),drop=F])
data<-reshape2::melt(data)
data<-data[!is.na(data$value),]
colnames(data)<-c('id','curve','pseudotime')
data$value<-FetchFeature(object,stem)[data$id]
data$cluster<-object@meta.data[[group.by]][data$id]
ggplot(data,aes(x=pseudotime,y=value))+ geom_point(aes(color=cluster),size=1,alpha=0.3) + geom_smooth(formula=y~x ,method='loess',aes(color=curve)) + scale_color_manual('fitted curve',values=topo.colors(2))  + theme_bw()
#3E


#3F
CRC1_RNA_SNV<-load('SNV data got from mutect2')
CRC1_RNA_SNV$DP[,Cells(cancer_list$CRC1)]->CRC1_SNV_DP
CRC1_RNA_SNV$ALT[,Cells(cancer_list$CRC1)]->CRC1_SNV_ALT
covered_cell_count<-apply(CRC1_SNV_DP,1,function(x){sum(x>0)})
kept.site<-which(covered_cell_count>=10)#at least cover 10 cells
covered_cell_count<-covered_cell_count[kept.site]
CRC1_SNV_DP<-CRC1_SNV_DP[kept.site,]
CRC1_SNV_ALT<-CRC1_SNV_ALT[kept.site,]
alert_cell_count<-apply(CRC1_SNV_ALT,1,function(x){sum(x>0)})
kept.site<-which(alert_cell_count>=covered_cell_count*0.05)#at least mutate 5% cells among covered cells
CRC1_SNV_DP<-CRC1_SNV_DP[kept.site,]
CRC1_SNV_ALT<-CRC1_SNV_ALT[kept.site,]
CRC1_SNV<-CRC1_SNV_ALT/CRC1_SNV_DP
CRC1_SNV<-as.matrix(CRC1_SNV)
CRC1_SNV[is.nan(CRC1_SNV)]<-NA
CRC1_SNV[CRC1_SNV<0.01]<-0 #at least 1% reads in a cell are alert reads 
CRC1_SNV[CRC1_SNV>=0.01]<-1
temp<-lapply(paste('subclone',1:5,sep='_'),function(x){which(CRC1_cell_anno==x)})
CRC1_sub_pval<-as.data.frame(t(as.data.frame(apply(CRC1_SNV,1,function(x){
  CRC1_sub_p<-unlist(lapply(temp,function(y){
    judged<-x[y]
    others<-x[-y]
    judged<-judged[!is.na(judged)]
    others<-others[!is.na(others)]
    if(length(judged)>=3&length(others)>=3){
      return(wilcox.test(judged,others)$p.value)
    }
    return(NA)
  }))
  if(all(is.na(CRC1_sub_p))){return(c(NaN,NaN))}
  return(c(min(CRC1_sub_p,na.rm = T),which.min(CRC1_sub_p)))
}))))
CRC1_sub_pval<-CRC1_sub_pval[which(is.finite(CRC1_sub_pval$V1)),]
colnames(CRC1_sub_pval)<-c('p_val','sig.subclone')
CRC1_sub_pval$p.adj<-p.adjust(CRC1_sub_pval$V1)
CRC1_sub_p5<-rownames(CRC1_sub_pval[which(CRC1_sub_pval$p.adj<0.05),])

CRC1_sub_val<-as.data.frame(t(as.data.frame(apply(CRC1_SNV[rownames(CRC1_sub_pval[which(CRC1_sub_pval$p.adj<1),]),],1,function(x){
  CRC1_sub_value<-unlist(lapply(temp,function(y){
    judged<-x[y]
    judged<-judged[!is.na(judged)]
    return(ifelse(length(judged)>=5,sum(judged)/length(judged),NA))
  }))
  return(CRC1_sub_value)
}))))
colnames(CRC1_sub_val)<-paste('subclone',1:5,sep='_')
CRC1_sub_val<-CRC1_sub_val[order(unlist(apply(CRC1_sub_val,1,which.max))),]
temp<-CRC1_sub_val
rownames(temp)<-gsub('_[ATCG]_[ATCG]_\\.','',rownames(temp))
temp.anno<-data.frame(subclone=colnames(CRC1_sub_val))
rownames(temp.anno)<-colnames(CRC1_sub_val)
temp.anno.color<-list(subclone=ggsci::pal_startrek("uniform")(5))
names(temp.anno.color$subclone)<-temp.anno$subclone
pheatmap::pheatmap(replace(temp,is.na(temp),0),scale = 'row',cluster_rows = F,cluster_cols = F,show_rownames = T,show_colnames = F,annotation_col=temp.anno,fontsize = 16, fontsize_row = 12,annotation_colors = temp.anno.color)
#3G
PlotSNV<-function(object,SNV,site,reduction='CNVumap'){
  SNV_site<-SNV[site,Cells(object)]
  SNV_site[SNV_site==0]<-'Ref'
  SNV_site[SNV_site==1]<-'Alert'
  plot_obj<-data.frame(Dim1=object@reductions[[reduction]]@cell.embeddings[,1],
                       Dim2=object@reductions[[reduction]]@cell.embeddings[,2],
                       SNV=SNV_site)
  plot_obj<-plot_obj[rev(order(plot_obj$SNV)),]
  fig<-ggplot(data=plot_obj,aes(x = Dim1, y = Dim2))+
    geom_point(aes(color=SNV, size = SNV, alpha = SNV))+
    scale_size_manual(values = c(Alert=3,Ref=3),na.value=1)+
    scale_alpha_manual(values = c(Alert=0.9,Ref=0.9),na.value=0.1)+ 
    scale_color_manual(name='SNV state',values=c(Alert='springgreen3',Ref='steelblue3'))+
    theme_bw()+xlab(paste0(reduction,'_1'))+ylab(paste0(reduction,'_2'))
  return(fig)
}
temp<-apply(CRC1_SNV[rownames(CRC1_sub_val),],1,function(x){sum(x,na.rm = T)/length(x[which(!is.na(x))])})
temp1<-apply(CRC1_sub_val,1,function(x){max(x)-min(x)})
temp<-names(sort(temp1[which(temp<0.9&temp>0.1)],decreasing = T))
temp1<-grep('^chrM',temp,value = T,invert = T)
cowplot::plot_grid(plotlist = lapply(temp1,function(x){
  return(PlotSNV(CRC1,CRC1_SNV,x))
}),ncol=4)
#S6B
pheatmap::pheatmap(CRC1_cancer_CNR$segment,annotation_row = as.data.frame(CRC1_cancer_CNR$CNV_cluster))
#S6C
temp1<-as.data.frame(lapply(paste('subclone',1:5,sep = '_'),function(x){rowMeans(as.matrix(cancer_list$CRC1@assays$inferCNV@counts[,which(cancer_list$CRC1$CNV_cluster==x)]))}))
colnames(temp1)<-paste('subclone',1:5,sep = '_')
nearest_subclone<-cor(CRC1_cancer_CNR$mapped,temp1)
CRC1_cancer_subclone<-paste('subclone',apply(nearest_subclone,1,which.max),sep = '_')
names(CRC1_cancer_subclone)<-rownames(nearest_subclone)
temp<-data.frame(old=CRC1_DNA_CNR$CNV_cluster[which(CRC1_DNA_CNR$CNV_cluster!='Epithelial')[-5]],CNV_subclone=CRC1_cancer_subclone)
pheatmap::pheatmap(ppp,annotation_row = temp)
#S6D
temp<-CRC1_cancer_DRNA@assays$inferCNV@meta.features
temp<-cbind(temp,lapply(c('DNA','RNA'),function(origin){lapply(unique(CRC1_cancer_DRNA$CNV_cluster),function(subclone){
  rowMeans(CRC1_cancer_DRNA@assays$inferCNV@data[,which(CRC1_cancer_DRNA$origin==origin&CRC1_cancer_DRNA$CNV_cluster==subclone)])
})}))
temp<-reshape2::melt(temp)
ggplot(temp) +geom_step(aes(x=start_pos,y=CNV_score,color=group))+facet_wrap(~CNV_cluster)
#S6E
temp<-data.frame(DNA=round(100*table(CRC1_DNA_CNR$CNV_cluster[-c(1,8)])/88),
                 RNA=round(100*c('Non-cancer'=5308,table(cancer_list$CRC1$CNV_cluster))/10778))
colnames(temp)<-c('subclone','DNA','RNA')
temp<-reshape2::melt(temp)
colnames(temp)<-c('subclone','Sequence','Freq')
ggplot(temp, aes(x="", y=Freq, fill=subclone)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  geom_text_repel(aes(label = paste0(Freq, "%")), position = position_stack(vjust=0.5),box.padding=0) +
  labs(x = NULL, y = NULL, fill = NULL) +
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  scale_fill_manual(values=c('grey',ggsci::pal_startrek("uniform")(5)))+facet_wrap('Sequence')
#S6F
temp<-subset(CRCKUL,cell=Cells(CRCKUL)[which(CRCKUL$celltype=='Cancer')])
temp<-data.frame(inferCNV_score=tapply(temp$inferCNV_score,factor(paste(temp$patient,temp$CNV_cluster)),mean),
                 entropy=tapply(temp$entropy,factor(paste(temp$patient,temp$CNV_cluster)),mean))
temp$CNV_cluster<-gsub('.*\\s','',rownames(temp))
temp$patient<-gsub('\\s.*','',rownames(temp))
lm_result<-lm(inferCNV_score~entropy, temp)
r_value<-as.numeric(sign(coef(lm_result)[2])*sqrt(summary(lm_result)$r.squared))
p_value<-signif(summary(lm_result)$coefficients[2,4],2)
ggplot(temp)+geom_point(aes(x=inferCNV_score,y=entropy,color=patient,shape=CNV_cluster),size=6)+theme_bw()+geom_smooth(method = "lm",inherit.aes = F,aes(inferCNV_score, entropy),formula = y~x)+ geom_text(x=0.013, y=0.87, label=paste0('R=',signif(r_value,2),"\np.val=",p_value),size=10)+theme(text = element_text(size = 25),legend.box='horizontal')

