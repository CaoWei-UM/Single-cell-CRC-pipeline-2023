#Author: Cao Wei
#update: 2023/03/09
#content: data & pipeline

CRC1_N<-AutoRead10X('/yourpath/CRC1_N')
CRC1_P<-AutoRead10X('/yourpath/CRC1_P')
CRC1_T<-AutoRead10X('/yourpath/CRC1_T')
temp<-mergeSeurat(CRC1_N,CRC1_P,CRC1_T)
temp <- NormalizeData(temp, verbose = FALSE)
temp <- FindVariableFeatures(temp)
temp <- ScaleData(temp, verbose = FALSE)
temp <- RunPCA(temp, verbose = FALSE)
temp <- FindNeighbors(temp, reduction = "pca", dims = 1:20,k.param = 20)
temp <- FindClusters(temp, resolution = 1)
temp <- RunUMAP(temp, dims = 1:20,return.model = T)
temp<-subset(temp,idents=c(16,21,78,50,32,43,53,65,40,72,45,56,58,30,48,42,71,35,18,20,22,79,54,69,29,60,23,36,25,63,61,55),
       subset= nFeature_RNA>500&nCount_RNA>1000)
temp <- NormalizeData(temp, verbose = FALSE)
temp <- FindVariableFeatures(temp)
temp <- ScaleData(temp, verbose = FALSE)
temp <- RunPCA(temp, verbose = FALSE)
temp <- FindNeighbors(temp, reduction = "pca", dims = 1:20,k.param = 20)
temp <- FindClusters(temp, resolution = 0.8)
temp <- RunUMAP(temp, dims = 1:20,return.model = T)
CRC1<-MarkCelltype(temp,list(
  'Cancer'=c(1,2,3),
  'T_cell'=c(4,7,15,20),
  'B_cell'=c(0,18,19,22),
  'Mast_cell'=c(24),
  'Macrophage'=c(8,12,16),
  'Epithelial'=c(5,9,11,25),
  'Stromal'=c(14,23,26,6,21,13,17),
  'Endothelial'=c(10)
))
CRC<-merge(CRC1,CRC2,CRC3)
CRCKUL<-merge(CRC,KUL)

cancer_list<-lapply(unique(CRCKUL$patient),function(patient){NormalizeData(subset(CRCKUL,patient=patient))})
names(cancer_list)<-unique(CRCKUL$patient)

temp<-subset(CRCKUL,celltype='Epithelial'&tissue='N')
temp <- NormalizeData(temp, verbose = FALSE)
temp <- FindVariableFeatures(temp)
temp <- ScaleData(temp, verbose = FALSE)
temp <- RunPCA(temp, verbose = FALSE)
temp <- FindNeighbors(temp, reduction = "pca", dims = 1:20,k.param = 20)
temp <- FindClusters(temp, resolution = 1)
temp <- RunUMAP(temp, dims = 1:20,return.model = T)
temp<-MarkCelltype(temp,list(
  'BEST4_OTOP2'=c(13,22),
  'Colonocytes'=c(0,1,9,6,2),
  'CT_colonocytes'=c(7,8),
  'Goblet'=c(10,20,21),
  'EECs'=c(25),
  'Undiff'=c(12,5,23,4,17,3,19),
  'Immuno_like'=c(14,15),
  'PLCG2_positive'=c(18),
  'CCL23_positive'=c(11),
  'SH2D6_positive'=c(16)
),slotname = 'epithelial_subtype',override = T)
CRCKUL$epithelial_subtype<-NA
CRCKUL$epithelial_subtype[Cells(temp)]<-temp$epithelial_subtype

chr.use=paste0('chr',c(seq(1,22),'X'))
CNR_list <- list.files(path = '/yourpath/', pattern = '.cnr$')
CNR_file_list <- lapply(CNR_list, function(x){data.table::fread(file = paste0(dir,'/',x))})
metadata <- CNR_file_list[[1]][,c('chromosome','start','end','gene')]
CRCSZ04_DNA_CNR <- list()
CRCSZ04_DNA_CNR$metadata <- metadata
CRCSZ04_DNA_CNR$log2 <- as.data.frame(lapply(CNR_file_list,function(x){subset(x,select = 'log2')}))
CRCSZ04_DNA_CNR$depth <- as.data.frame(lapply(CNR_file_list,function(x){subset(x,select = 'depth')}))
CRCSZ04_DNA_CNR$weight <- as.data.frame(lapply(CNR_file_list,function(x){subset(x,select = 'weight')}))
for (i in c('log2','depth','weight')){names(CRCSZ04_DNA_CNR[[i]])<-CNR_list}

CRCKUL_EpiCan<-subset(CRCKUL,cells=Cells(CRCKUL)[which(CRCKUL$celltype=='Cancer'|(CRCKUL$celltype=='Epithelial'&CRCKUL$tissue=='N'))])
CRCSZ04_RNA_CNR<-subset(CRCKUL_EpiCan,patient='CRC1')
InferCNR<-function(object,window.size=51,cut.off=3.5,ref.cell=NULL,ref.shift=0.2,exp.limit=3,verbose=T,
                          ref.quantile=0.95,chr.use=paste0('chr',c(seq(1,22),'X')), step.len=1, band=F,
                          gene_list=Human_hg38.bed){
  exp.mat <- as.matrix(GetAssayData(object, slot = "counts"))
  #1. Remove genes not in gene list
  raw_gene_count<-dim(exp.mat)[1]
  used_ensembl <- intersect(gene_list$gene_id, rownames(exp.mat))
  used_symbol <- intersect(gene_list$Gene, rownames(exp.mat))
  if(length(used_ensembl)>=length(used_symbol)){
    exp.mat<-exp.mat[match(used_ensembl,rownames(exp.mat)),]
    gene_list<-gene_list[match(used_ensembl,gene_list$gene_id),]
  }else{
    exp.mat<-exp.mat[match(used_symbol,rownames(exp.mat)),]
    gene_list<-gene_list[match(used_symbol,gene_list$Gene),]
  }
  if(verbose){print(paste0(raw_gene_count-dim(exp.mat)[1],' genes removed due to not in gene order file.'))}
  #2. Calculate TPM matrix
  total_counts<-apply(exp.mat,2,sum)
  TPM<-1000000*sweep(exp.mat, 2, total_counts, `/`)
  #3. Remove genes which has low global expression level
  gene_reserved<-sort(which(log2(apply(TPM,1,mean)+1)>=cut.off))
  TPM<-TPM[gene_reserved,]
  gene_list<-gene_list[gene_reserved,]
  if(verbose){print(paste0(dim(exp.mat)[1]-length(gene_reserved),' genes removed due to mean TPM less than cut off'))}
  if(verbose){print(paste0('Total ',dim(TPM)[1],' genes reserved after filter.'))}
  #4. Calculate infercnv relative expression matrix
  Exp_matrix<-log2((TPM/10)+1)
  Exp_gene_mean<-apply(Exp_matrix,1,mean)
  Exp_matrix<-apply(Exp_matrix,2,function(x){x-Exp_gene_mean})
  over_limit_ratio<-length(which(Exp_matrix>exp.limit|Exp_matrix<(-exp.limit)))/(dim(Exp_matrix)[1]*dim(Exp_matrix)[2])
  if(verbose){print(paste0('Expression limit is ',exp.limit,' and ',round(100*over_limit_ratio,2),'% of all elements are over limit.'))}
  Exp_matrix <- replace(Exp_matrix,Exp_matrix>exp.limit,exp.limit)
  Exp_matrix <- replace(Exp_matrix,Exp_matrix<(-exp.limit),(-exp.limit))
  if(verbose){print("Prepare infercnv input file over.")}
  #5. Split inferCNV into chromosome and band
  used_chr<-intersect(unique(gene_list$Chr),chr.use)
  gene_pos<-split(1:length(gene_list$Chr),factor(gene_list$Chr))[used_chr]
  if(band){
    bands<-split(1:length(gene_list$Chr),factor(substring(gene_list$band,0,1)))[c('p','q')]
    gene_pos<-unlist(lapply(gene_pos,function(x){list(p=intersect(x,bands$p),q=intersect(x,bands$q))}),recursive=F)
  }
  for(pos in names(gene_pos)){
    if(length(gene_pos[[pos]])<window.size){
      warning(paste0('Genes in ', pos, ' is less than window.size: ', window.size,', so skip it.'))
      gene_pos[pos]<-NULL
    }
  }
  #6. Calculate inferCNV
  CNV_list <- lapply(gene_pos,function(x){
    Pos_exp_matrix<-Exp_matrix[sort(x),]
    Pos_gene_list<-gene_list[sort(x),]
    CNV_start_pos<-rollapply(Pos_gene_list$GeneStart, window.size, min, by=step.len)
    CNV_end_pos<-rollapply(Pos_gene_list$GeneEnd, window.size, max, by=step.len)
    Pos_CNV<-apply(Pos_exp_matrix,2,function(y){rollapply(y, window.size, mean, by=step.len)})
    rownames(Pos_CNV)<-1:dim(Pos_CNV)[1]
    CNV_result<-list(CNV=t(Pos_CNV),start_pos=CNV_start_pos,end_pos=CNV_end_pos)
    return(CNV_result)
  })
  CNV<-as.data.frame(lapply(CNV_list,function(x){x$CNV}))
  start_pos<-unlist(lapply(CNV_list,function(x){x$start_pos}))
  end_pos<-unlist(lapply(CNV_list,function(x){x$end_pos}))
  CNV_cell_mean<-apply(CNV,1,mean)
  CNV<-t(apply(CNV,2,function(x){x-CNV_cell_mean}))
  if(!is.null(ref.cell)){
    ref_cells<-intersect(ref.cell,colnames(CNV))
    if(length(ref_cells)!=length(ref.cell)){
      warning('Some ref.cell not in seurat object!')
    }
    Max_baseline<-apply(CNV[,ref_cells,drop=F],1,function(x){quantile(x,ref.quantile)})
    Min_baseline<-apply(CNV[,ref_cells,drop=F],1,function(x){quantile(x,(1-ref.quantile))})
    CNV_pos<-rownames(CNV)
    CNV<-apply(CNV,2,function(x){
      cnvi<-rep(0,length(x))
      out_max<-which(x>(Max_baseline+ref.shift))
      out_min<-which(x<(Min_baseline-ref.shift))
      cnvi[out_max]<-(x-Max_baseline)[out_max]
      cnvi[out_min]<-(x-Min_baseline)[out_min]
      return(cnvi)
    })
    rownames(CNV)<-CNV_pos
  }
  if(verbose){print("Calculate infercnv over.")}
  chromosome<-gsub('\\..*$','',rownames(CNV),perl=T)
  meta_feature<-data.frame(chromosome)
  if(band){meta_feature$band<-gsub('^(.*\\..*)\\..*$','\\1',rownames(CNV),perl=T)}
  meta_feature$start_pos<-start_pos
  meta_feature$end_pos<-end_pos
  rownames(meta_feature)<-rownames(CNV)
  if("inferCNV" %in% names(exp.mat@assays)){exp.mat[["inferCNV"]] <- NULL}
  object[["inferCNV"]] <- CreateAssayObject(counts = CNV)
  object@assays$inferCNV@meta.features <- meta_feature
  object$inferCNV_score<-NA
  object$inferCNV_score<-apply(CNV,2,function(x){mean(x^2)})
  return(object)
}
NormalizeCNR<-function(object, window.size=30, mad.threshold=2.5, probe.min=0,gamma=40){
  raw_CNV <- as.matrix(GetAssayData(object, assay='inferCNV', slot='counts'))
  raw_CNV_feature <- object@assays$inferCNV@meta.features
  loc <- round((raw_CNV_feature$start_pos+raw_CNV_feature$end_pos)/2)
  chr <- raw_CNV_feature$chromosome
  raw_CNV <- cbind.data.frame(loc,raw_CNV)
  chr_index <- unique(chr)
  raw_CNV <- cbind.data.frame(match(chr,chr_index),raw_CNV)
  CNR_smooth <- winsorize(data=raw_CNV, tau=mad.threshold, k=floor(window.size/2),
                          assembly='hg38',verbose=FALSE)
  CNR_segment <- multipcf(data=CNR_smooth, Y=raw_CNV, gamma=gamma, assembly='hg38',verbose=FALSE)
  # probe.min not use in this version
  seg <- CNR_segment[,-c(1:5)]
  info <- CNR_segment[,c(1:5)]
  info$chrom <- chr_index[info$chrom]
  seg_list <- list(meta=info,segments=seg)
  object@assays$inferCNV@data <- as.matrix(seg_list[['segments']])
  seg <- seg_list[['meta']]$n.probes
  segment_region <- unlist(lapply(1:length(seg), function(i){rep(i,seg[i])}))
  object@assays$inferCNV@meta.features$segment <- segment_region
  return(object)
}
RunCNVPCA<-function(object,method='length',dims=50){
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
  new_obj <- t(object*segment_weight)
  pca.result<-prcomp(new_obj, rank. = dims, center = F)
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
ClusterCNV<-function(object, slot='pca', method='euclidean', k=5, dims=10){
  if(slot=='pca'){
    data.use <- Embeddings(object[['CNVpca']])[,1:dims]
  }else if(slot=='umap'){
    data.use <- Embeddings(object[['CNVumap']])
  }else if(slot=='raw'){
    data.use <- t(object@assays$inferCNV@data)
  }else{
    stop('not support yet')
  }
  clusters <- cutree(hclust(dist(data.use, method=method),method = 'ward.D2'), k=k)
  clusters<-paste('CNV',clusters,sep='_')
  names(clusters)<-colnames(object)
  object$CNV_cluster<-clusters
  return(object)
}
CRCSZ04_RNA_CNR<-InferCNR(CRCSZ04_RNA_CNR,cut.off = 0.1,window.size = 51,ref.shift = 0.1,step.len = 10,ref.cell = names(which(CRCSZ04_RNA_CNR$celltype=='Epithelial')))
CRCSZ04_RNA_CNR<-NormalizeCNR(CRCSZ04_RNA_CNR,ploidy = F)
CRCSZ04_RNA_CNR<-RunCNVPCA(CRCSZ04_RNA_CNR)
CRCSZ04_RNA_CNR<-ClusterCNV(CRCSZ04_RNA_CNR)

library(CellChat)
library(patchwork)
used.patient<-"CRCSZ1","CRCSZ2","KUL01","KUL21","KUL24","KUL27","KUL30","KUL31")
intra_NT_Chat_list<-list()
for(patient in used.patient){
  c_list<-list()
  for(t in c('N','T')){
    cellchat<-subset(CRCKUL,patient=patient&tissue='N')
    cellchat<-subset(cellchat,celltype=names(which(table(cellchat$celltype)>50)))
    cellchat <- NormalizeData(cellchat, verbose = FALSE)
    if(t=='T'){cellchat$celltype[which(cellchat$celltype=='Cancer')]<-'Epithelial'}
    options(stringsAsFactors = FALSE)
    cellchat <- createCellChat(cellchat,group.by='celltype')
    cellchat@DB<-CellChatDB.human
    cellchat <- subsetData(cellchat)
    future::plan("multiprocess", workers = 8) 
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    cellchat <- projectData(cellchat, PPI.human)
    cellchat <- computeCommunProb(cellchat,type='triMean',trim=NULL)
    cellchat <- filterCommunication(cellchat, min.cells = 10)
    cellchat <- computeCommunProbPathway(cellchat)
    c_list[t] <- aggregateNet(cellchat)
  }
  group.new<-unique(c(levels(c_list[['N']]@idents),levels(c_list[['T']]@idents)))
  intra_NT_Chat_list[[patient]] <- mergeCellChat(list(liftCellChat(c_list[['N']], group.new),liftCellChat(c_list[['T']], group.new)), add.names=c('Normal','Tumor'))
}

intra_Chat_list_CP<-list()
for(patient in used.patient){
  cellchat<-subset(CRCKUL,patient=patient&tissue='T')
  cellchat<-subset(cellchat,celltype=names(which(table(cellchat$celltype)>50)))
  c_list<-list()
  for(t in unique(cellchat$CNV_cluster)){
    cellchat<-subset(cellchat,cell=Cells(cellchat)[which(cellchat$celltype=='Cancer'&cellchat$CNV_cluster!=t)],invert=T)
    cellchat <- NormalizeData(cellchat, verbose = FALSE)
    options(stringsAsFactors = FALSE)
    cellchat <- createCellChat(cellchat,group.by='celltype')
    cellchat@DB<-CellChatDB.human
    cellchat <- subsetData(cellchat)
    future::plan("multiprocess", workers = 8) 
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    cellchat <- projectData(cellchat, PPI.human)
    cellchat <- computeCommunProb(cellchat,type='triMean',trim=NULL)
    cellchat <- filterCommunication(cellchat, min.cells = 10)
    cellchat <- computeCommunProbPathway(cellchat)
    c_list[t] <- aggregateNet(cellchat)
  }
  group.new<-unique(unlist(lapply(unique(cellchat$CNV_cluster),function(subclone){levels(c_list[[subclone]]@idents)})))
  merge_list<-lapply(unique(cellchat$CNV_cluster),function(subclone){liftCellChat(c_list[[subclone]], group.new)})
  intra_Chat_list_CP[[patient]] <- mergeCellChat(merge_list, add.names=unique(cellchat$CNV_cluster))
}
