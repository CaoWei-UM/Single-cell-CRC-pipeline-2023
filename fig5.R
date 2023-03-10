#Author: Cao Wei
#update: 2023/03/09
#content: figure5

#5A
temp.used<-c("Endothelial",'Epithelial',"Macrophage","Plasma","Stromal","T_cell")
netVisual_circle(get_total_weight_commu(intra_NT_Chat_list,'Normal')[temp.used,temp.used],edge.width.max=2,alpha.edge = 1,color.use = scales::hue_pal()(9)[c(6,8,2,5,7,1)],vertex.label.cex = 2)
dev.off()
temp.used<-get_total_weight_commu(intra_NT_Chat_list,'Tumor')[temp.used,temp.used]
colnames(temp.used)[2]<-'Cancer'
rownames(temp.used)[2]<-'Cancer'
netVisual_circle(temp.used,edge.width.max=2,alpha.edge = 1,color.use = scales::hue_pal()(9)[c(6,9,2,5,7,1)],vertex.label.cex = 2)
#5B
source <- unique(unlist(lapply(intra_NT_Chat_list,function(x){return(levels(x@idents[[1]]))})))
target <- unique(unlist(lapply(intra_NT_Chat_list,function(x){return(levels(x@idents[[1]]))})))
df<-lapply(intra_NT_Chat_list,function(x){
  samples<-levels(x@meta$datasets)
  data<-x@net[[samples[1]]]$weight
  group1<-intersect(source,rownames(data))
  group2<-intersect(target,rownames(data))
  data<-reshape2::melt(data[group1,group2])
  type.pairs<-paste0(data$Var1,'->',data$Var2)
  data1<-x@net[[samples[2]]]$weight
  data1<-reshape2::melt(data1[group1,group2])
  data<-data.frame(c1=data$value, c2=data1$value, row.names=type.pairs)
  sums<-data$c1+data$c2
  data$c1<-100*data$c1/sums
  data$c2<-100*data$c2/sums
  colnames(data)<-samples
  data$pair<-rownames(data)
  data<-t(reshape2::melt(data))
  return(data)
})
names(df)<-names(intra_NT_Chat_list)
df<-as.data.frame(t(as.data.frame(df)))
colnames(df)<-c('pair','sample','interaction_score')
df<-df[which(!is.na(df$interaction_score)),]
df$interaction_score<-as.numeric(df$interaction_score)
samples<-unique(df$sample)
used.pair<-unique(df$pair)
comp<-as.data.frame(t(as.data.frame(lapply(used.pair,function(x){
  pos.1<-which(df$pair==x&df$sample==samples[1])
  pos.2<-which(df$pair==x&df$sample==samples[2])
  p.val<-wilcox.test(x=df$interaction_score[pos.1],y=df$interaction_score[pos.2])$p.value
  mean<-abs(mean(df$interaction_score[pos.1])-mean(df$interaction_score[pos.2]))
  return(c(p.val,mean))
}))))
colnames(comp)<-c('pval','mean')
rownames(comp)<-used.pair
comp<-comp[which(comp$pval<0.05),]
final.pairs<-rownames(comp)[order(comp$mean,decreasing = T)][1:min(10,dim(comp)[1])]
df<-df[which(df$pair %in% final.pairs),]
df$pair<-factor(df$pair,levels=final.pairs)
ggplot(df,aes(x=pair,y=interaction_score))+geom_boxplot(aes(fill=sample),outlier.size=0.5)+geom_point(aes(fill=sample),position = position_jitterdodge(jitter.width = 0.1,dodge.width = 0.7),size=0.5)+theme_bw()+xlab('Source to Target communication')+ylab('Normalized interaction score')+theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5))+ggprism::add_pvalue(df %>% rstatix::group_by(pair) %>% rstatix::pairwise_wilcox_test(interaction_score ~ sample) %>% rstatix::add_xy_position(x = "pair", step.increase = 0,dodge = 0.8),xmin = "xmin", xmax = "xmax",label = "p.adj.signif",y.position = 105,tip.length = 0,label.size = 5)
#5C
#GO:0005126:cytokine receptor binding
#GO:0005125:cytokine activity
temp<-get('GO:0005126', org.Hs.egGO2ALLEGS) %>% mget(org.Hs.egSYMBOL) %>% unlist()
temp1<-get('GO:0005125', org.Hs.egGO2ALLEGS) %>% mget(org.Hs.egSYMBOL) %>% unlist()
temp<-unique(c(temp,temp1))
temp1<-CellChatDB.human
cytokine_pathways<-unique(temp1$interaction$pathway_name[unlist(lapply(temp1$interaction$interaction_name,function(x){any(strsplit(x,'_')[[1]] %in% temp)}))])
cytokine_LR<-unique(temp1$interaction$interaction_name[unlist(lapply(temp1$interaction$interaction_name,function(x){any(strsplit(x,'_')[[1]] %in% temp)}))])
PlotWeightBox(intra_NT_Chat_list,cytokine_LR)+ylab('Interaction weight of cytokine')
dat<-as.data.frame(t(as.data.frame(lapply(names(intra_NT_Chat_list),function(patient){
  patient.weight<-lapply(names(intra_NT_Chat_list[[patient]]@net),function(sample){
    sample.net<-intra_NT_Chat_list[[patient]]@net[[sample]]
    used.LR<-intersect(cytokine_LR,dimnames(sample.net$prob)[[3]])
    LR.weight<-sum(unlist(lapply(used.LR,function(x){
      temp.value<-as.numeric(sample.net$prob[,,x])
      temp.pval<-as.numeric(sample.net$pval[,,x])
      weight<-sum(temp.value[which(temp.pval<0.05)])
      return(weight)
    })))
    return(c(sample,LR.weight,patient))
  })
}))))
colnames(dat)<-c('sample','weight','patient')
dat$weight<-as.numeric(dat$weight)
ggplot(dat,aes(x=sample,y=weight))+geom_boxplot(aes(fill=sample),outlier.size=0.5)+geom_line(aes(group=patient),linetype='dashed',color='grey50')+theme_bw()+ylab('Weight of inferred interactions')+theme(legend.position = 'none')+ggprism::add_pvalue(dat %>% rstatix::wilcox_test(weight ~ sample,paired = T) %>% rstatix::add_xy_position(x = "sample",dodge = 0.8),xmin = "xmin", xmax = "xmax",label = "p={p}",tip.length = 0,label.size = 5)
#5D
dat<-as.data.frame(t(as.data.frame(lapply(names(intra_NT_Chat_list),function(patient){
  lapply(names(intra_NT_Chat_list[[patient]]@net),function(sample){
    sample.netP<-intra_NT_Chat_list[[patient]]@netP[[sample]]
    lapply(cytokine_pathways,function(path){
      if(path %in% sample.netP$pathways){
        path.weight<-sum(as.numeric(sample.netP$prob[,,path]))
      }else{path.weight<-0}
      return(c(sample,path.weight,patient,path))
    })
  })
}))))
colnames(dat)<-c('sample','weight','patient','pathway')
dat$weight<-as.numeric(dat$weight)
p.dat<-dat %>% rstatix::group_by(pathway) %>% rstatix::pairwise_wilcox_test(weight ~ sample, paired=T)
p.dat<-p.dat[which(p.dat$p.adj<0.05),]
pathway.order<-order(p.dat$p.adj)
dat<-dat[which(dat$pathway %in% p.dat$pathway),]
ggplot(dat,aes(x=pathway,y=weight))+geom_boxplot(aes(fill=sample),outlier.size=0.5)+geom_point(aes(fill=sample),position = position_jitterdodge(jitter.width = 0.1,dodge.width = 0.7),size=0.5)+theme_bw()+ylab('Interaction weight')+theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5))+ggprism::add_pvalue(p.dat %>% rstatix::add_xy_position(x = "pathway",dodge = 0.8), x = 'x',label = "p.adj.signif",label.size = 5)
#5E,5F
for(pathway in c('MIF','SPP1')){
  dat<-bind_rows(unlist(lapply(names(intra_NT_Chat_list),function(patient){
    lapply(names(intra_NT_Chat_list[[patient]]@net),function(sample){
      df<-intra_NT_Chat_list[[patient]]@netP[[sample]]
      if(!(pathway %in% df$pathways)){return(NULL)}
      sample.source<-df$prob['Epithelial',,pathway]
      sample.target<-df$prob[,'Epithelial',pathway]
      df<-data.frame(patient=patient,sample=sample,celltype=rownames(df$prob),source=sample.source,target=sample.target)
      return(df)
    })
  }),recursive = F))
  dat<-reshape2::melt(dat,id.vars=c('patient', 'sample', 'celltype'),value.name = 'prob',variable.name='direction')
  dat$id<-paste(dat$patient,dat$sample,sep='.')
  dat$direction<-paste('Epithelial','as',dat$direction)
  ggplot(dat,aes(x=id,y=celltype))+geom_tile(aes(fill=prob),color='black')+facet_wrap('direction')+ylab(paste('Interaction score between cell types with','Epithelial','in',pathway))+xlab('Patients and samples')+theme_bw()+theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5))+scale_fill_gradient(low = 'white',high = 'red',na.value = 'white')
}
#5G,5H
temp<-temp_get_subclone_interaction_strength(intra_Chat_list_CP,'MIF',source = 'measured',target = c('Macrophage','T_cell'))
temp$prob<-as.numeric(temp$prob)
temp$inferCNV_score<-temp_subclone_inferCNV_score[paste(temp$patient,temp$subclone)]
ggplot(temp,aes(x=inferCNV_score,y=prob))+geom_point(aes(color=patient),size=5)+scale_color_manual('Patient',values=temp.color)+geom_line(aes(group=patient),linetype='dashed',color='grey50')
#S9A
ppp<-intra_NT_Chat_list
for(i in names(ppp)){
  temp.pos<-which(colnames(ppp[[i]]@net$Normal$weight)=='Epithelial')
  colnames(ppp[[i]]@net$Tumor$weight)[temp.pos]<-'Epithelial/Cancer'
  rownames(ppp[[i]]@net$Tumor$weight)[temp.pos]<-'Epithelial/Cancer'
  colnames(ppp[[i]]@net$Normal$weight)[temp.pos]<-'Epithelial/Cancer'
  rownames(ppp[[i]]@net$Normal$weight)[temp.pos]<-'Epithelial/Cancer'
  temp.pos<-which(levels(ppp[[i]]@idents[[1]])=='Epithelial')
  levels(ppp[[i]]@idents[[1]])[temp.pos]<-'Epithelial/Cancer'
}
PlotDiffInteraction<-function (cellchat_list, comparison = c(1, 2), slot.name="net", measure, color.use = NULL, color.edge = c("#b2182b", "#2166ac"), title.name = NULL, sources.use = NULL, targets.use = NULL, remove.isolate = FALSE, top = 1, weight.scale = FALSE, vertex.weight = 20, vertex.weight.max = NULL, 
                               vertex.size.max = 15, vertex.label.cex = 1, vertex.label.color = "black", 
                               edge.weight.max = NULL, edge.width.max = 8, alpha.edge = 0.6, 
                               label.edge = FALSE, edge.label.color = "black", edge.label.cex = 0.8, 
                               edge.curved = 0.2, shape = "circle", layout = in_circle(), 
                               margin = 0.2, arrow.width = 1, arrow.size = 0.2) {
  options(warn = -1)
  #measure <- match.arg(measure)
  if(measure %in% c("count", "weight", "count.merged", "weight.merged")){
    measure.type<-'normal'
  }else{
    measure.type<-'single'
  }
  all.celltype<-lapply(cellchat_list,function(x){levels(x@idents[[1]])})
  names(all.celltype)<-names(cellchat_list)
  used.celltype<-names(which(table(unlist(all.celltype))>=length(cellchat_list)/2))
  net.diff <- matrix(0,ncol = length(used.celltype), nrow = length(used.celltype))
  colnames(net.diff)<-used.celltype
  rownames(net.diff)<-used.celltype
  for(i in names(cellchat_list)){
    kept.celltype<-intersect(all.celltype[[i]],used.celltype)
    if(measure.type=='normal'){
      obj1 <- methods::slot(cellchat_list[[i]], slot.name)[[comparison[1]]][[measure]][kept.celltype,kept.celltype]
      obj2 <- methods::slot(cellchat_list[[i]], slot.name)[[comparison[2]]][[measure]][kept.celltype,kept.celltype]
    }else if(measure.type=='single'){
      temp1<-methods::slot(cellchat_list[[i]], slot.name)[[comparison[1]]]$prob
      temp2<-methods::slot(cellchat_list[[i]], slot.name)[[comparison[2]]]$prob
      if(measure %in% dimnames(temp1)[[3]]){obj1 <- temp1[kept.celltype,kept.celltype,measure]}else{obj1 <- 0}
      if(measure %in% dimnames(temp2)[[3]]){obj2 <- temp2[kept.celltype,kept.celltype,measure]}else{obj2 <- 0}
    }else{stop('check para measure.type')}
    net.diff[kept.celltype,kept.celltype]<-net.diff[kept.celltype,kept.celltype]+obj2-obj1
  }
  if (measure %in% c("count", "count.merged")) {
    if (is.null(title.name)) {
      title.name = "Differential number of interactions"
    }
  }
  else if (measure %in% c("weight", "weight.merged")) {
    if (is.null(title.name)) {
      title.name = "Differential interaction strength"
    }
  }else{
    if (is.null(title.name)) {
      title.name = paste0("Differential interaction of ",measure," pathway")
    }
  }
  net <- net.diff
  if ((!is.null(sources.use)) | (!is.null(targets.use))) {
    df.net <- reshape2::melt(net, value.name = "value")
    colnames(df.net)[1:2] <- c("source", "target")
    if (!is.null(sources.use)) {
      if (is.numeric(sources.use)) {
        sources.use <- rownames(net.diff)[sources.use]
      }
      df.net <- subset(df.net, source %in% sources.use)
    }
    if (!is.null(targets.use)) {
      if (is.numeric(targets.use)) {
        targets.use <- rownames(net.diff)[targets.use]
      }
      df.net <- subset(df.net, target %in% targets.use)
    }
    cells.level <- rownames(net.diff)
    df.net$source <- factor(df.net$source, levels = cells.level)
    df.net$target <- factor(df.net$target, levels = cells.level)
    df.net$value[is.na(df.net$value)] <- 0
    net <- tapply(df.net[["value"]], list(df.net[["source"]], 
                                          df.net[["target"]]), sum)
  }
  if (remove.isolate) {
    idx1 <- which(Matrix::rowSums(net) == 0)
    idx2 <- which(Matrix::colSums(net) == 0)
    idx <- intersect(idx1, idx2)
    net <- net[-idx, ]
    net <- net[, -idx]
  }
  net[abs(net) < stats::quantile(abs(net), probs = 1 - top)] <- 0
  g <- graph_from_adjacency_matrix(net, mode = "directed", 
                                   weighted = T)
  edge.start <- igraph::ends(g, es = igraph::E(g), names = FALSE)
  coords <- layout_(g, layout)
  if (nrow(coords) != 1) {
    coords_scale = scale(coords)
  }
  else {
    coords_scale <- coords
  }
  if (is.null(color.use)) {
    color.use = scPalette(length(igraph::V(g)))
  }
  if (is.null(vertex.weight.max)) {
    vertex.weight.max <- max(vertex.weight)
  }
  vertex.weight <- vertex.weight/vertex.weight.max * vertex.size.max + 
    5
  loop.angle <- ifelse(coords_scale[igraph::V(g), 1] > 0, -atan(coords_scale[igraph::V(g), 
                                                                             2]/coords_scale[igraph::V(g), 1]), pi - atan(coords_scale[igraph::V(g), 
                                                                                                                                       2]/coords_scale[igraph::V(g), 1]))
  igraph::V(g)$size <- vertex.weight
  igraph::V(g)$color <- color.use[igraph::V(g)]
  igraph::V(g)$frame.color <- color.use[igraph::V(g)]
  igraph::V(g)$label.color <- vertex.label.color
  igraph::V(g)$label.cex <- vertex.label.cex
  if (label.edge) {
    igraph::E(g)$label <- igraph::E(g)$weight
    igraph::E(g)$label <- round(igraph::E(g)$label, digits = 1)
  }
  igraph::E(g)$arrow.width <- arrow.width
  igraph::E(g)$arrow.size <- arrow.size
  igraph::E(g)$label.color <- edge.label.color
  igraph::E(g)$label.cex <- edge.label.cex
  igraph::E(g)$color <- ifelse(igraph::E(g)$weight > 0, color.edge[1], 
                               color.edge[2])
  igraph::E(g)$color <- grDevices::adjustcolor(igraph::E(g)$color, 
                                               alpha.edge)
  igraph::E(g)$weight <- abs(igraph::E(g)$weight)
  if (is.null(edge.weight.max)) {
    edge.weight.max <- max(igraph::E(g)$weight)
  }
  if (weight.scale == TRUE) {
    igraph::E(g)$width <- 0.3 + igraph::E(g)$weight/edge.weight.max * 
      edge.width.max
  }
  else {
    igraph::E(g)$width <- 0.3 + edge.width.max * igraph::E(g)$weight
  }
  if (sum(edge.start[, 2] == edge.start[, 1]) != 0) {
    igraph::E(g)$loop.angle[which(edge.start[, 2] == edge.start[, 
                                                                1])] <- loop.angle[edge.start[which(edge.start[, 
                                                                                                               2] == edge.start[, 1]), 1]]
  }
  radian.rescale <- function(x, start = 0, direction = 1) {
    c.rotate <- function(x) (x + start)%%(2 * pi) * direction
    c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
  }
  label.locs <- radian.rescale(x = 1:length(igraph::V(g)), 
                               direction = -1, start = 0)
  label.dist <- vertex.weight/max(vertex.weight) + 2
  plot(g, edge.curved = edge.curved, vertex.shape = shape, 
       layout = coords_scale, margin = margin, vertex.label.dist = label.dist, 
       vertex.label.degree = label.locs, vertex.label.family = "Helvetica", 
       edge.label.family = "Helvetica")
  if (!is.null(title.name)) {
    text(0, 1.5, title.name, cex = 1.1)
  }
  gg <- recordPlot()
  return(gg)
}
PlotDiffInteraction(ppp, weight.scale = T, measure = "weight",color.edge=rev(scales::hue_pal()(2)),alpha.edge = 1)
#S9C,S9D
netVisual_bubble(intra_NT_Chat_list,source.use = c('Macrophage'),target.use = c('Epithelial'),signaling = c('SPP1'))
netVisual_bubble(intra_NT_Chat_list,source.use = c('Epithelial'),target.use = c('Macrophage','T_cell'),signaling = c('MIF'))
#S9E,S9F,S9G,S9H
PlotRankNetPairwise<-function(cellchat_list,group.1=NULL,group.2=NULL,slot.name='netP',topn=10,duplex=T,p.val=0.05,tol=0.5){
  if(is.null(group.1)){group.1 <- unique(unlist(lapply(cellchat_list,function(x){return(levels(x@idents[[1]]))})))}
  if(is.null(group.2)){group.2 <- unique(unlist(lapply(cellchat_list,function(x){return(levels(x@idents[[1]]))})))}
  getdata<-function(object,slot.name,group.1,group.2,do.stat,duplex){
    object.names <- names(methods::slot(object, slot.name))
    prob.list <- list()
    pSum <- list()
    pSum.original <- list()
    pair.name <- list()
    idx <- list()
    pSum.original.all <- c()
    object.names.comparison <- c()
    comparison <- seq(1,length(unique(object@meta$datasets)))
    for (i in 1:length(comparison)) {
      object.list <- methods::slot(object, slot.name)[[comparison[i]]]
      prob <- object.list$prob
      prob[object.list$pval > 0.05] <- 0
      prob.list[[i]] <- prob
      rawtype<-rownames(prob)
      g1<-intersect(rawtype,group.1)
      g2<-intersect(rawtype,group.2)
      g.twice<-intersect(g1,g2)
      data1<-prob
      data1[setdiff(rawtype,g1),,]<-0
      data1[,setdiff(rawtype,g2),]<-0
      if(duplex){
        data2<-prob[,,]
        data2[setdiff(rawtype,g2),,]<-0
        data2[,setdiff(rawtype,g1),]<-0
        if(length(g.twice)>0){
          data2[g.twice,g.twice,]<-0
        }
      }else{data2<-0}
      prob<-data1+data2
      if (sum(prob) == 0) {
        warning("No inferred communications for the input!")
      }
      pSum.original[[i]] <- apply(prob, 3, sum)
      pSum[[i]] <- -1/log(pSum.original[[i]])
      pair.name[[i]] <- names(pSum.original[[i]])
      pSum[[i]][is.na(pSum[[i]])] <- 0
      idx[[i]] <- which(is.infinite(pSum[[i]]) | pSum[[i]] <  0)
      pSum.original.all <- c(pSum.original.all, pSum.original[[i]][idx[[i]]])
      object.names.comparison <- c(object.names.comparison, 
                                   object.names[comparison[i]])
    }
    values.assign <- seq(max(unlist(pSum)) * 1.1, max(unlist(pSum)) * 
                           1.5, length.out = length(unlist(idx)))
    position <- sort(pSum.original.all, index.return = TRUE)$ix
    for (i in 1:length(comparison)) {
      if (i == 1) {
        pSum[[i]][idx[[i]]] <- values.assign[match(1:length(idx[[i]]), position)]
      }
      else {
        pSum[[i]][idx[[i]]] <- values.assign[match(length(unlist(idx[1:i - 1])) + 1:length(unlist(idx[1:i])), position)]
      }
    }
    pair.name.all <- as.character(unique(unlist(pair.name)))
    df <- list()
    for (i in 1:length(comparison)) {
      df[[i]] <- data.frame(name = pair.name.all, contribution = 0, 
                            contribution.scaled = 0, group = object.names[comparison[i]], 
                            row.names = pair.name.all)
      df[[i]][pair.name[[i]], 3] <- pSum[[i]]
      df[[i]][pair.name[[i]], 2] <- pSum.original[[i]]
    }
    contribution.relative <- list()
    for (i in 1:(length(comparison) - 1)) {
      contribution.relative[[i]] <- as.numeric(format(df[[length(comparison) - 
                                                            i + 1]]$contribution/df[[1]]$contribution, digits = 1))
      contribution.relative[[i]][is.na(contribution.relative[[i]])] <- 0
    }
    names(contribution.relative) <- paste0("contribution.relative.", 
                                           1:length(contribution.relative))
    for (i in 1:length(comparison)) {
      for (j in 1:length(contribution.relative)) {
        df[[i]][[names(contribution.relative)[j]]] <- contribution.relative[[j]]
      }
    }
    df[[1]]$contribution.data2 <- df[[length(comparison)]]$contribution
    if (length(comparison) == 2) {
      idx <- with(df[[1]], order(-contribution.relative.1, 
                                 contribution, -contribution.data2))
    }
    else if (length(comparison) == 3) {
      idx <- with(df[[1]], order(-contribution.relative.1, 
                                 -contribution.relative.2, contribution, -contribution.data2))
    }
    else if (length(comparison) == 4) {
      idx <- with(df[[1]], order(-contribution.relative.1, 
                                 -contribution.relative.2, -contribution.relative.3, 
                                 contribution, -contribution.data2))
    }
    else {
      idx <- with(df[[1]], order(-contribution.relative.1, 
                                 -contribution.relative.2, -contribution.relative.3, 
                                 -contribution.relative.4, contribution, -contribution.data2))
    }
    for (i in 1:length(comparison)) {
      df[[i]] <- df[[i]][idx, ]
      df[[i]]$name <- factor(df[[i]]$name, levels = as.character(df[[i]]$name))
    }
    df[[1]]$contribution.data2 <- NULL
    df <- do.call(rbind, df)
    if (do.stat & length(comparison) == 2) {
      for (i in 1:length(pair.name.all)) {
        if (nrow(prob.list[[j]]) != nrow(prob.list[[1]])) {
          stop("Statistical test is not applicable to datasets with different cellular compositions! Please set `do.stat = FALSE`")
        }
        prob.values <- matrix(0, nrow = nrow(prob.list[[1]]) * 
                                nrow(prob.list[[1]]), ncol = length(comparison))
        for (j in 1:length(comparison)) {
          if (pair.name.all[i] %in% pair.name[[j]]) {
            prob.values[, j] <- as.vector(prob.list[[j]][, 
                                                         , pair.name.all[i]])
          }
          else {
            prob.values[, j] <- NA
          }
        }
        prob.values <- prob.values[rowSums(prob.values, 
                                           na.rm = TRUE) != 0, , drop = FALSE]
        if (nrow(prob.values) > 3 & sum(is.na(prob.values)) == 
            0) {
          pvalues <- wilcox.test(prob.values[, 1], prob.values[, 
                                                               2], paired = TRUE)$p.value
        }
        else {
          pvalues <- 0
        }
        pvalues[is.na(pvalues)] <- 0
        df$pvalues[df$name == pair.name.all[i]] <- pvalues
      }
    }
    for (i in 1:length(pair.name.all)) {
      df.t <- df[df$name == pair.name.all[i], "contribution"]
      if (sum(df.t) == 0) {
        df <- df[-which(df$name == pair.name.all[i]), 
        ]
      }
    }
    df$contribution <- abs(df$contribution)
    df$contribution.scaled <- abs(df$contribution.scaled)
    return(df)
  }
  temp<-lapply(cellchat_list,function(x){return(getdata(x,slot.name,group.1,group.2,T,duplex))})
  groups<-unique(unlist(lapply(temp,function(x){x$group})))
  temp<-lapply(temp,function(x){return(x[which(x$pvalues<p.val),])})
  paths<-lapply(temp,function(x){
    strong<-list()
    if(dim(x)[1]==0){
      length(strong)<-length(groups)
      names(strong)<-groups
      return(strong)
    }
    contri<-reshape2::dcast(x[,c('name','group','contribution')],name~group,value.var = 'contribution')
    strong[[colnames(contri)[2]]]<-as.character(contri$name[which((contri[,2]-contri[,3])/(contri[,2]+contri[,3])>tol)])
    strong[[colnames(contri)[3]]]<-as.character(contri$name[which((contri[,3]-contri[,2])/(contri[,2]+contri[,3])>tol)])
    return(strong)
  })
  df.1<-100*table(unlist(lapply(paths,function(x){x[[groups[1]]]})))/length(cellchat_list)
  df.2<-100*table(unlist(lapply(paths,function(x){x[[groups[2]]]})))/length(cellchat_list)
  order.1<-order(as.numeric(df.1),decreasing = T)[1:min(topn,length(df.1))]
  order.2<-rev(order(as.numeric(df.2),decreasing = T)[1:min(topn,length(df.2))])
  df.1<-df.1[order.1]
  df.2<-df.2[order.2]
  shared.path<-intersect(names(df.1),names(df.2))
  if(length(shared.path)>0){
    names(df.1)[match(shared.path,names(df.1))]<-paste0(names(df.1)[match(shared.path,names(df.1))],'.1')
    names(df.2)[match(shared.path,names(df.2))]<-paste0(names(df.2)[match(shared.path,names(df.2))],'.2')
  }
  df<-data.frame(signal=factor(names(c(df.1,df.2)),levels=names(c(df.1,df.2))),
                 info_flow=as.numeric(c(df.1,df.2)),group=c(rep(groups[1],length(df.1)),rep(groups[2],length(df.2))))
  ggplot(df) + geom_bar(aes(x=signal,y=info_flow,fill=group),stat='identity',position = position_dodge2()) +
    theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5,size = 10)) +ylab('Pathway percent')+xlab('Singal')+ coord_cartesian(ylim=c(0,100))
}
PlotRankNetPairwise(intra_NT_Chat_list,slot.name = 'netP',topn = 10,duplex = T)
