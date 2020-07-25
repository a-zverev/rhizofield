library(phyloseq)
library(ggplot2)
library(plyr)
library(ggpubr)
library(picante)
library(reshape2)
library(psych)

setwd('/home/alexey/Analysis/RW/Field_RW/All_data_R/ver2/')

rarefy_collapse <- function(ps){
  #set.seed(7493)
  ps <- rarefy_even_depth(ps, verbose = FALSE)
  ps <- tax_glom(ps, "Genus")
  return(ps)
}

plot_correlation <- function(ps.1, ps.2, metric){
  #read data
  names <- c(paste0(ps.1), paste0(ps.2))
  ps.1 <- get(paste0(ps.1))
  ps.2 <- get(paste0(ps.2))
  require(psych)
  alpha1 <- alpha_div(ps.1)
  alpha2 <- alpha_div(ps.2)
  data <- merge(x=subset(alpha1, select = -Source), 
                y=subset(alpha2, select = -Spot), by=0, all=T)
  x = data[,paste0(metric, '.x')]
  y = data[,paste0(metric, '.y')]
  #correlation
  fit <- corr.test(x=x, y=y, method="spearman")
  header <- paste0('R: ', round(fit$r, 3), ' (p-value - ', round(fit$p, 3), ')')
  #plotting
  data$Spot <- as.factor(data$Spot)
  ggplot(data = data, aes(x = x, y = y)) + 
    geom_point(aes(color = Source, shape = Spot)) +
    geom_smooth(method = 'lm') + theme_light() +
    theme(plot.title = element_text(color = ifelse(fit$p < 0.05, "darkgreen", "brown2"))) +
    labs(title=paste(header), x = paste(metric, '-', names[1]), y = paste(metric, '-', names[2]))
}

#Draw just OTUs/abundance
plot_abundance <- function(ps_object){
  draw <- psmelt(ps_object)
  #draw <- draw['Abundance' != 0]
  #draw <- draw[order(draw$Abundance),]
  ggplot(data = draw, aes(x=Genus, y=Abundance)) + 
    geom_point(aes(colour=Source, shape = as.factor(Spot))) +
    # geom_line(aes(colour=Source)) +
    theme_light() +
    #theme(axis.text.x = element_text(angle = 90)) +
    theme(axis.text.x = element_blank()) +
    ggtitle(paste("Abundance of reads in different sources -", paste(deparse(substitute(ps_object))),"reads"))
}

#Calculate weighted PD index from tree and OTU-table
weighted_pd <- function(tree, table){
  res <- c()
  #calculate lenghts for every tip for tree
  n <- length(tree$tip.label)
  ee <- setNames(tree$edge.length[sapply(1:n,function(x,y) which(y==x),y=tree$edge[,2])],tree$tip.label)
  lenghts <- as.data.frame(ee) #df with lenght of every tip
  
  
  #for every sample in a table, multiply lenght of tips by abundance in table
  weights <- as.data.frame(table)
  for (sample in rownames(weights)) {
    row.weights <- t(weights[sample,])
    lenghts.n.weights <- merge(lenghts, row.weights, by = 0, all = T)
    #print(lenghts.n.weights)
    lenghts.n.weights$w.lenght <- as.numeric(lenghts.n.weights$ee) * as.numeric(lenghts.n.weights[,sample])
    res <- c(res, sum(lenghts.n.weights$w.lenght))
  }
  res <- data.frame(row.names = rownames(weights), w.pd = res)
  return(res)
}

#Calculate several alpha-diversity indexes, return one dataframe
alpha_div <- function(ps){
  ps <- prune_taxa(taxa_sums(ps) > 0, ps)
  pd <- pd(ps@otu_table@.Data, ps@phy_tree, include.root = FALSE)
  obs_sim <- estimate_richness(ps, split = TRUE, measures = c("Observed", "Simpson", "Shannon"))
  #w.pd <- weighted_pd(ps@phy_tree, ps@otu_table)
  Source <- ps@sam_data$Source
  Spot <- ps@sam_data$Spot
  mpd <- mpd(ps@otu_table@.Data, cophenetic(ps@phy_tree), 
             abundance.weighted = TRUE)
  alpha <- cbind(obs_sim, 
                 pd, 
                 #w.pd, 
                 mpd, 
                 Source, Spot)
  return(alpha)
}

#Plot alpha-diversity by selected metric
plot_alpha <- function(alpha_table, metric){
  ggplot(data = alpha_table, aes(x=rownames(alpha_table), y=alpha_table[,metric])) + 
    geom_point(aes(colour=Source), size=3, alpha=0.7) + 
    theme_light() + theme(axis.text.x = element_text(angle = 90)) +
    #theme(axis.text.x = element_blank()) +
    labs(x = element_blank(), y = metric, size = 6) +
    ggtitle(paste(deparse(substitute(alpha_table)))) +
    theme(plot.title = element_text(size = 10, face = "italic"))
  #header as metric, remove axis names and maybe x labels
}

#Plot correlation-like graph from two group and metric
plot_correlation <- function(alpha1, alpha2, metric){
  data <- merge(x=subset(alpha1, select = -Source), 
                y=subset(alpha2, select = -Spot), by=0, all=T)
  data$Spot <- as.factor(data$Spot)
  ggplot(data = data, aes(x = data[,paste(metric, '.x', sep = '')], 
                          y = data[,paste(metric, '.y', sep = '')], 
                          color = Source, shape = Spot)) + 
    geom_point() + ggtitle(paste(metric, '-', deparse(substitute(alpha1)), 'vs', deparse(substitute(alpha2)))) +
    labs(x = paste(metric, deparse(substitute(alpha1))), y = paste(metric, deparse(substitute(alpha2)))) + theme_light()
}