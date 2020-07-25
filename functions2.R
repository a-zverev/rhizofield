library(phyloseq)
library(ggplot2)
library(plyr)
library(ggpubr)
library(picante)
library(reshape2)
library(psych)

rarefy_collapse <- function(ps){
  set.seed(7493)
  ps <- rarefy_even_depth(ps, 0.7*min(sample_sums(ps)), verbose = FALSE)
  ps <- tax_glom(ps, "Genus")
  return(ps)
}

pattern <- function(soil, community){
  text = paste(soil, '\n', 'Roots -vs- ', community, '\ncommunities')
  ggplot() + 
    annotate("text", x = 4, y = 25, size=8, label = text) + 
    theme_light() +
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank()) +
    labs(title=element_blank(), x = element_blank(), y = element_blank())
}

bootstrapped_alpha_div <- function(ps){
  rarefy_collapse <- function(ps){
    #set.seed(7493)
    ps <- rarefy_even_depth(ps, 0.7*min(sample_sums(ps)), verbose = FALSE)
    ps <- tax_glom(ps, "Genus")
    return(ps)
  }
  data <- rdply(10, alpha_div(rarefy_collapse(ps)))
  mean <- stats::aggregate(. ~ Source + Spot, data, "mean")
  sd <- stats::aggregate(. ~ Source + Spot, data, "sd")
  res <- merge(x = mean, y = sd, by = c("Source", "Spot"), suffixes = c('.mean', '.sd'))
  return(res)
}

plot_correlation <- function(ps.1, ps.2, metric){
  require(psych)
  require(RColorBrewer)
  require(stats)
  #calculate rarefactions as bootstraps
  alpha1 <- bootstrapped_alpha_div(ps.1)
  alpha2 <- bootstrapped_alpha_div(ps.2)
  data <- merge(alpha1, alpha2,
                by = c("Source", "Spot"), 
                suffixes = c('.x', '.y'), all=T)
  #data <- merge(x=subset(alpha1, select = -Source), 
  #              y=subset(alpha2, select = -Spot), by=0, all=T)
  x <- data[,paste0(metric, '.mean.x')]
  y <- data[,paste0(metric, '.mean.y')]
  err.x <- data[,paste0(metric, '.sd.x')]
  err.y <- data[,paste0(metric, '.sd.x')]
  #correlation
  fit <- corr.test(x=x, y=y, method="spearman")
  header <- paste0('R: ', round(fit$r, 3), ' (p-value - ', round(fit$p, 3), ')')
  #plotting
  data$Spot <- as.factor(data$Spot)
  ggplot(data = data, aes(x = x, y = y)) + 
    geom_point(aes(color = Source)) +
    geom_errorbar(aes(ymin=(y-err.y), 
                      ymax=(y+err.y)),
                  size = 0.1) +
    geom_errorbarh(aes(xmin=(x-err.x), 
                       xmax=(x+err.x)),
                   size = 0.1) +
    geom_smooth(method = 'lm') + theme_light() +
    theme(plot.title = element_text(color = ifelse(fit$p < 0.05, "darkgreen", "brown2"))) +
    scale_color_brewer(palette = ifelse (deparse(substitute(ps.1)) == 'ng.roots.BS', "Set1", "Accent")) +
    labs(title=paste(header), x = paste(metric, '-', deparse(substitute(ps.1))), 
         y = paste(metric, '-', deparse(substitute(ps.2))))
}

plot_p_distance <- function(data, group, soil){
  require(psych)
  require(RColorBrewer)
  require(stats)
  
  x <- data[,"roots.mean"]
  y <- data[,paste0(group, '.mean')]
  err.x <- data[,'roots.sd']
  err.y <- data[,paste0(group, '.sd')]
  #correlation
  fit <- corr.test(x=x, y=y, method="spearman")
  header <- paste0('R: ', round(fit$r, 3), ' (p-value - ', round(fit$p, 3), ')')
  #plotting
  #data$Spot <- as.factor(data$Spot)
  ggplot(data = data, aes(x = x, y = y)) + 
    geom_point(aes(color = Source)) +
    geom_errorbar(aes(ymin=(y-err.y), 
                      ymax=(y+err.y)),
                  size = 0.1) +
    geom_errorbarh(aes(xmin=(x-err.x), 
                       xmax=(x+err.x)),
                   size = 0.1) +
    geom_smooth(method = 'lm') + theme_light() +
    theme(plot.title = element_text(color = ifelse(fit$p < 0.05, "darkgreen", "brown2"))) +
    scale_color_brewer(palette = ifelse (soil == 'BS', "Set1", "Accent")) +
    labs(title=paste(header), x = paste0('p-dist - ng.roots.', soil), 
         y = paste0('p-dist - ng.', group, '.', soil))
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
  alpha <- cbind(obs_sim, pd, mpd, Source, Spot)
  return(alpha)
}

#Calculate bootstrapped alpha-div
bootstrapped_alpha_div <- function(ps){
  rarefy_collapse <- function(ps){
    #set.seed(7493)
    ps <- rarefy_even_depth(ps, 0.7*min(sample_sums(ps)), verbose = FALSE)
    ps <- tax_glom(ps, "Genus")
    return(ps)
  }
  data <- rdply(2, alpha_div(rarefy_collapse(ps)))
  mean <- stats::aggregate(. ~ Source + Spot, data, "mean")
  sd <- stats::aggregate(. ~ Source + Spot, data, "sd")
  res <- merge(x = mean, y = sd, by = c("Source", "Spot"), suffixes = c('.mean', '.sd'))
  return(res)
}

beta_plot <- function(ps, metric){
  ps.prop <- transform_sample_counts(ps, function(x) x/sum(x))
  ord.nmds.bray <- ordinate(ps.prop, method='PCoA', distance=metric)
  plot_ordination(ps.prop, ord.nmds.bray, color = "Source", title=paste(metric, '-', deparse(substitute(ps)))) +
    geom_point(size=3, alpha=0.7) + labs() +
    theme_light()
}

beta_plot_drawer <- function(metric){
  p1 <- beta_plot(ng.roots.PS, metric)
  p2 <- beta_plot(ng.bact.PS, metric)
  p3 <- beta_plot(ng.roots.BS, metric)
  p4 <- beta_plot(ng.bact.BS, metric)
  
  ggarrange(ggarrange(p1, p2, common.legend = T, legend = 'bottom'),
            ggarrange(p3, p4, common.legend = T, legend = 'bottom'),
            nrow = 2)
}

beta_corr <- function(ps.1, ps.2, method){
  
  dist <- function(ps, method){
    beta <- phyloseq::distance(ps, method=method)
    mat <- as.matrix(beta)
    diag(mat) <- NA
    mat[lower.tri(mat)] <- NA
    melt <- melt(mat, na.rm = T)
    melt$Var3 <- gsub('.{2}$', '', melt$Var1)
    melt$Var4 <- gsub('.{2}$', '', melt$Var2)
    
    #intra <- melt[melt$Var3 == melt$Var4,][,1:3]
    extra <- melt[melt$Var3 != melt$Var4,][,1:3]
    return(extra)
    #if (what == 'intra') {return(intra)}
    #else if (what == 'extra') {return(extra)}
  }
  
  require(psych)
  a <- dist(ps.1, method)
  b <- dist(ps.2, method)
  
  data <- merge(a, b, by = c('Var1', 'Var2'))
  x = data$value.x
  y = data$value.y
  data$Color <- paste(gsub('.{2}$', '', data$Var1), '-', gsub('.{2}$', '', data$Var2))
  #correlation
  fit <- corr.test(x=x, y=y, method="spearman")
  header <- paste0('R: ', round(fit$r, 3), ' (p-value - ', round(fit$p, 3), ')')
  #plotting
  ggplot(data = data, aes(x = x, y = y)) + 
    #geom_point(aes(color = Color)) +
    geom_point(aes()) +
    geom_smooth(method = 'lm', se = F) + 
    theme_light() +
    theme(plot.title = element_text(color = ifelse(fit$p < 0.05, "darkgreen", "brown2"))) +
    labs(title=paste(header), 
         x = paste(method, '-', deparse(substitute(ps.1))), 
         y = paste(method, '-', deparse(substitute(ps.2))))
  
}

sig_table <- function(ps_object){
  ds <- phyloseq_to_deseq2(ps_object, ~Source)
  ds = estimateSizeFactors(ds, type="poscounts")
  ds = estimateDispersions(ds, fitType = "local")
  ds = DESeq(ds)
  #mcols(ds, use.names=TRUE)
  res = results(ds, cooksCutoff = FALSE)
  alpha = 0.05
  sigtab = res[which(res$padj < alpha), ]
  if (nrow(sigtab) == 0) {
    return(NA)
  }
  sigtab = cbind(as(sigtab, "data.frame"), 
                 as(tax_table(ps_object)[rownames(sigtab), ], "matrix")
  )
  return(sigtab)
}

bargraph <- function(ps, rank, threshold){
  physeq2 = filter_taxa(ps, function(x) mean(x) > 0.1, TRUE)
  physeq3 = transform_sample_counts(physeq2, function(x) x / sum(x) )
  ph_ps <- tax_glom(physeq3, taxrank = rank)
  data <- psmelt(ph_ps) # create dataframe from phyloseq object
  data$Plot <- as.character(data[,rank]) #convert to character
  data$Plot[data$Abundance < threshold] <- paste0("<", threshold, " abund.")
  medians <- ddply(data, ~Plot, function(x) c(median=median(x$Abundance)))
  remainder <- medians[medians$median <= threshold,]$Plot
  p <- ggplot(data=data, aes(x=Sample, y=Abundance, fill=Plot))
  p + geom_bar(aes(), stat="identity", position="stack") + theme_light() +
    scale_fill_manual(values = c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue",
                                 "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey")) +
    theme(legend.position="bottom") + guides() +
    theme(axis.text.x = element_text(angle = 90))
}

ps_from_deseq <- function(ps, soil){
  require(phyloseq)
  require(DESeq2)
  if (soil == "BS") {
    ps.1 <- subset_samples(ps, Source %in% c('Graminae', 'Wheat'))
    ps.2 <- subset_samples(ps, Source %in% c('Melilot', 'Wheat'))
    ps.3 <- subset_samples(ps, Source %in% c('Melilot', 'Graminae'))
  }
  if (soil == "PS") {
    ps.1 <- subset_samples(ps, Source %in% c('Graminae', 'Rye'))
    ps.2 <- subset_samples(ps, Source %in% c('Bedstraw', 'Rye'))
    ps.3 <- subset_samples(ps, Source %in% c('Bedstraw', 'Graminae'))
  }
  a <- sig_table(ps.1)
  b <- sig_table(ps.2)
  c <- sig_table(ps.3)
  names <- unique(c(rownames(a), rownames(b), rownames(c)))
  ps <- prune_taxa(names, ps)
  return(ps)
}