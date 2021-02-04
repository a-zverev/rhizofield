library(plyr)
library(psych)
library(phangorn)
library(phyloseq)
library(DECIPHER)
library(Biostrings)
library(tidyverse)
library(reshape2)
library(seqinr)
library(picante)
library(corrplot)
library(plotly)
library(ggpubr)

rarefy_collapse <- function(ps, collapse = FALSE){
  set.seed(7493)
  ps <- rarefy_even_depth(ps, min(sample_sums(ps)), verbose = FALSE)
  if (collapse == T) {
    ps <- tax_glom(ps, "Genus")
  }
  return(ps)
}

############################# ALPHA - DIV #############################

# pattern <- function(soil, community){
#   text = paste(soil, '\n', 'Roots -vs- ', community, '\ncommunities')
#   ggplot() + 
#     annotate("text", x = 4, y = 25, size=8, label = text) + 
#     theme_light() +
#     theme(panel.grid.major=element_blank(),
#           panel.grid.minor=element_blank(),
#           axis.text.x = element_blank(),
#           axis.text.y = element_blank()) +
#     labs(title=element_blank(), x = element_blank(), y = element_blank())
# }

# Calculate pairwise distance for OTUs in phyloseq-object
pairwise_distances_from_ps <- function(ps){
  # Multiple alignment
  alignment <- AlignSeqs(refseq(ps), gapOpening=c(-16, -12), gapExtension=c(-12, -8), processors = 25, verbose = FALSE)
  
  # Transform alignment to seqinr, calculate distance matrix
  seqinr.alignment <- as.alignment(nb = length(alignment), 
                                   nam = names(alignment), 
                                   seq = as.character(alignment), 
                                   com = width(alignment))
  am <- as.matrix(dist.alignment(seqinr.alignment, matrix = "identity"))
  
  #form lower triangle pairs data
  am[lower.tri(am, diag = TRUE)] <- NA
  dm <- melt(am, varnames = c("V1", "V2")) %>% 
    filter(!is.na(value)) %>% 
    mutate(pid = -(value^2 - 1)*100,
           pseudo.p.dist = 1 + (value^2 - 1)) %>% 
    select(-value)
  
  return(dm)
}

# Calculate alpha- metrics (including weighted p-distance) from ps and distance matrix
alpha_div <- function(ps, distance_matrix){
  ps <- prune_taxa(taxa_sums(ps) > 0, ps)
  df <- t(otu_table(ps)) %>% as.data.frame()
  
  multiple_p_dist <- function(column){
    df <- data.frame(column)
    df  %>% filter(.[,1] != 0)
    colnames(df) <- "Abd"
    df$V <- rownames(df)
    distance_matrix %>%
      inner_join(df, by = c("V1" = "V")) %>%
      inner_join(df, by = c("V2" = "V")) %>%
      mutate(weight = (Abd.y / sum(df$Abd) * (Abd.x / sum(df$Abd))),
             weighted.pseudo.p.dist = pseudo.p.dist * weight) %>%
      select(weighted.pseudo.p.dist) %>% sum()
  }
  
  p.dist <- apply(df, 2, multiple_p_dist)
  pd.tree <- pd(ps@otu_table@.Data, ps@phy_tree, include.root = FALSE)
  obs_sim <- estimate_richness(ps, split = TRUE, measures = c("Observed", "Simpson", "Shannon"))
  mpd <- mpd(ps@otu_table@.Data, cophenetic(ps@phy_tree), 
             abundance.weighted = TRUE)
  Source <- ps@sam_data$Source
  Spot <- ps@sam_data$Spot
  alpha <- cbind(p.dist, obs_sim, pd.tree, mpd, Source, Spot)
  return(alpha)
}

#Calculate bootstrapped alpha-div
bootstrapped_alpha_div <- function(ps){
  dm <- pairwise_distances_from_ps(ps)
  data <- rdply(3, alpha_div(rarefy_even_depth(ps, 0.7*min(sample_sums(ps)), verbose = FALSE), dm))
  mean <- stats::aggregate(. ~ Source + Spot, data, "mean")
  sd <- stats::aggregate(. ~ Source + Spot, data, "sd")
  res <- merge(x = mean, y = sd, by = c("Source", "Spot"), suffixes = c('.mean', '.sd'))
  return(res)
}

# Plot correlation for two alpha-metrics tables
plot_correlation <- function(alpha1, alpha2, metric){
  require(psych)
  require(RColorBrewer)
  require(stats)
  #calculate rarefactions as bootstraps
  # alpha1 <- bootstrapped_alpha_div(ps.1)
  # alpha2 <- bootstrapped_alpha_div(ps.2)
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
    scale_color_brewer(palette = ifelse (deparse(substitute(alpha1)) == 'ng.roots.BS.alpha', "Set1", "Accent")) +
    labs(title=paste(header), x = paste(metric, '-', deparse(substitute(alpha1))), 
         y = paste(metric, '-', deparse(substitute(alpha2))))
}

# Plot correlations in alpha-metrics
plot_internal_correlation <- function(df){
  col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
  cor.mtest <- function(mat, ...) {
    mat <- as.matrix(mat)
    n <- ncol(mat)
    p.mat<- matrix(NA, n, n)
    diag(p.mat) <- 0
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        tmp <- cor.test(mat[, i], mat[, j], ...)
        p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
      }
    }
    colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
    p.mat
  }
  p.mat <- cor.mtest(df)
  corrplot(cor(df), method="color", col=col(200),  
           type="upper", order="hclust", 
           addCoef.col = "black", # Add coefficient of correlation
           tl.col="black", tl.srt=45, #Text label color and rotation
           # Combine with significance
           p.mat = p.mat, sig.level = 0.01, insig = "blank", 
           # hide correlation coefficient on the principal diagonal
           diag=FALSE 
  )
  
}

############################# BETA - DIV #############################

# Plot beta-diversity
beta_plot <- function(ps, metric){
  ps.prop <- transform_sample_counts(ps, function(x) x/sum(x))
  ord.nmds.bray <- ordinate(ps.prop, method='PCoA', distance=metric)
  plot_ordination(ps.prop, ord.nmds.bray, color = "Source", title=paste(metric, '-', deparse(substitute(ps)))) +
    geom_point(size=3, alpha=0.7) + labs() +
    theme_light()
}

# Plot beta-distance correlation
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

############################# DESEQ 2 #############################

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
  ps2 <- tax_glom(ps, taxrank = rank)
  ps3 = transform_sample_counts(ps2, function(x) x / sum(x) )
  data <- psmelt(ps3) # create dataframe from phyloseq object
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

ps_from_deseq <- function(ps, soil, taxa){
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
  ps <- ps <- tax_glom(ps, taxa)
  return(ps)
}