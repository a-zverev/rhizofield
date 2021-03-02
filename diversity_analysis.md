---
title: "Diversity in soils"
output:
  html_document: 
    keep_md: yes
    number_section: yes
    toc: yes
    toc_float: yes
---


```r
library(phyloseq)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(ggpubr)
library(picante)
library(reshape2)
library(tidyr)
library(psych)
library(DESeq2)
library(fpc)

setwd('/home/alexey/Analysis/RW/3rd publication/rhizofield/')
source('functions2.R')
```


```r
#Black soil
roots.BS <- readRDS('roots.BS.RData')
bact.BS <- readRDS('bact.BS.RData')  %>% subset_taxa(Family != "Mitochondria" & Class != "Chloroplast")

bact.BS <- subset_samples(bact.BS, Spot %in% c('1', '2', '3'))

n.roots.BS <- rarefy_collapse(roots.BS)
n.bact.BS <- rarefy_collapse(bact.BS)
```


```r
#Podzol soil
roots.PS <- readRDS('roots.PS.RData')
bact.PS <- readRDS('bact.PS.RData') %>% subset_taxa(Family != "Mitochondria" & Class != "Chloroplast")

bact.PS <- subset_samples(bact.PS, Spot %in% c('1', '2', '3'))

n.roots.PS <- rarefy_collapse(roots.PS)
n.bact.PS <- rarefy_collapse(bact.PS)
```

## Представленность таксонов

### Состав корневых сообществ

Как задать корректный порог? Или все же отдельные графики?


```r
a <- phyloseq(otu_table(n.roots.BS),
              tax_table(n.roots.BS),
              sam_data(n.roots.BS))
b <- phyloseq(otu_table(n.roots.PS),
              tax_table(n.roots.PS),
              sam_data(n.roots.PS))
c <- merge_phyloseq(a, b)

bargraph(c, "Genus", 0.1) + facet_grid(~Soil, scale = 'free_x')
```

![](diversity_analysis_files/figure-html/unnamed-chunk-4-1.png)<!-- -->

```r
ggarrange(bargraph(n.roots.BS, "Genus", 0.03),
          bargraph(n.roots.PS, "Genus", 0.03),
          ncol = 2)
```

![](diversity_analysis_files/figure-html/unnamed-chunk-4-2.png)<!-- -->

### Состав бактериальных сообществ



```r
a <- phyloseq(otu_table(n.bact.BS),
              tax_table(n.bact.BS),
              sam_data(n.bact.BS))
b <- phyloseq(otu_table(n.bact.PS),
              tax_table(n.bact.PS),
              sam_data(n.bact.PS))
c <- merge_phyloseq(a, b)

bargraph(c, "Phylum", 0.03) + facet_grid(~Soil, scale = 'free_x')
```

![](diversity_analysis_files/figure-html/unnamed-chunk-5-1.png)<!-- -->


## Альфа-разнообразие

**Гипотеза: богатство корневых сообществ определяет богатство сообществ бактерий**

Четыре метрики. Индексы **richness** у нас представляют количество видов **Observed** и сумма длин веток дерева **PD**. За индексы **eveness** будут выступать индекс Симпсона **Simpson** и взвешенная средняя дистанция между последовательностями **mpd**.
Каким образом эти экологические метрики сочетаются с более "сырыми" дистанциями **p-value**, и можем ли мы увидеть "отпечатки" растительных сообществ в сообществах микроорганизмов?
Давайте посчитаем корреляцию между этими индексами, рассчитанными для разных почв в парах "корни - бактерии"


### Графики корреляции метрик альфа-разнообразия


```r
plot_internal_correlation(alpha_div(n.roots.BS, pairwise_distances_from_ps(n.roots.BS)) %>% select(-Source, -Spot))
```

![](diversity_analysis_files/figure-html/unnamed-chunk-6-1.png)<!-- -->

```r
plot_internal_correlation(alpha_div(n.roots.PS, pairwise_distances_from_ps(n.roots.PS)) %>% select(-Source, -Spot))
```

![](diversity_analysis_files/figure-html/unnamed-chunk-6-2.png)<!-- -->

```r
plot_internal_correlation(alpha_div(n.bact.BS, pairwise_distances_from_ps(n.bact.BS)) %>% select(-Source, -Spot))
```

![](diversity_analysis_files/figure-html/unnamed-chunk-6-3.png)<!-- -->

```r
plot_internal_correlation(alpha_div(n.bact.PS, pairwise_distances_from_ps(n.bact.PS)) %>% select(-Source, -Spot))
```

![](diversity_analysis_files/figure-html/unnamed-chunk-6-4.png)<!-- -->

### 


```r
# Metrics
n.roots.BS.alpha <- bootstrapped_alpha_div(n.roots.BS)
n.roots.PS.alpha <- bootstrapped_alpha_div(n.roots.PS)
n.bact.BS.alpha <- bootstrapped_alpha_div(n.bact.BS)
n.bact.PS.alpha <- bootstrapped_alpha_div(n.bact.PS)

drawer <- function(metric){
  p1 <- plot_correlation(n.roots.BS.alpha, n.bact.BS.alpha, metric)
  p3 <- plot_correlation(n.roots.PS.alpha, n.bact.PS.alpha, metric)
  ggarrange(p1, p3, ncol = 2)
}

drawer('Observed')
```

![](diversity_analysis_files/figure-html/unnamed-chunk-7-1.png)<!-- -->

```r
drawer('PD')
```

![](diversity_analysis_files/figure-html/unnamed-chunk-7-2.png)<!-- -->

```r
drawer('Simpson')
```

![](diversity_analysis_files/figure-html/unnamed-chunk-7-3.png)<!-- -->

```r
drawer('mpd')
```

![](diversity_analysis_files/figure-html/unnamed-chunk-7-4.png)<!-- -->

```r
drawer('p.dist')
```

![](diversity_analysis_files/figure-html/unnamed-chunk-7-5.png)<!-- -->


Вероятно, мы оказываемся на границе чувствительности, так как при очень схожем паттерне достоверность то есть, то нет. 

---


## Бета-разнообразие

**Гипотеза: богатство корневых сообществ определяет богатство сообществ бактерий и грибов**

### Бета-разнообразие


```r
beta_plot_drawer <- function(metric){
  p1 <- beta_plot(n.roots.PS, metric)
  p2 <- beta_plot(n.bact.PS, metric)
  p3 <- beta_plot(n.roots.BS, metric)
  p4 <- beta_plot(n.bact.BS, metric)
  
  ggarrange(ggarrange(p1, p2, common.legend = T, legend = 'bottom'),
            ggarrange(p3, p4, common.legend = T, legend = 'bottom'),
            nrow = 2)
}


beta_plot_drawer("unifrac")
```

![](diversity_analysis_files/figure-html/unnamed-chunk-8-1.png)<!-- -->

### Статистика для кластеров бета-разнообразия

Гипотеза: структура кластеров похожа. Для этого используем статистики из пакета кластерного анализа. Использованы следующие метрики:

Для кластера:

* Dunn - индекс Дунна. Рассчитывается один на весь плот, мера силы кластеризации. Отношение минимальной меры сепарации к максимальному диаметру кластера.
* Median_dist - медианная дистанция. Медианная дистанция между точками в выбранном классе.
* Min_separation - минимальная величина сепарации. Минимальное расстояние от между точкой, принадлежащей этому кластеру, и точкой, принадлежащей любому другому кластеру.

Для пары кластеров:

* Clusters_Separation - минимальная дистанция между точкой в одном кластере и точкой в другом кластере.
* Distance_Between - средняя дистанция между кластерами.



```r
# https://www.rdocumentation.org/packages/fpc/versions/2.2-9/topics/cluster.stats

# because every run of unifrac is unique, we should use bootstrap as multiple calculation
# (or, perhaps, root the tree?)
bootstrapped_wunifrac_cluster_stat <- function(ps, repeats){
  # dm <- as.matrix(UniFrac(ps, weighted=TRUE))
  dm.matrices <- replicate(repeats, as.matrix(UniFrac(ps, weighted=TRUE)), simplify = FALSE)
  true.dm <- Reduce("+", dm.matrices) / length(dm.matrices)
  info <- sample_data(ps)$Source
  clusters <- levels(info)
  levels(info) <- c(1, 2, 3)
  info <- as.numeric(as.vector(info))
  res <- cluster.stats(d = true.dm, clustering = info, silhouette = FALSE)
  common <- data.frame(Cluster = clusters,
             Median_Dist = res$median.distance,
             Min_separation = res$separation,
             # Sepr_m = res$separation.matrix,
             # AvBtw_m = res$ave.between.matrix,
             Dunn = res$dunn) %>% 
    mutate_if(is.numeric, round, digits = 3)
  matrix_to_df <- function(X, names){
      colnames(X) <- clusters
      rownames(X) <- clusters
      ind <- which(upper.tri(X, diag = F), arr.ind = TRUE)
      nn <- dimnames(X)
      data.frame(Cluster_1 = nn[[1]][ind[, 1]],
           Cluster_2 = nn[[2]][ind[, 2]],
           val = X[ind])
  }
  extra2 <- matrix_to_df(res$ave.between.matrix, clusters) %>% rename(Distance_Between = val)
  extra <- matrix_to_df(res$separation.matrix, clusters) %>% rename(Clusters_Separation = val)
  extra <- merge(extra, extra2, by = c('Cluster_1', 'Cluster_2'))
  list(common, extra)
}
```

#### Чернозем


```r
p1 <- beta_plot(n.roots.BS, 'wunifrac')
d1 <- bootstrapped_wunifrac_cluster_stat(n.roots.BS, 100)

p2 <- beta_plot(n.bact.BS, 'wunifrac')
d2 <- bootstrapped_wunifrac_cluster_stat(n.bact.BS, 100)


ggarrange(p1, p2, common.legend = T)
```

![](diversity_analysis_files/figure-html/unnamed-chunk-10-1.png)<!-- -->

```r
merge(d1[1] %>% as.data.frame(), d2[1] %>% as.data.frame, 
      by = 'Cluster', 
      suffixes = c(":Roots",":Bact")) %>%  select(order(colnames(.)))
```

```
##    Cluster Dunn:Bact Dunn:Roots Median_Dist:Bact Median_Dist:Roots
## 1 Graminae     0.422      0.605            0.049             0.501
## 2  Melilot     0.422      0.605            0.034             0.206
## 3    Wheat     0.422      0.605            0.055             0.096
##   Min_separation:Bact Min_separation:Roots
## 1               0.026                0.358
## 2               0.026                0.358
## 3               0.041                0.538
```

```r
merge(d1[2] %>% as.data.frame(), d2[2] %>% as.data.frame, 
      by = c('Cluster_1', 'Cluster_2'), 
      suffixes = c(":Roots",":Bact")) %>%  select(order(colnames(.)))
```

```
##   Cluster_1 Cluster_2 Clusters_Separation:Bact Clusters_Separation:Roots
## 1  Graminae   Melilot               0.02619387                 0.3581302
## 2  Graminae     Wheat               0.04128057                 0.5592192
## 3   Melilot     Wheat               0.04252493                 0.5375461
##   Distance_Between:Bact Distance_Between:Roots
## 1            0.04567953              0.4635901
## 2            0.05189913              0.6289813
## 3            0.05700721              0.6301929
```

#### Подзол


```r
p1 <- beta_plot(n.roots.PS, 'wunifrac')
d1 <- bootstrapped_wunifrac_cluster_stat(n.roots.PS, 100)

p2 <- beta_plot(n.bact.PS, 'wunifrac')
d2 <- bootstrapped_wunifrac_cluster_stat(n.bact.PS, 100)


ggarrange(p1, p2, common.legend = T)
```

![](diversity_analysis_files/figure-html/unnamed-chunk-11-1.png)<!-- -->

```r
merge(d1[1] %>% as.data.frame(), d2[1] %>% as.data.frame, 
      by = 'Cluster', 
      suffixes = c(":Roots",":Bact")) %>%  select(order(colnames(.)))
```

```
##    Cluster Dunn:Bact Dunn:Roots Median_Dist:Bact Median_Dist:Roots
## 1 Bedstraw     0.805      0.709            0.027             0.308
## 2 Graminae     0.805      0.709            0.025             0.476
## 3      Rye     0.805      0.709            0.027             0.197
##   Min_separation:Bact Min_separation:Roots
## 1               0.022                0.390
## 2               0.022                0.390
## 3               0.035                0.867
```

```r
merge(d1[2] %>% as.data.frame(), d2[2] %>% as.data.frame, 
      by = c('Cluster_1', 'Cluster_2'), 
      suffixes = c(":Roots",":Bact")) %>%  select(order(colnames(.)))
```

```
##   Cluster_1 Cluster_2 Clusters_Separation:Bact Clusters_Separation:Roots
## 1  Bedstraw  Graminae               0.02241464                 0.3896401
## 2  Bedstraw       Rye               0.03694616                 0.9305281
## 3  Graminae       Rye               0.03462622                 0.8671479
##   Distance_Between:Bact Distance_Between:Roots
## 1            0.02616957              0.5903996
## 2            0.04757406              0.9735252
## 3            0.04204957              0.9443660
```

<!-- ### Корреляця дистанций бета-разнообразия -->

<!-- Каждая точка на графике - численное выражение дистанция между двумя сообществами (например, Graminae.1 - Rye.2). По оси Х эта дистанция взята для растительных сообществ, по оси Y - для сообществ микроорганизмов. Насколько распределение дистанций в корневых сообществах повторяет его в сообществах микроорганизмов? -->


```r
# beta_drawer <- function(metric){
#   p1 <- beta_corr(n.roots.BS, n.bact.BS, metric)
#   p3 <- beta_corr(n.roots.PS, n.bact.PS, metric)
#   ggarrange(p1, p3, ncol = 2)
# }

# beta_drawer('unifrac')
# beta_drawer('wunifrac')
# beta_drawer('bray')
```




```r
# # Strange beta p-distance
# 
# all.beta.pd <- read.csv("beta_p-distance.csv", comment.char = "#", sep = "\t")
# 
# beta.pd.roots.BS <- all.beta.pd %>% filter(norm == 1, Source == "Roots", Soil == "BS")
# beta.pd.roots.PS <- all.beta.pd %>% filter(norm == 1, Source == "Roots", Soil == "PS")
# beta.pd.bact.BS <- all.beta.pd %>% filter(norm == 1, Source == "Bact", Soil == "BS")
# beta.pd.bact.PS <- all.beta.pd %>% filter(norm == 1, Source == "Bact", Soil == "PS")
# 
# convert_longer <- function(table_df){
#   table_df <- table_df  %>% 
#     select(-norm, -Source, -Soil)
#   colnames(table_df) <- c("X", as.character(table_df$X)) 
#   table_df %>%
#     pivot_longer(!X, names_to = "Y", values_to = "value")
# }
# 
# beta_p_distance_correlation <- function(root, bact){
#   a <- convert_longer(root)
#   b <- convert_longer(bact)
#   
#   data <- inner_join(a, b, by = c('X', 'Y'))
#   x = data$value.x
#   y = data$value.y
#   data$Color <- paste(gsub('.{2}$', '', data$X), '-', gsub('.{2}$', '', data$Y))
#   #correlation
#   fit <- corr.test(x=x, y=y)
#   print(data)
#   print(fit$p)
#   header <- paste0('R: ', round(fit$r, 3), ' (p-value - ', round(fit$p, 3), ')')
#   #plotting
#   ggplot(data = data, aes(x = x, y = y)) + 
#     geom_point(aes(color = Color)) +
#     geom_smooth(method = 'lm', se = F) + 
#     theme_light() +
#     theme(plot.title = element_text(color = ifelse(fit$p < 0.05, "darkgreen", "brown2"))) +
#     labs(title=paste(header), 
#          x = paste("beta_p-dist", '-', deparse(substitute(root))), 
#          y = paste("beta_p-dist", '-', deparse(substitute(bact))))
#   
# }
# 
# 
# 
# print(beta_p_distance_correlation(beta.pd.roots.BS, beta.pd.bact.BS))
# print(beta_p_distance_correlation(beta.pd.roots.PS, beta.pd.bact.PS))
```



## DeSEQ анализ

### Достоверно изменяющие численность таксоны


```r
plot_heatmap(ps_from_deseq(n.bact.PS, "PS", "Family"), taxa.label = "Family", na.value = "white",
             high = "#000099", low = "grey")
```

![](diversity_analysis_files/figure-html/unnamed-chunk-14-1.png)<!-- -->


```r
plot_heatmap(ps_from_deseq(n.bact.BS, "BS", "Family"), taxa.label = "Family", na.value = "white",
             high = "#000099", low = "grey")
```

![](diversity_analysis_files/figure-html/unnamed-chunk-15-1.png)<!-- -->


### Ternary plots

See [here](https://github.com/a-zverev/rhizofield/blob/master/ggtern_diagrams.md)







