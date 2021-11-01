# Landscape sctructure affects individual flower visitation network
# Oliveira R, Souza CS, Schneiberg I and Varassin IG
# 

#### Required packages ####
install.packages("bipartite")
install.packages("corrplot")
install.packages("ggplot2")
install.packages("matrixStats")
install.packages("broom")
install.packages("AICcmodavg")
install.packages("cowplot")

library (tidyverse)
library (ggplot2) # for graphics
library (matrixStats) # to aplly functions in matrix
library (corrplot) # correlation graphics
library (ggcorrplot) # better visualization for ggplot
library (bipartite) # network analysis
library (broom) # to visualize better linear regression
library (AICcmodavg) # package to models selections analysis
library (cowplot) # complement of graphics visualization

#### Network analysis ####

# 1) General network matrix
complete <- read.table("complete.txt")
plotweb(complete)
nulls <- nullmodel(complete)

# 2) Modularity Measure an CZ Values
# 2.1) General Modularity measure

mod.complete<-metaComputeModules(complete, N=100)
mod.complete@likelihood # Modulatiry value
png(filename = "mod.png", width = 7, height = 8, units = "in", res = 300)
plotModuleWeb(mod.complete,rank=TRUE, labsize=1)
dev.off()
printoutModuleInformation(mod.complete)

png(filename = "net.png", width = 8, height = 5, units = "in", res = 300)
plotweb(complete, col.interaction = "aquamarine", bor.col.interaction = "aquamarine4", col.high="coral2", col.low="lightskyblue4", text.high.col="red3", text.rot=60, labsize=2, high.lab.dis = 0.15,low.lablength = 0, high.y = 1.4)
dev.off()

czlocals<-czvalues(mod.complete, level="higher", weighted = T)
czbirds<-czvalues(mod.complete, level="lower", weighted = T)
# The following functions are just to replace NA Values (considering one species modules as 0)
czlocals<-rapply(czlocals, function(x) ifelse(is.na(x),0,x), how="replace" )
czbirds<-rapply(czbirds, function(x) ifelse(is.na(x),0,x), how="replace" )

# 2.2) Z-score and comparison to null models
modfun<-function(x) {metaComputeModules(x, N=100)} 
mod.nulls<-lapply(nulls, modfun) # Do not run again, it takes 8 hours
# use load ("C;/.../mod.nulls.RData")
mod.nulls1<-unlist(lapply(mod.nulls, function(x) x@likelihood)) # Converting to a list
zscore_mod <- (mod.complete@likelihood - mean(mod.nulls1))/sd(mod.nulls1) # Z-score

# 2.3) Species and Local roles + Null models
# Call landscape metrics matrix to plot the graphics in next section

# Locals
nullczlocals<-lapply(mod.nulls, function (x) czvalues(x, level="higher", weighted = T))

quantile(sapply(nullczlocals, function(x) x$c),0.95, na.rm = T)
quantile(sapply(nullczlocals, function(x) x$z),0.95, na.rm = T)

czlocals1<-cbind(czlocals[["c"]],czlocals[["z"]])
colnames(czlocals1)<-c("c","z")
czlocals1<-as.data.frame(czlocals1)

czlocals_graphic <- ggplot(data = czlocals1, aes(c, z)) +
  geom_point(color = "red", alpha = .2, aes(size = landvar$urban_500)) +
  geom_text(aes(label=row.names(czlocals1)),hjust=-.8, vjust=-.8,
            colour = "red", alpha = .5, size = 5) +
  geom_hline (yintercept = 1.4227, alpha = .2, size = 0.9) +
  geom_vline (xintercept = 0.7168, alpha = .2, size = 0.9) +
  labs (x = "Among-module connectivity (c)", y = "Within-module degree (z)") +
  theme(text=element_text(size=16),
        legend.position = "none",
        axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black"),
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"))
czlocals_graphic

ggsave("czvalue.png", plot = czlocals_graphic, dpi = 300)

# Birds
nullczbirds<-lapply(mod.nulls, function (x) czvalues(x, level="lower", weighted = T))

quantile(sapply(nullczbirds, function(x) x$c),0.95, na.rm = T)
quantile(sapply(nullczbirds, function(x) x$z),0.95, na.rm = T)

czbirds1<-cbind(czbirds[["c"]],czbirds[["z"]])
colnames(czbirds1)<-c("c","z")
czbirds1<-as.data.frame(czbirds1)

czbirds_graphic <- ggplot(data = czbirds1, aes(c, z)) +
  geom_point(color = "blue", size = 5, alpha = .2) +
  geom_hline (yintercept = 1.7873, alpha = .2, size = .9) +
  geom_vline (xintercept = 0.7292, alpha = .2, size = .9) +
  labs (x = "Among-module connectivity (c)", y = "Within-module degree (z)") +
  theme (text=element_text (size = 16),
         axis.line = element_line(colour = "black"),
         axis.text = element_text(colour = "black"),
         panel.grid = element_blank(),
         plot.background = element_rect(fill = "transparent"),
         panel.background = element_rect (fill = "transparent"))

# 3) Node metrics

# 3,1) Caculation of normalised degree, closeness centrality, betweenness centrality and interactions selectivity index


nodelocals <- cbind(specieslevel(complete, index = c("normalised degree",
                                               "closeness", "betweenness"),
                                 level = "higher")) %>% select(c(1,2,4))


# correlation between "node level" metrics for local
nodelocalscor<-cor(nodelocals)
row.names(nodelocalscor) <- c("ND", "BC", "CC")
colnames(nodelocalscor) <- c("ND", "BC", "CC")

teste <- ggcorrplot(nodelocalscor, type = "lower", show.diag = FALSE,
                    lab = TRUE, show.legend = FALSE)
teste

png(filename = "corresp.png", width = 5, height = 5, units = "in", res = 300)
corrplot(nodelocalscor, method="shade", shade.col=NA, tl.col="black", tl.srt=45, addCoef.col="black", addcolorlabel="no", bg = "transparent")
dev.off()

# correlation between species level metrics for bird species
nodebirdscor<-cor(nodebirds)
row.names(nodebirdscor) <- c("ND", "BC", "wBC", "CC", "wCC", "d'", "c", "z")
colnames(nodebirdscor) <- c("ND", "BC", "wBC", "CC", "wCC", "d'", "c", "z")
corrplot(nodebirdscor, method="shade", shade.col=NA, tl.col="black", tl.srt=45, addCoef.col="black", addcolorlabel="no", bg = "transparent")

# 3.2) null models for centrality measures and z-score

# applying functions in the nulls matrices
nodelocals_nulls<-lapply(nulls, function(x) specieslevel(x, level = "higher",
                                              index = c("normalised degree",
                                                        "closeness",
                                                        "betweenness")))

nodebirds_nulls<-lapply(nulls, function(x) specieslevel(x, level = "lower",
                                                         index = c("normalised degree",
                                                                   "closeness",
                                                                   "betweenness")))

ndlocals_null <- sapply(nodelocals_nulls, function(x) x[,1])
bclocals_null <- sapply(nodelocals_nulls, function(x) x[,2])
cclocals_null <- sapply(nodelocals_nulls, function(x) x[,4])


zscore_nodelocals <- cbind((nodelocals[,1] - rowMeans(ndlocals_null))/ rowSds(ndlocals_null),
                           (nodelocals[,2] - rowMeans(bclocals_null))/ rowSds(bclocals_null),
                           (nodelocals[,3] - rowMeans(cclocals_null))/ rowSds(cclocals_null))
colnames(zscore_nodelocals) <- c("ND", "BC", "CC")
row.names(zscore_nodelocals) <- colnames(complete)



ndbirds_null <- sapply(nodebirds_nulls, function(x) x[,1])
bcbirds_null <- sapply(nodebirds_nulls, function(x) x[,2])
ccbirds_null <- sapply(nodebirds_nulls, function(x) x[,4])

zscore_nodebirds <- cbind((nodebirds[,1] - rowMeans(ndbirds_null))/ rowSds(ndbirds_null),
                           (nodebirds[,2] - rowMeans(bcbirds_null))/ rowSds(bcbirds_null),
                           (nodebirds[,4] - rowMeans(ccbirds_null))/ rowSds(ccbirds_null))
colnames(zscore_nodebirds) <- c("ND", "BC", "CC")
row.names(zscore_nodebirds) <- rownames(complete)

# Summary of node metrics
# nodelocals - all of calculated metrics for each local
# nodebirds - all of calculated metrics for each bird species
# nodelocals_nulls - node metrics applied to all nulls matrices for locals
# nodebirds_nulls -  node metrics applied to all nulls matrices for birds
# zscore_nodelocals - difference values of observed and expected by null models for node level metrics for each local
# zscore_nodebirds - difference values of observed and expected by null models for node level metrics for each bird species





#### Landscape Metrics ####

# Correlation between predictible variables

landvar<-read.table("variables.txt")
landvar_cor<-cor(landvar)

rownames(landvar_cor)<-c("FC 0.5", "FC 1.0", "FC 2.0", "GA 0.5", "GA 1.0", "GA 2.0", "NF 0.5", "NF 1.0", "NF 2.0", "ISO 0.5", "ISO 1.0", "ISO 2.0", "UC 0.5", "UC 1.0", "UC 2.0")
colnames(landvar_cor)<-c("FC 0.5", "FC 1.0", "FC 2.0", "GA 0.5", "GA 1.0", "GA 2.0", "NF 0.5", "NF 1.0", "NF 2.0", "ISO 0.5", "ISO 1.0", "ISO 2.0", "UC 0.5", "UC 1.0", "UC 2.0")


plot_cor <- ggcorrplot(landvar_cor, type = "lower", show.diag = FALSE,
                      lab = TRUE)

ggsave(plot = plot_cor, filename = "corpred.png",
       dpi = 700)


col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
png(filename = "corpred.png", width = 8, height = 8, units = 'in', res = 300)
corrplot(landvar_cor, method="shade", shade.col=NA, tl.col="black", tl.srt=45,
         col=col(200), addCoef.col="black", addcolorlabel="no", order="AOE")
dev.off()

# making linear regressions for 


nodelocals1<-subset.data.frame(nodelocals, select = -c(weighted.closeness, weighted.betweenness))

regressions_net_land <- NULL

for (i in 1:ncol(nodelocals)){
  for (j in 1:ncol(landvar)){
    modelsnode<-lm(nodelocals[,i]~landvar[,j])
    
    regressions_net_land<-rbind(regressions_net_land,
                       cbind(network_metrics = i, 
                                 landscape_metrics = j,
                                 glance(modelsnode),
                                 intercept=model_$coefficients[1],
                                 slope=model_$coefficients[2]
                       )
    )
  }
}



### GLMS
### 
### In this we are going to use the following data frames
attach(nodelocals)
attach(landvar)

### normalised.degree
Cand.models <- list( )
Cand.models[[1]]<-glm(normalised.degree~floresta_500,family = gaussian())
Cand.models[[2]]<-glm(normalised.degree~vegetacao_500,family = gaussian())
Cand.models[[3]]<-glm(normalised.degree~fragmentos_500,family = gaussian())
Cand.models[[4]]<-glm(normalised.degree~iso_2000,family = gaussian())
Cand.models[[5]]<-glm(normalised.degree~urban_500,family = gaussian())
Cand.models[[6]]<-glm(normalised.degree~floresta_500 + I(floresta_500^2),family = gaussian())
Cand.models[[7]]<-glm(normalised.degree~urban_500 + I(urban_500^2),family = gaussian())
Cand.models[[8]]<-glm(normalised.degree~floresta_500 + vegetacao_500,family = gaussian())
Cand.models[[9]]<-glm(normalised.degree~vegetacao_500 + iso_2000,family = gaussian())
Cand.models[[10]]<-glm(normalised.degree~1)

Modnames<-paste("mod",1:length(Cand.models),sep=" ")
glm_nd<-aictab(cand.set = Cand.models, modnames = Modnames,sort = TRUE)


## betweenness
Cand.models <- list( )
Cand.models[[1]]<-glm(betweenness~floresta_500,family = gaussian())
Cand.models[[2]]<-glm(betweenness~vegetacao_500,family = gaussian())
Cand.models[[3]]<-glm(betweenness~fragmentos_500,family = gaussian())
Cand.models[[4]]<-glm(betweenness~iso_2000,family = gaussian())
Cand.models[[5]]<-glm(betweenness~urban_500,family = gaussian())
Cand.models[[6]]<-glm(betweenness~urban_500 + I(urban_500^2),family = gaussian())
Cand.models[[7]]<-glm(betweenness~floresta_500 + I(floresta_500^2),family = gaussian())
Cand.models[[8]]<-glm(betweenness~floresta_500 + vegetacao_500,family = gaussian())
Cand.models[[9]]<-glm(betweenness~vegetacao_500 + iso_2000,family = gaussian())
Cand.models[[10]]<-glm(betweenness~1)

Modnames<-paste("mod",1:length(Cand.models),sep=" ")
glm_bc<-aictab(cand.set = Cand.models, modnames = Modnames,sort = TRUE)

## closeness
Cand.models <- list( )
Cand.models[[1]]<-glm(closeness~floresta_1000,family = gaussian())
Cand.models[[2]]<-glm(closeness~vegetacao_500,family = gaussian())
Cand.models[[3]]<-glm(closeness~fragmentos_500,family = gaussian())
Cand.models[[4]]<-glm(closeness~iso_1000,family = gaussian())
Cand.models[[5]]<-glm(closeness~urban_1000,family = gaussian())
Cand.models[[6]]<-glm(closeness~urban_1000 + I(urban_1000^2),family = gaussian())
Cand.models[[7]]<-glm(closeness~floresta_1000 + I(floresta_1000^2),family = gaussian())
Cand.models[[8]]<-glm(closeness~floresta_1000 + vegetacao_500,family = gaussian())
Cand.models[[9]]<-glm(closeness~vegetacao_500 + iso_1000,family = gaussian())
Cand.models[[10]]<-glm(closeness~vegetacao_500 + urban_1000,family = gaussian())
Cand.models[[11]]<-glm(closeness~1)

Modnames<-paste("mod",1:length(Cand.models),sep=" ")
glm_cc<-aictab(cand.set = Cand.models, modnames = Modnames,sort = TRUE)

##


glm_result1 <- list(glm_nd, glm_bc, glm_cc)
names(glm_result1) <- c("ND", "BC", "CC")

install.packages("openxlsx")
library(openxlsx)
glm_result1 <- createWorkbook()
addWorksheet(glm_result1, "ND")
addWorksheet(glm_result1, "BC")
addWorksheet(glm_result1, "CC")
writeData(glm_result1, 1, glm_result$ND)
writeData(glm_result1, 2, glm_result$BC)
writeData(glm_result1, 3, glm_result$CC)

saveWorkbook(glm_result1, file = "glm_result.xlsx", overwrite = TRUE)


# most parsimonious models
#
final_graphics <- cbind(nodelocals, landvar)


# Normalised degree
library(png)

net <- readPNG("net.png")
net2 <- as.raster(net)

netplot <- ggplot() + theme(plot.background = element_blank(),
                            panel.background = element_blank()) +
  annotation_raster(net2, xmin = 0, xmax = 1, ymin = 0, ymax = 1)

netplot

col_site = c("deeppink", "deeppink","deeppink", "orange1",
             "deeppink", "green2", "green2", "deeppink",
             "green2", "orange1", "blue")


# green areas
ndveg <- ggplot(data = final_graphics, aes(x = vegetacao_500, y = normalised.degree)) +
  geom_point (size = 1.5, alpha=.6, colour = col_site) +
  stat_smooth(method = lm, alpha = .1, colour = "darkgray") +
  geom_text(aes(label=ifelse(normalised.degree>0.001,
                             as.character(rownames(final_graphics)), "")),
            alpha = .5, size = 3, vjust=-0.8, hjust = -0.3) +
  labs (x = "Green areas (500 m)", y = "Normalised degree") +
  theme_minimal() + theme(axis.title = element_text(size = 8))

ndveg


ndfl <- ggplot(data = final_graphics, aes(x = floresta_500, y = normalised.degree)) +
  geom_point (size = 1.5, alpha=.6, colour = col_site) +
  stat_smooth(method = lm, formula = y~x+I(x^2), alpha = .1, colour = "darkgray") +
  geom_text(aes(label=ifelse(normalised.degree>0.001,
                             as.character(rownames(final_graphics)), "")),
            alpha = .5, size = 3, vjust=1.5, hjust = 0.5) +
  labs (x = "Forest Cover (500 m)", y = "") +
  theme_minimal() + theme(axis.title = element_text(size = 8))

ndfl


ndfr <- ggplot(data = final_graphics, aes(x = fragmentos_500, y = normalised.degree)) +
  geom_point (size = 1.5, alpha=.6, colour = col_site) +
  stat_smooth(method = lm, alpha = .1, colour = "darkgray") +
  geom_text(aes(label=ifelse(normalised.degree>0.001,
                             as.character(rownames(final_graphics)), "")),
            alpha = .5, size = 3, vjust=1.5, hjust = 0.5) +
  labs (x = "Number of Fragments (500 m)", y="") +
  theme_minimal() + theme(axis.title = element_text(size = 8))

ndfr

ndgrap <- ggarrange(ndveg, ndfl, ndfr, ncol=3, labels = c("b","c","d"))

####

#Betwenness

bcfr <- ggplot(data = final_graphics, aes(x = fragmentos_500, y = betweenness)) +
  geom_point (size = 3, alpha=.6, colour = col_site) +
  stat_smooth(method = lm, alpha = .1, colour = "darkgray") +
  geom_text(aes(label=ifelse(betweenness>0.001,
                             as.character(rownames(final_graphics)), "")),
            alpha = .5, size = 3, vjust=1.2, hjust = -.5) +
  labs (x = "Number of Fragments (500 m)", y = "Betweenness") +
  theme_minimal()

bcfr


bcvg <- ggplot(data = final_graphics, aes(x = vegetacao_500, y = betweenness)) +
  geom_point (size = 3, alpha=.6, colour = col_site) +
  stat_smooth(method = lm, alpha = .1, colour = "darkgray") +
  geom_text(aes(label=ifelse(betweenness>0.001,
                             as.character(rownames(final_graphics)), "")),
            alpha = .5, size = 3, vjust=-1.5) +
  labs (x = "Green Areas (500 m)", y = "") +
  theme_minimal()

bcvg

bcgraphs <- ggarrange(bcfr, bcvg, ncol=2, labels = c("e","f"))



finalplot <- ggarrange(netplot, ndgrap, bcgraphs, nrow = 3, labels = c("a","",""),
                   heights = c(2.25, 1, 1.5))

ggsave(plot = finalplot, filename = "net-rel.png", dpi = 700, width = 8,
       height = 9.5)

