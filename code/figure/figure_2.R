library(readr)
library(readxl)
library(tidyverse)





# set the factor level by the mean value of occur number in the graph
Factor <-Result0 %>% group_by(pathway) %>% summarise(median=mean(Occur_num))
Factor <- Factor[order(Factor$median),]
Result0$pathway <-factor(Result0$pathway, levels=Factor$pathway)

ggplot(data=Result0, aes(x=pathway, y=Occur_num)) + geom_boxplot(fill='grey93') +
  xlab('') + ylab('') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 75, hjust = 1)) +
  theme(legend.position = c(0.85, 0.2)) +
  theme(axis.text=element_text(size=10, family="Arial"),
        axis.title=element_text(size=12,family="Arial"),
        legend.text = element_text(size=10, family="Arial")) +
  ggtitle('') +
  theme(panel.background = element_rect(fill = "white", color="black", size = 1)) #+







##batch plot for file
for (input in files){
  output <- basename(file.path(input,fsep=".txt")) 
  graph  <- read.table(input,header=T)
  pdf(paste0("output for ", input,".pdf"))
  plot(graph$x,graph$y)
  dev.off()
} 


##batch plot for dataframe
enrGSANit <- read_csv("enrGSANit.csv")
data <- enrGSANit
data <- as.data.frame(data)
for ( name in colnames(data)[-1]) {
  pdf(paste(name, '.pdf', sep = ""))
  lbls <- c("A", "B", "C")
  pie(data[, name], labels = lbls, main=name)
  dev.off()
}



## extract the shortest path between two nodes
require(igraph)
require(tidyverse)
require(stringr)
c1 <- 'A -B'
edges = graph.formula(c1)



edges= graph.formula(A - B,A - C,B - C,B - C,C - D,B - D,D - E,f)
components(edges)

plot(edges, edge.arrow.size=.4)
all_shortest_paths(edges,from="A", to="D")
all_simple_paths(edges,from="A", to="D")

all_simple_paths(edges)




nodes <- read.csv("Dataset1-Media-Example-NODES.csv", sep = ";", header=T, as.is=T)
links <- read.csv("Dataset1-Media-Example-EDGES.csv", sep = ";", header=T, as.is=T)

# adjust the weight of each link
links$weight <- 10

#links <- aggregate(links[,3], links[,-3], sum) # remove the duplicated connect

net <- graph_from_data_frame(d=links, vertices=nodes, directed=T) 
net <- graph_from_data_frame(d=links,directed=FALSE) 

class(net)
ss <- E(net)       # The edges of the "net" object
V(net)       # The vertices of the "net" object
E(net)$type  # Edge attribute "type"
V(net)$media # Vertex attribute "media"

#net <- simplify(net, remove.multiple = F, remove.loops = T) 
plot(net, edge.arrow.size=.4)

connectionAll <- all_simple_paths(net,from="s13", to="s06")
connectionShort <- all_shortest_paths(net,from="s13", to="s06")

#Breadth-first search
ss <- bfs(net, root= "s13", "out",
    order=TRUE, rank=TRUE, father=TRUE, pred=TRUE,
    succ=TRUE, dist=TRUE)


node_parameter <- select(nodes, id)
node_parameter$dist <- ss$dist

distance_target <- node_parameter$dist[which(node_parameter$id=="s06")]

## choose all the medium node between start node and the target node
medium_node <- filter(node_parameter, dist <= distance_target)


## simplify all the connection based on the medium nodes
ss1 <- which(nodes$id %in% medium_node$id ==TRUE)
nodes0 <- nodes[ss1, ]

ss2 <- which(links$from %in% medium_node$id ==TRUE)
ss3 <- which(links$to %in% medium_node$id ==TRUE)
ss4 <- intersect(ss2,ss3)

links0 <- links[ss4,]

net_new0 <- graph_from_data_frame(d=links00, vertices=nodes00, directed=T) 
plot(net_new0, edge.arrow.size=.4)

all_simple_paths(net_new,from="s13", to="s06")
shortest_path <- all_shortest_paths(net_new,from="s13", to="s06")

shortest_path$res[[1]]

###new way
ss <- bfs(net, root= "s13", "out",
          order=TRUE, rank=TRUE, father=TRUE, pred=TRUE,
          succ=TRUE, dist=TRUE)

node_parameter <- select(nodes, id)
node_parameter$dist <- ss$dist
nodes <- nodes[,1]




# this two functions could find all the children node based on input node
getNextLayer <- function(NET=net, NODE = nodes, root ){
  #root <- "s13"
  ss <- bfs(NET, root, "out",
            order=TRUE, rank=TRUE, father=TRUE, pred=TRUE,
            succ=TRUE, dist=TRUE)
  
  node_parameter <- select(NODE, id)
  node_parameter$dist <- ss$dist
  node1 <- filter(node_parameter, dist==1)
  node2 <- node1$id
  return(node2)
}

getMultiNodes <- function (multiRoot){
  #mutiRoot = tt[[5]]
  ss <- unlist(multiRoot)
  nn <- list()
  mm <- length(ss)
  for (i in 1:mm){
    nn[[i]] <- getNextLayer(NET = net, NODE = nodes, ss[i] )
    
  }
  
  tt <- unique(unlist(nn))
  return(tt)
}

# first step
tt <- list()
tt[[1]] <- "s13"
tt[[2]] <- getMultiNodes(multiRoot = tt[[1]])
tt[[3]] <- getMultiNodes(multiRoot = tt[[2]])
tt[[4]] <- getMultiNodes(multiRoot = tt[[3]])
tt[[5]] <- getMultiNodes(multiRoot = tt[[4]])
tt[[6]] <- getMultiNodes(multiRoot = tt[[5]])
tt[[7]] <- getMultiNodes(multiRoot = tt[[6]])

##########################example
# Second step
cp <- list()
cp[[1]] <- "s13"

#second
tt2 <- list()
tt2[[1]] <- c("s13","s17")
tt2[[2]] <- c("s13","s12")

cp[[2]] <- tt2


#thid
ss1 <- getNextLayer(root = "s17")
ss2 <- getNextLayer(root = "s12")

tt3 <- list()
tt3[[1]] <- c(tt2[[1]],ss1[1]) 
tt3[[2]] <- c(tt2[[2]],ss2[1]) 
tt3[[3]] <- c(tt2[[2]],ss2[2]) 
tt3[[4]] <- c(tt2[[2]],ss2[3]) 

cp[[3]] <- tt3
###########################################example


#function get the connection

getPath <- function (t){

 #t<- list()
 #t[[1]] <- "s13"
tt2 <- t
nn <- length(tt2)
length0 <- vector()
n0 <- list()
for (i in 1:nn){
  n0[[i]] <- getNextLayer(root = tt2[[i]][length(tt2[[i]])])
  length0[i] <- length(n0[[i]])
}

#next layer
nn0 <- unlist(n0)

tt20 <- list()
for (i in 1:nn){
  tt20[[i]] <- paste0(tt2[[i]], collapse = "@")
}

mm0 <- list()

for (i in 1:nn){
  mm0[[i]] <- replicate(length0[i], tt20[[i]])
}

mm1 <- unlist(mm0)

# combine the current layer and next layer
total <- sum(length0)

mm1 <- paste(mm1,nn0, sep = "@")  


mm2 <- str_split(mm1,"@")

return(mm2)

}

tt1 <- list()
tt1[[1]] <- "s13"
tt2 <- getPath(tt1)
tt3 <- getPath(tt2)
tt4 <- getPath(tt3)
tt5 <- getPath(tt4)
tt6 <- getPath(tt5)
tt7 <- getPath(tt6)
tt8 <- getPath(tt7)
tt9 <- getPath(tt8)
tt10 <- getPath(tt9)
tt11 <- getPath(tt10)
tt12 <- getPath(tt11)
tt13 <- getPath(tt12)


#clust the graph to obtain the connected subgraph
g <- simplify(
graph.compose(
  graph.ring(10), 
  graph.star(5, mode = "undirected")
)
) + edge("7", "8")


clusters(g)
dg <- decompose.graph(g)
plot(g, edge.arrow.size=.4)
plot(dg[[1]])
plot(dg[[2]])
plot(dg[[3]])


g <- sample_gnp(10, 1/5)
plot(g, edge.arrow.size=.4)
clu  <- components(g)


g <- make_ring(10)
ei <- get.edge.ids(g, c(1,2, 4,5))
E(g)[ei]



##mfuzz
source("https://bioconductor.org/biocLite.R")
biocLite("Mfuzz")
library(Mfuzz)
library(stringr)
library(tidyverse)
Mfuzzgui()

FPKM <- read.delim2('FPKM.txt', header = TRUE, sep = "\t", stringsAsFactors = FALSE)

FPKM$X24h[478]

#two missing 0 0 1 1
FPKM1 <- filter(FPKM, X16h ==""& X24h=="" & X42h !=""& X66h !="")
FPKM1$X42h <- as.numeric(FPKM1$X42h)
FPRM10 <- filter(FPKM1,X42h <= 10 )
#two missing 1 1 0 0
FPKM2 <- filter(FPKM, X16h !=""& X24h!="" & X42h ==""& X66h =="")

FPKM2$X16h <- as.numeric(FPKM2$X16h)
FPRM20 <- filter(FPKM2,X16h <= 10 )

## cut-off 0.25 remove 3782 genes
## cut-off 0.50 remove 2193 genes

library(readxl)
old <- read_excel("不同发酵阶段-聚类分析结果.xlsx")
cluster12 <- filter(old, cluster.12 >= 0.5)
cluster12 <- cluster12$gene.id


index0 <- which(FPKM$GENE.ID %in%cluster12 ==TRUE )
cluster12_d <- FPKM[index0,]


vs <- list(
  c(1, 3, 5, 6, 10),
  c(1, 2, 3, 7, 8, 10),
  c(1, 2, 3, 4, 8, 9, 10)
)

mu <- list(5,10,-3)
seq_along(vs)

ss <- c('a','b','c')
ss1 <- c(1,2,3)
test <-data.frame(name=ss, num=ss1)
for (i in seq_along(ss)){
  print(ss[i])
}

for (i in seq_along(test$name)){
  print(test$num[i])
}


#plot
library(readxl)
library(tidyverse)
library(ggplot2)
score_of_different_yeastGEM <- read_excel("score of different yeastGEM.xlsx", sheet = "model test memote")
s0 <- c("EC_missing","Stoichiometric Consistency")
score_of_different_yeastGEM <- score_of_different_yeastGEM[which(score_of_different_yeastGEM$Paramter %in% s0==FALSE),]
score_of_different_yeastGEM$Paramter <- as.factor(score_of_different_yeastGEM$Paramter)


p <- ggplot(score_of_different_yeastGEM, aes(x=Model_source, y=Score, shape=Paramter, colour=Paramter)) + geom_line() +
  geom_point(size=4) + 
  xlab('\nVersion of yeastGEM') +
  ylab('Score from memote test\n') +
  theme_bw() +
  theme(legend.position = c(0.70, 0.8))

p+theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text = element_text(size=10))




score_of_different_yeastGEM <- read_excel("score of different yeastGEM.xlsx", sheet = "new")
score_of_different_yeastGEM$Paramter <- as.factor(score_of_different_yeastGEM$Paramter)




p <- ggplot(score_of_different_yeastGEM, aes(x=Model, y=Score, shape=Paramter, colour=Paramter)) + geom_line(size=1.5) +
  geom_point(size=6) + 
  xlab('\nVersion of yeastGEM') +
  ylab('Score from memote test\n') +
  theme_bw() +
  theme(legend.position = c(0.85, 0.2)) +
  theme(axis.text=element_text(size=20,face="bold", family="Arial"),
        axis.title=element_text(size=24,face="bold", family="Arial"),
        legend.text = element_text(size=13,face="bold", family="Arial")) +
  ggtitle('') +
  theme(panel.background = element_rect(fill = "white", color="black", size = 1)) +
  ggsave(out <- paste('result/','memote test score','.eps', sep = ""), width=8, height=6, dpi=300)





#------------------------------------------------------
#heatmap
#-----------------------------------------------------

# heastmap-superheat
library(tidyverse)
library(readxl)
library(superheat)
growth <- read_excel("maxGrowth_data.xlsx")
rownames0 <- growth$carbon
growth <- growth[,-1]
rownames(growth) <- rownames0
superheat(growth,
          # scale the matrix columns
          scale =FALSE,
          # add row dendrogram
          left.label.size = 0.2,
          bottom.label.size = 0.4,
          left.label.text.size = 3,
          bottom.label.text.size = 3,
          grid.hline.col = "white",
          grid.vline.col = "white",
          bottom.label.text.angle = 90,
          legend = TRUE)



#fcc
fcc <- read_excel("fcc.xlsx")
rownames0 <- fcc$Carbon
fcc <- fcc[,-1]
rownames(fcc) <- rownames0
superheat(fcc,
          # scale the matrix columns
          scale =FALSE,
          # add row dendrogram
          row.dendrogram = FALSE,
          left.label.size = 0.4,
          bottom.label.size = 0.2,
          left.label.text.size = 4,
          bottom.label.text.size = 3,
          grid.hline.col = "white",
          grid.vline.col = "white",
          # change the color (#b35806 = brown and #542788 = purple)
          # change the color (#b35806 = brown and #542788 = purple)
          heat.pal = c("navy", "white", "firebrick3"),
          heat.pal.values = c(0, 0.5, 1),
          bottom.label.text.angle = 45,
          legend = TRUE)

  



#ggplot2-heatmap
library(ggplot2)

library(plyr)

require(reshape2)

require(scales)

nba <- read.csv( "http://datasets.flowingdata.com/ppg2008.csv")

nba$Name <- with(nba, reorder(Name, PTS))

nba.m <- melt(nba) #对数据进行融合

nba.m <- ddply(nba.m, .(variable), transform, rescale = rescale(value))

ggplot(nba.m, aes(variable, Name)) + 
  geom_tile(aes(fill = rescale), colour = "white") + 
  scale_fill_gradient(low = "white", high = "steelblue")



#pheatmap
library(pheatmap)
library(RColorBrewer)

library(grid)
draw_colnames_45 <- function (coln, ...) {
  m = length(coln)
  x = (1:m)/m - 1/2/m
  grid.text(coln, x = x, y = unit(0.96, "npc"), vjust = .5, 
            hjust = 1, rot = 90, gp = gpar(...)) ## 注意缺省值为 'hjust=0' 和'rot=270'
}
#然后将default draw_colnames覆盖
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                  ns=asNamespace("pheatmap"))



pheatmap(fcc,
         scale = 'none', cluster_rows=FALSE, cluster_cols=FALSE)
pheatmap(fcc,
         scale = 'column')

pheatmap(fcc,
         scale = 'none', #纵轴的zscore进行标准化一下
         color = colorRampPalette(colors = c("blue","white","red"))(100)
)

bk = unique(c(seq(0,0.1, length=100)))
pheatmap(fcc,breaks = bk
         ,color = colorRampPalette(c("white", "blue", "firebrick3"))(100))

pheatmap(fcc,
         method = c("pearson"),
         clustering_method = "complete",
         treeheight_row = 40,
         treeheight_col = 40,
         cluster_row = FALSE,
         cluster_col = FALSE,
         show_rownames = T,
         show_colnames = T,
         legend = T,
         fontsize = 15,
         color = colorRampPalette(brewer.pal(9,"Reds"))(400))


pheatmap(fcc,
         method = c("pearson"),
         clustering_method = "complete",
         treeheight_row = 40,
         treeheight_col = 40,
         cluster_row = FALSE,
         cluster_col = FALSE,
         show_rownames = T,
         show_colnames = T,
         legend = T,
         fontsize = 18,
         color = colorRampPalette(c("white", "SandyBrown", "firebrick3"))(100))







pheatmap(growth,
         scale = 'column',
         method = c("pearson"),
         clustering_method = "complete",
         treeheight_row = 40,
         treeheight_col = 40,
         cluster_row = TRUE,
         cluster_col = TRUE,
         show_rownames = T,
         show_colnames = T,
         legend = T,
         fontsize = 15,
         color = colorRampPalette(brewer.pal(9,"Reds"))(400))



pheatmap(growth,
         scale = 'none',
         method = c("pearson"),
         clustering_method = "complete",
         treeheight_row = 40,
         treeheight_col = 40,
         cluster_row = TRUE,
         cluster_col = FALSE,
         show_rownames = T,
         show_colnames = T,
         legend = T,
         fontsize = 15,
         color = colorRampPalette(c("white", "SandyBrown", "firebrick3"))(100))

#compare the growth of simulated and calculated
ecGEM_growth <- read_excel("calculate growth.xlsx", 
                               sheet = "ecGEM")
GEM_growth <- read_excel("calculate growth.xlsx", 
                           sheet = "GEM")
scale_grotwh <- read_excel("calculate growth.xlsx", 
                           sheet = "scale_measured")

growth_all <- data.frame(ecGEM_growth=unlist(ecGEM_growth), GEM_growth=unlist(GEM_growth), scale_grotwh=unlist(scale_grotwh), stringsAsFactors = FALSE)

growth_all$scale_grotwh <- as.numeric(growth_all$scale_grotwh)
#growth_all$scale_grotwh[is.na(growth_all$scale_grotw)] <- 0

growth_all <- filter(growth_all, scale_grotwh <0.5)


par(mgp=c(2,0.5,0),mar=c(2.5,2.5,1,1)+0.1,tcl=-0)
par(mfrow=c(1,2))
x1 <- growth_all$scale_grotwh
x2 <- growth_all$GEM_growth
df <- data.frame(x1,x2)

## Use densCols() output to get density at each point
x <- densCols(x1,x2, colramp=colorRampPalette(c("black", "white")))
df$dens <- col2rgb(x)[1,] + 1L

## Map densities to colors
cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", 
                            "#FCFF00", "#FF9400", "#FF3100"))(256)
df$col <- cols[df$dens]
df <- filter(df, x2 >0.5)
## Plot it, reordering rows so that densest points are plotted on top
plot(x2~x1, data=df[order(df$dens),], pch=20, col=col, cex=2, xlim=c(0, 0.5),ylim=c(0.5,4.5),xlab='',ylab='',
     cex.lab=1.2)


#title(xlab="Rescaled maximum growth rate(1/h)", line=1.6, cex.lab=1.2, family="Arial")
#title(ylab="Predicted maximum growth rate(1/h)", line=1.6, cex.lab=1.2, family="Arial")


x1 <- growth_all$scale_grotwh
x2 <- growth_all$ecGEM_growth
df <- data.frame(x1,x2)

## Use densCols() output to get density at each point
x <- densCols(x1,x2, colramp=colorRampPalette(c("black", "white")))
df$dens <- col2rgb(x)[1,] + 1L

## Map densities to colors
cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", 
                            "#FCFF00", "#FF9400", "#FF3100"))(256)
df$col <- cols[df$dens]
df <- filter(df, x2 >0.5)
## Plot it, reordering rows so that densest points are plotted on top
plot(x2~x1, data=df[order(df$dens),], pch=20, col=col, cex=2, xlim=c(0, 0.5),ylim=c(0.5, 4.5),xlab='',ylab='',
     cex.lab=1.2)
#title(xlab="Rescaled maximum growth rate(1/h)", line=1.6, cex.lab=1.2, family="Arial")
#title(ylab="Predicted maximum growth rate(1/h)", line=1.6, cex.lab=1.2, family="Arial")




#graph2
#par(mfrow=c(1,2))
#par(mgp=c(2,0.5,0),mar=c(2.5,2.5,1,1)+0.1,tcl=-0.2) #mar(左右上下)
x1 <- growth_all$scale_grotwh
x2 <- growth_all$GEM_growth
df <- data.frame(x1,x2)

## Use densCols() output to get density at each point
x <- densCols(x1,x2, colramp=colorRampPalette(c("black", "white")))
df$dens <- col2rgb(x)[1,] + 1L

## Map densities to colors
cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", 
                            "#FCFF00", "#FF9400", "#FF3100"))(256)
df$col <- cols[df$dens]
#df <- filter(df, x2 >0.5)
## Plot it, reordering rows so that densest points are plotted on top
plot(x2~x1, data=df[order(df$dens),], pch=20, col=col, cex=2, xlim=c(0, 0.5),ylim=c(0,0.5),xlab='',ylab='',
     cex.lab=1.2)
title(xlab="Rescaled maximum growth rate(1/h)", line=1.6, cex.lab=1.2, family="Arial")
title(ylab="Predicted maximum growth rate(1/h)", line=1.6, cex.lab=1.2, family="Arial")



x1 <- growth_all$scale_grotwh
x2 <- growth_all$ecGEM_growth
df <- data.frame(x1,x2)

## Use densCols() output to get density at each point
x <- densCols(x1,x2, colramp=colorRampPalette(c("black", "white")))
df$dens <- col2rgb(x)[1,] + 1L

## Map densities to colors
cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", 
                            "#FCFF00", "#FF9400", "#FF3100"))(256)
df$col <- cols[df$dens]

## Plot it, reordering rows so that densest points are plotted on top
plot(x2~x1, data=df[order(df$dens),], pch=20, col=col, cex=2, xlim=c(0, 0.5),ylim=c(0, 0.5),xlab='',ylab='',
     cex.lab=1.2)
title(xlab="Rescaled maximum growth rate(1/h)", line=1.6, cex.lab=1.2, family="Arial")
title(ylab="Predicted maximum growth rate(1/h)", line=1.6, cex.lab=1.2, family="Arial")


#plot the uncontinuous
#############
library(plotrix)
x <- c(1:5, 6.9, 7)
y <- 2^x
from <- 33
to <- 110
plot(x, y, type="b", xlab="index", ylab="value")
gap.plot(x, y, gap=c(from,to), type="b", xlab="index", ylab="value")
axis.break(2, from, breakcol="snow", style="gap")
axis.break(2, from*(1+0.02), breakcol="black", style="slash")
axis.break(4, from*(1+0.02), breakcol="black", style="slash")
axis(2, at=from)



# easy barplot
rotate_x <- function(data, column_to_plot, labels_vec, rot_angle) {
    plt <- barplot(data[[column_to_plot]], col='steelblue', xaxt="n")
    text(plt, par("usr")[3], labels = labels_vec, srt = rot_angle, adj = c(1.1,1.1), xpd = TRUE, cex=0.6) 
}
rotate_x(mtcars, 'mpg', row.names(mtcars), 45)







# plot graph for MM equation
kcat2=10
E2=8
KM3=15
M3=0:200
v=kcat2*E2*M3/(KM3+M3)
V <- c()
for (i in 1:length(M3)){
  S <- M3[i]
  #print(S)
  v0 <- kcat2*E2*S/(KM3+S)
  print(v0)
  V[i] <- v0
}

plot(M3, V, type="l", col="black", pch="o", lty=1, ylim=c(0,90), ylab="reaction rate (mmol/min)", xlab = "Substrate concentration (mMol)" )



# enzyme
kcat2=10
E2=6
KM3=15
M3=0:200
v=kcat2*E2*M3/(KM3+M3)
EE1 <- c()
for (i in 1:length(M3)){
  S <- M3[i]
  #print(S)
  v0 <- kcat2*E2*S/(KM3+S)
  print(v0)
  EE1[i] <- v0
}
plot(M3, EE1)

kcat2=10
E2=9
KM3=15
M3=0:200
v=kcat2*E2*M3/(KM3+M3)
EE2 <- c()
for (i in 1:length(M3)){
  S <- M3[i]
  #print(S)
  v0 <- kcat2*E2*S/(KM3+S)
  print(v0)
  EE2[i] <- v0
}
plot(M3, EE2)

plot(M3, EE1, type="l", lwd = 2, col="black", pch="o", lty=1, ylim=c(0,90), ylab="reaction rate (mmol/min)", xlab = "Substrate concentration (mMol)" )
lines(M3, V, type="l",  lwd = 2, col="red", pch="o",lty=2)
lines(M3, EE2, type="l", lwd = 2, col="blue", pch="o",lty=2)
 

?lines



########### KM3
kcat2=10
E2=8
KM3=3
M3=0:200
v=kcat2*E2*M3/(KM3+M3)
V2 <- c()
for (i in 1:length(M3)){
  S <- M3[i]
  #print(S)
  v0 <- kcat2*E2*S/(KM3+S)
  print(v0)
  V2[i] <- v0
}
plot(M3, V2)


kcat2=10
E2=8
KM3=9
M3=0:200
v=kcat2*E2*M3/(KM3+M3)
V3 <- c()
for (i in 1:length(M3)){
  S <- M3[i]
  #print(S)
  v0 <- kcat2*E2*S/(KM3+S)
  print(v0)
  V3[i] <- v0
}
plot(M3, V3)


plot(M3, V, type="l", col="black", pch="o", lty=1, ylim=c(0,90), ylab="reaction rate (mmol/min)", xlab = "Substrate concentration (mMol)" )
lines(M3, V2, col="blue",lty=2)
lines(M3, V3, col="red",lty=2)








p <- ggplot(score_of_different_yeastGEM, aes(x=Model, y=Score, shape=Paramter, colour=Paramter)) + geom_line(size=1.5) +
  geom_point(size=6) + 
  xlab('\nVersion of yeastGEM') +
  ylab('Score from memote test\n') +
  theme_bw() +
  theme(legend.position = c(0.85, 0.2)) +
  theme(axis.text=element_text(size=20,face="bold", family="Arial"),
        axis.title=element_text(size=24,face="bold", family="Arial"),
        legend.text = element_text(size=13,face="bold", family="Arial")) +
  ggtitle('') +
  theme(panel.background = element_rect(fill = "white", color="black", size = 1))




# small task
library(readxl)
library(hongR)
library(stringr)
parameter_m <- read_excel("parameter_sce.xlsx", 
                            sheet = "ec_with_measured_km")

parameter_gem <- read_excel("parameter_sce.xlsx", 
                            sheet = "ec_yeast_GEM")
parameter_gem$rr <- NA
parameter_gem$rr <- paste("r", 1:length(parameter_gem$`EC Number`),sep = "_")
parameter_gem$`EC Number`[67] <- "2.3"
parameter_gem0 <- splitAndCombine(parameter_gem$`EC Number`, parameter_gem$rr, ";")
unique_ec_gem <- unique(parameter_gem0$v1)

measure_gem <- str_trim(unique(parameter_m$EC), side="both")

length(intersect(unique_ec_gem, measure_gem))



# small task for the review paper
models_table <- read_excel("models_table.xlsx")
organism_year <- models_table[, c("Organism","Year")]
organism_year <- organism_year[!is.na(organism_year$Year),]
Factor <- as.data.frame(table(organism_year$Organism))
colnames(Factor) <- c("organism","GEM_versions")
Factor <- Factor[order(Factor$GEM_versions),]

Factor$organism <-factor(Factor$organism, levels=Factor$organism)
#Result0$pathway <-factor(Result0$pathway, levels=Factor$pathway)
ggplot(data=Factor, aes(x=organism, y=GEM_versions)) + geom_bar(stat="identity",fill = "#FF6666") +
  xlab('') + ylab('Number of GEM version') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = c(0.85, 0.2)) +
  theme(axis.text=element_text(size=16, family="Arial"),
        axis.title=element_text(size=18,family="Arial"),
        legend.text = element_text(size=10, family="Arial")) +
  ggtitle('') +
  theme(panel.background = element_rect(fill = "white", color="black", size = 1)) #+

# second samll figure in review paper
sce_model <- models_table[models_table$Organism=="S. cerevisiae", ]
sce_model <- sce_model[sce_model$Genes !="NA", ]
sce_model <- sce_model[sce_model$Rxns !="NA", ]
sce_model$Genes <- as.numeric(sce_model$Genes)
sce_model$Rxns <- as.numeric(sce_model$Rxns)
sce_model$Mets <- as.numeric(sce_model$Mets)
Factor2 <- sce_model[,c("Model_ID","Genes")]
Factor2 <- Factor2[order(Factor2$Genes),]
sce_model$Model_ID1 <- sce_model$Model_ID
sce_model$Model_ID <-factor(sce_model$Model_ID, levels=Factor2$Genes)
sce_model$model_order <- 1: nrow(sce_model)
model_gene <- sce_model[,c("Genes","model_order")]
model_gene$type <- "Gene"
colnames(model_gene) <- c("Number", "model_order", "type")

model_rxn <- sce_model[,c("Rxns","model_order")]
model_rxn$type <- "Reaction"
colnames(model_rxn) <- c("Number", "model_order", "type")

model_met <- sce_model[,c("Mets","model_order")]
model_met$type <- "Metabolite"
colnames(model_met) <- c("Number", "model_order", "type")

# combine the three dataframe
model_inf <- rbind.data.frame(model_gene, model_rxn, model_met)
model_inf %>%
  ggplot(aes(x=model_order, y=Number, group=type, shape=type, color=type)) +geom_line(size=1) +
  geom_point(size=4) + 
  xlab('Version of yeast GEM') +
  ylab('Number') +
  theme_bw() +
  theme(legend.position = c(0.2, 0.8)) +
  theme(axis.text=element_text(size=20, family="Arial"),
        axis.title=element_text(size=24, family="Arial"),
        legend.text = element_text(size=13, family="Arial")) +
  ggtitle('') +
  theme(panel.background = element_rect(fill = "white", color="black", size = 1))

# stacked and percent stacked barplot
model_inf %>%
  ggplot(aes(fill=type, y=Number, x=model_order)) + 
  geom_bar(position="stack", stat="identity") +
  xlab('') +
  ylab('Number') +
  theme_bw() +
  theme(legend.position = c(0.2, 0.8)) +
  theme(axis.text=element_text(size=20, family="Arial"),
        axis.title=element_text(size=24, family="Arial"),
        legend.text = element_text(size=13, family="Arial")) +
  ggtitle('') +
  theme(panel.background = element_rect(fill = "white", color="black", size = 1))



# a small example
df <- data.frame(x=1:30, y=NA)
ss <- seq(from = 0, to = 3, by = 0.2)
df$y <- sample(ss, 30, replace = TRUE)
df %>%
  ggplot(aes(y= y, x= x)) + 
  geom_bar(position="stack", stat="identity", color="blue", fill="grey") +
  xlab('Residue sites') +
  ylab('dN/dS') +
  #theme_bw() +
  theme(legend.position = c(0.2, 0.8)) +
  theme(axis.text=element_text(size=20, family="Arial"),
        axis.title=element_text(size=24, family="Arial"),
        legend.text = element_text(size=13, family="Arial")) +
  ggtitle('') +
  theme(panel.background = element_rect(fill = "white", color="black", size = 1))


df <- data.frame(x=1:30, y=NA)
ss <- seq(from = 0, to = 3, by = 0.2)
df$y <- sample(ss, 30, replace = TRUE)
df %>%
  ggplot(aes(y= y, x= x)) + 
  geom_bar(position="stack", stat="identity", color="blue", fill="grey") +
  xlab('Site coordinate') +
  ylab('dN/dS') +
  #theme_bw() +
  theme(legend.position = c(0.2, 0.8)) +
  theme(axis.text=element_text(size=20, family="Arial"),
        axis.title=element_text(size=24, family="Arial"),
        legend.text = element_text(size=13, family="Arial")) +
  ggtitle('') +
  theme(panel.background = element_rect(fill = "white", color="black", size = 1))




###
#pheatmap
library(pheatmap)
library(RColorBrewer)

library(grid)
draw_colnames_45 <- function (coln, ...) {
  m = length(coln)
  x = (1:m)/m - 1/2/m
  grid.text(coln, x = x, y = unit(0.96, "npc"), vjust = .5, 
            hjust = 1, rot = 90, gp = gpar(...)) ## 注意缺省值为 'hjust=0' 和'rot=270'
}
#然后将default draw_colnames覆盖
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                  ns=asNamespace("pheatmap"))


library(readxl)
organelle_ratio_combine <- read_excel("organelle_ratio_combine.xlsx")
rownames(organelle_ratio_combine) <- organelle_ratio_combine$organelle
organelle_ratio_combine0 <-  subset(organelle_ratio_combine, select = -c(organelle, C_N_30))
rownames(organelle_ratio_combine0) <- organelle_ratio_combine$organelle




dat <- cbind(matrix(rnorm(120), 30, 40), matrix(sample(15, 120, T), 30))
my.breaks <- c(seq(-1.5, 0, by=0.1), seq(0.1, 1.5, by=0.1)) 
my.colors <- c(colorRampPalette(colors = c("blue", "white"))(length(my.breaks)/2), colorRampPalette(colors = c("white", "orange", "red", "purple"))(length(my.breaks)/2))
pheatmap(organelle_ratio_combine0,
         method = c("pearson"),
         clustering_method = "complete",
         treeheight_row = 40,
         treeheight_col = 40,
         cluster_row = TRUE,
         cluster_col = TRUE,
         show_rownames = T,
         show_colnames = T,
         legend = T,
         fontsize = 3,
         color = my.colors,
         breaks = my.breaks)

