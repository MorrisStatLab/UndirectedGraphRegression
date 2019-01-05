library('coda')
load('fdrm_2_2.rda')
load('Path.rda')
path = path[[6]]
heatmapR = pr.heatmapR; row.names(heatmapR) = 1:nrow(heatmapR)
nG = length(path); ncount = 0
for(m in 1:(nG-1)){
  for(n in (m+1):nG){
    ncount = ncount + 1
    row.names(heatmapR)[ncount] = paste(m,n,sep='_')
  }
}
ug.symbol = as.character(read.table('protein_ab.name.txt', header = F)[path,1])
#normal:1 ; tumor:2; 0.5 tumor: 3
library(igraph)
## prepare for a tumor graph
tumor.igraph = heatmapR[,2]; summary(tumor.igraph)
tumor.igraph.edge = names(tumor.igraph)
tumor.igraph.weight = abs(tumor.igraph)
tumor.igraph.color = tumor.igraph
tumor.igraph.color = rep('white', nrow(pr.heatmapR))
tumor.igraph.color[(tumor.igraph > 0)] = 'green'; 
tumor.igraph.color[(tumor.igraph < 0)] = 'red'; 
tumor.igraph.di = c()

for (i in 1 : length(tumor.igraph.edge)){
  tmp = as.numeric(strsplit(tumor.igraph.edge[i], split='_' )[[1]])
  tumor.igraph.di = c(tumor.igraph.di, tmp)  
}
tumor.igraph.di1 = tumor.igraph.di
for (i in 1 : length(ug.symbol)) tumor.igraph.di1[tumor.igraph.di1 == i] = ug.symbol[i]
tumor.igraph.di2 = graph(tumor.igraph.di1)
tumor.igraph.un = as.undirected(tumor.igraph.di2, mode = 'each')

## prepare for a normal graph
normal.igraph = heatmapR[,1]; summary(normal.igraph)
normal.igraph.edge = names(normal.igraph)
normal.igraph.weight = abs(normal.igraph)
normal.igraph.color = rep('white', nrow(pr.heatmapR))
normal.igraph.color[(normal.igraph > 0)] = 'green'; 
normal.igraph.color[(normal.igraph < 0)] = 'red'; 
normal.igraph.di = c()
for (i in 1 : length(normal.igraph.edge)){
  tmp = as.numeric(strsplit(normal.igraph.edge[i], split='_' )[[1]])
  normal.igraph.di = c(normal.igraph.di, tmp)  
}
normal.igraph.di1 = normal.igraph.di
for (i in 1 : length(ug.symbol)) normal.igraph.di1[normal.igraph.di1 == i] = ug.symbol[i]
normal.igraph.di2 = graph(normal.igraph.di1)
normal.igraph.un = as.undirected(normal.igraph.di2, mode = 'each')

## prepare for a 0.5 tumor (half) graph
half.igraph = heatmapR[,3]; summary(half.igraph)
half.igraph.edge = names(half.igraph)
half.igraph.weight = abs(half.igraph)
half.igraph.color = rep('white', nrow(pr.heatmapR))
half.igraph.color[(half.igraph > 0)] = 'green'; 
half.igraph.color[(half.igraph < 0)] = 'red'; 
half.igraph.di = c()
for (i in 1 : length(half.igraph.edge)){
  tmp = as.numeric(strsplit(half.igraph.edge[i], split='_' )[[1]])
  half.igraph.di = c(half.igraph.di, tmp)  
}
half.igraph.di1 = half.igraph.di
for (i in 1 : length(ug.symbol)) half.igraph.di1[half.igraph.di1 == i] = ug.symbol[i]
half.igraph.di2 = graph(half.igraph.di1)
half.igraph.un = as.undirected(half.igraph.di2, mode = 'each')
ecount(tumor.igraph.un); ecount(normal.igraph.un);ecount(half.igraph.un)
#############################################################################################################################
#find the commone plot for tumor and normal graph and half graph
#use another procedure for finding separate graph
normal.igraph.edgeC = names(normal.igraph)[normal.igraph != 0]; normal.igraph.diC = c()
for (i in 1 : length(normal.igraph.edgeC)) normal.igraph.diC = c(normal.igraph.diC, as.numeric(strsplit(normal.igraph.edgeC[i], split='_' )[[1]]))  
for (i in 1 : length(ug.symbol)) normal.igraph.diC[normal.igraph.diC == i] = ug.symbol[i]
normal.igraph.unC = as.undirected(graph(normal.igraph.diC), mode = 'each')
tumor.igraph.edgeC = names(tumor.igraph)[tumor.igraph != 0]; tumor.igraph.diC = c()
for (i in 1 : length(tumor.igraph.edgeC)) tumor.igraph.diC = c(tumor.igraph.diC, as.numeric(strsplit(tumor.igraph.edgeC[i], split='_' )[[1]]))  
for (i in 1 : length(ug.symbol)) tumor.igraph.diC[tumor.igraph.diC == i] = ug.symbol[i]
tumor.igraph.unC = as.undirected(graph(tumor.igraph.diC), mode = 'each')
half.igraph.edgeC = names(half.igraph)[half.igraph != 0]; half.igraph.diC = c()
for (i in 1 : length(half.igraph.edgeC)) half.igraph.diC = c(half.igraph.diC, as.numeric(strsplit(half.igraph.edgeC[i], split='_' )[[1]]))  
for (i in 1 : length(ug.symbol)) half.igraph.diC[half.igraph.diC == i] = ug.symbol[i]
half.igraph.unC = as.undirected(graph(half.igraph.diC), mode = 'each')

common.igraph.unC = graph.intersection(normal.igraph.unC, tumor.igraph.unC, half.igraph.unC, keep.all.vertices = F)
##################################################################do the commone graph by considering the same sign###################################################################
rule = which(!is.na(normal.igraph)&(abs(normal.igraph) > 0)&
               !is.na(tumor.igraph)&(abs(tumor.igraph) > 0)&
               !is.na(half.igraph)&(abs(half.igraph) > 0)&
               ((tumor.igraph*normal.igraph > 0)&
                  (tumor.igraph*half.igraph > 0)&
                  (normal.igraph*half.igraph > 0)))
rule1 = which(!is.na(normal.igraph)&(abs(normal.igraph) > 0)&
                !is.na(tumor.igraph)&(abs(tumor.igraph) > 0)&
                !is.na(half.igraph)&(abs(half.igraph) > 0)&
                !((tumor.igraph*normal.igraph > 0)&
                (tumor.igraph*half.igraph > 0)&
                (normal.igraph*half.igraph > 0)))
                
common.igraph.edge = row.names(heatmapR)
common.igraph.color = rep('white', nrow(pr.heatmapR))
tmp = rep(F, nrow(pr.heatmapR)); tmp1 = rep(F, nrow(pr.heatmapR)); tmp[rule] = T; tmp1[rule1] = T; rule = tmp; rule1 = tmp1
common.igraph.color[(normal.igraph > 0)&rule] = 'green'; common.igraph.color[(normal.igraph < 0)&rule] = 'red'
common.igraph.color[rule1] = 'black'; 
common.igraph.di = c()
for (i in 1 : length(common.igraph.edge)){
  tmp = as.numeric(strsplit(common.igraph.edge[i], split='_' )[[1]])
  common.igraph.di = c(common.igraph.di, tmp)  
}
common.igraph.di1 = common.igraph.di
for (i in 1 : length(ug.symbol)) common.igraph.di1[common.igraph.di1 == i] = ug.symbol[i]
common.igraph.di2 = graph(common.igraph.di1)
common.igraph.un = as.undirected(common.igraph.di2, mode = 'each')
common.igraph.weight = rep(1, nrow(pr.heatmapR)); common.igraph.weight[!(rule|rule1)] = 0
#
normal.igraph.color[rule] = 'blue'; tumor.igraph.color[rule] = 'blue'; half.igraph.color[rule] = 'blue'

#
tumor.igraph.pi <- tumor.igraph
tumor.igraph.pi.un <- tumor.igraph.un
tumor.igraph.pi.edge <- tumor.igraph.edge
normal.igraph.pi <- normal.igraph
normal.igraph.pi.un <- normal.igraph.un
normal.igraph.pi.edge <- normal.igraph.edge
half.igraph.pi <- half.igraph
half.igraph.pi.un <- half.igraph.un
half.igraph.pi.edge <- half.igraph.edge

##########finally for group union, we derive node subgroup for each node
load('path_afp_immune.rda')
p1 = path[[1]]; p2 = path[[3]]; p4 = path[[4]]
path = path[[6]]
t1 = (path %in% p1)&(!(path %in% p2))&(!(path %in% p4))
t2 = (!(path %in% p1))&((path %in% p2))&(!(path %in% p4))
t3 = (!(path %in% p1))&(!(path %in% p2))&((path %in% p4))
t4 = (!(path %in% p1))&((path %in% p2))&((path %in% p4))
t5 = ((path %in% p1))&(!(path %in% p2))&((path %in% p4))
t6 = ((path %in% p1))&((path %in% p2))&(!(path %in% p4))
t7 = ((path %in% p1))&((path %in% p2))&((path %in% p4))


ptype = rep("0", length(path))
ptype[t1] = "1"; ptype[t2] = "3"; ptype[t3] = "5"; ptype[t6] = "2"; ptype[t4] = "4";ptype[2] = "6"

###################################################################################################################################################################
##Convert to Rcytoscape object
library(RCytoscape)
#high
tumor.cyto.un <- igraph.to.graphNEL(tumor.igraph.un)
tumor.cyto.un <- initEdgeAttribute(graph=tumor.cyto.un, attribute.name='weight', attribute.type='numeric', default.value=1.0)
tumor.cyto.un <- initEdgeAttribute(graph=tumor.cyto.un, attribute.name='colors', attribute.type='char', default.value='red')
for (i in 1:length(tumor.igraph.color)){
  colorname <- tumor.igraph.edge[i]
  colorvalue <- tumor.igraph.color[i]
  nodeA <- ug.symbol[as.numeric(strsplit(colorname, split='_' )[[1]])[1]]
  nodeB <- ug.symbol[as.numeric(strsplit(colorname, split='_' )[[1]])[2]]
  edgeData(tumor.cyto.un, nodeA, nodeB, 'colors') <- colorvalue
}
color.values = c('red', 'blue', 'green', 'white'); colorhex = c('#FF0000', '#0000FF', '#00FF00', '#FFFFFF')
cw <- new.CytoscapeWindow('tumor.cyto', graph=tumor.cyto.un) 
edge.names = as.character(cy2.edge.names (cw@graph))
tumor.width = unname(tumor.igraph.weight); tumor.width1 = tumor.width
tumor.width[tumor.width1 == 0] = 0
tumor.width[(tumor.width1 > 0)&(tumor.width1 <= 0.1)] = 1;tumor.width[(tumor.width1 > 0.1)&(tumor.width1 <= 0.2)] = 2;
tumor.width[(tumor.width1 > 0.2)&(tumor.width1 <= 0.3)] = 3;tumor.width[(tumor.width1 > 0.3)&(tumor.width1 <= 0.4)] = 4;
tumor.width[(tumor.width1 > 0.4)&(tumor.width1 <= 0.5)] = 5;tumor.width[(tumor.width1 > 0.5)&(tumor.width1 <= 0.6)] = 6;
tumor.width[(tumor.width1 > 0.6)&(tumor.width1 <= 0.7)] = 7;tumor.width[(tumor.width1 > 0.7)&(tumor.width1 <= 0.8)] = 8;
tumor.width[(tumor.width1 > 0.8)&(tumor.width1 <= 0.9)] = 9;tumor.width[(tumor.width1 > 0.9)] = 10;

for (i in 1:length(edge.names)) {
  setEdgeLineWidthDirect (cw, edge.names[i], tumor.width[i])
  if(tumor.width[i] == 0) setEdgeOpacityDirect(cw, edge.names[i], 0)
}
setEdgeColorRule(cw, 'colors', color.values, colorhex, mode = 'lookup')
##set node size for different hub genes with hub size
node.names <- ug.symbol
ntumor.width <- igraph::degree(tumor.igraph.unC)
nw.min <- min(ntumor.width); nw.max <- max(ntumor.width)
#node size from 45 to 90
for(i in 1:length(node.names)){
  tmp = (ntumor.width[names(ntumor.width) == node.names[i]] - nw.min)/(nw.max - nw.min)*55 + 50
  if(sum(names(ntumor.width) == node.names[i]) == 0) tmp = 35
  setNodeSizeDirect(cw, node.names[i], tmp)
  setNodeFontSizeDirect(cw, node.names[i], tmp/80*25)
}
g <- cw@graph
nodeDataDefaults(g, attr = "moleculeType") <- "undefined"
attr(nodeDataDefaults(g, attr = "moleculeType"), "class") <- "string"
for(i in 1:length(node.names)){
  nodeData(g, node.names[i], "moleculeType") <- ptype[i]
}
cw@graph <- g
displayGraph(cw)
layoutNetwork(cw, layout.name = "attribute-circle" )
redraw(cw)

#low
normal.cyto.un <- igraph.to.graphNEL(normal.igraph.un)
normal.cyto.un <- initEdgeAttribute(graph=normal.cyto.un, attribute.name='weight', attribute.type='numeric', default.value=1.0)
normal.cyto.un <- initEdgeAttribute(graph=normal.cyto.un, attribute.name='colors', attribute.type='char', default.value='red')
for (i in 1:length(normal.igraph.color)){
  colorname <- normal.igraph.edge[i]
  colorvalue <- normal.igraph.color[i]
  nodeA <- ug.symbol[as.numeric(strsplit(colorname, split='_' )[[1]])[1]]
  nodeB <- ug.symbol[as.numeric(strsplit(colorname, split='_' )[[1]])[2]]
  edgeData(normal.cyto.un, nodeA, nodeB, 'colors') <- colorvalue
}
color.values = c('red', 'blue', 'green', 'white'); colorhex = c('#FF0000', '#0000FF', '#00FF00', '#FFFFFF')
cw <- new.CytoscapeWindow('normal.cyto', graph=normal.cyto.un) 
edge.names = as.character(cy2.edge.names (cw@graph))
normal.width = unname(normal.igraph.weight); normal.width1 = normal.width
normal.width[normal.width1 == 0] = 0
normal.width[(normal.width1 > 0)&(normal.width1 <= 0.1)] = 1;normal.width[(normal.width1 > 0.1)&(normal.width1 <= 0.2)] = 2;
normal.width[(normal.width1 > 0.2)&(normal.width1 <= 0.3)] = 3;normal.width[(normal.width1 > 0.3)&(normal.width1 <= 0.4)] = 4;
normal.width[(normal.width1 > 0.4)&(normal.width1 <= 0.5)] = 5;normal.width[(normal.width1 > 0.5)&(normal.width1 <= 0.6)] = 6;
normal.width[(normal.width1 > 0.6)&(normal.width1 <= 0.7)] = 7;normal.width[(normal.width1 > 0.7)&(normal.width1 <= 0.8)] = 8;
normal.width[(normal.width1 > 0.8)&(normal.width1 <= 0.9)] = 9;normal.width[(normal.width1 > 0.9)] = 10;
for (i in 1:length(edge.names)) {
  setEdgeLineWidthDirect (cw, edge.names[i], normal.width[i])
  if(normal.width[i] == 0) setEdgeOpacityDirect(cw, edge.names[i], 0)
}
setEdgeColorRule(cw, 'colors', color.values, colorhex, mode = 'lookup')
##set node size for different hub genes with hub size
node.names <- ug.symbol
nnormal.width <- igraph::degree(normal.igraph.unC)
nw.min <- min(nnormal.width); nw.max <- max(nnormal.width)
#node size from 45 to 90
for(i in 1:length(node.names)){
  tmp = (nnormal.width[names(nnormal.width) == node.names[i]] - nw.min)/(nw.max - nw.min)*55 + 50
  if(sum(names(nnormal.width) == node.names[i]) == 0) tmp = 35
  setNodeSizeDirect(cw, node.names[i], tmp)
  setNodeFontSizeDirect(cw, node.names[i], tmp/80*25)
}
g <- cw@graph
nodeDataDefaults(g, attr = "moleculeType") <- "undefined"
attr(nodeDataDefaults(g, attr = "moleculeType"), "class") <- "string"
for(i in 1:length(node.names)){
  nodeData(g, node.names[i], "moleculeType") <- ptype[i]
}
cw@graph <- g
#setDefaultNodeSize (cw, 80)
#setDefaultNodeFontSize (cw, 25)
displayGraph(cw)
layoutNetwork(cw, layout.name = "attribute-circle" )
#redraw(cw)

#medium
half.cyto.un <- igraph.to.graphNEL(half.igraph.un)
half.cyto.un <- initEdgeAttribute(graph=half.cyto.un, attribute.name='weight', attribute.type='numeric', default.value=1.0)
half.cyto.un <- initEdgeAttribute(graph=half.cyto.un, attribute.name='colors', attribute.type='char', default.value='red')
for (i in 1:length(half.igraph.color)){
  colorname <- half.igraph.edge[i]
  colorvalue <- half.igraph.color[i]
  nodeA <- ug.symbol[as.numeric(strsplit(colorname, split='_' )[[1]])[1]]
  nodeB <- ug.symbol[as.numeric(strsplit(colorname, split='_' )[[1]])[2]]
  edgeData(half.cyto.un, nodeA, nodeB, 'colors') <- colorvalue
}
color.values = c('red', 'blue', 'green', 'white'); colorhex = c('#FF0000', '#0000FF', '#00FF00', '#FFFFFF')
cw <- new.CytoscapeWindow('half.cyto', graph=half.cyto.un) 
edge.names = as.character(cy2.edge.names (cw@graph))
half.width = unname(half.igraph.weight); half.width1 = half.width
half.width[half.width1 == 0] = 0
half.width[(half.width1 > 0)&(half.width1 <= 0.1)] = 1;half.width[(half.width1 > 0.1)&(half.width1 <= 0.2)] = 2;
half.width[(half.width1 > 0.2)&(half.width1 <= 0.3)] = 3;half.width[(half.width1 > 0.3)&(half.width1 <= 0.4)] = 4;
half.width[(half.width1 > 0.4)&(half.width1 <= 0.5)] = 5;half.width[(half.width1 > 0.5)&(half.width1 <= 0.6)] = 6;
half.width[(half.width1 > 0.6)&(half.width1 <= 0.7)] = 7;half.width[(half.width1 > 0.7)&(half.width1 <= 0.8)] = 8;
half.width[(half.width1 > 0.8)&(half.width1 <= 0.9)] = 9;half.width[(half.width1 > 0.9)] = 10;

for (i in 1:length(edge.names)) {
  setEdgeLineWidthDirect (cw, edge.names[i], half.width[i])
  if(half.width[i] == 0) setEdgeOpacityDirect(cw, edge.names[i], 0)
}
setEdgeColorRule(cw, 'colors', color.values, colorhex, mode = 'lookup')
##set node size for different hub genes with hub size
node.names <- ug.symbol
nhalf.width <- igraph::degree(half.igraph.unC)
nw.min <- min(nhalf.width); nw.max <- max(nhalf.width)
#node size from 45 to 90
for(i in 1:length(node.names)){
  tmp = (nhalf.width[names(nhalf.width) == node.names[i]] - nw.min)/(nw.max - nw.min)*55 + 50
  if(sum(names(nhalf.width) == node.names[i]) == 0) tmp = 35
  setNodeSizeDirect(cw, node.names[i], tmp)
  setNodeFontSizeDirect(cw, node.names[i], tmp/80*25)
}
g <- cw@graph
nodeDataDefaults(g, attr = "moleculeType") <- "undefined"
attr(nodeDataDefaults(g, attr = "moleculeType"), "class") <- "string"
for(i in 1:length(node.names)){
  nodeData(g, node.names[i], "moleculeType") <- ptype[i]
}
cw@graph <- g
#setDefaultNodeSize (cw, 80)
#setDefaultNodeFontSize (cw, 25)
displayGraph(cw)
layoutNetwork(cw, layout.name = "attribute-circle" )
#redraw(cw)



#common
common.cyto.un <- igraph.to.graphNEL(common.igraph.un)
common.cyto.un <- initEdgeAttribute(graph=common.cyto.un, attribute.name='weight', attribute.type='numeric', default.value=1.0)
common.cyto.un <- initEdgeAttribute(graph=common.cyto.un, attribute.name='colors', attribute.type='char', default.value='red')
for (i in 1:length(common.igraph.color)){
  colorname <- common.igraph.edge[i]
  colorvalue <- common.igraph.color[i]
  nodeA <- ug.symbol[as.numeric(strsplit(colorname, split='_' )[[1]])[1]]
  nodeB <- ug.symbol[as.numeric(strsplit(colorname, split='_' )[[1]])[2]]
  edgeData(common.cyto.un, nodeA, nodeB, 'colors') <- colorvalue
}
color.values = c('red', 'black', 'green', 'white'); colorhex = c('#FF0000', '#000000', '#00FF00', '#FFFFFF')
cw <- new.CytoscapeWindow('common.cyto', graph=common.cyto.un) 
edge.names = as.character(cy2.edge.names (cw@graph))
common.width = unname(common.igraph.weight); common.width1 = common.width
common.width[common.width1 == 0] = 0; common.width[(common.width1 > 0)] = 3
for (i in 1:length(edge.names)) {
  setEdgeLineWidthDirect (cw, edge.names[i], common.width[i])
  if(common.width[i] == 0) setEdgeOpacityDirect(cw, edge.names[i], 0)
}
setEdgeColorRule(cw, 'colors', color.values, colorhex, mode = 'lookup')
##set node size for different hub genes with hub size
node.names <- ug.symbol
ncommon.width <- igraph::degree(common.igraph.unC)
nw.min <- min(ncommon.width); nw.max <- max(ncommon.width)
#node size from 45 to 90
for(i in 1:length(node.names)){
  tmp = (ncommon.width[names(ncommon.width) == node.names[i]] - nw.min)/(nw.max - nw.min)*55 + 50
  if(sum(names(ncommon.width) == node.names[i]) == 0) tmp = 35
  setNodeSizeDirect(cw, node.names[i], tmp)
  setNodeFontSizeDirect(cw, node.names[i], tmp/80*25)
}
#setDefaultNodeSize (cw, 80)
#setDefaultNodeFontSize (cw, 25)
g <- cw@graph
nodeDataDefaults(g, attr = "moleculeType") <- "undefined"
attr(nodeDataDefaults(g, attr = "moleculeType"), "class") <- "string"
for(i in 1:length(node.names)){
  nodeData(g, node.names[i], "moleculeType") <- ptype[i]
}
cw@graph <- g
displayGraph(cw)
layoutNetwork(cw, layout.name = "attribute-circle" )
#redraw(cw)
