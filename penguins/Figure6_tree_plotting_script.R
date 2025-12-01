library(cli)
library(treeio)
library(ggplot2)
library(ggtree)
library(ggrepel)
library(ape)
library(tidytree)
library(deeptime)

# Load all paths from external file
source("penguin_paths.R")


# ### The path to tree file containing only the tree maximizing the posterior
tree_path<- file.path(base_dir, "morph_at_start/penguins_inf_morph_at_start.map.tree.txt")

sRangesPath

################################################################################
################################################################################
############################ Overriding the functions   ########################
################################################################################
################################################################################

fortify.phylo.range <- function (model, data, layout = "rectangular", ladderize = TRUE, 
                                right = FALSE, branch.length = "branch.length", mrsd = NULL, 
                                as.Date = FALSE, yscale = "none", root.position = 0, ...) 
{
  x <- as.phylo(model)
  label<-x$tip.label[x$edge[,2][x$edge[,2]<=length(x$tip.label)]]
  if (ladderize == TRUE) {
    x <- ladderize(x, right = right)
  }
  if (!is.null(x$edge.length)) {
    if (anyNA(x$edge.length)) {
      warning("'edge.length' contains NA values...\n## setting 'edge.length' to NULL automatically when plotting the tree...")
      x$edge.length <- NULL
    }
  }
  if (layout %in% c("equal_angle", "daylight", "ape")) {
    res <- layout.unrooted(model, layout.method = layout, 
                           branch.length = branch.length, ...)
  }
  else {
    dd <- get.data(model) %>% arrange(node)
    
    ypos <- getYcoord.range(x, node.orientation=dd$orientation, match="ancestor",
                            tip.order=label)
    N <- Nnode(x, internal.only = FALSE)
    if (is.null(x$edge.length) || branch.length == "none") {
      if (layout == "slanted") {
        sbp <- .convert_tips2ancestors_sbp(x, include.root = TRUE)
        xpos <- getXcoord_no_length_slanted(sbp)
        ypos <- getYcoord_no_length_slanted(sbp)
      }
      else {
        xpos <- getXcoord_no_length(x)
      }
    }
    else {
      xpos <- getXcoord(x)
    }
    xypos <- tibble::tibble(node = 1:N, x = xpos + root.position, 
                            y = ypos)
    df <- as_tibble(model) %>% mutate(isTip = !.data$node %in% 
                                        .data$parent)
    res <- full_join(df, xypos, by = "node")
  }
  res <- calculate_branch_mid(res, layout = layout)
  if (!is.null(mrsd)) {
    res <- scaleX_by_time_from_mrsd(res, mrsd, as.Date)
  }
  if (layout == "slanted") {
    res <- add_angle_slanted(res)
  }
  else {
    res <- calculate_angle(res)
  }
  res <- scaleY(as.phylo(model), res, yscale, layout, ...)
  res <- adjust_hclust_tip.edge.len(res, x)
  class(res) <- c("tbl_tree", class(res))
  attr(res, "layout") <- layout
  return(res)
}

getYcoord.range <- function(tr, step=5, tip.order = NULL, node.orientation = NULL, match = NULL) {
  Ntip <- length(tr[["tip.label"]])
  N <- getNodeNum(tr)
  
  edge <- tr[["edge"]]
  edge_length <- tr[["edge.length"]]
  parent <- edge[,1]
  child <- edge[,2]
  
  label<-tr$tip.label[edge[,2][edge[,2]<=Ntip]]
  
  cl <- split(child, parent)
  child_list <- list()
  child_list[as.numeric(names(cl))] <- cl
  
  y <- numeric(N)
  if (is.null(tip.order)) {
    tip.idx <- child[child <= Ntip]
    y[tip.idx] <- 1:Ntip * step
  } else {
    tip.idx <- 1:Ntip
    y[tip.idx] <- match(tr$tip.label, tip.order) * step
  }
  
  for (t in 1:Ntip){
    t_ <- tip.order[1]
    i <-1
    s = 0
    while(tr$tip.label[t]!=t_){
      tip_n = which(tr$tip.label==t_)
      row_id = which(child==tip_n)
      if (edge_length[row_id]==0){
        s=s+1
      }
      i <- i+1
      t_<-tip.order[i]
    }
    y[t] <- y[t]-step*s
  }
  
  
  
  y[-tip.idx] <- NA
  
  
  pvec <- edge2vec(tr)
  
  currentNode <- 1:Ntip
  while(anyNA(y)) {
    ## pNode <- unique(parent[child %in% currentNode])
    pNode <- unique(pvec[currentNode])
    
    ## piping of magrittr is slower than nested function call.
    ## pipeR is fastest, may consider to use pipeR
    ##
    ## child %in% currentNode %>% which %>% parent[.] %>% unique
    ## idx <- sapply(pNode, function(i) all(child[parent == i] %in% currentNode))
    idx <- sapply(pNode, function(i) all(child_list[[i]] %in% currentNode))
    newNode <- pNode[idx]
    
    if (!is.null(node.orientation)){
      for (n in newNode){
        if (n==31){
          print("")
        }
        ch <- child_list[[n]]
        if (node.orientation[ch[1]]==node.orientation[ch[2]]){
          # choose the one with longer branch length
          l1 <- edge_length[which(edge[,1]==n & edge[,2]==ch[1])]
          if (l1!=0){
            y[n] <- y[ch[1]]
          } else{
            y[n] <- y[ch[2]]
          }
        } else{
          l1 <- edge_length[which(edge[,1]==n & edge[,2]==ch[1])]
          l2 <- edge_length[which(edge[,1]==n & edge[,2]==ch[2])]
          if ((node.orientation[ch[1]]==match & l1!=0) | l2==0){
            y[n] <- y[ch[1]]
          } else {
            y[n] <- y[ch[2]]
          }
        }
      }
      
    } else {
      y[newNode] <- sapply(newNode, function(i) {
        mean(y[child_list[[i]]], na.rm=TRUE)
        ##child[parent == i] %>% y[.] %>% mean(na.rm=TRUE)
      })
    }
    
    for (t in 1:Ntip){
      row_id = which(child==t)
      if (edge_length[row_id]==0){
        p = parent[row_id]
        ch <- child_list[[p]]
        ch <- ch[which(ch!=t)]
        y[t] <- y[ch]
      }
    }
    
    currentNode <- c(currentNode[!currentNode %in% unlist(child_list[newNode])], newNode)
    ## currentNode <- c(currentNode[!currentNode %in% child[parent %in% newNode]], newNode)
    ## parent %in% newNode %>% child[.] %>%
    ##     `%in%`(currentNode, .) %>% `!` %>%
    ##         currentNode[.] %>% c(., newNode)
  }
  
  return(y)
}

### override functions in namespace
unlockBinding("getYcoord", getNamespace("ggtree"))
assign("getYcoord", getYcoord.range, getNamespace("ggtree"))

unlockBinding("fortify.phylo", getNamespace("ggtree"))
assign("fortify.phylo", fortify.phylo.range, getNamespace("ggtree"))

################################################################################
################################################################################
remove_last_substring <- function(str) {
  gsub("_[^_]*$", "", str)
}

t<-read.beast(tree_path)
dd <-get.data(t)
data_range<-dd$range

labels <- as.phylo(t)$tip.label
species <- remove_last_substring(labels)
species_count <- as.data.frame(table(species))

## to avoid labeling the range twice, we make amty labels for first occurence
new_labels <- character(length(labels))
for (l in 1:length(labels)){
  id <- which(species_count[,1]==species[l])
  if (species_count[id, 2]==1 | grepl("last", labels[l])){
    new_labels[l] <- gsub("_", " ", species[l])
  } else {
    new_labels[l] <- ""
  }
}

## create a dataframe where each occurence label maps to newly created label
d <- data.frame(label = labels,
                new_labels = new_labels)

## construct a tidytree object
dd$range_mod<- data_range
tt <- treedata(phylo=as.phylo(t), data=dd)

zero_edge<-which(as.phylo(t)$edge.length==0)
zero_edge_end_node_names<-species[as.phylo(t)$edge[zero_edge, 2]]
neg_id<- which(zero_edge_end_node_names %in% unique(data_range))
sa_nodes <- as.phylo(t)$edge[zero_edge, 2][-neg_id]

## plot simple tree, legend outside, ranges coloured
p1 <- ggplot(tt, aes(x, y, color=range_mod)) + geom_tree(linewidth=2.5) + geom_rootedge(colour="azure4",rootedge = 3, linewidth=2.5) +
  geom_point2(aes(subset=(node %in% dd$node[-which(is.na(dd$range_mod))]), color=range_mod), size=2) + 
  geom_point2(aes(subset=(node %in% sa_nodes)),
              size=5, shape=15)+
  geom_point2(aes(subset=(node %in% which(species %in% c("Marplesornis_novaezealandiae", "Madrynornis_mirandus")))),
              color="red",size=10, shape=17)
  
## if you want simple tip labels
# p1 <- p1 + geom_tiplab(align=F, linetype=NA, size=2)

## smartly spaced out tip/range labels
p3 <-p1 +
  coord_geo(neg = T, xlim = c(-65, 15), ylim = c(-2, NA),
            abbrv = list(T,F), size = list(12, 14),
            height = list(unit(4, "lines"), unit(4, "line")),
            pos = list("bottom","bottom"),
            skip=c("Quaternary", "Holocene", "Meghalayan","Northgrippian","Greenlandian", "Late Pleistocene" ),
            rot=list(0, 0),
            dat=list("epochs","periods"),
            center_end_labels = TRUE) +
  scale_x_continuous(breaks = -rev(sort(c(round(epochs$max_age,1), 0, 10, 20, 30, 40, 50, 60))), labels = rev(sort(c(round(epochs$max_age,1), 0, 10, 20, 30, 40, 50, 60))))


p4 <- p3 %<+% d + 
  geom_text_repel(aes(label=new_labels, fontface=3), size=9,
                  xlim = c(NA, Inf), min.segment.length = 0, force = 0.3,
                  nudge_x = 0.5,direction = "both",hjust = 0,segment.size = 0.2) +  
  theme_tree2(legend.position="None")



## need to reverse time
p5 <- revts(p4+geom_vline(alpha=0.2, xintercept = -rev(sort(unique(round(stages$max_age,1))))))+theme(text=element_text(size=50))
ggsave(paste0(figure6_dir,"/penguin_map_plot.pdf"),p5, width=42, height=45, limitsize = FALSE)


