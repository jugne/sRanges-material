library(cli)
library(treeio)
library(ggplot2)
library(ggtree)
library(ggrepel)
library(ape)
library(tidytree)
library(deeptime)
options(scipen = 999)


# Load all paths from external file
source("canid_helper.R")

###############################################################################
## LABEL PLACEMENT CONTROLS -- adjust these to fine-tune label positioning   ##
###############################################################################

## -- SA square-point node labels (sampled ancestors inside the tree) --
sa_label_size     <- 12    # text size
sa_nudge_x        <- -1  # horizontal nudge (negative = left)
sa_nudge_y        <- 0.9   # vertical nudge   (positive = above)
sa_force          <- 2     # repulsion between labels
sa_force_pull     <- 0.5   # pull back toward anchor point
sa_hjust          <- 0.5   # 0 = left-aligned, 0.5 = centered, 1 = right-aligned

## -- Internal range labels (ranges fully inside the tree) --
irange_label_size <- 12    # text size
irange_nudge_x    <- -0.3    # horizontal nudge (positive = right)
irange_nudge_y    <- 2.2   # vertical nudge magnitude (sign chosen automatically)
irange_force      <- 4     # repulsion between labels
irange_force_pull <- 0.3   # pull back toward anchor point
irange_hjust      <- 0.5   # horizontal justification

## -- Tip labels (single tips + tip-ending ranges) --
tip_label_size    <- 12    # text size
tip_nudge_x       <- 0.5   # horizontal nudge
tip_nudge_y       <- -0.3   # vertical nudge
tip_force         <- 4     # repulsion between labels
tip_force_pull    <- 0.5   # pull back toward anchor point
tip_hjust         <- 0     # left-aligned (label starts at node)

## -- Node marker sizes --
sa_square_size    <- 8     # size of SA square markers
highlight_tri_size <- 14   # size of red triangle markers

###############################################################################

### Automatically find the MAP tree (tree with highest posterior) from the
### combined log and trees files in the sRanges output directory.

find_map_tree <- function(run_dir, force = FALSE) {
  ## Output path for the extracted MAP tree
  log_file  <- list.files(run_dir, pattern = "\\.combined\\.log$",  full.names = TRUE)[1]
  tree_file <- list.files(run_dir, pattern = "\\.combined\\.trees$", full.names = TRUE)[1]
  out_file  <- sub("\\.combined\\.trees$", ".map.tree.txt", tree_file)

  if (!force && file.exists(out_file)) {
    message("MAP tree already exists: ", out_file)
    return(out_file)
  }

  ## 1. Read the log, find the sample with the highest posterior
  log <- read.table(log_file, header = TRUE, sep = "\t", comment.char = "#")
  map_row    <- which.max(log$posterior)
  map_sample <- as.character(as.numeric(log$Sample[map_row]))
  message("MAP sample: STATE_", map_sample,
          "  (posterior = ", log$posterior[map_row], ")")

  ## 2. Read the trees file and extract the header + matching tree line
  tree_con    <- file(tree_file, open = "r")
  header      <- character()
  in_header   <- TRUE
  map_tree    <- NULL
  while (length(line <- readLines(tree_con, n = 1)) > 0) {
    if (in_header && grepl("^\\s*tree\\s+", line)) {
      in_header <- FALSE              # first tree line: header is complete
    }
    if (in_header) {
      header <- c(header, line)
    }
    if (!in_header && grepl(paste0("^\\s*tree\\s+STATE_", map_sample, "\\s*="), line)) {
      map_tree <- line
      break
    }
  }
  close(tree_con)

  if (is.null(map_tree)) stop("Could not find tree STATE_", map_sample,
                               " in ", tree_file)

  ## 3. Write a minimal NEXUS file with only the MAP tree
  writeLines(c(header, map_tree, "End;"), out_file)
  message("Wrote MAP tree to: ", out_file)
  return(out_file)
}

tree_path <- find_map_tree(sRangesPath)

### Import internal ggtree functions needed by the overrides below.
### This avoids having to download and source the entire ggtree package locally.
getXcoord              <- ggtree:::getXcoord
getXcoord_no_length    <- ggtree:::getXcoord_no_length
getNodeNum             <- ggtree:::getNodeNum
edge2vec               <- ggtree:::edge2vec
calculate_branch_mid   <- ggtree:::calculate_branch_mid
calculate_angle        <- ggtree:::calculate_angle
scaleY                 <- ggtree:::scaleY
adjust_hclust_tip.edge.len <- ggtree:::adjust_hclust_tip_edge_len
layout.unrooted        <- ggtree:::layout.unrooted
scaleX_by_time_from_mrsd   <- ggtree:::scaleX_by_time_from_mrsd
add_angle_slanted      <- ggtree:::add_angle_slanted
getXcoord_no_length_slanted   <- ggtree:::getXcoord_no_length_slanted
getYcoord_no_length_slanted   <- ggtree:::getYcoord_no_length_slanted
.convert_tips2ancestors_sbp   <- ggtree:::.convert_tips2ancestors_sbp
## get.data is exported by tidytree (loaded above), no import needed

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

## construct a tidytree object
dd$range_mod<- data_range
tt <- treedata(phylo=as.phylo(t), data=dd)

## Identify zero-length edges (sampled ancestors + internal range endpoints)
zero_edge<-which(as.phylo(t)$edge.length==0)
zero_edge_end_node_names<-species[as.phylo(t)$edge[zero_edge, 2]]
neg_id<- which(zero_edge_end_node_names %in% unique(data_range))
sa_nodes <- as.phylo(t)$edge[zero_edge, 2][-neg_id]

## Flag range labels (species with both _first and _last, label on _last)
is_range <- species_count$Freq[match(species, species_count$species)] > 1 & grepl("last", labels)

## Identify which range _last nodes are internal to the tree (sampled-ancestor
## ranges) vs ending at a true tip.  Internal range endpoints sit on zero-length
## edges that were filtered into neg_id.
internal_range_nodes <- as.phylo(t)$edge[zero_edge[neg_id], 2]
internal_range_labels <- as.phylo(t)$tip.label[internal_range_nodes]
## is_range_internal: range label AND internal (needs above/below placement)
## is_range_tip:      range label AND terminal (label to the right, like tips)
is_range_internal <- is_range & labels %in% internal_range_labels
is_range_tip      <- is_range & !labels %in% internal_range_labels

## Flag SA tip labels (the square-point sampled ancestors inside the tree)
sa_labels <- as.phylo(t)$tip.label[sa_nodes]

## create a dataframe where each occurence label maps to newly created label
d <- data.frame(label = labels,
                new_labels = new_labels,
                is_range = is_range,
                is_range_internal = is_range_internal,
                is_range_tip = is_range_tip,
                is_sa = labels %in% sa_labels)

## plot simple tree, legend outside, ranges coloured
p1 <- ggplot(tt, aes(x, y, color=range_mod)) + geom_tree(linewidth=2.5) + geom_rootedge(colour="azure4",rootedge = 3, linewidth=2.5) +
  geom_point2(aes(subset=(node %in% dd$node[-which(is.na(dd$range_mod))]), color=range_mod), size=2) + 
  geom_point2(aes(subset=(node %in% sa_nodes)),
              size=sa_square_size, shape=15)
  
## if you want simple tip labels
# p1 <- p1 + geom_tiplab(align=F, linetype=NA, size=2)

## Compute per-label nudge direction for internal range labels.
## For each internal-range _last node, check whether the nearest tree edge is
## above or below, and nudge the label in the direction with more free space.
plot_data <- p1$data
all_y <- sort(unique(plot_data$y))

range_nudge <- rep(0, nrow(d))
for (i in which(is_range_internal)) {
  node_idx <- which(plot_data$label == labels[i])
  if (length(node_idx) == 0) next
  yval <- plot_data$y[node_idx[1]]
  # distance to nearest neighbour above and below
  above <- all_y[all_y > yval]
  below <- all_y[all_y < yval]
  gap_above <- if (length(above) > 0) min(above) - yval else Inf
  gap_below <- if (length(below) > 0) yval - min(below) else Inf
  # nudge towards the side with more room
  range_nudge[i] <- if (gap_below >= gap_above) -0.8 else 0.8
}
d$range_nudge <- range_nudge

## ggplot() does not set the "ggtree" subclass (unlike ggtree()), so
## %<+% has no applicable method.  Add the class manually before calling it.
class(p1) <- c("ggtree", class(p1))
p1 <- p1 %<+% d

## smartly spaced out tip/range labels
p3 <-p1 +
    coord_geo(neg = T, xlim = c(-37, 10), ylim = c(-2, NA),
              abbrv = list(T,F), size = list(8, 10),
              height = list(unit(4, "lines"), unit(4, "line")),
              pos = list("bottom", "bottom"),
              skip=c("Quaternary", "Holocene", "Meghalayan","Northgrippian","Greenlandian", "Late Pleistocene" ),
              rot=list(0, 0),
              dat=list("epochs","periods"),
              center_end_labels = TRUE) +
    scale_x_continuous(breaks = -rev(sort(c(round(epochs$max_age,1), 0, 5, 10, 15, 20, 25, 30))), labels = rev(sort(c(round(epochs$max_age,1), 0, 5, 10, 15, 20, 25, 30))))


p4 <- p3 +
  ## Labels for SA square-point nodes: centered, nudged above the point
  geom_text_repel(aes(label = ifelse(is_sa & new_labels != "", new_labels, NA),
                      fontface = 3), size = sa_label_size,
                  min.segment.length = 0, force = sa_force,
                  force_pull = sa_force_pull, max.overlaps = Inf,
                  nudge_x = sa_nudge_x, nudge_y = sa_nudge_y,
                  direction = "both", hjust = sa_hjust, vjust = 0,
                  segment.size = 0.2) +
  ## Labels for internal ranges nudged above
  geom_text_repel(aes(label = ifelse(is_range_internal & range_nudge > 0 & new_labels != "",
                                     new_labels, NA),
                      fontface = 3), size = irange_label_size,
                  min.segment.length = 0, force = irange_force,
                  force_pull = irange_force_pull, max.overlaps = Inf,
                  nudge_x = irange_nudge_x, nudge_y = irange_nudge_y,
                  direction = "y", hjust = irange_hjust,
                  segment.size = 0.2) +
  ## Labels for internal ranges nudged below
  geom_text_repel(aes(label = ifelse(is_range_internal & range_nudge < 0 & new_labels != "",
                                     new_labels, NA),
                      fontface = 3), size = irange_label_size,
                  min.segment.length = 0, force = irange_force,
                  force_pull = irange_force_pull, max.overlaps = Inf,
                  nudge_x = irange_nudge_x, nudge_y = -irange_nudge_y,
                  direction = "y", hjust = irange_hjust,
                  segment.size = 0.2) +
  ## Labels for tip-ending ranges and other single tips: to the right
  geom_text_repel(aes(label = ifelse(!is_sa & !is_range_internal & new_labels != "",
                                     new_labels, NA),
                      fontface = 3), size = tip_label_size,
                  min.segment.length = 0, force = tip_force,
                  force_pull = tip_force_pull, max.overlaps = Inf,
                  nudge_x = tip_nudge_x, nudge_y = tip_nudge_y,
                  direction = "both", hjust = tip_hjust, segment.size = 0.2) +
  theme_tree2(legend.position="None") +
  theme(axis.line.x = element_blank())



## need to reverse time
p5 <- revts(p4)+theme(text=element_text(size=50))
ggsave(paste0(figure9_dir,"/canid_map_plot.pdf"),p5, width=52, height=60, limitsize = FALSE)


