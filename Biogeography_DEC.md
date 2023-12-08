
# Biogeography with DEC models

**Author**: Léo-Paul Dagallier  
**Last update**: 2023-08-08

------------------------------------------------------------------------

## BioGeoBEARS package

<http://phylo.wikidot.com/biogeobears#script>

``` r
install.packages("rexpokit")
install.packages("cladoRcpp")
library(devtools)
devtools::install_github(repo="nmatzke/BioGeoBEARS")
```

## Input data

Prepare the output folder:

``` bash
path_to_output="/outputs/Monodoreae_3";
cd $path_to_output
mkdir Biogeo_DEC
```

``` r
path_to_tree = c("data/name_MCC_monodoreae3_monod_pruned.newick")
path_to_tree_beast = c("data/name_MCC_monodoreae3_monod_pruned.tree")
path_to_output = c("outputs/Biogeo_DEC/")
data_suffix <- "Monodoreae_3"
data_prefix <- "Monodoreae_3"
wd = "" # insert here your working directory, if necessary
path_to_treefile <- paste0(wd, path_to_tree)
path_to_treefile_beast <- paste0(wd, path_to_tree_beast)
path_to_geotf <- "data/Biogeo_DEC/monodoreae_3_DEC_ranges.txt"
data_text = "Monodoreae 3"
path_to_disp_multip_fn <- "data/Biogeo_DEC/dispersal_multipliers"
path_to_areas_adj_fn <- "data/Biogeo_DEC/areas_adjacency"
```

⚠️ For the tree used as input in BioGeoBEARS, it is better if the tree has always been manipulated in R as a "phylo" object. Especially if subsetting the tree at a node with `treeio::tree_subset`, make sure to always manipulate the tree as a "phylo" object and not as a "treedata" object (with e.g. `as.phylo()` function, or subsetting the treedata with `tree@phylo`). From my experience, it seems `treeio::tree_subset` messes the tree structure in such a way that `BioGeoBEARS::prt()` (used to represent the tree as a tabular object in BioGeoBEARS) do not correctly associate the node number with the correct tip.

## Read the tree

``` r
library(ape)
library(picante)
library(treeio)
tree <-  read.tree(file = path_to_treefile)
tree_beast <-  read.beast(file = path_to_treefile_beast)
tot_time <- max(node.age(tree)$ages)
```

## Packages load

``` r
library(rexpokit)
library(cladoRcpp)
library(BioGeoBEARS)
library(parallel)
setwd(dir = paste0(wd, path_to_output))
```

Prepare a file with the ranges of each species.

``` r
# Look at the raw geography text file:
moref(path_to_geotf)
# Look at your geographic range data:
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=path_to_geotf)
tipranges
# Maximum range size observed:
max(rowSums(dfnums_to_numeric(tipranges@df)))
# Set the maximum number of areas any species may occupy; this cannot be larger 
# than the number of areas you set up, but it can be smaller.
max_range_size = 3
```

Modify the list of possible ranges:

``` r
# Get your states list (assuming, say, 4-area analysis, with max. rangesize=4)
max_range_size = 3
areas = getareas_from_tipranges_object(tipranges)
#areas = c("A", "B", "C", "D")

# This is the list of states/ranges, where each state/range
# is a list of areas, counting from 0
states_list_0based = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=TRUE)

# How many states/ranges, by default: 163
length(states_list_0based)

# Make the list of ranges
ranges_list = NULL
for (i in 1:length(states_list_0based))
    {    
    if ( (length(states_list_0based[[i]]) == 1) && (is.na(states_list_0based[[i]])) )
        {
        tmprange = "_"
        } else {
        tmprange = paste(areas[states_list_0based[[i]]+1], collapse="")
        }
    ranges_list = c(ranges_list, tmprange)
    }

# Look at the ranges list
ranges_list

# How many states/ranges, by default: 163
length(ranges_list)

# Let's remove some non-adjacent ranges
nonadjacent=c("WE","WM","CM","WEM","CEM", "WCM")
keepTF = ranges_list %in% nonadjacent == FALSE

ranges_list_NEW = ranges_list[keepTF]
length(ranges_list_NEW)     # now 148

states_list_0based_NEW = states_list_0based[keepTF]
length(states_list_0based_NEW)     # now 148
```

## DEC analysis

Run DEC analysis (see http://phylo.wikidot.com/biogeobears#script).

Export the result table to .csv file:

``` r
write.csv(resDEC[["outputs"]]@params_table, file = paste0("DEC_parameters_results_final_", data_suffix, ".csv"))
```

## Plot the ancestral ranges

Prepare the input.

``` r
states_prob = data.frame(resDEC[["ML_marginal_prob_each_state_at_branch_top_AT_node"]])
colnames(states_prob) <- ranges_list_NEW
# In this table:
# - columns are states/ranges
# - rows are nodes, in APE order (tips, then root, then internal)

#  You can see the node numbers in the same APE order with:
trtable = prt(tree, printflag=FALSE)
head(trtable)
tail(trtable)

states_prob$node <-  trtable$node

library(ggtree)
library(tidytree)
library(ggplot2)
library(colorspace)

tree_data <-  full_join(tree, states_prob)
```

Prepare the pie plots.

``` r
colors = c("W" = "#e41a1c", "C" = "#008cff", "E" = "#1da819", "M" = "#ffff33", "WC" = "#984ea3", "CE" = "#14c9a9", "EM" = "#ff7f00", "WCE" = "#a65628") # Color mixing set perso 6 (= clearer set5)
demoplot(colors, "pie")

pies=nodepie(states_prob, cols = 2:9, outline.color = "black",outline.size = 0.2)
pies <- lapply(pies, function(g) g+scale_fill_manual(values = colors))

legend_df = data.frame(name = names(x = colors), position_x = rep(x = 1,8),  position_y = 1:8, color = colors, label = c("West", "Centre", "East", "Madagascar", "West-Centre", "Centre-East", "East-Madagascar", "West-Centre-East"))
```

Custom geological timescale

``` r
library(deeptime)
GTS <- force(epochs)
lmio <-  c("Late Miocene", 11.6300, 5.3330, "L.Mio", "#FFFF66")
mmio <-  c("Middle Miocene", 15.9700, 11.6300, "M.Mio", "#FFFF4D")
emio <-  c("Early Miocene", 23.0300, 15.9700, "E.Mio", "#FFFF33")
GTS_perso <- rbind(GTS[c(1:3),], lmio, mmio, emio, GTS[5,])
GTS_perso$max_age <- as.numeric(GTS_perso$max_age)
GTS_perso$min_age <- as.numeric(GTS_perso$min_age)
```

Final plot

``` r
library(ggpp)
gg =(ggtree(tree_data) +
    geom_point(data = legend_df, aes(x = position_x-25, y = position_y+75, color = name))+
    scale_color_manual(values = legend_df$color, name = "Geographical range", labels = legend_df$label) +
      # uncomment the following line to plot the pies at the nodes (works only when GTS is not plot)
    # geom_inset(pies, width = 0.05, hjust = max+0.065,  vjust = 0.10, reverse_x = F, x = "node")+
      # uncomment the following line to plot node numbers along with the pies (useful for identifying vicariance and sympatry events)
    # geom_nodelab(aes(x=x, label=node), hjust=-1, size=2, color = "black") +
    geom_tiplab(offset = 0.5, size = 4)+
   
    theme_tree2() +
    theme(axis.line.x.bottom = element_line("#bdbdbd"),
          panel.grid.major.x = element_line("#bdbdbd"),
          panel.grid.minor.x = element_line("#f0f0f0"),
          legend.position=c(0.2, 0.92),
          legend.text=element_text(size=15))) %>% revts()+ scale_x_continuous(labels=abs, breaks = c(0,-10,-20,-30), limits = c(-27, 12))+ guides(color = (guide_legend(override.aes = list(size=8))))
gg

df <- tibble::tibble(node=as.numeric(states_prob$node), pies=pies)
gg2 = gg %<+% df
gg3 = gg2 + geom_plot(data = td_filter(node %in% 1:500), mapping=aes(x=x,y=y, label=pies), vp.width=0.027, hjust=0.58, vjust=0.502) + coord_geo(neg = T, pos = "b", dat = GTS_perso, abbrv = F, height = unit(1, "line"), size = 2.7, bord = c(), skip = c("Holocene"), expand = T, center_end_labels = T)

pdf(file = paste0("figures/",data_prefix,"_DEC3_GTS_ok_for_publi",".pdf"), width = 11.8 , height = 19.7)
gg3
dev.off()
```

Check the node numbers equivalency between newick tree and beast tree

``` r
ggtree(tree_data)+
  geom_tiplab(offset = 0.5, size = 4)+
  geom_nodelab(aes(x=x, label=node), hjust=-1, size=2, color = "black")
ggtree(tree_beast)+
  geom_tiplab(offset = 0.5, size = 4)+
  geom_nodelab(aes(x=x, label=node), hjust=-1, size=2, color = "red")
```

Merge the data from tree_beast into tree_data

``` r
tree_beast@data$node <- as.numeric(tree_beast@data$node)
tree_data_full <-  full_join(tree_data@data, tree_beast@data, by = "node")
```

## Extract the vicariance and sympatry splits events

Get the nodes number with vicariance events. The anagenetic events
(range expansion or contraction) are considered along the branches and
not at the nodes.

``` r
tree_data_events <- tree_data_full[, c("node", "height_0.95_HPD", "height")]
tree_data_events$event <-  NA
tree_data_events$direction <- NA
# vicariance 
centre_east <- c(119, 102, 97, 129, 174)
east_west <- c(147)
centre_west <- c(109, 162, 156)
for (n in centre_east){
  tree_data_events$event[which(tree_data_events$node == n)] <- "vicariance"
  tree_data_events$direction[which(tree_data_events$node == n)] <- "centre_east"
}
for (n in east_west){
  tree_data_events$event[which(tree_data_events$node == n)] <- "vicariance"
  tree_data_events$direction[which(tree_data_events$node == n)] <- "east_west"
}
for (n in centre_west){
  tree_data_events$event[which(tree_data_events$node == n)] <- "vicariance"
  tree_data_events$direction[which(tree_data_events$node == n)] <- "centre_west"
}

# founder event
east_mada <- c(112)
for (n in east_mada ){
  tree_data_events$event[which(tree_data_events$node == n)] <- "founder_event"
  tree_data_events$direction[which(tree_data_events$node == n)] <- "east_mada"
}

library(dplyr)
a = bind_cols(tree_data_events$height_0.95_HPD)
tree_data_events$height_0.95_HPD_lower <- as.numeric(a[1,])
tree_data_events$height_0.95_HPD_upper <- as.numeric(a[2,])

# sympatry (subset) = range contraction
# Here for each event, we need to retrieve the 2 nodes around the corresponding branch

widespread_east_from <- c(111, 105, 131, 169, 146, 136)
widespread_east_to <- c(112, 34, 41, 79, 60, 137)

CE_WCE_centre_from <- c(107, 106, 92, 130, 168, 139, 172, 125)
CE_WCE_centre_to <- c(108, 120, 93, 132, 80, 140, 175, 126)

WC_centre_from <-  c(145)
WC_centre_to <-  c(52)

widespread_west_from <- c(117, 154, 144, 173)
widespread_west_to <- c(118, 170, 50, 85)

for (n in widespread_east_from ){
  tree_data_events$event[which(tree_data_events$node == n)] <- "sympatry_subset"
  tree_data_events$direction[which(tree_data_events$node == n)] <- "widespread_east"
  tree_data_events$height_0.95_HPD_upper[which(tree_data_events$node == n)] <- tree_data_events$height[which(tree_data_events$node == n)]
  tree_data_events$height_0.95_HPD_lower[which(tree_data_events$node == n)]  <- tree_data_events$height[which(tree_data_events$node == widespread_east_to[which(widespread_east_from == n)])]
}
for (n in CE_WCE_centre_from ){
  tree_data_events$event[which(tree_data_events$node == n)] <- "sympatry_subset"
  tree_data_events$direction[which(tree_data_events$node == n)] <- "CE_WCE_centre"
  tree_data_events$height_0.95_HPD_upper[which(tree_data_events$node == n)] <- tree_data_events$height[which(tree_data_events$node == n)]
  tree_data_events$height_0.95_HPD_lower[which(tree_data_events$node == n)]  <- tree_data_events$height[which(tree_data_events$node == CE_WCE_centre_to[which(CE_WCE_centre_from == n)])]
}
for (n in WC_centre_from ){
  tree_data_events$event[which(tree_data_events$node == n)] <- "sympatry_subset"
  tree_data_events$direction[which(tree_data_events$node == n)] <- "WC_centre"
  tree_data_events$height_0.95_HPD_upper[which(tree_data_events$node == n)] <- tree_data_events$height[which(tree_data_events$node == n)]
  tree_data_events$height_0.95_HPD_lower[which(tree_data_events$node == n)]  <- tree_data_events$height[which(tree_data_events$node == WC_centre_to[which(WC_centre_from == n)])]
}
for (n in widespread_west_from ){
  tree_data_events$event[which(tree_data_events$node == n)] <- "sympatry_subset"
  tree_data_events$direction[which(tree_data_events$node == n)] <- "widespread_west"
  tree_data_events$height_0.95_HPD_upper[which(tree_data_events$node == n)] <- tree_data_events$height[which(tree_data_events$node == n)]
  tree_data_events$height_0.95_HPD_lower[which(tree_data_events$node == n)]  <- tree_data_events$height[which(tree_data_events$node == widespread_west_to[which(widespread_west_from == n)])]
}

# range expansion
east_CE_WCE_from <- c(171, 138, 124, 91, 135)
east_CE_WCE_to <- c(172, 139, 125, 92, 136)

centre_WC_from <- c(143, 161, 94, 93, 96, 108)
centre_WC_to <- c(144, 162, 1, 3, 4, 109)

centre_CE_from <-  c(127)
centre_CE_to <-  c(128)

west_WC_from <- c(170, 118)
west_WC_to <- c(81, 25)

for (n in east_CE_WCE_from ){
  tree_data_events$event[which(tree_data_events$node == n)] <- "range_expansion"
  tree_data_events$direction[which(tree_data_events$node == n)] <- "east_CE_WCE"
  tree_data_events$height_0.95_HPD_upper[which(tree_data_events$node == n)] <- tree_data_events$height[which(tree_data_events$node == n)]
  tree_data_events$height_0.95_HPD_lower[which(tree_data_events$node == n)]  <- tree_data_events$height[which(tree_data_events$node == east_CE_WCE_to[which(east_CE_WCE_from == n)])]
}
for (n in centre_WC_from ){
    tree_data_events$event[which(tree_data_events$node == n)] <- "range_expansion"
  tree_data_events$direction[which(tree_data_events$node == n)] <- "centre_WC"
  tree_data_events$height_0.95_HPD_upper[which(tree_data_events$node == n)] <- tree_data_events$height[which(tree_data_events$node == n)]
  tree_data_events$height_0.95_HPD_lower[which(tree_data_events$node == n)]  <- tree_data_events$height[which(tree_data_events$node == centre_WC_to[which(centre_WC_from == n)])]
}
for (n in centre_CE_from ){
    tree_data_events$event[which(tree_data_events$node == n)] <- "range_expansion"
  tree_data_events$direction[which(tree_data_events$node == n)] <- "centre_CE"
  tree_data_events$height_0.95_HPD_upper[which(tree_data_events$node == n)] <- tree_data_events$height[which(tree_data_events$node == n)]
  tree_data_events$height_0.95_HPD_lower[which(tree_data_events$node == n)]  <- tree_data_events$height[which(tree_data_events$node == centre_CE_to[which(centre_CE_from == n)])]
}
for (n in west_WC_from ){
    tree_data_events$event[which(tree_data_events$node == n)] <- "range_expansion"
  tree_data_events$direction[which(tree_data_events$node == n)] <- "west_WC"
  tree_data_events$height_0.95_HPD_upper[which(tree_data_events$node == n)] <- tree_data_events$height[which(tree_data_events$node == n)]
  tree_data_events$height_0.95_HPD_lower[which(tree_data_events$node == n)]  <- tree_data_events$height[which(tree_data_events$node == west_WC_to[which(west_WC_from == n)])]
}

tree_data_events$height_0.95_HPD_lower <- as.numeric(tree_data_events$height_0.95_HPD_lower)
tree_data_events$height_0.95_HPD_upper <- as.numeric(tree_data_events$height_0.95_HPD_upper)
```

## Plot the vicariance and sympatry events

``` r
library(ggrepel)

colors <- c("centre_east" = "#14c9a9", "east_west" = "#a65628", "centre_west" = "#984ea3", "east_mada" = "#ffff33", "widespread_east" = "#1da819", "CE_WCE_centre" = "#008cff", "WC_centre" = "#008cff", "widespread_west" = "#e41a1c", "east_CE_WCE" = "#1da819", "centre_WC" = "#008cff", "centre_CE" = "#008cff", "west_WC" = "#e41a1c")
demoplot(colors, "pie")

labels <- c("centre_east" = "Centre-East", "east_west" = "East-West", "centre_west" = "Centre-West", "east_mada" = "East to Madagascar", "widespread_east" = "Widespread to East", "CE_WCE_centre" = "Centre-East or West-Centre-East to Centre", "WC_centre" = "West-Centre to Centre", "widespread_west" = "Widespread to West", "east_CE_WCE" = "East to Widespread", "centre_WC" = "Centre to West-Centre", "centre_CE" = "Centre to Centre-East", "west_WC" = "West to West-Centre")

tree_data_events$event[which(tree_data_events$event =="vicariance")] <- "Vicariance"
tree_data_events$event[which(tree_data_events$event =="founder_event")] <- "Founder event"
tree_data_events$event[which(tree_data_events$event =="range_expansion")] <- "Range expansion (dispersion)"
tree_data_events$event[which(tree_data_events$event =="sympatry_subset")] <- "Range contraction (sympatry subset)"

labels_facet = c("vicariance" = "Vicariance", "founder_event" = "Founder event", "sympatry_subset" = "Range contraction (sympatry subset)", "range_expansion" = "Range expansion")

data <- tree_data_events[which(!is.na(tree_data_events$event)),]
data <- arrange(data, event, desc(direction), (as.numeric(height)))
data$nudge_y <- NA
data$levels <-  paste0(data$event, "_", data$direction)
table_lev <-  table(data$levels)

for (lev in unique(data$levels)){
  by <-  0.15
  length <-  length(which(data$levels == lev))
  if (length < 9){
      data$nudge_y[which(data$levels == lev)] <-  round(seq(from = -(by*length/2)+by/2, length.out = length, by = by), 2)    
  } else {
     data$nudge_y[which(data$levels == lev)] <-  c(round(seq(from = -(by*(length-1)/2)+by/2, length.out = length-1, by = by), 2), 0)   
  }
}

g <- ggplot(data)+
  geom_point(aes(x = as.numeric(height), y = direction, color = direction), size = 2.5, position = position_nudge(y = data$nudge_y))+
  geom_segment(aes(xend = height_0.95_HPD_lower, x = height_0.95_HPD_upper, y = direction, yend = direction, color = direction), size = 1.7, alpha = 1, position = position_nudge(y = data$nudge_y))+
  geom_text(aes(x = as.numeric(height_0.95_HPD_lower), y = direction, label = paste0(round(height_0.95_HPD_upper, 2)," - ", round(height_0.95_HPD_lower, 2))), size = 2, position = position_nudge(y = data$nudge_y, x = 0.1), hjust = "left")+
  geom_text(aes(label = event), x = -25, y = Inf, hjust = "left", vjust = 1.7, check_overlap = T)+
  scale_fill_discrete(labels = labels_facet)+
  scale_color_manual(values = colors, labels = labels)+
  scale_y_discrete(labels = labels)+
scale_x_reverse(breaks = seq(from = 0, to = 27, by = 1), limits = c(25,-1), minor_breaks = seq(from = 0.5, to = 26.5, by = 1))+
  facet_grid(rows = vars(event), scales = "free_y", space = "free_y", labeller = labeller(event = labels_facet))+
    theme_bw() +
  theme(legend.position = "none",
        panel.grid.major.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank())+
  labs(x = "My before present", y ="")
g
# ggsave(filename = paste0("FIGURE_DEC_Splitting_events_Monodoreae_3_test.pdf"), units = "cm", width = 30, height = 20)
```

With geological timescale

``` r
library(deeptime)
GTS <- force(epochs)
lmio <-  c("Late Miocene", 11.6300, 5.3330, "L.Mio", "#FFFF66")
mmio <-  c("Middle Miocene", 15.9700, 11.6300, "M.Mio", "#FFFF4D")
emio <-  c("Early Miocene", 23.0300, 15.9700, "E.Mio", "#FFFF33")
GTS_perso <- rbind(GTS[c(1:3),], lmio, mmio, emio, GTS[5,])
GTS_perso$max_age <- as.numeric(GTS_perso$max_age)
GTS_perso$min_age <- as.numeric(GTS_perso$min_age)

g2 = g+ coord_geo(neg = F, pos = "b", dat = GTS_perso, abbrv = F, height = unit(1, "line"), size = 3, bord = c(), skip = c("Holocene"), expand = T, center_end_labels = T)

pdf(file = paste0("figures/","FIGURE_DEC_Splitting_events_Monodoreae_3_GTS_ok_for_publi.pdf"), width = 11.8 , height = 7.9)
g2
dev.off()
```
