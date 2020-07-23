if (!requireNamespace("RBGL", quietly = TRUE))
  BiocManager::install("RBGL")
library(RBGL)
#Command args
args <- commandArgs(TRUE)
fid_intersection_file <- as.character(args[1])
fid_clusters_file <- as.character(args[2])

#FIDs=read.table(file="/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/test_R_graph_cluster/DIPG.all_edges.tsv", header = TRUE, stringsAsFactors = F)
FIDs=read.table(fid_intersection_file, header = TRUE, stringsAsFactors = F)

# For some reason pairToPair output has duplicates, need to remove them
FIDs.nodup <- FIDs[!duplicated(FIDs ), ]

# need to trick the function by specifying edgemode to be directional
g <- ftM2graphNEL(as.matrix(FIDs.nodup), W=NULL, V=NULL, edgemode="directed")
# convert to undirected graph
edgemode(g) <- "undirected"

# "connections" is the object containing FID clusters
connections <- connectedComp(g)
connections <- lapply(connections, function(x){paste0(x, collapse=",")})
connections <- data.frame(FIDs=matrix(unlist(connections)),stringsAsFactors=FALSE)

#write tsv
#write.table(connections, file="/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/test_R_graph_cluster/FID.clusters.tsv-TMP", quote=FALSE, sep='\t', row.names=FALSE)
write.table(connections, file=fid_clusters_file, quote=FALSE, sep='\t', row.names=FALSE)
