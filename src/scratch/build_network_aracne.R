require(minet)
require(igraph)
require(org.Hs.eg.db)

script_dir      = "/Users/xxx/Documents/workspace/phospho_network/script_files"
mapk_dir        = '/Users/xxx/Documents/workspace/phospho_network/RAWDATA/mapk_info'
interaction_dir = '/Users/xxx/Documents/workspace/phospho_network/script_files/interaction_data_processed'

literal_file    = 'network_all_gene.csv'
kegg_file       = 'mapk_network.csv'

kegg_network    <- read.csv(paste(mapk_dir,kegg_file,sep = '/'),as.is = T,row.names = 1)
annotation.col1 <- select(org.Hs.eg.db, keys=colnames(kegg_network), columns = c('SYMBOL'), keytype="UNIPROT")
mapk_genes      <- annotation.col1$SYMBOL
kegg_network2   <- kegg_network
colnames(kegg_network2) <- mapk_genes
rownames(kegg_network2) <- mapk_genes

literal_interact  <- read.csv(paste(interaction_dir,literal_file,sep = '/'),header = T,as.is = T)
#remove self interaction
literal_interact2 <- literal_interact[literal_interact$geneA != literal_interact$geneB,]
literature_mapk   <- literal_interact2[literal_interact2$geneA %in% mapk_genes & literal_interact2$geneB %in% mapk_genes,]

g_literal <- graph_from_edgelist(as.matrix(literature_mapk), directed = TRUE)
g_kegg    <- graph_from_adjacency_matrix(as.matrix(kegg_network2))

literal_network  <- as.matrix(as_adjacency_matrix(g_literal))
literal_network2 <- literal_network[mapk_genes,mapk_genes]
overlap_edges    <- sum(literal_network2 == 1 & kegg_network2 == 1)/sum(kegg_network2 == 1)
unique_edges     <- sum(literal_network2 == 1 & kegg_network2 == 0)/sum(literal_network2 == 1)
print(c(overlap_edges,unique_edges))

