library('biomaRt')
PATH_ROOT <- "D:/Desktop/Northeastern_University/Research/Proteomics/ProteinProteinAssociation/Development"
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

pQTL_protein_path <- paste(PATH_ROOT,"data_sources","pQTL","pQTL_protein.csv",sep="/")
# proteomeHD_df_path <- paste(PATH_ROOT,"data_sources","ProteomeHD","ProteomeHD_v1_1.csv",sep="/")
# human_ensp_df_path <- paste(PATH_ROOT,"data_sources","StringDB","human","9606.protein.links.full.v11.0.txt",sep="/")
df <- read.table(pQTL_protein_path,sep=",",header = TRUE)
# df[c("protein1", "protein2")] <- lapply(df[c("protein1", "protein2")], function(x) gsub("9606.","",x))
ensg_col <- df$ENSG
#protein1_col <- df$protein1
#protein2_col <- df$protein2

protein1_converted_list <- getBM(filters = "ensembl_gene_id", 
                attributes = c("ensembl_gene_id", "uniprotswissprot"),
                values = ensg_col, mart = mart)

protein1_complete <- protein1_converted_list[!(protein1_converted_list$uniprotswissprot==""), ]

protein1_complete_no_dup <- protein1_complete[!duplicated(protein1_complete$ensembl_gene_id), ]

df_joined <- merge(df,protein1_complete_no_dup,by.x="ENSG",by.y="ensembl_gene_id")

df_joined_no_na <- df_joined[!(df_joined$uniprotswissprot==""), ]

df_joined_reordered <- df_joined[c(1,ncol(df_joined),2:(ncol(df_joined)-4))]

write.csv(df_joined_reordered,paste(PATH_ROOT,"data_sources","pQTL","pQTL_converted.csv",sep="/"), row.names=FALSE)

joined_df <- merge(df,protein1_converted_list,by.x="protein1",by.y="ensembl_peptide_id")

protein2_converted_list <- getBM(filters = "ensembl_peptide_id", 
                                 attributes = c("ensembl_peptide_id", "uniprotswissprot"),
                                 values = protein2_col, mart = mart)

joined_df <- merge(joined_df,protein2_converted_list,by.x="protein2",by.y="ensembl_peptide_id")

joined_df <- joined_df[!(joined_df$uniprotswissprot.x=="" | joined_df$uniprotswissprot.y==""), ]
joined_df_reordered <- joined_df[c(2,1,17,18,3:ncol(joined_df)-2)]
drops <- c("protein2.1","protein1.1")
joined_df_dropped <- joined_df_reordered[ , !(names(joined_df_reordered) %in% drops)]
colnames(joined_df_dropped)[which(names(joined_df_dropped) == "uniprotswissprot.x")] <- "protein1_uniprot"
colnames(joined_df_dropped)[which(names(joined_df_dropped) == "uniprotswissprot.y")] <- "protein2_uniprot"

write.csv(joined_df_dropped,paste(PATH_ROOT,"data_sources","StringDB","human","stringDB_protein_interactions.csv",sep="/"), row.names = FALSE)
