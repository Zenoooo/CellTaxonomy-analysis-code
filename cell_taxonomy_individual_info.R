
library(plyr)
library(stringr)
#setwd("D:/R/cell taxonomy/13_cell_type_database/01_data/03_all/")
options(stringsAsFactors = F)

cell_taxonomy_cell_tissue_v5 = read.table("cell_taxonomy_final.txt",sep="\t",header=T,quote="",fill=T)
cell_taxonomy_cell_tissue_v5[cell_taxonomy_cell_tissue_v5 == ""] <- "-"
cell_taxonomy_cell_tissue_v5[cell_taxonomy_cell_tissue_v5 == "-"] <- "-"
cell_taxonomy_cell_tissue_v5[is.na(cell_taxonomy_cell_tissue_v5)] <- "-"


##########################################
### individual cell information ##########
##########################################
#CL_info3 = read.table("CL_info3.txt",header=T,sep="\t",fill=T,quote="")
CL_info3 = read.table("CL_info_final_tree.txt",header=T,sep="\t",fill=T,quote="")

CL_info3[CL_info3 == ""] <- "-"
CL_info3[CL_info3 == "-"] <- "-"
CL_info3[is.na(CL_info3)] <- "-"

#CL_info3$Class.ID = str_replace_all(CL_info3$Class.ID,"CL_","CL:")

#Definitions
cell_description = unique(cell_taxonomy_cell_tissue_v5[,c("CT_ID","Cell_standard","Specific_Cell_Ontology_ID","Is_New_Cell_Type","Cell_alias_change")])
#1262

cell_description2.1 = merge(cell_description,CL_info3,by.x="Specific_Cell_Ontology_ID",by.y="Class.ID")
#782


cell_description2.1$New_Cell_Type_Annotation = cell_description2.1$Definitions

cell_description2.1 = cell_description2.1[,c("Specific_Cell_Ontology_ID","CT_ID","Cell_standard","New_Cell_Type_Annotation","Is_New_Cell_Type","Cell_alias_change")]

cell_description2.1[cell_description2.1$New_Cell_Type_Annotation == "","New_Cell_Type_Annotation"]="-"

#cell_description2.2 = cell_description[!cell_description$Specific_Cell_Ontology_ID %in% CL_info3$Class.ID,]


cell_description_result = cell_description2.1

write.table(cell_description_result,"cell_description_result.txt",sep="\t",quote=F,col.names = T,row.names = F)
#used for database, individual cell description

cell_table = unique(cell_taxonomy_cell_tissue_v5[,c("CT_ID","Specific_Cell_Ontology_ID","Species","Tissue_standard","Tissue_UberonOntology_ID","Tissue_UberonOntology_ID2","Disease_Type","Disease_Ontology_ID","Disease_Ontology_ID2","Additional_Information2","Cell_Marker","Gene_ENTREZID","Gene_ENTREZID2","PMID")])

#cell_table2 = adply(cell_table,1,function(df){df$vocabulary = paste(df$Specific_Cell_Ontology_ID,df$Species,df$Tissue_UberonOntology_ID2,df$Disease_Ontology_ID2,sep = "_");return(df)})

#resource_info2 = ddply(cell_table2,"vocabulary",function(df){

#df$Source = paste(unique(df$Additional_Information2),collapse=",")
#return(df)

#})
#resource_info_result = unique(resource_info2[,c("vocabulary","Source")])

#cell_table3 = merge(cell_table2,resource_info_result,by.x="vocabulary",by.y="vocabulary")
#3925

cell_table[cell_table$Disease_Type=="","Disease_Type"]="-"

cell_table2 = adply(cell_table,1,function(df){

df$n_NA = length(df[df=="-"])
return(df)

})


write.table(cell_table2,"individual_cell_table.txt",sep="\t",quote=F,col.names = T,row.names = F)
#used for database, individual cell table



cell_marker_info = unique(cell_taxonomy_cell_tissue_v5[,c("CT_ID","Specific_Cell_Ontology_ID","Cell_Marker","Gene_ENTREZID","Gene_ENTREZID2","Species","Tissue_standard","Tissue_UberonOntology_ID","Tissue_UberonOntology_ID2","Disease_Type","Disease_Ontology_ID","Disease_Ontology_ID2","PMID","Additional_Information2","Marker_Resource2")])
#86086

cell_marker_info2 = adply(cell_marker_info,1,function(df){df$vocabulary = paste(df$Gene_ENTREZID2,df$Specific_Cell_Ontology_ID,df$Species,df$Tissue_UberonOntology_ID2,df$Disease_Ontology_ID2,df$PMID,df$Marker_Resource2,sep = "_");return(df)})


resource_info = ddply(cell_marker_info2,"vocabulary",function(df){

df$Source = paste(unique(df$Additional_Information2),collapse=",")
return(df)

})
resource_info_result = unique(resource_info[,c("vocabulary","Source")])

cell_marker_info3 = merge(cell_marker_info2,resource_info_result,by.x="vocabulary",by.y="vocabulary")
#3925


cell_marker_info3[cell_marker_info3$Disease_Type=="","Disease_Type"]="-"
cell_marker_info3[cell_marker_info3$PMID=="","PMID"]="-"

cell_marker_info4 = cell_marker_info3[cell_marker_info3$Specie != "-",]


cell_marker_info4 = adply(cell_marker_info4,1,function(df){

df$n_NA = length(df[df=="-"])
return(df)

})


write.table(cell_marker_info4,"individual_cell_cell_marker.txt",sep="\t",quote=F,col.names = T,row.names = F)
#used for database, individual cell's cell marker table

#################################################
### individual cell marker information ##########
#################################################
marker_description = unique(cell_taxonomy_cell_tissue_v5[,c("Cell_Marker","Gene_ENTREZID","Gene_ENTREZID2","Gene_Alias","Gene_Ensembl_ID","Uniprot","PFAM2","GO2")])
#29952

#n_occur <- data.frame(table(marker_description$Gene_ENTREZID2))
#length(n_occur[n_occur$Freq > 1,"Var1"])
#5294

#marker_description_duplicate = marker_description[marker_description$Gene_ENTREZID2 %in% n_occur[n_occur$Freq > 1,"Var1"],]
#marker_description_no_duplicate = marker_description[marker_description$Gene_ENTREZID2 %in% n_occur[n_occur$Freq == 1,"Var1"],]


#marker_description_duplicate_tmp = adply(marker_description_duplicate,1,function(df){

#df$len = nchar(df$Gene_Alias)
#return(df)

#})

#marker_description_duplicate_tmp = marker_description_duplicate_tmp[with(marker_description_duplicate_tmp, order(-len)), ]
#marker_description_duplicate_remove = marker_description_duplicate_tmp[!duplicated(marker_description_duplicate_tmp$Gene_ENTREZID2),]
#5334 <- 10903
#marker_description_duplicate_remove = marker_description_duplicate_remove[,-which(colnames(marker_description_duplicate_remove) %in% "len") ,]

#marker_description_result = unique(rbind(marker_description_duplicate_remove,marker_description_no_duplicate))
#24308

#length(unique(marker_description_result$Gene_ENTREZID2))
#24138

marker_description_result = marker_description

HomoloGene = read.table("HomoloGene.txt",sep="\t",header=T,quote="",fill=T)

marker_description_result2 = merge(marker_description_result,HomoloGene,by.y="Gene.ID",by.x="Gene_ENTREZID2")

marker_description_result3 = marker_description_result[!marker_description_result$Gene_ENTREZID2 %in% HomoloGene$Gene.ID,]
marker_description_result3$HomoloGene="-"
marker_description_result3$HomoloGene_num="0"

marker_description_result_final = rbind(marker_description_result2,marker_description_result3)


write.table(marker_description_result_final,"marker_description_result.txt",sep="\t",quote=F,col.names = T,row.names = F)
#used for database, individual marker's description

#cell_taxonomy_cell_tissue_v6 = cell_taxonomy_cell_tissue_v5[ , -which(names(cell_taxonomy_cell_tissue_v5) %in% c("Gene_Alias","Gene_Ensembl_ID","Uniprot","PFAM","GO2"))]

#cell_taxonomy_cell_tissue_v7 = unique(merge(cell_taxonomy_cell_tissue_v6,marker_description_result[,c(3:8)],by.x="Gene_ENTREZID2",by.y="Gene_ENTREZID2"))

####### marker table

cell_marker_table = unique(cell_taxonomy_cell_tissue_v5[,c("CT_ID","Cell_Marker","Gene_ENTREZID","Gene_ENTREZID2","Species","Tissue_standard","Tissue_UberonOntology_ID","Tissue_UberonOntology_ID2","Disease_Type","Disease_Ontology_ID","Disease_Ontology_ID2","Cell_standard","Specific_Cell_Ontology_ID","PMID","Additional_Information2","Marker_Resource2")])
#87002

cell_marker_table2 = adply(cell_marker_table,1,function(df){

df$n_NA = length(df[df=="-"])
return(df)

})


write.table(cell_marker_table2,"cell_marker_table.txt",sep="\t",quote=F,col.names = T,row.names = F)
#used for database, individual marker's table

#cell_taxonomy_cell_tissue_v5 = read.table("20210605_cell_taxonomy_final.txt",sep="\t",header=T,quote="",fill=T)

#all_search_table = unique(cell_taxonomy_cell_tissue_v5[,c("Cell_Marker","Gene_ENTREZID","Gene_ENTREZID2","Gene_Alias","Gene_Ensembl_ID","Uniprot","PFAM","GO2","Species","Tissue_standard","Tissue_UberonOntology_ID","Tissue_UberonOntology_ID2","Disease_Type","Disease_Ontology_ID","Disease_Ontology_ID2","Cell_standard","CT_ID","Specific_Cell_Ontology_ID","Is_New_Cell_Type","PMID","Additional_Information2","Marker_Resource2")])

all_search_table = unique(cell_taxonomy_cell_tissue_v5[,c("Cell_Marker","Gene_ENTREZID","Gene_ENTREZID2","Gene_Alias","Gene_Ensembl_ID","Uniprot","PFAM","GO2","Species","Tissue_standard","Tissue_UberonOntology_ID","Tissue_UberonOntology_ID2","Disease_Type","Disease_Ontology_ID","Disease_Ontology_ID2","Cell_standard","CT_ID","Specific_Cell_Ontology_ID","Is_New_Cell_Type","PMID","Additional_Information2","Marker_Resource2","Species_alias","Species_tax_ID","Cell_alias_change")])

all_search_table = adply(all_search_table,1,function(df){

df$n_NA = length(df[df=="-"])
return(df)

})

write.table(all_search_table,"all_search_table.txt",sep="\t",quote=F,col.names = T,row.names = F)



#################################################
### individual Tissue information ###############
#################################################

tissue_description = unique(cell_taxonomy_cell_tissue_v5[,c("Tissue_standard","Tissue_UberonOntology_ID")])
#326

tissue_info = read.table("tissue_UBERON_info_all",header=T,sep="\t",fill=T,quote="")

tissue_info$Class.ID = str_replace_all(tissue_info$Class.ID,"UBERON_","UBERON:")

tissue_tmp = merge(tissue_description,tissue_info,by.x="Tissue_UberonOntology_ID",by.y="Class.ID")
tissue_tmp = tissue_tmp[,c("Tissue_UberonOntology_ID","Tissue_standard","Synonyms","Definitions")]


capFirst <- function(s) {
	if(!is.na(s)) return(paste(toupper(substring(s, 1, 1)), substring(s, 2), sep = ""))
	if(is.na(s)) return(NA)
}

tissue_tmp$Synonyms = capFirst(tissue_tmp$Synonyms)
tissue_tmp$Definitions = capFirst(tissue_tmp$Definitions)

tissue_tmp2 = tissue_description[!tissue_description$Tissue_UberonOntology_ID %in% tissue_info$Class.ID,]
tissue_tmp2$Synonyms = "-"
tissue_tmp2$Definitions = "-"

tissue_description_result = rbind(tissue_tmp,tissue_tmp2)

#n_occur <- data.frame(table(tissue_description_result$Tissue_standard))
#length(n_occur[n_occur$Freq > 1,"Var1"])


#tissue_description_result = tissue_description_result[!(tissue_description_result$Tissue_standard == "Lymph" & tissue_description_result$Tissue_UberonOntology_ID =="-"),]
#tissue_description_result = tissue_description_result[!(tissue_description_result$Tissue_standard == "Placenta" & tissue_description_result$Tissue_UberonOntology_ID =="-"),]

tissue_description_result[tissue_description_result$Synonyms == "","Synonyms"]="-"
tissue_description_result[tissue_description_result$Definitions == "","Definitions"]="-"

tissue_description_result[is.na(tissue_description_result$Synonyms),"Synonyms"]="-"
tissue_description_result[is.na(tissue_description_result$Definitions),"Definitions"]="-"


write.table(tissue_description_result,"tissue_description.txt",sep="\t",quote=F,col.names = T,row.names = F)
#used for database, individual tissue's description

tissue_table = unique(cell_taxonomy_cell_tissue_v5[,c("Tissue_standard","Tissue_UberonOntology_ID","Tissue_UberonOntology_ID2","Cell_standard","CT_ID","Specific_Cell_Ontology_ID","Species","Disease_Type","Disease_Ontology_ID","Disease_Ontology_ID2","Additional_Information2")])
#3907


tissue_table2 = adply(tissue_table,1,function(df){df$vocabulary = paste(df$Specific_Cell_Ontology_ID,df$Species,df$Tissue_UberonOntology_ID2,df$Disease_Ontology_ID2,sep = "_");return(df)})

resource_info2 = ddply(tissue_table2,"vocabulary",function(df){

df$Source = paste(unique(df$Additional_Information2),collapse=",")
return(df)

})
resource_info_result = unique(resource_info2[,c("vocabulary","Source")])

tissue_table3 = merge(tissue_table2,resource_info_result,by.x="vocabulary",by.y="vocabulary")
#4411

tissue_table3 = tissue_table3[tissue_table3$Species!="-",]

tissue_table3 = adply(tissue_table3,1,function(df){

df$n_NA = length(df[df=="-"])
return(df)

})


write.table(tissue_table3,"tissue_table.txt",sep="\t",quote=F,col.names = T,row.names = F)
#used for database, individual tissue's table

