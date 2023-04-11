library(plyr)
library(stringr)
#setwd("D:/R/cell taxonomy/13_cell_type_database/01_data/03_all/")
options(stringsAsFactors = F)

cell_taxonomy_cell_tissue_v5 = read.table("cell_taxonomy_final.txt",sep="\t",header=T,quote="",fill=T)

cell_taxonomy_3 = cell_taxonomy_cell_tissue_v5


cell_similarity_species = read.table("cell_similarity_species.txt",sep="\t",header=T,quote="",fill=T)
cell_similarity_all_species = read.table("cell_similarity_all_species.txt",sep="\t",header=T,quote="",fill=T)
cell_similarity_Tissue_UberonOntology_ID2 = read.table("cell_similarity_Tissue_UberonOntology_ID2.txt",sep="\t",header=T,quote="",fill=T)
cell_similarity_Tissue_species = read.table("cell_similarity_Tissue_species.txt",sep="\t",header=T,quote="",fill=T)


tissue_ID = unique(cell_taxonomy_3[,c("Tissue_standard","Tissue_UberonOntology_ID2")])
cell_ID = unique(cell_taxonomy_3[,c("CT_ID","Specific_Cell_Ontology_ID","Cell_standard")])

cell_similarity_species2 = merge(cell_similarity_species,cell_ID,by.x="Specific_Cell_Ontology_ID",by.y="Specific_Cell_Ontology_ID")
cell_similarity_species2$Tissue_UberonOntology_ID2 = NA
colnames(cell_similarity_species2)[5] = "Tissue_standard"

cell_similarity_all_species2 = merge(cell_similarity_all_species,cell_ID,by.x="Specific_Cell_Ontology_ID",by.y="Specific_Cell_Ontology_ID")
cell_similarity_all_species2$Tissue_UberonOntology_ID2 = NA
colnames(cell_similarity_all_species2)[5] = "Tissue_standard"

cell_similarity_Tissue_UberonOntology_ID2.1 = merge(cell_similarity_Tissue_UberonOntology_ID2,cell_ID,by.x="Specific_Cell_Ontology_ID",by.y="Specific_Cell_Ontology_ID")
cell_similarity_Tissue_UberonOntology_ID2.2 = merge(cell_similarity_Tissue_UberonOntology_ID2.1,tissue_ID,by.x="Tissue_UberonOntology_ID2",by.y="Tissue_UberonOntology_ID2")


cell_similarity_Tissue_species2 = merge(cell_similarity_Tissue_species,cell_ID,by.x="Specific_Cell_Ontology_ID",by.y="Specific_Cell_Ontology_ID")
cell_similarity_Tissue_species3 = merge(cell_similarity_Tissue_species2,tissue_ID,by.x="Tissue_UberonOntology_ID2",by.y="Tissue_UberonOntology_ID2")

result = rbind(cell_similarity_species2,cell_similarity_all_species2,cell_similarity_Tissue_UberonOntology_ID2.2,cell_similarity_Tissue_species3)

result2 = merge(result,cell_ID,by.x="other_cell",by.y="Specific_Cell_Ontology_ID")
colnames(result2)[1]=c("other_cell_ID")
colnames(result2)[7]="Cell_standard"
colnames(result2)[10]="Other_Cell_standard"

write.table(result2,"cell_cell_similarity.txt",sep="\t",quote=F,col.names = T,row.names = F)

