library(plyr)
library(stringr)
#setwd("D:/R/cell taxonomy/13_cell_type_database/01_data/03_all/")
options(stringsAsFactors = F)


human_surface_ID_final2  = read.table("human_surface_ID_final2.txt",sep="\t",header=T,quote="",fill=T)

cell_marker_fisher_test_result = read.table("cell_marker_fisher_test_result.txt",sep="\t",header=T,quote="",fill=T)

cell_taxonomy_5_23_surface = merge(cell_marker_fisher_test_result,human_surface_ID_final2,by.x="Gene_ENTREZID2",by.y="ENTREZID")
#32496

cell_taxonomy_5_23_surface = cell_taxonomy_5_23_surface[cell_taxonomy_5_23_surface$Species == "Homo sapiens",]

length(unique(cell_taxonomy_5_23_surface$Gene_ENTREZID2))
#[1] 3489 surface cell markers

human_surface_marker = unique(cell_taxonomy_5_23_surface[,c("CT_ID","Specific_Cell_Ontology_ID","Gene_ENTREZID2","marker_cell_number","FDR","Protein.class", "Evidence.summary", "HPA.evidence", "UniProt.evidence", "MS.evidence","Antibody", "Reliability..IH.", "Reliability..IF.", "Subcellular.location")])

cell_taxonomy_5_27 = read.table("20210627_cell_taxonomy_final_NA.txt",sep="\t",header=T,quote="",fill=T)

gene_name = unique(cell_taxonomy_5_27[,c("Gene_ENTREZID2","Cell_Marker")])

human_surface_marker2 = merge(human_surface_marker,gene_name,by.x="Gene_ENTREZID2",by.y="Gene_ENTREZID2")

human_surface_marker2$Species = "Homo sapiens"

write.table(human_surface_marker2,"human_surface_marker.txt",sep="\t",quote=F,col.names = T,row.names = F)
