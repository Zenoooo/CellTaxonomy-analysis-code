#cell similarity
library(plyr)
library(stringr)
#setwd("D:/R/cell taxonomy/13_cell_type_database/01_data/03_all/")
options(stringsAsFactors = F)

cell_taxonomy_cell_tissue_v5 = read.table("cell_taxonomy_final.txt",sep="\t",header=T,quote="",fill=T)

cell_taxonomy_3 = cell_taxonomy_cell_tissue_v5

tissue_ID = unique(cell_taxonomy_3[,c("Tissue_standard","Tissue_UberonOntology_ID2")])
gene_ID = unique(cell_taxonomy_3[,c("Cell_Marker","Gene_ENTREZID2")])
cell_ID = unique(cell_taxonomy_3[,c("Specific_Cell_Ontology_ID","Cell_standard")])

all_cell = unique(cell_taxonomy_3[,c("Specific_Cell_Ontology_ID")])


#cell_taxonomy_3 = cell_taxonomy_3[cell_taxonomy_3$Species %in% c("homo sapiens","mus musculus"),]

all_cell_Gene = unique(cell_taxonomy_3[,c("CT_ID","Specific_Cell_Ontology_ID","Gene_ENTREZID2","Species")])


cell_similarity_species = ddply(all_cell_Gene,"Specific_Cell_Ontology_ID",function(df){

result_tmp1 = ddply(df,"Species",function(df2){
	
	Cell = unique(df2$Specific_Cell_Ontology_ID)
	Species = unique(df2$Species)
	Cell = Cell[!is.na(Cell)]
	Species = Species[!is.na(Species)]
	other_cell = all_cell[!all_cell %in% Cell]
	
	other_cell = as.data.frame(other_cell)
	other_cell$score=NA
	
	this_one_cell_to_other_cell = adply(other_cell,1,function(x){	
	cell_this =x[,1]
	Cell_marker = unique(df2$Gene_ENTREZID2)
	other_cell_marker = unique(all_cell_Gene[all_cell_Gene$Specific_Cell_Ontology_ID==cell_this & all_cell_Gene$Species == Species,"Gene_ENTREZID2"])
	
	Cell_marker=Cell_marker[!is.na(Cell_marker)]
	other_cell_marker=other_cell_marker[!is.na(other_cell_marker)]
	
	all_Cell_marker_num = length(unique(c(Cell_marker,other_cell_marker)))
	common_num = length(intersect(Cell_marker, other_cell_marker))
	score = common_num/all_Cell_marker_num
	
	x$score = score
	return(x)
	})
	
	return(this_one_cell_to_other_cell)
	
	})

return(result_tmp1)
})

cell_similarity_species = cell_similarity_species[cell_similarity_species$score>0,]

cell_similarity_species$score = signif(cell_similarity_species$score,2)

cell_similarity_species$Tissue = "All"

write.table(cell_similarity_species,"cell_similarity_species.txt",sep="\t",quote=F,col.names = T,row.names = F)



