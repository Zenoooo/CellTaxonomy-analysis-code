#cell similarity
library(plyr)
library(stringr)
#setwd("D:/R/cell taxonomy/13_cell_type_database/01_data/03_all/")
options(stringsAsFactors = F)

cell_taxonomy_3 = read.table("cell_taxonomy_final.txt",sep="\t",header=T,quote="",fill=T)

HomoloGene = read.table("HomoloGene_raw.txt",sep="\t",header=T,quote="",fill=T)

cell_ID = unique(cell_taxonomy_3[,c("Specific_Cell_Ontology_ID","Cell_standard")])

all_cell = unique(cell_taxonomy_3[,c("Specific_Cell_Ontology_ID")])


all_cell_Gene = unique(cell_taxonomy_3[,c("CT_ID","Specific_Cell_Ontology_ID","Gene_ENTREZID2","Species","Cell_Marker")])

all_cell_Gene_HomoloGene = merge(all_cell_Gene,HomoloGene,by.x="Gene_ENTREZID2",by.y="original_gene")
table(all_cell_Gene_HomoloGene$Species)


all_cell_Gene_HomoloGene2 = all_cell_Gene[!all_cell_Gene$Gene_ENTREZID2 %in% HomoloGene$original_gene,]
all_cell_Gene_HomoloGene2$Gene.ID=NA
all_cell_Gene_HomoloGene2$Gene.Symbol=NA
all_cell_Gene_HomoloGene2$species=NA

all_cell_Gene_HomoloGene_final = rbind(all_cell_Gene_HomoloGene,all_cell_Gene_HomoloGene2)


cell_type_process = unique(all_cell_Gene_HomoloGene_final[,c("Specific_Cell_Ontology_ID","Species")])

test2 = data.frame(table(cell_type_process$Specific_Cell_Ontology_ID))
cell_type_with_multiple_species = unique(test2[test2$Freq!=1,"Var1"])

all_cell_Gene_HomoloGene2 = all_cell_Gene_HomoloGene_final[all_cell_Gene_HomoloGene_final$Specific_Cell_Ontology_ID %in% cell_type_with_multiple_species, ]



species1_vs_species2_cell_type <- function(data,species1_name,species2_name) {

result = ddply(data,"Specific_Cell_Ontology_ID",function(df){

species1_marker = unique(df[df$Species==species1_name,"Gene_ENTREZID2"])
species1_marker_info = df[df$Species==species1_name,]
species2_marker = unique(df[df$Species==species2_name,"Gene_ENTREZID2"])

#compare species1_marker_orthologous_in_species2 vs species2_marker

#species is the species of orthologous gene
species1_marker_orthologous_in_species2 = unique(df[df$Species==species1_name & df$species==species2_name,"Gene.ID"])

species1_marker_orthologous_common_in_species2 = unique(species1_marker_orthologous_in_species2[species1_marker_orthologous_in_species2 %in% species2_marker])

#check species1_marker_orthologous_common_in_species2 correspond to #number of species1 markers

common_num = length(species1_marker_orthologous_common_in_species2)


#one species1 gene ~ multiple orthologous genes in species2
common_species1_marker_with_species2 = species1_marker_info[species1_marker_info$Gene.ID %in% species1_marker_orthologous_common_in_species2,]

if(dim(common_species1_marker_with_species2)[1]>0){
#output common gene info
tmp = adply(common_species1_marker_with_species2,1,function(df){

each_gene=df[,"Gene_ENTREZID2"]
gene_orthologous_info = common_species1_marker_with_species2[common_species1_marker_with_species2$Gene_ENTREZID2==each_gene,c("Gene_ENTREZID2","Cell_Marker","species","Gene.Symbol","Gene.ID")]

gene_orthologous_info$HomoloGene = paste(gene_orthologous_info$Cell_Marker, gene_orthologous_info$Gene_ENTREZID2, gene_orthologous_info$Gene.Symbol,gene_orthologous_info$Gene.ID, gene_orthologous_info$species,sep="+JS+")

return(gene_orthologous_info)

})

common_marker = unique(tmp[,c("HomoloGene")])
common_marker_number = length(common_marker)

species1_specific_marker = unique(species1_marker_info[!species1_marker_info$Gene.ID %in% species1_marker_orthologous_common_in_species2,"Gene_ENTREZID2"])
species1_specific_marker_num = length(species1_specific_marker)

species2_specific_marker_ID = unique(species2_marker[!species2_marker %in% species1_marker_orthologous_in_species2])
species2_specific_marker = unique(df[df$Gene_ENTREZID2 %in% species2_specific_marker_ID,"Gene_ENTREZID2"])
species2_specific_marker_num  = length(species2_specific_marker)


}

if(dim(common_species1_marker_with_species2)[1]==0){

common_marker = NA
common_marker_number = 0

species1_specific_marker = unique(unique(df[df$Species==species1_name,"Gene_ENTREZID2"]))
species1_specific_marker_num = length(species1_specific_marker)

species2_specific_marker = unique(unique(df[df$Species==species2_name,"Gene_ENTREZID2"]))
species2_specific_marker_num  = length(species2_specific_marker)

}



all_Cell_marker_num = common_marker_number + species1_specific_marker_num + species2_specific_marker_num

score = common_marker_number/all_Cell_marker_num


df$score = score
df$common_marker = paste(common_marker,collapse=", ")
df$common_marker_number = common_marker_number
df$species1_specific_marker = paste(species1_specific_marker,collapse=", ")
df$species1_specific_marker_num = species1_specific_marker_num
df$species2_specific_marker = paste(species2_specific_marker,collapse=", ")
df$species2_specific_marker_num = species2_specific_marker_num


return(df)
})

return(result)
}

cell_similarity_species = ddply(all_cell_Gene_HomoloGene2,"Specific_Cell_Ontology_ID",function(df){

#combination without repetition
#combination = as.data.frame(t(combn(sort(unique(df$Species)), 2)))

#combination with repetition
combination = as.data.frame(expand.grid(rep(list(unique(df$Species)), 2)))

colnames(combination) = c("Species1","Species2")

combination = combination[combination$Species1 != combination$Species2,]

result2 = adply(combination,1,function(species){

species1_name = species[,"Species1"]
species2_name = species[,"Species2"]

result = species1_vs_species2_cell_type(df,species1_name,species2_name)
result = unique(result[,c("CT_ID","Specific_Cell_Ontology_ID","score","common_marker","common_marker_number","species1_specific_marker","species1_specific_marker_num","species2_specific_marker","species2_specific_marker_num")])

return(result)

})

return(result2)

})


cell_similarity_species$score = signif(cell_similarity_species$score,2)

cell_similarity_species[cell_similarity_species$species2_specific_marker=="","species2_specific_marker"]="-"
cell_similarity_species[cell_similarity_species$species1_specific_marker=="","species1_specific_marker"]="-"
cell_similarity_species[cell_similarity_species$common_marker=="NA","common_marker"]="-"
cell_similarity_species[is.na(cell_similarity_species$common_marker),"common_marker"]="-"

cell_similarity_species$Tissue_standard = "All"
cell_similarity_species$Tissue_UberonOntology_ID2 = "-"
#ready

######################################
########## consider tissue ###########
######################################


all_cell_Gene = unique(cell_taxonomy_3[,c("CT_ID","Specific_Cell_Ontology_ID","Gene_ENTREZID2","Species","Cell_Marker","Tissue_standard","Tissue_UberonOntology_ID2")])
all_cell_Gene = all_cell_Gene[!is.na(all_cell_Gene$Tissue_UberonOntology_ID2),]

all_cell_Gene_HomoloGene = merge(all_cell_Gene,HomoloGene,by.x="Gene_ENTREZID2",by.y="original_gene")
#table(all_cell_Gene_HomoloGene$Species)
#            Danio rerio Drosophila melanogaster           Gallus gallus
#                     76                      36                      18
#           Homo sapiens            Mus musculus       Rattus norvegicus
#                  39202                   21434                      24


all_cell_Gene_HomoloGene2 = all_cell_Gene[!all_cell_Gene$Gene_ENTREZID2 %in% HomoloGene$original_gene,]
all_cell_Gene_HomoloGene2$Gene.ID=NA
all_cell_Gene_HomoloGene2$Gene.Symbol=NA
all_cell_Gene_HomoloGene2$species=NA

all_cell_Gene_HomoloGene_final = rbind(all_cell_Gene_HomoloGene,all_cell_Gene_HomoloGene2)
all_cell_Gene_HomoloGene_final$keyword = paste(all_cell_Gene_HomoloGene_final$Specific_Cell_Ontology_ID, all_cell_Gene_HomoloGene_final$Tissue_UberonOntology_ID2,sep="_")

cell_type_process = unique(all_cell_Gene_HomoloGene_final[,c("Specific_Cell_Ontology_ID","Species","Tissue_UberonOntology_ID2")])

cell_type_process2 = ddply(cell_type_process,"Specific_Cell_Ontology_ID",function(df){

tmp = ddply(df,"Tissue_UberonOntology_ID2",function(df2){

df2$n_species = length(unique(df2$Species))

return(df2)
})

return(tmp)
})

cell_type_process2$keyword = paste(cell_type_process2$Specific_Cell_Ontology_ID, cell_type_process2$Tissue_UberonOntology_ID2,sep="_")

cell_tissue_with_multiple_species = cell_type_process2[cell_type_process2$n_species > 1,"keyword"]

all_cell_Gene_HomoloGene2 = all_cell_Gene_HomoloGene_final[all_cell_Gene_HomoloGene_final$keyword %in% cell_tissue_with_multiple_species, ]



cell_similarity_species_tissue = ddply(all_cell_Gene_HomoloGene2,"keyword",function(df){

#combination with repetition
combination = as.data.frame(expand.grid(rep(list(unique(df$Species)), 2)))

colnames(combination) = c("Species1","Species2")

combination = combination[combination$Species1 != combination$Species2,]

result2 = adply(combination,1,function(species){

species1_name = species[,"Species1"]
species2_name = species[,"Species2"]

result = species1_vs_species2_cell_type(df,species1_name,species2_name)
result = unique(result[,c("CT_ID","Specific_Cell_Ontology_ID","Tissue_standard","Tissue_UberonOntology_ID2","score","common_marker","common_marker_number","species1_specific_marker","species1_specific_marker_num","species2_specific_marker","species2_specific_marker_num")])

return(result)

})

return(result2)

})

cell_similarity_species_tissue = cell_similarity_species_tissue[,-which(names(cell_similarity_species_tissue) %in% "keyword")]

cell_similarity_species_tissue$score = signif(cell_similarity_species_tissue$score,2)

cell_similarity_species_final = rbind(cell_similarity_species,cell_similarity_species_tissue)

cell_similarity_species_final[is.na(cell_similarity_species_final$Tissue_standard),"Tissue_standard"]="-"
cell_similarity_species_final[is.na(cell_similarity_species_final$Tissue_UberonOntology_ID2),"Tissue_UberonOntology_ID2"]="-"

write.table(cell_similarity_species_final,"same_cell_similarity_species.txt",sep="\t",quote=F,col.names = T,row.names = F)

#specific table of common and specific markers
library(stringr)

#common_marker = unique(cell_similarity_species_final[,c("CT_ID","Tissue_UberonOntology_ID2","Species1","Species2","common_marker")])          
common_marker = unique(cell_similarity_species_final[,c("CT_ID","Tissue_standard","Tissue_UberonOntology_ID2","Species1","Species2","common_marker","common_marker_number","score")])          


common_marker_table = adply(common_marker,1,function(df){

string = unlist(strsplit(df$common_marker,split=',', fixed=TRUE))
string = str_replace_all(string,"^ ","")

l = strsplit(string,split='+JS+', fixed=TRUE)

result = data.frame(matrix(unlist(l), ncol=5, byrow=TRUE),stringsAsFactors=FALSE)
colnames(result) = c("Cell_Marker","Gene_ENTREZID2","Ortholog","Ortholog_ID","Ortholog_species")

return(result)
})

common_marker_table = unique(common_marker_table[, -which(names(common_marker_table) %in% c("common_marker"))])

common_marker_table2 = unique(merge(common_marker_table,unique(cell_taxonomy_3[,c("Gene_ENTREZID","Gene_ENTREZID2")]),by.x="Gene_ENTREZID2",by.y="Gene_ENTREZID2"))

common_marker_table2[is.na(common_marker_table2)]="-"

write.table(common_marker_table2,"common_marker_table.txt",sep="\t",quote=F,col.names = T,row.names = F)

species1_specific_marker = unique(cell_similarity_species_final[,c("CT_ID","Tissue_standard","Tissue_UberonOntology_ID2","Species1","Species2","species1_specific_marker")])          

species1_specific_marker_table = adply(species1_specific_marker,1,function(df){

string = unlist(strsplit(df$species1_specific_marker,split=', ', fixed=TRUE))
string = str_replace_all(string,"^ ","")


result = as.data.frame(string)
colnames(result) = c("species1_specific_marker_ID")

return(result)
})

species1_specific_marker_table = unique(species1_specific_marker_table[, -which(names(species1_specific_marker_table) %in% c("species1_specific_marker"))])

species1_specific_marker_table2 = unique(merge(species1_specific_marker_table,unique(cell_taxonomy_3[,c("Cell_Marker","Gene_ENTREZID","Gene_ENTREZID2","Gene_Alias","Uniprot","PFAM")]),by.x="species1_specific_marker_ID",by.y="Gene_ENTREZID2"))

species1_specific_marker_table2[is.na(species1_specific_marker_table2)]="-"

write.table(species1_specific_marker_table2,"species1_specific_marker_table.txt",sep="\t",quote=F,col.names = T,row.names = F)


species2_specific_marker = unique(cell_similarity_species_final[,c("CT_ID","Tissue_standard","Tissue_UberonOntology_ID2","Species1","Species2","species2_specific_marker")])          

species2_specific_marker_table = adply(species2_specific_marker,1,function(df){

string = unlist(strsplit(df$species2_specific_marker,split=', ', fixed=TRUE))
string = str_replace_all(string,"^ ","")


result = as.data.frame(string)
colnames(result) = c("species2_specific_marker_ID")

return(result)
})

species2_specific_marker_table = unique(species2_specific_marker_table[, -which(names(species2_specific_marker_table) %in% c("species2_specific_marker"))])

species2_specific_marker_table2 = unique(merge(species2_specific_marker_table,unique(cell_taxonomy_3[,c("Cell_Marker","Gene_ENTREZID","Gene_ENTREZID2","Gene_Alias","Uniprot","PFAM")]),by.x="species2_specific_marker_ID",by.y="Gene_ENTREZID2"))

species2_specific_marker_table2[is.na(species2_specific_marker_table2)]="-"

write.table(species2_specific_marker_table2,"species2_specific_marker_table.txt",sep="\t",quote=F,col.names = T,row.names = F)


