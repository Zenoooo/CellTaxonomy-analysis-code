#### homologene #######

library(plyr)
library(stringr)
#setwd("D:/R/cell taxonomy/13_cell_type_database/01_data/03_all/")
options(stringsAsFactors = F)

cell_taxonomy_5_23 = read.table("cell_taxonomy_final.txt",sep="\t",header=T,quote="",fill=T)

data = read.table("NCBI_homologene_v68.data",header=T,sep="\t",quote="")

species = unique(cell_taxonomy_5_23[,c("Species","Species_tax_ID")])

colnames(species)=c("species","taxID")
data_species = merge(data,species,by.y="taxID",by.x="Taxonomy.ID")
#table(data_species$species)

#Caenorhabditis elegans             Danio rerio Drosophila melanogaster           Gallus gallus            Homo sapiens 
#                   7575                   20897                    8438                   14600                   19129 
#         Macaca mulatta            Mus musculus            Oryza sativa       Rattus norvegicus 
#                  16843                   21207                   16112                   20616


gene = unique(cell_taxonomy_5_23$Gene_ENTREZID2)
gene_ID_name= unique(cell_taxonomy_5_23[,c("Gene_ENTREZID2","Cell_Marker")])

data_species_gene = data_species[data_species$Gene.ID %in% gene,]

library(plyr)

tmp = adply(data_species_gene,1,function(df){

each_gene=df[,"Gene.ID"]
each_gene_species = data_species_gene[data_species_gene$Gene.ID==each_gene,"species"]

HomoloGene_group = data_species_gene[data_species_gene$Gene.ID==each_gene,"HID.HomoloGene.group.id."]

#get orthologous in different species
HomoloGene = data_species_gene[data_species_gene$"HID.HomoloGene.group.id."==HomoloGene_group & data_species_gene$species != each_gene_species,c("Gene.ID","Gene.Symbol","species")]
HomoloGene = HomoloGene[HomoloGene$Gene.ID != each_gene,]

if(dim(HomoloGene)[1]!=0){
HomoloGene$gene_species = paste(HomoloGene$Gene.Symbol," (",HomoloGene$Gene.ID,", ",HomoloGene$species,")",sep="")

HomoloGene2 = paste(HomoloGene$gene_species,collapse=", ")
HomoloGene_num = length(unique(HomoloGene$gene_species))
}

if(dim(HomoloGene)[1]==0){
HomoloGene2=NA
HomoloGene_num=0
}

df$HomoloGene = HomoloGene2
df$HomoloGene_num = HomoloGene_num

return(df)

})


result = tmp[,c("Gene.ID","HomoloGene","HomoloGene_num")]

result = result[!is.na(result$HomoloGene),]

write.table(result,"HomoloGene.txt",sep="\t",quote=F,col.names = T,row.names = F)

length(unique(result$Gene.ID))
#[1] found Orthologous genes for 15408 cell markers



HomoloGene_table = adply(data_species_gene,1,function(df){

each_gene=df[,"Gene.ID"]
each_gene_species = data_species_gene[data_species_gene$Gene.ID==each_gene,"species"]

HomoloGene_group = data_species_gene[data_species_gene$Gene.ID==each_gene,"HID.HomoloGene.group.id."]

#get orthologous in different species
HomoloGene = data_species_gene[data_species_gene$"HID.HomoloGene.group.id."==HomoloGene_group & data_species_gene$species != each_gene_species,c("Gene.ID","Gene.Symbol","species")]
HomoloGene = HomoloGene[HomoloGene$Gene.ID != each_gene,]

if(dim(HomoloGene)[1]!=0){
HomoloGene$Original_gene = each_gene
}

if(dim(HomoloGene)[1]==0){
HomoloGene = data.frame(Gene.ID=NA,Gene.Symbol=NA,species=NA,Original_gene=each_gene)
}


return(HomoloGene)

})

HomoloGene_table2 = unique(HomoloGene_table[,c("Protein","Protein.accession","Gene.ID","Gene.Symbol","species","Original_gene")])

HomoloGene_table_final = unique(merge(HomoloGene_table2,gene_ID_name,by.x="Original_gene",by.y="Gene_ENTREZID2"))
colnames(HomoloGene_table_final)=c("Gene_ID","Orthologous_protein","Orthologous_protein_ID","Orthologous_gene_ID","Orthologous_gene","Orthologous_species","Cell_Marker")

write.table(HomoloGene_table_final,"HomoloGene_table_final_for_web.txt",sep="\t",quote=F,col.names = T,row.names = F)


tmp2 = adply(data_species_gene,1,function(df){

each_gene=df[,"Gene.ID"]

HomoloGene_group = data_species_gene[data_species_gene$Gene.ID==each_gene,"HID.HomoloGene.group.id."]

Homolo_gene = data_species_gene[data_species_gene$"HID.HomoloGene.group.id."==HomoloGene_group,c("Gene.ID","Gene.Symbol","species")]
Homolo_gene$original_gene = each_gene

Homolo_gene = Homolo_gene[Homolo_gene$Gene.ID!=Homolo_gene$original_gene,]
return(Homolo_gene)

})

result2 = unique(tmp2[,c("Gene.ID","Gene.Symbol","species","original_gene")])

write.table(result2,"HomoloGene_raw.txt",sep="\t",quote=F,col.names = T,row.names = F)
