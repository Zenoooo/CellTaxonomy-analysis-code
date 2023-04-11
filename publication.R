#publication.R
library(plyr)
library(stringr)
#setwd("D:/R/cell taxonomy/13_cell_type_database/01_data/03_all/")
options(stringsAsFactors = F)

cell_taxonomy_cell_tissue_v5 = read.table("cell_taxonomy_final.txt",sep="\t",header=T,quote="",fill=T)

cell_taxonomy_5_2 = cell_taxonomy_cell_tissue_v5

#write.table(as.data.frame(unique(cell_taxonomy_5_2$PMID)),"Publication.txt",sep="\t",quote=F,col.names = T,row.names = F)

publication_infomation = unique(cell_taxonomy_5_2[,c("PMID","Cell_standard","CT_ID","Specific_Cell_Ontology_ID","Species","Cell_Marker","Gene_ENTREZID","Gene_ENTREZID2")])

#resource_num
tmp1 = ddply(publication_infomation,"PMID",function(df){

df$cell_num = length(unique(df[!is.na(df$Specific_Cell_Ontology_ID),"Specific_Cell_Ontology_ID"]))
df$marker_num = length(unique(df[!is.na(df$Gene_ENTREZID2),"Gene_ENTREZID2"]))
return(df)

})

resource_num = unique(tmp1[,c("PMID","cell_num","marker_num")])
resource_num = resource_num[!is.na(resource_num$PMID),]

#### pre
#publication_info = read.table("publication_info.txt",sep="\t",header=T,quote="",fill=T)
#publication_info_add  = read.table("publication_info_add.txt",sep="\t",header=F,quote="",fill=T)
#colnames(publication_info_add) = colnames(publication_info)

#publication_info2 = rbind(publication_info,publication_info_add)
#

publication_info = read.table("PUBLICATION_20220420.txt",sep="\t",header=F,quote="",fill=T)
colnames(publication_info) = c("PMID","title","author","journalvolume","issue","pageInfo","Null","pubYear","pubDate", "type","citations","abstract")

publication_info1 = read.table("pub-ct-1.txt",sep="\t",header=F,quote="",fill=T)
colnames(publication_info1) = c("PMID","doi","title","author","journalvolume","pubYear","citations","abstract","Null")


publication_info2 = rbind(publication_info[,c("PMID","title","journalvolume","pubYear","citations")],publication_info1[,c("PMID","title","journalvolume","pubYear","citations")])

publication_info2 = publication_info2[publication_info2$pubYear !=0,]

publication_info2$title = str_replace_all(publication_info2$title,".$","")

publication_info3 = adply(publication_info2,1,function(df){

df$journalvolume = unlist(strsplit(df$journalvolume, split = "\\("))[1]
df$journalvolume = unlist(strsplit(df$journalvolume, split = ":"))[1]
df$journalvolume = unlist(strsplit(df$journalvolume, split = "\\."))[1]
df$journalvolume = unlist(strsplit(df$journalvolume, split = "\\="))[1]
df$journalvolume = str_replace_all(df$journalvolume," $","")
return(df)

})

publication_info2 = publication_info3

#add 32214235
#add_publication = data.frame(PMID=32214235,title="Construction of a human cell landscape at single-cell level",journalvolume="Nature",pubYear=2020,citations=92)

#publication_info3 = rbind(publication_info2,add_publication)

result = merge(resource_num,publication_info2,by.x="PMID",by.y="PMID")


write.table(result,"Publication_for_browse.txt",sep="\t",quote=F,col.names = T,row.names = F)

publication_infomation = publication_infomation[!is.na(publication_infomation$PMID),]

write.table(publication_infomation,"Publication_detail_for_browse.txt",sep="\t",quote=F,col.names = T,row.names = F)


#write.table(as.data.frame(unique(resource_num[!resource_num$PMID %in% publication_info2$PMID,"PMID"])),"no_PMID.txt",sep="\t",quote=F,col.names = T,row.names = F)


