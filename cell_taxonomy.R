library(plyr)
library(stringr)
#setwd("D:/R/cell taxonomy/13_cell_type_database/01_data/03_all/")
options(stringsAsFactors = F)

cell_taxonomy_cell_tissue_v5 = read.table("cell_taxonomy_final.txt",sep="\t",header=T,quote="",fill=T)

cell_taxonomy_3 = cell_taxonomy_cell_tissue_v5

tissue_ID = unique(cell_taxonomy_3[,c("Tissue_standard","Tissue_UberonOntology_ID2")])
gene_ID = unique(cell_taxonomy_3[,c("Cell_Marker","Gene_ENTREZID2")])
cell_ID = unique(cell_taxonomy_3[,c("CT_ID","Specific_Cell_Ontology_ID","Cell_standard")])

########### paper supported number #######

cell_taxonomy_3_paper_support = ddply(cell_taxonomy_3,"Specific_Cell_Ontology_ID",function(df){

result_tmp = ddply(df,"Gene_ENTREZID2",function(df2){
df2$supported_paper_number = length(unique(df2[!is.na(df2$PMID),"PMID"]))
return(df2)
})

return(result_tmp)
})



paper_support = unique(cell_taxonomy_3_paper_support[,c("CT_ID","Specific_Cell_Ontology_ID","Cell_standard","Gene_ENTREZID2","Cell_Marker","supported_paper_number")])
#64091

paper_support$Species="All"


cell_taxonomy_3_paper_support_species = ddply(cell_taxonomy_3,"Specific_Cell_Ontology_ID",function(df){

result_tmp1 = ddply(df,"Species",function(df2){
	result_tmp2 = ddply(df2,"Gene_ENTREZID2",function(df3){
	df3$supported_paper_number = length(unique(df3[!is.na(df3$PMID),"PMID"]))
	return(df3)
	})
	return(result_tmp2)
})

return(result_tmp1)
})

paper_support2 = unique(cell_taxonomy_3_paper_support_species[,c("CT_ID","Specific_Cell_Ontology_ID","Cell_standard","Gene_ENTREZID2","Cell_Marker","supported_paper_number","Species")])

paper_support3 = rbind(paper_support,paper_support2)

paper_support3[paper_support3 == ""] <- "-"
paper_support3[paper_support3 == "-"] <- "-"
paper_support3[is.na(paper_support3)] <- "-"

write.table(paper_support3,"cell_taxonomy_paper_support.txt",sep="\t",quote=F,col.names = T,row.names = F)
#overall marker evaluation for each cell type
#result used in database

#Specific_Cell_Ontology_ID Gene_ENTREZID2 supported_paper_number



cell_taxonomy_3_paper_support_tissue = ddply(cell_taxonomy_3,"Specific_Cell_Ontology_ID",function(df){

result_tmp1 = ddply(df,"Tissue_UberonOntology_ID2",function(df2){
	result_tmp2 = ddply(df2,"Gene_ENTREZID2",function(df3){
	df3$supported_paper_number = length(unique(df3[!is.na(df3$PMID),"PMID"]))
	return(df3)
	})
	return(result_tmp2)
})

return(result_tmp1)
})


paper_support_tissue = unique(cell_taxonomy_3_paper_support_tissue[,c("CT_ID","Specific_Cell_Ontology_ID","Cell_standard","Tissue_UberonOntology_ID2","Tissue_standard","Gene_ENTREZID2","Cell_Marker","supported_paper_number")])
#73207

paper_support_tissue = paper_support_tissue[!is.na(paper_support_tissue$Tissue_UberonOntology_ID2),]
#66335

paper_support_tissue$Species="All"


cell_taxonomy_3_paper_support_tissue_species = ddply(cell_taxonomy_3,"Specific_Cell_Ontology_ID",function(df){

result_tmp0 = ddply(df,"Species",function(dff){

	result_tmp1 = ddply(dff,"Tissue_UberonOntology_ID2",function(df2){
		result_tmp2 = ddply(df2,"Gene_ENTREZID2",function(df3){
		df3$supported_paper_number = length(unique(df3[!is.na(df3$PMID),"PMID"]))
		return(df3)
		})
		return(result_tmp2)
	})
	return(result_tmp1)
	})
return(result_tmp0)

})

paper_support_tissue2 = unique(cell_taxonomy_3_paper_support_tissue_species[,c("CT_ID","Specific_Cell_Ontology_ID","Cell_standard","Tissue_UberonOntology_ID2","Tissue_standard","Gene_ENTREZID2","Cell_Marker","supported_paper_number","Species")])
# 

paper_support_tissue2 = paper_support_tissue2[!is.na(paper_support_tissue2$Tissue_UberonOntology_ID2),]
# 


paper_support_tissue_3 = rbind(paper_support_tissue,paper_support_tissue2)
#138677

paper_support_tissue_3[paper_support_tissue_3 == ""] <- "-"
paper_support_tissue_3[paper_support_tissue_3 == "-"] <- "-"
paper_support_tissue_3[is.na(paper_support_tissue_3)] <- "-"

write.table(paper_support_tissue_3,"cell_taxonomy_paper_support_tissue.txt",sep="\t",quote=F,col.names = T,row.names = F)

cell_taxonomy_tissue_paper_support = ddply(cell_taxonomy_3,"Tissue_UberonOntology_ID2",function(df){

result_tmp = ddply(df,"Specific_Cell_Ontology_ID",function(df2){
df2$supported_paper_number = length(unique(df2[!is.na(df2$PMID),"PMID"]))
return(df2)
})

return(result_tmp)
})


tissue_paper_support = unique(cell_taxonomy_tissue_paper_support[,c("CT_ID","Specific_Cell_Ontology_ID","Cell_standard","Tissue_standard","Tissue_UberonOntology_ID","Tissue_UberonOntology_ID2","supported_paper_number")])
#5865

tissue_paper_support$Species="All"


cell_taxonomy_3_tissue_paper_support_species = ddply(cell_taxonomy_3,"Tissue_UberonOntology_ID2",function(df){

result_tmp1 = ddply(df,"Species",function(df2){
	result_tmp2 = ddply(df2,"Specific_Cell_Ontology_ID",function(df3){
	df3$supported_paper_number = length(unique(df3[!is.na(df3$PMID),"PMID"]))
	return(df3)
	})
	return(result_tmp2)
})

return(result_tmp1)
})

tissue_paper_support2 = unique(cell_taxonomy_3_tissue_paper_support_species[,c("CT_ID","Specific_Cell_Ontology_ID","Cell_standard","Tissue_standard","Tissue_UberonOntology_ID","Tissue_UberonOntology_ID2","supported_paper_number","Species")])

tissue_paper_support3 = rbind(tissue_paper_support,tissue_paper_support2)

tissue_paper_support3[tissue_paper_support3 == ""] <- "-"
tissue_paper_support3[tissue_paper_support3 == "-"] <- "-"
tissue_paper_support3[is.na(tissue_paper_support3)] <- "-"

write.table(tissue_paper_support3,"cell_taxonomy_tissue_paper_support.txt",sep="\t",quote=F,col.names = T,row.names = F)


############ example ##############
#				tissue    other-tissue    
#marker           7       2             sum = marker_all_tissue
#other-marker     0       5
#         tissue_all_marker_number     all_sum

#dat <- data.frame(
#  tissue_yes = c(7, 0),
#  tissue_no = c(2, 5),
#  row.names = c("marker", "other-marker"),
#  stringsAsFactors = FALSE
#)
#colnames(dat) <- c("tissue", "other-tissue")

#fisher.test(dat,alternative ="greater")
#######################################


marker_tissue_data = unique(cell_taxonomy_3[,c("Specific_Cell_Ontology_ID","Gene_ENTREZID2","Tissue_UberonOntology_ID2","PMID")])

marker_tissue_data = marker_tissue_data[!is.na(marker_tissue_data$Tissue_UberonOntology_ID2) & marker_tissue_data$Tissue_UberonOntology_ID2!="" & !is.na(marker_tissue_data$PMID) & marker_tissue_data$PMID!="" & !is.na(marker_tissue_data$Gene_ENTREZID2),]
#78723

marker_tissue_num = ddply(marker_tissue_data,"Tissue_UberonOntology_ID2",function(df){

df2 = df

tmp1 = ddply(df2,"Specific_Cell_Ontology_ID",function(df3){

cell = unique(df3$Specific_Cell_Ontology_ID)

tmp0 = ddply(df3,"Gene_ENTREZID2",function(df4){

gene = unique(df4$Gene_ENTREZID2)

df4$marker_cell_number = length(unique(df4$PMID))
df4$other_marker_cell_number = length(unique(df2[df2$Specific_Cell_Ontology_ID==cell & df2$Gene_ENTREZID2 != gene,"PMID"]))
df4$marker_other_cell_number = length(unique(df2[df2$Specific_Cell_Ontology_ID!=cell & df2$Gene_ENTREZID2 == gene,"PMID"]))
df4$other_marker_other_cell = length(unique(df2[df2$Specific_Cell_Ontology_ID!=cell & df2$Gene_ENTREZID2 != gene,"PMID"]))

return(df4)

})
return(tmp0)
})

return(tmp1)
})

marker_tissue_num2 = unique(marker_tissue_num[,c("Specific_Cell_Ontology_ID","Tissue_UberonOntology_ID2","Gene_ENTREZID2","marker_cell_number","other_marker_cell_number","marker_other_cell_number","other_marker_other_cell")])

data_for_fisher_test2 = marker_tissue_num2

fisher_test_result2 = adply(data_for_fisher_test2,1,function(df){

dat <- data.frame(
  cell_yes = c(df$marker_cell_number, df$other_marker_cell_number),
  cell_no = c(df$marker_other_cell_number, df$other_marker_other_cell),
  row.names = c("marker", "other-marker"),
  stringsAsFactors = FALSE
)
colnames(dat) <- c("cell", "other-cell")

result = fisher.test(dat,alternative ="greater")
df$P_value = result$p.value
df$OR = result$estimate
return(df)

})

fisher_test_result2$FDR = signif(p.adjust(fisher_test_result2$P_value, "BH"), 3)

#dim(fisher_test_result2[fisher_test_result2$FDR < 0.05 & ((fisher_test_result2$OR > 0) | (fisher_test_result2$OR=="Inf")),])
#270  10

fisher_test_result2$P_value = signif(fisher_test_result2$P_value, 3)
#dim(fisher_test_result2[fisher_test_result2$P_value < 0.05 & ((fisher_test_result2$OR > 0) | (fisher_test_result2$OR=="Inf")),])
#22952   10


fisher_test_result3 = merge(fisher_test_result2,tissue_ID,by.x="Tissue_UberonOntology_ID2",by.y="Tissue_UberonOntology_ID2")
fisher_test_result4 = merge(fisher_test_result3,cell_ID,by.x="Specific_Cell_Ontology_ID",by.y="Specific_Cell_Ontology_ID")

fisher_test_result4$Species="All"

fisher_test_result_all_species = fisher_test_result4


############ 考虑species

marker_tissue_data_species = unique(cell_taxonomy_3[,c("Specific_Cell_Ontology_ID","Gene_ENTREZID2","Tissue_UberonOntology_ID2","PMID","Species")])

marker_tissue_data_species = marker_tissue_data_species[!is.na(marker_tissue_data_species$Tissue_UberonOntology_ID2) & marker_tissue_data_species$Tissue_UberonOntology_ID2!="" & !is.na(marker_tissue_data_species$PMID) & marker_tissue_data_species$PMID!="" & !is.na(marker_tissue_data_species$Gene_ENTREZID2) & !is.na(marker_tissue_data_species$Species),]
#70069

marker_tissue_num = ddply(marker_tissue_data_species,"Species",function(df){

tmp2 = ddply(df,"Tissue_UberonOntology_ID2",function(df2){

tmp1 = ddply(df2,"Specific_Cell_Ontology_ID",function(df3){

cell = unique(df3$Specific_Cell_Ontology_ID)

tmp0 = ddply(df3,"Gene_ENTREZID2",function(df4){

gene = unique(df4$Gene_ENTREZID2)

df4$marker_cell_number = length(unique(df4$PMID))
df4$other_marker_cell_number = length(unique(df2[df2$Specific_Cell_Ontology_ID==cell & df2$Gene_ENTREZID2 != gene,"PMID"]))
df4$marker_other_cell_number = length(unique(df2[df2$Specific_Cell_Ontology_ID!=cell & df2$Gene_ENTREZID2 == gene,"PMID"]))
df4$other_marker_other_cell = length(unique(df2[df2$Specific_Cell_Ontology_ID!=cell & df2$Gene_ENTREZID2 != gene,"PMID"]))

return(df4)

})
return(tmp0)
})

return(tmp1)

})
return(tmp2)
})

marker_tissue_num2 = unique(marker_tissue_num[,c("Species","Specific_Cell_Ontology_ID","Tissue_UberonOntology_ID2","Gene_ENTREZID2","marker_cell_number","other_marker_cell_number","marker_other_cell_number","other_marker_other_cell")])

data_for_fisher_test2 = marker_tissue_num2

fisher_test_result2 = adply(data_for_fisher_test2,1,function(df){

dat <- data.frame(
  cell_yes = c(df$marker_cell_number, df$other_marker_cell_number),
  cell_no = c(df$marker_other_cell_number, df$other_marker_other_cell),
  row.names = c("marker", "other-marker"),
  stringsAsFactors = FALSE
)
colnames(dat) <- c("cell", "other-cell")

result = fisher.test(dat,alternative ="greater")
df$P_value = result$p.value
df$OR = result$estimate
return(df)

})

fisher_test_result2$FDR = signif(p.adjust(fisher_test_result2$P_value, "BH"), 3)

dim(fisher_test_result2[fisher_test_result2$FDR < 0.05 & ((fisher_test_result2$OR > 0) | (fisher_test_result2$OR=="Inf")),])
#196 11

fisher_test_result2$P_value = signif(fisher_test_result2$P_value, 3)
dim(fisher_test_result2[fisher_test_result2$P_value < 0.05 & ((fisher_test_result2$OR > 0) | (fisher_test_result2$OR=="Inf")),])
#18639   11

fisher_test_result3 = merge(fisher_test_result2,tissue_ID,by.x="Tissue_UberonOntology_ID2",by.y="Tissue_UberonOntology_ID2")
fisher_test_result4 = merge(fisher_test_result3,cell_ID,by.x="Specific_Cell_Ontology_ID",by.y="Specific_Cell_Ontology_ID")

fisher_test_result_individual_species = fisher_test_result4

fisher_test_result5 = rbind(fisher_test_result_individual_species,fisher_test_result_all_species)

max_OR = max(fisher_test_result5[fisher_test_result5$OR!="Inf","OR"])

fisher_test_result5$OR2 = fisher_test_result5$OR

fisher_test_result5[fisher_test_result5$OR2=="Inf","OR2"]=max_OR*1.1

fisher_test_result5$OR2 = log(fisher_test_result5$OR2+1)

fisher_test_result5[fisher_test_result5 == ""] <- "-"
fisher_test_result5[fisher_test_result5 == "-"] <- "-"
fisher_test_result5[is.na(fisher_test_result5)] <- "-"

fisher_test_result5$OR = signif(fisher_test_result5$OR,2)
fisher_test_result5$FDR = signif(fisher_test_result5$FDR,2)
fisher_test_result5$P_value = signif(fisher_test_result5$P_value,2)

write.table(fisher_test_result5,"tissue_marker_fisher_test_result_for_each_cell.txt",sep="\t",quote=F,col.names = T,row.names = F)

marker_cell_data = unique(cell_taxonomy_3[,c("Gene_ENTREZID2","Specific_Cell_Ontology_ID","PMID")])
#72363

marker_cell_data = marker_cell_data[ !is.na(marker_cell_data$PMID) & marker_cell_data$PMID!="" & !is.na(marker_cell_data$Gene_ENTREZID2),]
#78723

result = ddply(marker_cell_data,"Specific_Cell_Ontology_ID",function(df){

cell = unique(df$Specific_Cell_Ontology_ID)

tmp0 = ddply(df,"Gene_ENTREZID2",function(df4){

gene = unique(df4$Gene_ENTREZID2)

df4$marker_cell_number = length(unique(df4$PMID))
df4$other_marker_cell_number = length(unique(marker_cell_data[marker_cell_data$Specific_Cell_Ontology_ID==cell & marker_cell_data$Gene_ENTREZID2 != gene,"PMID"]))
df4$marker_other_cell_number = length(unique(marker_cell_data[marker_cell_data$Specific_Cell_Ontology_ID!=cell & marker_cell_data$Gene_ENTREZID2 == gene,"PMID"]))
df4$other_marker_other_cell = length(unique(marker_cell_data[marker_cell_data$Specific_Cell_Ontology_ID!=cell & marker_cell_data$Gene_ENTREZID2 != gene,"PMID"]))

return(df4)

})
return(tmp0)
})

data_for_fisher_test2 = unique(result[,c("Specific_Cell_Ontology_ID","Gene_ENTREZID2","marker_cell_number","other_marker_cell_number","marker_other_cell_number","other_marker_other_cell")])

fisher_test_result2 = adply(data_for_fisher_test2,1,function(df){

dat <- data.frame(
  cell_yes = c(df$marker_cell_number, df$other_marker_cell_number),
  cell_no = c(df$marker_other_cell_number, df$other_marker_other_cell),
  row.names = c("marker", "other-marker"),
  stringsAsFactors = FALSE
)
colnames(dat) <- c("cell", "other-cell")

result = fisher.test(dat,alternative ="greater")
df$P_value = result$p.value
df$OR = result$estimate
return(df)

})

fisher_test_result2$FDR = signif(p.adjust(fisher_test_result2$P_value, "BH"), 3)

#dim(fisher_test_result2[fisher_test_result2$FDR < 0.05 & ((fisher_test_result2$OR > 0) | (fisher_test_result2$OR=="Inf")),])
#270  10

fisher_test_result2$P_value = signif(fisher_test_result2$P_value, 3)
#dim(fisher_test_result2[fisher_test_result2$P_value < 0.05 & ((fisher_test_result2$OR > 0) | (fisher_test_result2$OR=="Inf")),])
#22952   10


fisher_test_result_no_tissue_all_species = merge(fisher_test_result2,cell_ID,by.x="Specific_Cell_Ontology_ID",by.y="Specific_Cell_Ontology_ID")

fisher_test_result_no_tissue_all_species$Species="All"


marker_cell_data_species = unique(cell_taxonomy_3[,c("Gene_ENTREZID2","Specific_Cell_Ontology_ID","PMID","Species")])
#72363

marker_cell_data_species = marker_cell_data_species[ !is.na(marker_cell_data_species$PMID) & marker_cell_data_species$PMID!="" & !is.na(marker_cell_data_species$Gene_ENTREZID2) & !is.na(marker_cell_data_species$Species),]
#78723

result = ddply(marker_cell_data_species,"Species",function(df){

tmp1 = ddply(df,"Specific_Cell_Ontology_ID",function(df2){

cell = unique(df2$Specific_Cell_Ontology_ID)

tmp0 = ddply(df2,"Gene_ENTREZID2",function(df4){

gene = unique(df4$Gene_ENTREZID2)

df4$marker_cell_number = length(unique(df4$PMID))
df4$other_marker_cell_number = length(unique(df[df$Specific_Cell_Ontology_ID==cell & df$Gene_ENTREZID2 != gene,"PMID"]))
df4$marker_other_cell_number = length(unique(df[df$Specific_Cell_Ontology_ID!=cell & df$Gene_ENTREZID2 == gene,"PMID"]))
df4$other_marker_other_cell = length(unique(df[df$Specific_Cell_Ontology_ID!=cell & df$Gene_ENTREZID2 != gene,"PMID"]))

return(df4)

})
return(tmp0)
})

return(tmp1)
})


data_for_fisher_test2 = unique(result[,c("Species","Specific_Cell_Ontology_ID","Gene_ENTREZID2","marker_cell_number","other_marker_cell_number","marker_other_cell_number","other_marker_other_cell")])

fisher_test_result2 = adply(data_for_fisher_test2,1,function(df){

dat <- data.frame(
  cell_yes = c(df$marker_cell_number, df$other_marker_cell_number),
  cell_no = c(df$marker_other_cell_number, df$other_marker_other_cell),
  row.names = c("marker", "other-marker"),
  stringsAsFactors = FALSE
)
colnames(dat) <- c("cell", "other-cell")

result = fisher.test(dat,alternative ="greater")
df$P_value = result$p.value
df$OR = result$estimate
return(df)

})

fisher_test_result2$FDR = signif(p.adjust(fisher_test_result2$P_value, "BH"), 3)

#dim(fisher_test_result2[fisher_test_result2$FDR < 0.05 & ((fisher_test_result2$OR > 0) | (fisher_test_result2$OR=="Inf")),])
#270  10

fisher_test_result2$P_value = signif(fisher_test_result2$P_value, 3)
#dim(fisher_test_result2[fisher_test_result2$P_value < 0.05 & ((fisher_test_result2$OR > 0) | (fisher_test_result2$OR=="Inf")),])
#22952   10


fisher_test_result_no_tissue_individual_species = merge(fisher_test_result2,cell_ID,by.x="Specific_Cell_Ontology_ID",by.y="Specific_Cell_Ontology_ID")


fisher_test_result3 = rbind(fisher_test_result_no_tissue_individual_species,fisher_test_result_no_tissue_all_species)

max_OR = max(fisher_test_result3[fisher_test_result3$OR!="Inf","OR"])

fisher_test_result3$OR2 = fisher_test_result3$OR

fisher_test_result3[fisher_test_result3$OR2=="Inf","OR2"]=max_OR*1.1
fisher_test_result3$OR2 = log(fisher_test_result3$OR2+1)


fisher_test_result3[fisher_test_result3 == ""] <- "-"
fisher_test_result3[fisher_test_result3 == "-"] <- "-"
fisher_test_result3[is.na(fisher_test_result3)] <- "-"

fisher_test_result3$OR = signif(fisher_test_result3$OR,2)
fisher_test_result3$FDR = signif(fisher_test_result3$FDR,2)
fisher_test_result3$P_value = signif(fisher_test_result3$P_value,2)

write.table(fisher_test_result3,"cell_marker_fisher_test_result.txt",sep="\t",quote=F,col.names = T,row.names = F)
#可在marker type整体和cell页面展示
#result used in database



cell_taxonomy_3 = cell_taxonomy_3[cell_taxonomy_3$Additional_Information2!="",]

#############################################
####### cell-type level statistics ##########
#############################################
#Species_num				
tmp1 = ddply(cell_taxonomy_3,"Specific_Cell_Ontology_ID",function(df){

df$Species_num = length(unique(df[!is.na(df$Species) & df$Species!="","Species"]))
return(df)

})

Species_num = unique(tmp1[,c("Specific_Cell_Ontology_ID","Species_num")])
#1262

#Tissue_num
tmp1 = ddply(cell_taxonomy_3,"Specific_Cell_Ontology_ID",function(df){

df$Tissue_num = length(unique(df[!is.na(df$Tissue_UberonOntology_ID2) & df$Tissue_UberonOntology_ID2!="","Tissue_UberonOntology_ID2"]))
return(df)

})

Tissue_num = unique(tmp1[,c("Specific_Cell_Ontology_ID","Tissue_num")])
#1262

#Disease_num

tmp1 = ddply(cell_taxonomy_3,"Specific_Cell_Ontology_ID",function(df){

df$Disease_num = length(unique(df[!is.na(df$Disease_Ontology_ID2) & df$Disease_Ontology_ID2!="","Disease_Ontology_ID2"]))
return(df)

})

Disease_num = unique(tmp1[,c("Specific_Cell_Ontology_ID","Disease_num")])
#1262

#Marker_num
tmp1 = ddply(cell_taxonomy_3,"Specific_Cell_Ontology_ID",function(df){

df$Marker_num = length(unique(df[!is.na(df$Gene_ENTREZID2) & df$Gene_ENTREZID2!="","Gene_ENTREZID2"]))
return(df)

})

Marker_num = unique(tmp1[,c("Specific_Cell_Ontology_ID","Marker_num")])
#1262


#PMID_num
tmp1 = ddply(cell_taxonomy_3,"Specific_Cell_Ontology_ID",function(df){

df$PMID_num = length(unique(df[!is.na(df$PMID) & df$PMID!="","PMID"]))
return(df)

})

PMID_num = unique(tmp1[,c("Specific_Cell_Ontology_ID","PMID_num")])
#1262


#resource_num
tmp1 = ddply(cell_taxonomy_3,"Specific_Cell_Ontology_ID",function(df){

df$resource_num = length(unique(df[!is.na(df$Additional_Information2) & df$Additional_Information2!="","Additional_Information2"]))
return(df)

})

resource_num = unique(tmp1[,c("Specific_Cell_Ontology_ID","resource_num")])
#1262


cell_statistics = Reduce(function(x, y) merge(x, y, all=TRUE), list(Species_num,Tissue_num,Disease_num,Marker_num,PMID_num,resource_num))
cell_statistics2 = merge(cell_statistics,cell_ID,by.x="Specific_Cell_Ontology_ID",by.y="Specific_Cell_Ontology_ID")

Is_new_cell = unique(cell_taxonomy_3[,c("Is_New_Cell_Type","Specific_Cell_Ontology_ID")])

cell_statistics3 = merge(cell_statistics2,Is_new_cell,by.x="Specific_Cell_Ontology_ID",by.y="Specific_Cell_Ontology_ID")

cell_statistics3[cell_statistics3 == ""] <- "-"
cell_statistics3[cell_statistics3 == "-"] <- "-"
cell_statistics3[is.na(cell_statistics3)] <- "-"

write.table(cell_statistics3,"cell_browse_statistics.txt",sep="\t",quote=F,col.names = T,row.names = F)
#used in database cell browser
#

#############################################
####### marker level statistics #############
#############################################

#cell_num
tmp1 = ddply(cell_taxonomy_3,"Gene_ENTREZID2",function(df){

df$cell_num = length(unique(df[!is.na(df$Specific_Cell_Ontology_ID) & df$Specific_Cell_Ontology_ID!="","Specific_Cell_Ontology_ID"]))
return(df)

})

cell_num = unique(tmp1[,c("Gene_ENTREZID2","cell_num")])
#29575

#Tissue_num
tmp1 = ddply(cell_taxonomy_3,"Gene_ENTREZID2",function(df){

df$Tissue_num = length(unique(df[!is.na(df$Tissue_UberonOntology_ID2) & df$Tissue_UberonOntology_ID2!="","Tissue_UberonOntology_ID2"]))
return(df)

})

Tissue_num = unique(tmp1[,c("Gene_ENTREZID2","Tissue_num")])
#29575

#Disease_num

cell_taxonomy_3[is.na(cell_taxonomy_3$Disease_Ontology_ID2) | cell_taxonomy_3$Disease_Ontology_ID2=="","Disease_Ontology_ID2"] = cell_taxonomy_3[is.na(cell_taxonomy_3$Disease_Ontology_ID2) | cell_taxonomy_3$Disease_Ontology_ID2=="","Disease_Type"]

tmp1 = ddply(cell_taxonomy_3,"Gene_ENTREZID2",function(df){

df$Disease_num = length(unique(df[!is.na(df$Disease_Ontology_ID2) & df$Disease_Ontology_ID2!="","Disease_Ontology_ID2"]))
return(df)

})

Disease_num = unique(tmp1[,c("Gene_ENTREZID2","Disease_num")])
#29575



#PMID_num
tmp1 = ddply(cell_taxonomy_3,"Gene_ENTREZID2",function(df){

df$PMID_num = length(unique(df[!is.na(df$PMID) & df$PMID!="","PMID"]))
return(df)

})

PMID_num = unique(tmp1[,c("Gene_ENTREZID2","PMID_num")])
#29575

#cell_num	Tissue_num	Disease_num	PMID_num


#resource_num
tmp1 = ddply(cell_taxonomy_3,"Gene_ENTREZID2",function(df){

df$resource_num = length(unique(df[!is.na(df$Additional_Information2) & df$Additional_Information2!="","Additional_Information2"]))
return(df)

})

resource_num = unique(tmp1[,c("Gene_ENTREZID2","resource_num")])
#1262


gene_ID = unique(cell_taxonomy_3[,c("Gene_ENTREZID2","Gene_ENTREZID","Cell_Marker")])

marker_statistics = Reduce(function(x, y) merge(x, y, all=TRUE), list(cell_num,Tissue_num,Disease_num,PMID_num,resource_num))
marker_statistics2 = merge(marker_statistics,gene_ID,by.x="Gene_ENTREZID2",by.y="Gene_ENTREZID2")

marker_statistics2[marker_statistics2 == ""] <- "-"
marker_statistics2[marker_statistics2 == "-"] <- "-"
marker_statistics2[is.na(marker_statistics2)] <- "-"

HomoloGene = read.table("HomoloGene.txt",sep="\t",header=T,quote="",fill=T)

marker_statistics3 = merge(marker_statistics2,HomoloGene,by.x="Gene_ENTREZID2", by.y="Gene.ID")

marker_statistics4 = marker_statistics2[!marker_statistics2$Gene_ENTREZID2 %in% HomoloGene$Gene.ID,]
marker_statistics4$HomoloGene="-"
marker_statistics4$HomoloGene_num="0"

marker_statistics_final = rbind(marker_statistics3,marker_statistics4)


write.table(marker_statistics_final,"marker_browse_statistics.txt",sep="\t",quote=F,col.names = T,row.names = F)
#used in database marker browser
#24383


#############################################
####### tissue level statistics #############
#############################################


#cell_num
tmp1 = ddply(cell_taxonomy_3,"Tissue_UberonOntology_ID2",function(df){

df$cell_num = length(unique(df[!is.na(df$Specific_Cell_Ontology_ID) & df$Specific_Cell_Ontology_ID!="","Specific_Cell_Ontology_ID"]))
return(df)

})

cell_num = unique(tmp1[,c("Tissue_UberonOntology_ID2","cell_num")])
cell_num = cell_num[cell_num$Tissue_UberonOntology_ID2 !="",]
#325



#Species_num
tmp1 = ddply(cell_taxonomy_3,"Tissue_UberonOntology_ID2",function(df){

df$Species_num = length(unique(df[!is.na(df$Species) & df$Species!="","Species"]))
return(df)

})

Species_num = unique(tmp1[,c("Tissue_UberonOntology_ID2","Species_num")])
Species_num = Species_num[Species_num$Tissue_UberonOntology_ID2 !="",]
#325

#Disease_num

cell_taxonomy_3[is.na(cell_taxonomy_3$Disease_Ontology_ID2) | cell_taxonomy_3$Disease_Ontology_ID2=="","Disease_Ontology_ID2"] = cell_taxonomy_3[is.na(cell_taxonomy_3$Disease_Ontology_ID2) | cell_taxonomy_3$Disease_Ontology_ID2=="","Disease_Type"]

tmp1 = ddply(cell_taxonomy_3,"Tissue_UberonOntology_ID2",function(df){

df$Disease_num = length(unique(df[!is.na(df$Disease_Ontology_ID2) & df$Disease_Ontology_ID2!="","Disease_Ontology_ID2"]))
return(df)

})

Disease_num = unique(tmp1[,c("Tissue_UberonOntology_ID2","Disease_num")])
Disease_num = Disease_num[Disease_num$Tissue_UberonOntology_ID2 !="",]

#Marker_num

tmp1 = ddply(cell_taxonomy_3,"Tissue_UberonOntology_ID2",function(df){

df$Marker_num = length(unique(df[!is.na(df$Gene_ENTREZID2) & df$Gene_ENTREZID2!="","Gene_ENTREZID2"]))
return(df)

})

Marker_num = unique(tmp1[,c("Tissue_UberonOntology_ID2","Marker_num")])
Marker_num = Marker_num[Marker_num$Tissue_UberonOntology_ID2 !="",]



#PMID_num
tmp1 = ddply(cell_taxonomy_3,"Tissue_UberonOntology_ID2",function(df){

df$PMID_num = length(unique(df[!is.na(df$PMID) & df$PMID!="","PMID"]))
return(df)

})

PMID_num = unique(tmp1[,c("Tissue_UberonOntology_ID2","PMID_num")])
PMID_num = PMID_num[PMID_num$Tissue_UberonOntology_ID2 !="",]


#resource_num
tmp1 = ddply(cell_taxonomy_3,"Tissue_UberonOntology_ID2",function(df){

df$resource_num = length(unique(df[!is.na(df$Additional_Information2) & df$Additional_Information2!="","Additional_Information2"]))
return(df)

})

resource_num = unique(tmp1[,c("Tissue_UberonOntology_ID2","resource_num")])
resource_num = resource_num[resource_num$Tissue_UberonOntology_ID2 !="",]
#1262


#Species_num	Disease_num	cell_num	Marker_num	PMID_num

tissue_statistics = Reduce(function(x, y) merge(x, y, all=TRUE), list(Species_num,Disease_num,cell_num,Marker_num,PMID_num,resource_num))

tissue_ID = unique(cell_taxonomy_3[,c("Tissue_standard","Tissue_UberonOntology_ID2","Tissue_UberonOntology_ID")])

tissue_statistics2 = merge(tissue_statistics,tissue_ID,by.x="Tissue_UberonOntology_ID2",by.y="Tissue_UberonOntology_ID2")

tissue_statistics2[tissue_statistics2 == ""] <- "-"
tissue_statistics2[tissue_statistics2 == "-"] <- "-"
tissue_statistics2[is.na(tissue_statistics2)] <- "-"

write.table(tissue_statistics2,"tissue_browse_statistics.txt",sep="\t",quote=F,col.names = T,row.names = F)
#used in database tissue browser


