library("clusterProfiler")
library("enrichplot")
library("goProfiles")
library("topGO")
library("ggplot2")
library("openxlsx") # 読み込み
library("org.Hs.eg.db")
library("KEGGREST")
# ディレクトリの指定
setwd("C:\\R_program\\GO")

# ファイルの読み込み
dataset <- read.xlsx("ECM_dataset.xlsx",) # excel_file

#クラス１の細分化
class1 <- dataset[dataset$class==1,]
class1_collagen <- class1[class1$CATEGORY=="E_collagen",]
class1_glycoprotein  <- class1[class1$CATEGORY=="E_glycoprotein",]
class1_proteoglycans  <- class1[class1$CATEGORY=="E_proteoglycans",]
class1_affiliated <- class1[class1$CATEGORY=="EA_affiliated",]
class1_regulater <- class1[class1$CATEGORY=="EA_regulater",]
class1_serected <- class1[class1$CATEGORY=="EA_serected",]
class1_ID <- class1$KEGGID
class1_collagen_ID <- class1_collagen$KEGGID
class1_glycoprotein_ID <- class1_glycoprotein$KEGGID
class1_proteoglycans_ID <- class1_proteoglycans$KEGGID
class1_affiliated_ID <- class1_affiliated$KEGGID
class1_regulater_ID <- class1_regulater$KEGGID
class1_serected_ID <- class1_serected$KEGGID

#クラス２の細分化
class2 <- dataset[dataset$class==2,]
class2_collagen <- class2[class2$CATEGORY=="E_collagen",]
class2_glycoprotein  <- class2[class2$CATEGORY=="E_glycoprotein",]
class2_proteoglycans  <- class2[class2$CATEGORY=="E_proteoglycans",]
class2_affiliated <- class2[class2$CATEGORY=="EA_affiliated",]
class2_regulater <- class2[class2$CATEGORY=="EA_regulater",]
class2_serected <- class2[class2$CATEGORY=="EA_serected",]
class2_ID <- class2$KEGGID
class2_collagen_ID <- class2_collagen$KEGGID
class2_glycoprotein_ID <- class2_glycoprotein$KEGGID
class2_proteoglycans_ID <- class2_proteoglycans$KEGGID
class2_affiliated_ID <- class2_affiliated$KEGGID
class2_regulater_ID <- class2_regulater$KEGGID
class2_serected_ID <- class2_serected$KEGGID

#クラス３の細分化
class3 <- dataset[dataset$class==3,]
class3_collagen <- class3[class3$CATEGORY=="E_collagen",]
class3_glycoprotein  <- class3[class3$CATEGORY=="E_glycoprotein",]
class3_proteoglycans  <- class3[class3$CATEGORY=="E_proteoglycans",]
class3_affiliated <- class3[class3$CATEGORY=="EA_affiliated",]
class3_regulater <- class3[class3$CATEGORY=="EA_regulater",]
class3_serected <- class3[class3$CATEGORY=="EA_serected",]
class3_ID <- class3$KEGGID
class3_collagen_ID <- class3_collagen$KEGGID
class3_glycoprotein_ID <- class3_glycoprotein$KEGGID
class3_proteoglycans_ID <- class3_proteoglycans$KEGGID
class3_affiliated_ID <- class3_affiliated$KEGGID
class3_regulater_ID <- class3_regulater$KEGGID
class3_serected_ID <- class3_serected$KEGGID

#クラス４の細分化
class4 <- dataset[dataset$class==4,]
class4_collagen <- class4[class4$CATEGORY=="E_collagen",]
class4_glycoprotein  <- class4[class4$CATEGORY=="E_glycoprotein",]
class4_proteoglycans  <- class4[class4$CATEGORY=="E_proteoglycans",]
class4_affiliated <- class4[class4$CATEGORY=="EA_affiliated",]
class4_regulater <- class4[class4$CATEGORY=="EA_regulater",]
class4_serected <- class4[class4$CATEGORY=="EA_serected",]
class4_ID <- class4$KEGGID
class4_collagen_ID <- class4_collagen$KEGGID
class4_glycoprotein_ID <- class4_glycoprotein$KEGGID
class4_proteoglycans_ID <- class4_proteoglycans$KEGGID
class4_affiliated_ID <- class4_affiliated$KEGGID
class4_regulater_ID <- class4_regulater$KEGGID
class4_serected_ID <- class4_serected$KEGGID

#class1
class1_uni <- keggConv("uniprot",class1_ID)
class1_uniID <- c()
for (i in class1_uni){
  j <- substring(i, 4,)
  class1_uniID <- append(class1_uniID,j)
}
go1 <- enrichGO(gene = class1_uniID,
               OrgDb = "org.Hs.eg.db",
               keyType = 'UNIPROT',
               ont = 'all',
               readable = TRUE, # regard as gene names
               pvalueCutoff = 0.01,
               qvalueCutoff  = 0.05,
               )
dotplot(go1, split="ONTOLOGY", font.size = 20, showCategory=5) + facet_grid(ONTOLOGY~., scale="free") #Dot plot all
barplot(go1, showCategory=8,x = "Count")
write.csv(go1, "dot plot all.csv")

egoMF1<- enrichGO(gene = class1_uniID,
                 OrgDb = "org.Hs.eg.db",
                 keyType = 'UNIPROT',
                 ont = 'MF',
                 )
egoBP1<- enrichGO(gene = class1_uniID,
                  OrgDb = "org.Hs.eg.db",
                  keyType = 'UNIPROT',
                  ont = 'BP',
)
egoCC1<- enrichGO(gene = class1_uniID,
                  OrgDb = "org.Hs.eg.db",
                  keyType = 'UNIPROT',
                  ont = 'CC',
)
goplot(egoMF1) #Induced GO DAG graph
goplot(egoBP1) 
goplot(egoCC1) 

ego21 <- simplify(egoBP1)
cnetplot(ego21, foldChange=class1_uniID) #Gene-Concept Network


#class1_collagen
class1_uni1 <- keggConv("uniprot",class1_collagen_ID)
class1_uniID1 <- c()
for (i in class1_uni1){
  j <- substring(i, 4,)
  class1_uniID1 <- append(class1_uniID1,j)
}
go1_1 <- enrichGO(gene = class1_uniID1,
                OrgDb = "org.Hs.eg.db",
                keyType = 'UNIPROT',
                ont = 'all',
                readable = TRUE, # regard as gene names
                pvalueCutoff = 0.01,
                qvalueCutoff  = 0.05,
)
dotplot(go1_1, split="ONTOLOGY", font.size = 20, showCategory=5) + facet_grid(ONTOLOGY~., scale="free") #Dot plot all
barplot(go1_1, showCategory=8,x = "Count")
write.csv(go1_1, "class1_collagen.csv")

egoMF1_1<- enrichGO(gene = class1_uniID1,
                  OrgDb = "org.Hs.eg.db",
                  keyType = 'UNIPROT',
                  ont = 'MF',
)
egoBP1_1<- enrichGO(gene = class1_uniID1,
                  OrgDb = "org.Hs.eg.db",
                  keyType = 'UNIPROT',
                  ont = 'BP',
)
egoCC1_1<- enrichGO(gene = class1_uniID1,
                  OrgDb = "org.Hs.eg.db",
                  keyType = 'UNIPROT',
                  ont = 'CC',
)
goplot(egoMF1_1) #Induced GO DAG graph
goplot(egoBP1_1) 
goplot(egoCC1_1) 

ego21_1 <- simplify(egoBP1_1)
cnetplot(ego21_1, foldChange=class1_uniID1) #Gene-Concept Network


#class1_glycoprotein
class1_uni2 <- keggConv("uniprot",class1_glycoprotein_ID)
class1_uniID2 <- c()
for (i in class1_uni2){
  j <- substring(i, 4,)
  class1_uniID2 <- append(class1_uniID2,j)
}
go1_2 <- enrichGO(gene = class1_uniID2,
                  OrgDb = "org.Hs.eg.db",
                  keyType = 'UNIPROT',
                  ont = 'all',
                  readable = TRUE, # regard as gene names
                  pvalueCutoff = 0.01,
                  qvalueCutoff  = 0.05,
)
dotplot(go1_2, split="ONTOLOGY", font.size = 20, showCategory=5) + facet_grid(ONTOLOGY~., scale="free") #Dot plot all
barplot(go1_2, showCategory=8,x = "Count")
write.csv(go1_2, "class1_glycoprotein.csv")

egoMF1_2<- enrichGO(gene = class1_uniID2,
                    OrgDb = "org.Hs.eg.db",
                    keyType = 'UNIPROT',
                    ont = 'MF',
)
egoBP1_2<- enrichGO(gene = class1_uniID2,
                    OrgDb = "org.Hs.eg.db",
                    keyType = 'UNIPROT',
                    ont = 'BP',
)
egoCC1_2<- enrichGO(gene = class1_uniID2,
                    OrgDb = "org.Hs.eg.db",
                    keyType = 'UNIPROT',
                    ont = 'CC',
)
goplot(egoMF1_2) #Induced GO DAG graph
goplot(egoBP1_2) 
goplot(egoCC1_2) 

ego21_2 <- simplify(egoBP1_2)
cnetplot(ego21_2, foldChange=class1_uniID2) #Gene-Concept Network


#class1_proteoglycans

class1_uni3 <- keggConv("uniprot",class1_proteoglycans_ID)
class1_uniID3 <- c()
for (i in class1_uni3){
  j <- substring(i, 4,)
  class1_uniID3 <- append(class1_uniID3,j)
}
go1_3 <- enrichGO(gene = class1_uniID3,
                  OrgDb = "org.Hs.eg.db",
                  keyType = 'UNIPROT',
                  ont = 'all',
                  readable = TRUE, # regard as gene names
                  pvalueCutoff = 0.01,
                  qvalueCutoff  = 0.05,
)
dotplot(go1_3, split="ONTOLOGY", font.size = 20, showCategory=5) + facet_grid(ONTOLOGY~., scale="free") #Dot plot all
barplot(go1_3, showCategory=8,x = "Count")
write.csv(go1_3, "class1_proteoglycans.csv")

egoMF1_3<- enrichGO(gene = class1_uniID3,
                    OrgDb = "org.Hs.eg.db",
                    keyType = 'UNIPROT',
                    ont = 'MF',
)
egoBP1_3<- enrichGO(gene = class1_uniID3,
                    OrgDb = "org.Hs.eg.db",
                    keyType = 'UNIPROT',
                    ont = 'BP',
)
egoCC1_3<- enrichGO(gene = class1_uniID3,
                    OrgDb = "org.Hs.eg.db",
                    keyType = 'UNIPROT',
                    ont = 'CC',
)
goplot(egoMF1_3) #Induced GO DAG graph
goplot(egoBP1_3) 
goplot(egoCC1_3) 

ego21_3 <- simplify(egoBP1_3)
cnetplot(ego21_3, foldChange=class1_uniID1) #Gene-Concept Network


#class1_affiated
class1_uni4 <- keggConv("uniprot",class1_affiliated_ID)
class1_uniID4 <- c()
for (i in class1_uni4){
  j <- substring(i, 4,)
  class1_uniID4 <- append(class1_uniID4,j)
}
go1_4 <- enrichGO(gene = class1_uniID4,
                  OrgDb = "org.Hs.eg.db",
                  keyType = 'UNIPROT',
                  ont = 'all',
                  readable = TRUE, # regard as gene names
                  pvalueCutoff = 0.01,
                  qvalueCutoff  = 0.05,
)
dotplot(go1_4, split="ONTOLOGY", font.size = 20, showCategory=5) + facet_grid(ONTOLOGY~., scale="free") #Dot plot all
barplot(go1_4, showCategory=8,x = "Count")
write.csv(go1_4, "class1_affiated.csv")

egoMF1_4<- enrichGO(gene = class1_uniID4,
                    OrgDb = "org.Hs.eg.db",
                    keyType = 'UNIPROT',
                    ont = 'MF',
)
egoBP1_4<- enrichGO(gene = class1_uniID4,
                    OrgDb = "org.Hs.eg.db",
                    keyType = 'UNIPROT',
                    ont = 'BP',
)
egoCC1_4<- enrichGO(gene = class1_uniID4,
                    OrgDb = "org.Hs.eg.db",
                    keyType = 'UNIPROT',
                    ont = 'CC',
)
goplot(egoMF1_4) #Induced GO DAG graph
goplot(egoBP1_4) 
goplot(egoCC1_4) 

ego21_4 <- simplify(egoBP1_4)
cnetplot(ego21_4, foldChange=class1_uniID4) #Gene-Concept Network

#class1_regulater
class1_uni5 <- keggConv("uniprot",class1_regulater_ID)
class1_uniID5 <- c()
for (i in class1_uni5){
  j <- substring(i, 4,)
  class1_uniID5 <- append(class1_uniID5,j)
}
go1_5 <- enrichGO(gene = class1_uniID5,
                  OrgDb = "org.Hs.eg.db",
                  keyType = 'UNIPROT',
                  ont = 'all',
                  readable = TRUE, # regard as gene names
                  pvalueCutoff = 0.01,
                  qvalueCutoff  = 0.05,
)
dotplot(go1_5, split="ONTOLOGY", font.size = 20, showCategory=5) + facet_grid(ONTOLOGY~., scale="free") #Dot plot all
barplot(go1_5, showCategory=8,x = "Count")
write.csv(go1_5, "class1_regulater.csv")

egoMF1_5<- enrichGO(gene = class1_uniID5,
                    OrgDb = "org.Hs.eg.db",
                    keyType = 'UNIPROT',
                    ont = 'MF',
)
egoBP1_5<- enrichGO(gene = class1_uniID5,
                    OrgDb = "org.Hs.eg.db",
                    keyType = 'UNIPROT',
                    ont = 'BP',
)
egoCC1_5<- enrichGO(gene = class1_uniID5,
                    OrgDb = "org.Hs.eg.db",
                    keyType = 'UNIPROT',
                    ont = 'CC',
)
goplot(egoMF1_5) #Induced GO DAG graph
goplot(egoBP1_5) 
goplot(egoCC1_5) 

ego21_5 <- simplify(egoBP1_5)
cnetplot(ego21_5, foldChange=class1_uniID5) #Gene-Concept Network


#class1_serected
class1_uni6 <- keggConv("uniprot",class1_serected_ID)
class1_uniID6 <- c()
for (i in class1_uni6){
  j <- substring(i, 4,)
  class1_uniID6 <- append(class1_uniID6,j)
}
go1_6 <- enrichGO(gene = class1_uniID6,
                  OrgDb = "org.Hs.eg.db",
                  keyType = 'UNIPROT',
                  ont = 'all',
                  readable = TRUE, # regard as gene names
                  pvalueCutoff = 0.01,
                  qvalueCutoff  = 0.05,
)
dotplot(go1_6, split="ONTOLOGY", font.size = 20, showCategory=5) + facet_grid(ONTOLOGY~., scale="free") #Dot plot all
barplot(go1_6, showCategory=8,x = "Count")
write.csv(go1_6, "class1_serected.csv")

egoMF1_6<- enrichGO(gene = class1_uniID6,
                    OrgDb = "org.Hs.eg.db",
                    keyType = 'UNIPROT',
                    ont = 'MF',
)
egoBP1_6<- enrichGO(gene = class1_uniID6,
                    OrgDb = "org.Hs.eg.db",
                    keyType = 'UNIPROT',
                    ont = 'BP',
)
egoCC1_6<- enrichGO(gene = class1_uniID6,
                    OrgDb = "org.Hs.eg.db",
                    keyType = 'UNIPROT',
                    ont = 'CC',
)
goplot(egoMF1_6) #Induced GO DAG graph
goplot(egoBP1_6) 
goplot(egoCC1_6) 

ego21_6 <- simplify(egoBP1_6)
cnetplot(ego21_6, foldChange=class1_uniID6) #Gene-Concept Network

#class2
class2_ID <- class2$KEGGID
class2_uni <- keggConv("uniprot",class2$KEGGID)
class2_uniID <- c()
for (i in class2_uni){
  j <- substring(i, 4,)
  class2_uniID <- append(class2_uniID,j)
}

go2 <- enrichGO(gene = class2_uniID,
               OrgDb = "org.Hs.eg.db",
               keyType = 'UNIPROT',
               ont = 'all',
               readable = TRUE, # regard as gene names
               pvalueCutoff = 0.01,
               qvalueCutoff  = 0.05,
               )
dotplot(go2, split="ONTOLOGY", font.size = 20, showCategory=5) + facet_grid(ONTOLOGY~., scale="free") #Dot plot all
barplot(go2, showCategory=8,x = "Count")
write.csv(go2, "go_class2.csv")

egoBP2<- enrichGO(gene = class2_uniID,
                 OrgDb = "org.Hs.eg.db",
                 keyType = 'UNIPROT',
                 ont = 'BP',
                 )
egoMF2<- enrichGO(gene = class2_uniID,
                  OrgDb = "org.Hs.eg.db",
                  keyType = 'UNIPROT',
                  ont = 'MF',
)
egoCC2<- enrichGO(gene = class2_uniID,
                  OrgDb = "org.Hs.eg.db",
                  keyType = 'UNIPROT',
                  ont = 'CC',
)
goplot(egoMF2) #Induced GO DAG graph
goplot(egoBP2)
goplot(egoCC2)
ego22 <- simplify(egoBP2)
cnetplot(ego22, foldChange=class2_uniID) #Gene-Concept Network


#class2_collagen
class2_uni1 <- keggConv("uniprot",class2_collagen_ID)
class2_uniID1 <- c()
for (i in class2_uni1){
  j <- substring(i, 4,)
  class2_uniID1 <- append(class2_uniID1,j)
}
go2_1 <- enrichGO(gene = class2_uniID1,
                  OrgDb = "org.Hs.eg.db",
                  keyType = 'UNIPROT',
                  ont = 'all',
                  readable = TRUE, # regard as gene names
                  pvalueCutoff = 0.01,
                  qvalueCutoff  = 0.05,
)
dotplot(go2_1, split="ONTOLOGY", font.size = 20, showCategory=5) + facet_grid(ONTOLOGY~., scale="free") #Dot plot all
barplot(go2_1, showCategory=8,x = "Count")
write.csv(go2_1, "class2_collagen.csv")

egoMF2_1<- enrichGO(gene = class2_uniID1,
                    OrgDb = "org.Hs.eg.db",
                    keyType = 'UNIPROT',
                    ont = 'MF',
)
egoBP2_1<- enrichGO(gene = class2_uniID1,
                    OrgDb = "org.Hs.eg.db",
                    keyType = 'UNIPROT',
                    ont = 'BP',
)
egoCC2_1<- enrichGO(gene = class2_uniID1,
                    OrgDb = "org.Hs.eg.db",
                    keyType = 'UNIPROT',
                    ont = 'CC',
)
goplot(egoMF2_1) #Induced GO DAG graph
goplot(egoBP2_1) 
goplot(egoCC2_1) 

ego22_1 <- simplify(egoBP2_1)
cnetplot(ego22_1, foldChange=class2_uniID1) #Gene-Concept Network


#class2_glycoprotein
class2_uni2 <- keggConv("uniprot",class2_glycoprotein_ID)
class2_uniID2 <- c()
for (i in class2_uni2){
  j <- substring(i, 4,)
  class2_uniID2 <- append(class2_uniID2,j)
}
go2_2 <- enrichGO(gene = class2_uniID2,
                  OrgDb = "org.Hs.eg.db",
                  keyType = 'UNIPROT',
                  ont = 'all',
                  readable = TRUE, # regard as gene names
                  pvalueCutoff = 0.01,
                  qvalueCutoff  = 0.05,
)
dotplot(go2_2, split="ONTOLOGY", font.size = 20, showCategory=5) + facet_grid(ONTOLOGY~., scale="free") #Dot plot all
barplot(go2_2, showCategory=8,x = "Count")
write.csv(go2_2, "class2_glycoprotein.csv")

egoMF2_2<- enrichGO(gene = class2_uniID2,
                    OrgDb = "org.Hs.eg.db",
                    keyType = 'UNIPROT',
                    ont = 'MF',
)
egoBP2_2<- enrichGO(gene = class2_uniID2,
                    OrgDb = "org.Hs.eg.db",
                    keyType = 'UNIPROT',
                    ont = 'BP',
)
egoCC2_2<- enrichGO(gene = class2_uniID2,
                    OrgDb = "org.Hs.eg.db",
                    keyType = 'UNIPROT',
                    ont = 'CC',
)
goplot(egoMF2_2) #Induced GO DAG graph
goplot(egoBP2_2) 
goplot(egoCC2_2) 

ego22_2 <- simplify(egoBP2_2)
cnetplot(ego22_2, foldChange=class2_uniID2) #Gene-Concept Network


#class2_proteoglycans

class2_uni3 <- keggConv("uniprot",class2_proteoglycans_ID)
class2_uniID3 <- c()
for (i in class2_uni3){
  j <- substring(i, 4,)
  class2_uniID3 <- append(class2_uniID3,j)
}
go2_3 <- enrichGO(gene = class2_uniID3,
                  OrgDb = "org.Hs.eg.db",
                  keyType = 'UNIPROT',
                  ont = 'all',
                  readable = TRUE, # regard as gene names
                  pvalueCutoff = 0.01,
                  qvalueCutoff  = 0.05,
)
dotplot(go2_3, split="ONTOLOGY", font.size = 20, showCategory=5) + facet_grid(ONTOLOGY~., scale="free") #Dot plot all
barplot(go2_3, showCategory=8,x = "Count")
write.csv(go2_3, "class2_proteoglycans.csv")

egoMF2_3<- enrichGO(gene = class2_uniID3,
                    OrgDb = "org.Hs.eg.db",
                    keyType = 'UNIPROT',
                    ont = 'MF',
)
egoBP2_3<- enrichGO(gene = class2_uniID3,
                    OrgDb = "org.Hs.eg.db",
                    keyType = 'UNIPROT',
                    ont = 'BP',
)
egoCC2_3<- enrichGO(gene = class2_uniID3,
                    OrgDb = "org.Hs.eg.db",
                    keyType = 'UNIPROT',
                    ont = 'CC',
)
goplot(egoMF2_3) #Induced GO DAG graph
goplot(egoBP2_3) 
goplot(egoCC2_3) 

ego22_3 <- simplify(egoBP2_3)
cnetplot(ego22_3, foldChange=class2_uniID3) #Gene-Concept Network


#class2_affiated
class2_uni4 <- keggConv("uniprot",class2_affiliated_ID)
class2_uniID4 <- c()
for (i in class2_uni4){
  j <- substring(i, 4,)
  class2_uniID4 <- append(class2_uniID4,j)
}
go2_4 <- enrichGO(gene = class2_uniID4,
                  OrgDb = "org.Hs.eg.db",
                  keyType = 'UNIPROT',
                  ont = 'all',
                  readable = TRUE, # regard as gene names
                  pvalueCutoff = 0.01,
                  qvalueCutoff  = 0.05,
)
dotplot(go2_4, split="ONTOLOGY", font.size = 20, showCategory=5) + facet_grid(ONTOLOGY~., scale="free") #Dot plot all
barplot(go2_4, showCategory=8,x = "Count")
write.csv(go2_4, "class2_affiated.csv")

egoMF2_4<- enrichGO(gene = class2_uniID4,
                    OrgDb = "org.Hs.eg.db",
                    keyType = 'UNIPROT',
                    ont = 'MF',
)
egoBP2_4<- enrichGO(gene = class2_uniID4,
                    OrgDb = "org.Hs.eg.db",
                    keyType = 'UNIPROT',
                    ont = 'BP',
)
egoCC2_4<- enrichGO(gene = class2_uniID4,
                    OrgDb = "org.Hs.eg.db",
                    keyType = 'UNIPROT',
                    ont = 'CC',
)
goplot(egoMF2_4) #Induced GO DAG graph
goplot(egoBP2_4) 
goplot(egoCC2_4) 

ego22_4 <- simplify(egoBP2_4)
cnetplot(ego22_4, foldChange=class2_uniID4) #Gene-Concept Network

#class2_regulater
class2_uni5 <- keggConv("uniprot",class2_regulater_ID)
class2_uniID5 <- c()
for (i in class2_uni5){
  j <- substring(i, 4,)
  class2_uniID5 <- append(class2_uniID5,j)
}
go2_5 <- enrichGO(gene = class2_uniID5,
                  OrgDb = "org.Hs.eg.db",
                  keyType = 'UNIPROT',
                  ont = 'all',
                  readable = TRUE, # regard as gene names
                  pvalueCutoff = 0.01,
                  qvalueCutoff  = 0.05,
)
dotplot(go2_5, split="ONTOLOGY", font.size = 20, showCategory=5) + facet_grid(ONTOLOGY~., scale="free") #Dot plot all
barplot(go2_5, showCategory=8,x = "Count")
write.csv(go2_5, "class2_regulater.csv")

egoMF2_5<- enrichGO(gene = class2_uniID5,
                    OrgDb = "org.Hs.eg.db",
                    keyType = 'UNIPROT',
                    ont = 'MF',
)
egoBP2_5<- enrichGO(gene = class2_uniID5,
                    OrgDb = "org.Hs.eg.db",
                    keyType = 'UNIPROT',
                    ont = 'BP',
)
egoCC2_5<- enrichGO(gene = class2_uniID5,
                    OrgDb = "org.Hs.eg.db",
                    keyType = 'UNIPROT',
                    ont = 'CC',
)
goplot(egoMF2_5) #Induced GO DAG graph
goplot(egoBP2_5) 
goplot(egoCC2_5) 

ego22_5 <- simplify(egoBP2_5)
cnetplot(ego22_5, foldChange=class2_uniID5) #Gene-Concept Network


#class1_serected
class2_uni6 <- keggConv("uniprot",class2_serected_ID)
class2_uniID6 <- c()
for (i in class2_uni6){
  j <- substring(i, 4,)
  class2_uniID6 <- append(class2_uniID6,j)
}
go2_6 <- enrichGO(gene = class2_uniID6,
                  OrgDb = "org.Hs.eg.db",
                  keyType = 'UNIPROT',
                  ont = 'all',
                  readable = TRUE, # regard as gene names
                  pvalueCutoff = 0.01,
                  qvalueCutoff  = 0.05,
)
dotplot(go2_6, split="ONTOLOGY", font.size = 20, showCategory=5) + facet_grid(ONTOLOGY~., scale="free") #Dot plot all
barplot(go2_6, showCategory=8,x = "Count")
write.csv(go2_6, "class2_serected.csv")

egoMF2_6<- enrichGO(gene = class2_uniID6,
                    OrgDb = "org.Hs.eg.db",
                    keyType = 'UNIPROT',
                    ont = 'MF',
)
egoBP2_6<- enrichGO(gene = class2_uniID6,
                    OrgDb = "org.Hs.eg.db",
                    keyType = 'UNIPROT',
                    ont = 'BP',
)
egoCC2_6<- enrichGO(gene = class2_uniID6,
                    OrgDb = "org.Hs.eg.db",
                    keyType = 'UNIPROT',
                    ont = 'CC',
)
goplot(egoMF2_6) #Induced GO DAG graph
goplot(egoBP2_6) 
goplot(egoCC2_6) 

ego22_6 <- simplify(egoBP2_6)
cnetplot(ego22_6, foldChange=class2_uniID6) #Gene-Concept Network


#class3
class3_ID <- class3$KEGGID
class3_uni <- keggConv("uniprot",class3$KEGGID)
class3_uniID <- c()
for (i in class3_uni){
  j <- substring(i, 4,)
  class3_uniID <- append(class3_uniID,j)
}

go3 <- enrichGO(gene = class3_uniID,
               OrgDb = "org.Hs.eg.db",
               keyType = 'UNIPROT',
               ont = 'all',
               readable = TRUE, # regard as gene names
               pvalueCutoff = 0.01,
               qvalueCutoff  = 0.05,
               )
dotplot(go3, split="ONTOLOGY", font.size = 20, showCategory=5) + facet_grid(ONTOLOGY~., scale="free") #Dot plot all
write.csv(go3, "class3.csv")

egoMF3<- enrichGO(gene = class3_uniID,
                 OrgDb = "org.Hs.eg.db",
                 keyType = 'UNIPROT',
                 ont = 'MF',
                 )
egoBP3<- enrichGO(gene = class3_uniID,
                  OrgDb = "org.Hs.eg.db",
                  keyType = 'UNIPROT',
                  ont = 'BP',
)
egoCC3<- enrichGO(gene = class3_uniID,
                  OrgDb = "org.Hs.eg.db",
                  keyType = 'UNIPROT',
                  ont = 'CC',
)
goplot(egoMF3) #Induced GO DAG graph
goplot(egoBP3)
goplot(egoCC3)

egoBP3<- enrichGO(gene = class3_uniID,
                 OrgDb = "org.Hs.eg.db",
                 keyType = 'UNIPROT',
                 ont = 'BP',
                 )
ego23 <- simplify(egoBP3)
cnetplot(ego23, foldChange=class3_uniID) #Gene-Concept Network

#class3_collagen
class3_uni1 <- keggConv("uniprot",class3_collagen_ID)
class3_uniID1 <- c()
for (i in class3_uni1){
  j <- substring(i, 4,)
  class3_uniID1 <- append(class3_uniID1,j)
}
go3_1 <- enrichGO(gene = class3_uniID1,
                  OrgDb = "org.Hs.eg.db",
                  keyType = 'UNIPROT',
                  ont = 'all',
                  readable = TRUE, # regard as gene names
                  pvalueCutoff = 0.01,
                  qvalueCutoff  = 0.05,
)
dotplot(go3_1, split="ONTOLOGY", font.size = 20, showCategory=5) + facet_grid(ONTOLOGY~., scale="free") #Dot plot all
barplot(go3_1, showCategory=8,x = "Count")
write.csv(go3_1, "class3_collagen.csv")

egoMF3_1<- enrichGO(gene = class3_uniID1,
                    OrgDb = "org.Hs.eg.db",
                    keyType = 'UNIPROT',
                    ont = 'MF',
)
egoBP3_1<- enrichGO(gene = class3_uniID1,
                    OrgDb = "org.Hs.eg.db",
                    keyType = 'UNIPROT',
                    ont = 'BP',
)
egoCC3_1<- enrichGO(gene = class3_uniID1,
                    OrgDb = "org.Hs.eg.db",
                    keyType = 'UNIPROT',
                    ont = 'CC',
)
goplot(egoMF3_1) #Induced GO DAG graph
goplot(egoBP3_1) 
goplot(egoCC3_1) 

ego23_1 <- simplify(egoBP3_1)
cnetplot(ego23_1, foldChange=class3_uniID1) #Gene-Concept Network


#class3_glycoprotein
class3_uni2 <- keggConv("uniprot",class3_glycoprotein_ID)
class3_uniID2 <- c()
for (i in class3_uni2){
  j <- substring(i, 4,)
  class3_uniID2 <- append(class3_uniID2,j)
}
go3_2 <- enrichGO(gene = class3_uniID2,
                  OrgDb = "org.Hs.eg.db",
                  keyType = 'UNIPROT',
                  ont = 'all',
                  readable = TRUE, # regard as gene names
                  pvalueCutoff = 0.01,
                  qvalueCutoff  = 0.05,
)
dotplot(go3_2, split="ONTOLOGY", font.size = 20, showCategory=5) + facet_grid(ONTOLOGY~., scale="free") #Dot plot all
barplot(go3_2, showCategory=8,x = "Count")
write.csv(go3_2, "class3_glycoprotein.csv")

egoMF3_2<- enrichGO(gene = class3_uniID2,
                    OrgDb = "org.Hs.eg.db",
                    keyType = 'UNIPROT',
                    ont = 'MF',
)
egoBP3_2<- enrichGO(gene = class3_uniID2,
                    OrgDb = "org.Hs.eg.db",
                    keyType = 'UNIPROT',
                    ont = 'BP',
)
egoCC3_2<- enrichGO(gene = class3_uniID2,
                    OrgDb = "org.Hs.eg.db",
                    keyType = 'UNIPROT',
                    ont = 'CC',
)
goplot(egoMF3_2) #Induced GO DAG graph
goplot(egoBP3_2) 
goplot(egoCC3_2) 

ego23_2 <- simplify(egoBP3_2)
cnetplot(ego23_2, foldChange=class3_uniID2) #Gene-Concept Network


#class3_proteoglycans

class3_uni3 <- keggConv("uniprot",class3_proteoglycans_ID)
class3_uniID3 <- c()
for (i in class3_uni3){
  j <- substring(i, 4,)
  class3_uniID3 <- append(class3_uniID3,j)
}
go3_3 <- enrichGO(gene = class3_uniID3,
                  OrgDb = "org.Hs.eg.db",
                  keyType = 'UNIPROT',
                  ont = 'all',
                  readable = TRUE, # regard as gene names
                  pvalueCutoff = 0.01,
                  qvalueCutoff  = 0.05,
)
dotplot(go3_3, split="ONTOLOGY", font.size = 20, showCategory=5) + facet_grid(ONTOLOGY~., scale="free") #Dot plot all
barplot(go3_3, showCategory=8,x = "Count")
write.csv(go3_3, "class3_proteoglycans.csv")

egoMF3_3<- enrichGO(gene = class3_uniID3,
                    OrgDb = "org.Hs.eg.db",
                    keyType = 'UNIPROT',
                    ont = 'MF',
)
egoBP3_3<- enrichGO(gene = class3_uniID3,
                    OrgDb = "org.Hs.eg.db",
                    keyType = 'UNIPROT',
                    ont = 'BP',
)
egoCC3_3<- enrichGO(gene = class3_uniID3,
                    OrgDb = "org.Hs.eg.db",
                    keyType = 'UNIPROT',
                    ont = 'CC',
)
goplot(egoMF3_3) #Induced GO DAG graph
goplot(egoBP3_3) 
goplot(egoCC3_3) 

ego23_3 <- simplify(egoBP3_3)
cnetplot(ego23_3, foldChange=class3_uniID3) #Gene-Concept Network


#class3_affiated
class3_uni4 <- keggConv("uniprot",class3_affiliated_ID)
class3_uniID4 <- c()
for (i in class3_uni4){
  j <- substring(i, 4,)
  class3_uniID4 <- append(class3_uniID4,j)
}
go3_4 <- enrichGO(gene = class3_uniID4,
                  OrgDb = "org.Hs.eg.db",
                  keyType = 'UNIPROT',
                  ont = 'all',
                  readable = TRUE, # regard as gene names
                  pvalueCutoff = 0.01,
                  qvalueCutoff  = 0.05,
)
dotplot(go3_4, split="ONTOLOGY", font.size = 20, showCategory=5) + facet_grid(ONTOLOGY~., scale="free") #Dot plot all
barplot(go3_4, showCategory=8,x = "Count")
write.csv(go3_4, "class3_affiated.csv")

egoMF3_4<- enrichGO(gene = class3_uniID4,
                    OrgDb = "org.Hs.eg.db",
                    keyType = 'UNIPROT',
                    ont = 'MF',
)
egoBP3_4<- enrichGO(gene = class3_uniID4,
                    OrgDb = "org.Hs.eg.db",
                    keyType = 'UNIPROT',
                    ont = 'BP',
)
egoCC3_4<- enrichGO(gene = class3_uniID4,
                    OrgDb = "org.Hs.eg.db",
                    keyType = 'UNIPROT',
                    ont = 'CC',
)
goplot(egoMF3_4) #Induced GO DAG graph
goplot(egoBP3_4) 
goplot(egoCC3_4) 

ego23_4 <- simplify(egoBP3_4)
cnetplot(ego23_4, foldChange=class3_uniID4) #Gene-Concept Network

#class3_regulater
class3_uni5 <- keggConv("uniprot",class3_regulater_ID)
class3_uniID5 <- c()
for (i in class3_uni5){
  j <- substring(i, 4,)
  class3_uniID5 <- append(class3_uniID5,j)
}
go3_5 <- enrichGO(gene = class3_uniID5,
                  OrgDb = "org.Hs.eg.db",
                  keyType = 'UNIPROT',
                  ont = 'all',
                  readable = TRUE, # regard as gene names
                  pvalueCutoff = 0.01,
                  qvalueCutoff  = 0.05,
)
dotplot(go3_5, split="ONTOLOGY", font.size = 20, showCategory=5) + facet_grid(ONTOLOGY~., scale="free") #Dot plot all
barplot(go3_5, showCategory=8,x = "Count")
write.csv(go3_5, "class3_regulater.csv")

egoMF3_5<- enrichGO(gene = class3_uniID5,
                    OrgDb = "org.Hs.eg.db",
                    keyType = 'UNIPROT',
                    ont = 'MF',
)
egoBP3_5<- enrichGO(gene = class3_uniID5,
                    OrgDb = "org.Hs.eg.db",
                    keyType = 'UNIPROT',
                    ont = 'BP',
)
egoCC3_5<- enrichGO(gene = class3_uniID5,
                    OrgDb = "org.Hs.eg.db",
                    keyType = 'UNIPROT',
                    ont = 'CC',
)
goplot(egoMF3_5) #Induced GO DAG graph
goplot(egoBP3_5) 
goplot(egoCC3_5) 

ego23_5 <- simplify(egoBP3_5)
cnetplot(ego23_5, foldChange=class3_uniID5) #Gene-Concept Network


#class3_serected
class3_uni6 <- keggConv("uniprot",class3_serected_ID)
class3_uniID6 <- c()
for (i in class3_uni6){
  j <- substring(i, 4,)
  class3_uniID6 <- append(class3_uniID6,j)
}
go3_6 <- enrichGO(gene = class3_uniID6,
                  OrgDb = "org.Hs.eg.db",
                  keyType = 'UNIPROT',
                  ont = 'all',
                  readable = TRUE, # regard as gene names
                  pvalueCutoff = 0.01,
                  qvalueCutoff  = 0.05,
)
dotplot(go3_6, split="ONTOLOGY", font.size = 20, showCategory=5) + facet_grid(ONTOLOGY~., scale="free") #Dot plot all
barplot(go3_6, showCategory=8,x = "Count")
write.csv(go3_6, "class3_serected.csv")

egoMF3_6<- enrichGO(gene = class3_uniID6,
                    OrgDb = "org.Hs.eg.db",
                    keyType = 'UNIPROT',
                    ont = 'MF',
)
egoBP3_6<- enrichGO(gene = class3_uniID6,
                    OrgDb = "org.Hs.eg.db",
                    keyType = 'UNIPROT',
                    ont = 'BP',
)
egoCC3_6<- enrichGO(gene = class3_uniID6,
                    OrgDb = "org.Hs.eg.db",
                    keyType = 'UNIPROT',
                    ont = 'CC',
)
goplot(egoMF3_6) #Induced GO DAG graph
goplot(egoBP3_6) 
goplot(egoCC3_6) 

ego23_6 <- simplify(egoBP3_6)
cnetplot(ego23_6, foldChange=class3_uniID6) #Gene-Concept Network


#class4
class4_ID <- class4$KEGGID
class4_uni <- keggConv("uniprot",class4$KEGGID)
class4_uniID <- c()
for (i in class4_uni){
  j <- substring(i, 4,)
  class4_uniID <- append(class4_uniID,j)
}

go4 <- enrichGO(gene = class4_uniID,
               OrgDb = "org.Hs.eg.db",
               keyType = 'UNIPROT',
               ont = 'all',
               readable = TRUE, # regard as gene names
               pvalueCutoff = 0.01,
               qvalueCutoff  = 0.05,
)
dotplot(go4, split="ONTOLOGY", font.size = 20, showCategory=5) + facet_grid(ONTOLOGY~., scale="free") #Dot plot all
write.csv(go4, "class4.csv")
egoMF4<- enrichGO(gene = class4_uniID,
                 OrgDb = "org.Hs.eg.db",
                 keyType = 'UNIPROT',
                 ont = 'MF',
)
egoBP4<- enrichGO(gene = class4_uniID,
                  OrgDb = "org.Hs.eg.db",
                  keyType = 'UNIPROT',
                  ont = 'BP',
)
egoCC4<- enrichGO(gene = class4_uniID,
                  OrgDb = "org.Hs.eg.db",
                  keyType = 'UNIPROT',
                  ont = 'CC',
)
goplot(egoMF4) #Induced GO DAG graph
goplot(egoBP4)
goplot(egoCC4)

ego24 <- simplify(egoBP4)
cnetplot(ego24, foldChange=class4_uniID) #Gene-Concept Network


#class4_collagen
class4_uni1 <- keggConv("uniprot",class4_collagen_ID)
class4_uniID1 <- c()
for (i in class4_uni1){
  j <- substring(i, 4,)
  class4_uniID1 <- append(class4_uniID1,j)
}
go4_1 <- enrichGO(gene = class4_uniID1,
                  OrgDb = "org.Hs.eg.db",
                  keyType = 'UNIPROT',
                  ont = 'all',
                  readable = TRUE, # regard as gene names
                  pvalueCutoff = 0.01,
                  qvalueCutoff  = 0.05,
)
dotplot(go4_1, split="ONTOLOGY", font.size = 20, showCategory=5) + facet_grid(ONTOLOGY~., scale="free") #Dot plot all
barplot(go4_1, showCategory=8,x = "Count")
write.csv(go4_1, "class4_collagen.csv")

egoMF4_1<- enrichGO(gene = class4_uniID1,
                    OrgDb = "org.Hs.eg.db",
                    keyType = 'UNIPROT',
                    ont = 'MF',
)
egoBP4_1<- enrichGO(gene = class4_uniID1,
                    OrgDb = "org.Hs.eg.db",
                    keyType = 'UNIPROT',
                    ont = 'BP',
)
egoCC4_1<- enrichGO(gene = class4_uniID1,
                    OrgDb = "org.Hs.eg.db",
                    keyType = 'UNIPROT',
                    ont = 'CC',
)
goplot(egoMF4_1) #Induced GO DAG graph
goplot(egoBP4_1) 
goplot(egoCC4_1) 

ego24_1 <- simplify(egoBP4_1)
cnetplot(ego24_1, foldChange=class4_uniID1) #Gene-Concept Network


#class4_glycoprotein
class4_uni2 <- keggConv("uniprot",class4_glycoprotein_ID)
class4_uniID2 <- c()
for (i in class4_uni2){
  j <- substring(i, 4,)
  class4_uniID2 <- append(class4_uniID2,j)
}
go4_2 <- enrichGO(gene = class4_uniID2,
                  OrgDb = "org.Hs.eg.db",
                  keyType = 'UNIPROT',
                  ont = 'all',
                  readable = TRUE, # regard as gene names
                  pvalueCutoff = 0.01,
                  qvalueCutoff  = 0.05,
)
dotplot(go4_2, split="ONTOLOGY", font.size = 20, showCategory=5) + facet_grid(ONTOLOGY~., scale="free") #Dot plot all
barplot(go4_2, showCategory=8,x = "Count")
write.csv(go4_2, "class4_glycoprotein.csv")

egoMF4_2<- enrichGO(gene = class4_uniID2,
                    OrgDb = "org.Hs.eg.db",
                    keyType = 'UNIPROT',
                    ont = 'MF',
)
egoBP4_2<- enrichGO(gene = class4_uniID2,
                    OrgDb = "org.Hs.eg.db",
                    keyType = 'UNIPROT',
                    ont = 'BP',
)
egoCC4_2<- enrichGO(gene = class4_uniID2,
                    OrgDb = "org.Hs.eg.db",
                    keyType = 'UNIPROT',
                    ont = 'CC',
)
goplot(egoMF4_2) #Induced GO DAG graph
goplot(egoBP4_2) 
goplot(egoCC4_2) 

ego24_2 <- simplify(egoBP4_2)
cnetplot(ego24_2, foldChange=class4_uniID2) #Gene-Concept Network


#class4_proteoglycans

class4_uni3 <- keggConv("uniprot",class4_proteoglycans_ID)
class4_uniID3 <- c()
for (i in class4_uni3){
  j <- substring(i, 4,)
  class4_uniID3 <- append(class4_uniID3,j)
}
go4_3 <- enrichGO(gene = class4_uniID3,
                  OrgDb = "org.Hs.eg.db",
                  keyType = 'UNIPROT',
                  ont = 'all',
                  readable = TRUE, # regard as gene names
                  pvalueCutoff = 0.01,
                  qvalueCutoff  = 0.05,
)
dotplot(go4_3, split="ONTOLOGY", font.size = 20, showCategory=5) + facet_grid(ONTOLOGY~., scale="free") #Dot plot all
barplot(go4_3, showCategory=8,x = "Count")
write.csv(go4_3, "class4_proteoglycans.csv")

egoMF4_3<- enrichGO(gene = class4_uniID3,
                    OrgDb = "org.Hs.eg.db",
                    keyType = 'UNIPROT',
                    ont = 'MF',
)
egoBP4_3<- enrichGO(gene = class4_uniID3,
                    OrgDb = "org.Hs.eg.db",
                    keyType = 'UNIPROT',
                    ont = 'BP',
)
egoCC4_3<- enrichGO(gene = class4_uniID3,
                    OrgDb = "org.Hs.eg.db",
                    keyType = 'UNIPROT',
                    ont = 'CC',
)
goplot(egoMF4_3) #Induced GO DAG graph
goplot(egoBP4_3) 
goplot(egoCC4_3) 

ego24_3 <- simplify(egoBP4_3)
cnetplot(ego24_3, foldChange=class4_uniID3) #Gene-Concept Network


#class3_affiated
class4_uni4 <- keggConv("uniprot",class4_affiliated_ID)
class4_uniID4 <- c()
for (i in class4_uni4){
  j <- substring(i, 4,)
  class4_uniID4 <- append(class4_uniID4,j)
}
go4_4 <- enrichGO(gene = class4_uniID4,
                  OrgDb = "org.Hs.eg.db",
                  keyType = 'UNIPROT',
                  ont = 'all',
                  readable = TRUE, # regard as gene names
                  pvalueCutoff = 0.01,
                  qvalueCutoff  = 0.05,
)
dotplot(go4_4, split="ONTOLOGY", font.size = 20, showCategory=5) + facet_grid(ONTOLOGY~., scale="free") #Dot plot all
barplot(go4_4, showCategory=8,x = "Count")
write.csv(go4_4, "class4_affiated.csv")

egoMF4_4<- enrichGO(gene = class4_uniID4,
                    OrgDb = "org.Hs.eg.db",
                    keyType = 'UNIPROT',
                    ont = 'MF',
)
egoBP4_4<- enrichGO(gene = class4_uniID4,
                    OrgDb = "org.Hs.eg.db",
                    keyType = 'UNIPROT',
                    ont = 'BP',
)
egoCC4_4<- enrichGO(gene = class4_uniID4,
                    OrgDb = "org.Hs.eg.db",
                    keyType = 'UNIPROT',
                    ont = 'CC',
)
goplot(egoMF4_4) #Induced GO DAG graph
goplot(egoBP4_4) 
goplot(egoCC4_4) 

ego24_4 <- simplify(egoBP4_4)
cnetplot(ego24_4, foldChange=class4_uniID4) #Gene-Concept Network

#class4_regulater
class4_uni5 <- keggConv("uniprot",class4_regulater_ID)
class4_uniID5 <- c()
for (i in class4_uni5){
  j <- substring(i, 4,)
  class4_uniID5 <- append(class4_uniID5,j)
}
go4_5 <- enrichGO(gene = class4_uniID5,
                  OrgDb = "org.Hs.eg.db",
                  keyType = 'UNIPROT',
                  ont = 'all',
                  readable = TRUE, # regard as gene names
                  pvalueCutoff = 0.01,
                  qvalueCutoff  = 0.05,
)
dotplot(go4_5, split="ONTOLOGY", font.size = 20, showCategory=5) + facet_grid(ONTOLOGY~., scale="free") #Dot plot all
barplot(go4_5, showCategory=8,x = "Count")
write.csv(go4_5, "class4_regulater.csv")

egoMF4_5<- enrichGO(gene = class4_uniID5,
                    OrgDb = "org.Hs.eg.db",
                    keyType = 'UNIPROT',
                    ont = 'MF',
)
egoBP4_5<- enrichGO(gene = class4_uniID5,
                    OrgDb = "org.Hs.eg.db",
                    keyType = 'UNIPROT',
                    ont = 'BP',
)
egoCC4_5<- enrichGO(gene = class4_uniID5,
                    OrgDb = "org.Hs.eg.db",
                    keyType = 'UNIPROT',
                    ont = 'CC',
)
goplot(egoMF4_5) #Induced GO DAG graph
goplot(egoBP4_5) 
goplot(egoCC4_5) 

ego24_5 <- simplify(egoBP4_5)
cnetplot(ego24_5, foldChange=class4_uniID5) #Gene-Concept Network


#class4_serected
class4_uni6 <- keggConv("uniprot",class4_serected_ID)
class4_uniID6 <- c()
for (i in class4_uni6){
  j <- substring(i, 4,)
  class4_uniID6 <- append(class4_uniID6,j)
}
go4_6 <- enrichGO(gene = class4_uniID6,
                  OrgDb = "org.Hs.eg.db",
                  keyType = 'UNIPROT',
                  ont = 'all',
                  readable = TRUE, # regard as gene names
                  pvalueCutoff = 0.01,
                  qvalueCutoff  = 0.05,
)
dotplot(go4_6, split="ONTOLOGY", font.size = 20, showCategory=5) + facet_grid(ONTOLOGY~., scale="free") #Dot plot all
barplot(go4_6, showCategory=8,x = "Count")
write.csv(go4_6, "class4_serected.csv")

egoMF4_6<- enrichGO(gene = class4_uniID6,
                    OrgDb = "org.Hs.eg.db",
                    keyType = 'UNIPROT',
                    ont = 'MF',
)
egoBP4_6<- enrichGO(gene = class4_uniID6,
                    OrgDb = "org.Hs.eg.db",
                    keyType = 'UNIPROT',
                    ont = 'BP',
)
egoCC4_6<- enrichGO(gene = class4_uniID6,
                    OrgDb = "org.Hs.eg.db",
                    keyType = 'UNIPROT',
                    ont = 'CC',
)
goplot(egoMF4_6) #Induced GO DAG graph
goplot(egoBP4_6) 
goplot(egoCC4_6) 

ego24_6 <- simplify(egoBP4_6)
cnetplot(ego24_6, foldChange=class4_uniID6) #Gene-Concept Network
# コメント
# class5_ID <- class5$KEGGID
# class5_uni <- keggConv("uniprot",class5$KEGGID)
# class5_uniID <- c()
# for (i in class5_uni){
#   j <- substring(i, 4,)
#   class5_uniID <- append(class5_uniID,j)
# }
# go <- enrichGO(gene = class5_uniID,
#                OrgDb = "org.Hs.eg.db",
#                keyType = 'UNIPROT',
#                ont = 'all',
#                readable = TRUE, # regard as gene names
#                pvalueCutoff = 0.01,
#                qvalueCutoff  = 0.05,
# )
# dotplot(go, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free") #Dot plot all
# 
# egoBP<- enrichGO(gene = class5_uniID,
#                  OrgDb = "org.Hs.eg.db",
#                  keyType = 'UNIPROT',
#                  ont = 'BP',
# )
# goplot(egoBP) #Induced GO DAG graph
# 
# egoBP<- enrichGO(gene = class5_uniID,
#                  OrgDb = "org.Hs.eg.db",
#                  keyType = 'UNIPROT',
#                  ont = 'BP',
# )
# ego2 <- simplify(egoBP)
# cnetplot(ego2, foldChange=class5_uniID) #Gene-Concept Network
# 
# 
# class6_ID <- class6$KEGGID
# class6_uni <- keggConv("uniprot",class6$KEGGID)
# class6_uniID <- c()
# for (i in class6_uni){
#   j <- substring(i, 4,)
#   class6_uniID <- append(class6_uniID,j)
# }
# 
# go <- enrichGO(gene = class6_uniID,
#                OrgDb = "org.Hs.eg.db",
#                keyType = 'UNIPROT',
#                ont = 'all',
#                readable = TRUE, # regard as gene names
#                pvalueCutoff = 0.01,
#                qvalueCutoff  = 0.05,
# )
# dotplot(go, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free") #Dot plot all
# 
# egoBP<- enrichGO(gene = class6_uniID,
#                  OrgDb = "org.Hs.eg.db",
#                  keyType = 'UNIPROT',
#                  ont = 'BP',
# )
# goplot(egoBP) #Induced GO DAG graph
# 
# egoBP<- enrichGO(gene = class6_uniID,
#                  OrgDb = "org.Hs.eg.db",
#                  keyType = 'UNIPROT',
#                  ont = 'BP',
# )
# ego2 <- simplify(egoBP)
# cnetplot(ego2, foldChange=class6_uniID) #Gene-Concept Network
# 
# 
# class7_ID <- class7$KEGGID
# class7_uni <- keggConv("uniprot",class7$KEGGID)
# class7_uniID <- c()
# for (i in class7_uni){
#   j <- substring(i, 4,)
#   class7_uniID <- append(class7_uniID,j)
# }
# 
# go <- enrichGO(gene = class7_uniID,
#                OrgDb = "org.Hs.eg.db",
#                keyType = 'UNIPROT',
#                ont = 'all',
#                readable = TRUE, # regard as gene names
#                pvalueCutoff = 0.01,
#                qvalueCutoff  = 0.05,
# )
# dotplot(go, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free") #Dot plot all
# 
# egoBP<- enrichGO(gene = class7_uniID,
#                  OrgDb = "org.Hs.eg.db",
#                  keyType = 'UNIPROT',
#                  ont = 'BP',
# # )
# goplot(egoBP) #Induced GO DAG graph
# 
# egoBP<- enrichGO(gene = class7_uniID,
#                  OrgDb = "org.Hs.eg.db",
#                  keyType = 'UNIPROT',
#                  ont = 'BP',
# )
# ego2 <- simplify(egoBP)
# cnetplot(ego2, foldChange=class7_uniID) #Gene-Concept Network
# 
# 
# class8_ID <- class8$KEGGID
# class8_uni <- keggConv("uniprot",class8$KEGGID)
# class8_uniID <- c()
# for (i in class8_uni){
#   j <- substring(i, 4,)
#   class8_uniID <- append(class8_uniID,j)
# }
# 
# go <- enrichGO(gene = class8_uniID,
#                OrgDb = "org.Hs.eg.db",
#                keyType = 'UNIPROT',
#                ont = 'all',
#                readable = TRUE, # regard as gene names
#                pvalueCutoff = 0.01,
#                qvalueCutoff  = 0.05,
# )
# dotplot(go, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free") #Dot plot all
# 
# egoBP<- enrichGO(gene = class8_uniID,
#                  OrgDb = "org.Hs.eg.db",
#                  keyType = 'UNIPROT',
#                  ont = 'BP',
# )
# goplot(egoBP) #Induced GO DAG graph
# 
# egoBP<- enrichGO(gene = class8_uniID,
#                  OrgDb = "org.Hs.eg.db",
#                  keyType = 'UNIPROT',
#                  ont = 'BP',
# )
#  ego2 <- simplify(egoBP)
#  cnetplot(ego2, foldChange=class8_uniID) #Gene-Concept Network

