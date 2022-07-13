#ライブラリの読み込み
library(ReactomePA)
library("openxlsx")
# ディレクトリの指定
setwd("C:\\R_practice\\Reactome_pathway")

# ファイルの読み込み 
dataset <- read.xlsx("ECM_Entrez.xlsx",) # excel_file

class1 <- dataset[dataset$class==1,]
class1_ID <- class1$Entrez_Gene_ID
class2 <- dataset[dataset$class==2,]
class2_ID <- class2$Entrez_Gene_ID
class3 <- dataset[dataset$class==3,]
class3_ID <- class3$Entrez_Gene_ID
class4 <- dataset[dataset$class==4,]
class4_ID <- class4$Entrez_Gene_ID


#Reactomeパスウェイ解析を実行
Reactome_enrichment_result1 <- enrichPathway(gene=class1_ID,pvalueCutoff=0.05, readable=T)
Reactome_enrichment_result2 <- enrichPathway(gene=class2_ID,pvalueCutoff=0.05, readable=T)
Reactome_enrichment_result3 <- enrichPathway(gene=class3_ID,pvalueCutoff=0.05, readable=T)
Reactome_enrichment_result4 <- enrichPathway(gene=class4_ID,pvalueCutoff=0.05, readable=T)

#関係図を表示
d1 <- GOSemSim::godata("org.Hs.eg.db", ont = "BP")    
Reactome_enrichment_result1 <- enrichplot::pairwise_termsim(Reactome_enrichment_result1, semData = d1)
emapplot(Reactome_enrichment_result1)
barplot(Reactome_enrichment_result1, showCategory=5,x = "Count")
dotplot(Reactome_enrichment_result1, showCategory=10, font.size = 20)

d2 <- GOSemSim::godata("org.Hs.eg.db", ont = "BP")    
Reactome_enrichment_result2 <- enrichplot::pairwise_termsim(Reactome_enrichment_result2, semData = d2)
emapplot(Reactome_enrichment_result2)
barplot(Reactome_enrichment_result2, showCategory=8,x = "Count")
dotplot(Reactome_enrichment_result2, showCategory=10, font.size = 20)

d3 <- GOSemSim::godata("org.Hs.eg.db", ont = "BP")    
Reactome_enrichment_result3 <- enrichplot::pairwise_termsim(Reactome_enrichment_result3, semData = d3)
emapplot(Reactome_enrichment_result3)
barplot(Reactome_enrichment_result3, showCategory=8,x = "Count")
dotplot(Reactome_enrichment_result3, showCategory=10, font.size = 20)

d4 <- GOSemSim::godata("org.Hs.eg.db", ont = "BP")    
Reactome_enrichment_result4 <- enrichplot::pairwise_termsim(Reactome_enrichment_result4, semData = d4)
emapplot(Reactome_enrichment_result4)
barplot(Reactome_enrichment_result4, showCategory=8,x = "Count")
dotplot(Reactome_enrichment_result4, showCategory=10, font.size = 20)

#結果の保存
write.csv(Reactome_enrichment_result1, "Reactome_pathway1.csv")
write.csv(Reactome_enrichment_result2, "Reactome_pathway2.csv")
write.csv(Reactome_enrichment_result3, "Reactome_pathway3.csv")
write.csv(Reactome_enrichment_result4, "Reactome_pathway4.csv")
