#internetのアクセス
import time
import requests
from bs4 import BeautifulSoup
import pandas as pd
import re

#データセットを読み込む
path = "C:/Users/user/OneDrive/デスクトップ/研究室/python_program/kegg_disease/"
dataset = "ECM_dataset.xlsx"
Gene_class_all = pd.read_excel(path + dataset)

#クラス分け
group = Gene_class_all.groupby('class')
Gene_class_1 = group.get_group(1)
Gene_class_2 = group.get_group(2)
Gene_class_3 = group.get_group(3)
Gene_class_4 = group.get_group(4)

#カテゴリ分け
group1 = Gene_class_1.groupby('CATEGORY')
#Gene_class_11 = group1.get_group('E_collagen')
Gene_class_12 = group1.get_group('E_glycoprotein')
Gene_class_13 = group1.get_group('E_proteoglycans')
Gene_class_14 = group1.get_group('EA_affiliated')
Gene_class_15 = group1.get_group('EA_regulater')
Gene_class_16 = group1.get_group('EA_serected')

group2 = Gene_class_2.groupby('CATEGORY')
Gene_class_21 = group2.get_group('E_collagen')
Gene_class_22 = group2.get_group('E_glycoprotein')
Gene_class_23 = group2.get_group('E_proteoglycans')
Gene_class_24 = group2.get_group('EA_affiliated')
Gene_class_25 = group2.get_group('EA_regulater')
Gene_class_26 = group2.get_group('EA_serected')

group3 = Gene_class_3.groupby('CATEGORY')
Gene_class_31 = group3.get_group('E_collagen')
Gene_class_32 = group3.get_group('E_glycoprotein')
Gene_class_33 = group3.get_group('E_proteoglycans')
Gene_class_34 = group3.get_group('EA_affiliated')
Gene_class_35 = group3.get_group('EA_regulater')
Gene_class_36 = group3.get_group('EA_serected')

group4 = Gene_class_4.groupby('CATEGORY')
#Gene_class_41 = group4.get_group('E_collagen')
Gene_class_42 = group4.get_group('E_glycoprotein')
#Gene_class_43 = group4.get_group('E_proteoglycans')
Gene_class_44 = group4.get_group('EA_affiliated')
Gene_class_45 = group4.get_group('EA_regulater')
#Gene_class_46 = group4.get_group('EA_serected')

#各クラスの遺伝子群のseries
Gene_class_1 = Gene_class_1['KEGGID']
#Gene_class_11 = Gene_class_11['KEGGID']
Gene_class_12 = Gene_class_12['KEGGID']
Gene_class_13 = Gene_class_13['KEGGID']
Gene_class_14 = Gene_class_14['KEGGID']
Gene_class_15 = Gene_class_15['KEGGID']
Gene_class_16 = Gene_class_16['KEGGID']

Gene_class_2 = Gene_class_2['KEGGID']
Gene_class_21 = Gene_class_21['KEGGID']
Gene_class_22 = Gene_class_22['KEGGID']
Gene_class_23 = Gene_class_23['KEGGID']
Gene_class_24 = Gene_class_24['KEGGID']
Gene_class_25 = Gene_class_25['KEGGID']
Gene_class_26 = Gene_class_26['KEGGID']

Gene_class_3 = Gene_class_3['KEGGID']
Gene_class_31 = Gene_class_31['KEGGID']
Gene_class_32 = Gene_class_32['KEGGID']
Gene_class_33 = Gene_class_33['KEGGID']
Gene_class_34 = Gene_class_34['KEGGID']
Gene_class_35 = Gene_class_35['KEGGID']
Gene_class_36 = Gene_class_36['KEGGID']

Gene_class_4 = Gene_class_4['KEGGID']
#Gene_class_41 = Gene_class_41['KEGGID']
Gene_class_42 = Gene_class_42['KEGGID']
#Gene_class_43 = Gene_class_43['KEGGID']
Gene_class_44 = Gene_class_44['KEGGID']
Gene_class_45 = Gene_class_45['KEGGID']
#Gene_class_46 = Gene_class_46['KEGGID']


#出力用のDataframe
Kegg_disease = pd.DataFrame({}, columns = ['KEGGID', 'DiseaseID', 'Diseaese_Name', 'Disease_category1', 'Disease_category2'])
NULL = ''

#　KEDD ID　→　疾患ID　

for gene in Gene_class_45:
    time.sleep(3)
    #HTMLテキスト取得
    url = "https://www.genome.jp/entry/" + gene
    r = requests.get(url)
    soup = BeautifulSoup(r.content, 'html.parser')
    #trタグのみ取得
    th = soup.find_all('tr')
    tmp = 0
    for i in th:
        j = i.text
        if tmp == 2:
            tex = j
        tmp += 1
    #疾患の情報を取得
    columns = tex.splitlines()
    tmp = 0
    loc = 0
    for i in columns:
        if i == 'Disease':
            loc = tmp
        tmp += 1
    colID = columns[loc+1]
    diseaseID = re.findall('[H]\d\d\d\d\d', colID)

    print(gene)
    #print(diseaseID)
    if diseaseID == []:
        print('なかった')
        now = pd.DataFrame({'KEGGID':[gene],'DiseaseID':[NULL],'Diseaese_Name':[NULL],'Disease_category1':[NULL],'Disease_category2':[NULL]})
        Kegg_disease = Kegg_disease.append(now)
    else:
        #疾患ID　→　疾患名、疾患カテゴリ
        for id in diseaseID:
            time.sleep(3)
            #HTMLテキスト取得
            url2 = "https://www.genome.jp/entry/" + id
            r2 = requests.get(url2)
            soup2 = BeautifulSoup(r2.content, 'html.parser')
            #trタグのみ取得
            th2 = soup2.find_all('tr')
            tmp = 0
            for i in th2:
                j = i.text
                m = re.search('Brite', str(i))
                if m:
                    tex_category = j
                tmp += 1


            #疾患カテゴリを取得
            list_category = tex_category.splitlines()
            tmp=0
            loc=0
            for i in list_category:
                m = re.search(id, i)
                if m:
                    loc = tmp
                tmp += 1
            if loc!=0:
                    Disease_name_id = list_category[4]#疾患名がリストの4番目にあった
                    Disease_category_small = list_category[3]#疾患名より一つ上の階層カテゴリ
                    Disease_category_big = list_category[2]#疾患名より二つ上の階層カテゴリ

                    #diseaseの階層構造のインデントを消去
                    Disease_category_big = Disease_category_big.strip()
                    Disease_category_small = Disease_category_small.strip()
                    Disease_name_id = Disease_name_id.strip()

                    #疾患名の前のid消去
                    Disease_name_id = Disease_name_id.lstrip('H0123456789')
                    Disease_name = Disease_name_id.strip()

                    print(id)
                    print(Disease_name)
                    print(Disease_category_small)
                    print(Disease_category_big)
                    now = pd.DataFrame({'KEGGID':[gene],'DiseaseID':[id],'Diseaese_Name':[Disease_name],'Disease_category1':[Disease_category_small],'Disease_category2':[Disease_category_big]})
                    Kegg_disease = Kegg_disease.append(now)
                    #rint(Kegg_disease)
Kegg_disease.to_excel(path + 'Kegg_disease.xlsx') 
