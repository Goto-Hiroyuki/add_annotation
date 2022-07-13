#internetのアクセス
from asyncio.windows_events import NULL
import string
import time
import requests
from bs4 import BeautifulSoup
import pandas as pd
import re

def STRING(path, name, name2):
    #出力用のデータを作成
    string_data = pd.DataFrame({})
    #データセットを読み込む
    name = "ECM_dataset.xlsx"
    dataset = pd.read_excel(path + "/" + name)

    #データセットからGene Symbolを抜き出す
    gene_symbol = dataset['Gene']

    #データセットを、joinしやすい形にしておく
    dataset = dataset.set_index('Gene')

    for gene in gene_symbol:
        #STRINGのAPIに接続
        time.sleep(3)
        url = "https://string-db.org/api/tsv-no-header/interaction_partners?identifiers=" + gene
        r = requests.get(url)
        
        #取得した情報をDataFrameに変換
        soup = BeautifulSoup(r.content, 'html.parser') 
        str_soup = str(soup)

        #APIで取得した情報が存在するかの判断
        if not str_soup:
            print('No data')
            no_data = pd.DataFrame({'gene':[gene],'associate gene':[NULL],'CATEGORY_gene':[NULL],'CATEGORY_associate':[NULL],'KEGGID_gene':[NULL], 'KEGGID_associate':[NULL], 'string_id_gene':[NULL], 'string_id_associate':[NULL]})
            string_data = string_data.append(no_data)
        else:
            list_soup = str_soup.split("\n") # 行区切りにリスト化
            dataframe = pd.DataFrame(list_soup)
            dataframe = dataframe[0].str.split('\t', expand=True) # 列区切り
            dataframe = dataframe.drop(range(4, len(dataframe.columns)), axis=1) # 必要のない要素の消去
            dataframe = dataframe.drop(len(dataframe)-1) # nullデータの消去
            
            #前のアノテーション情報の付加
            dataframe = dataframe.set_index(2)
            dataframe = dataframe.join(dataset, how = "left")

            #後ろのアノテーション情報の付加
            dataframe = dataframe.reset_index()#インデックスを振り直す前に、resetすることでデータの消失が防げる
            dataframe = dataframe.set_index(3)
            dataframe = dataframe.join(dataset, how = "left", lsuffix='_gene', rsuffix='_associate')
            dataframe = dataframe.reset_index()

            #見やすい表データに変更
            dataframe = dataframe.rename(columns = {'index':'gene', 3:'associate gene', 0:'string_id_gene', 1:'string_id_associate'})
            dataframe = dataframe.reindex(columns=['gene','associate gene','CATEGORY_gene','CATEGORY_associate','class_gene','class_associate','KEGGID_gene','KEGGID_associate','string_id_gene','string_id_associate'])
            string_data = string_data.append(dataframe)
            print(string_data)
    #エクセルファイルに保存
    string_data.to_excel(path + "/" + name2) 

if __name__ == '__main__':
    print("ファイル名を除いたpathを入力してください：")
    path = input()
    print("ファイル名を入力してください：")
    name = input()
    print("出力するファイル名を入力してください：")
    name2 = input()
    STRING(path, name, name2)