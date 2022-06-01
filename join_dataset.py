import pandas as pd

def join_dataset(path, path1, path2, path3, name):
    
    #pathの設定
    path1 = path + "/"+ path1
    path2 = path + "/"+ path2
    path3 = path + "/"+ path3

    #データセットを読み込む
    dataset1 = pd.read_excel(path1)
    dataset2 = pd.read_excel(path2)

    dataset1 = dataset1.set_index(name)
    dataset2 = dataset2.set_index(name)

    #共通の遺伝子を取り出す
    new_dataset = dataset1.join(dataset2, how='inner')

    print(new_dataset)
    new_dataset.to_excel(path3) 

if __name__ == '__main__':
    print("ファイル名を除いたpathを入力してください：")
    path = input()
    print("1つ目のファイル名を入力してください：")
    path1 = input()
    print("2つ目のファイル名を入力してください：")
    path2 = input()
    print("出力するファイル名を入力してください：")
    path3 = input()
    print("joinする行の名前を入力してください：")
    name = input()
    join_dataset(path, path1, path2, path3, name)