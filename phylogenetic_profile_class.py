import pandas as pd

path = "C:/Users/user/OneDrive/デスクトップ/研究室/python_program/phylogenetic_profile/"

profile = pd.read_csv(path + 'results/reslt_af_output_2.csv')
info_class = pd.read_excel(path + 'results/profile_result_fu_result.xlsx')

profile = profile.set_index("KEGGID")
info_class = info_class.set_index("KEGGID")

profile = profile.join(info_class, how = 'inner')

profile.to_csv(path + "results/reslt_af_output_3.csv")
