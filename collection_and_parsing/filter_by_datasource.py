import pandas as pd

inter=pd.read_csv('lncRNA_interaction.txt',sep='\t')
prot=pd.read_csv('combined_prot_table.tsv',sep='\t')
rna=pd.read_csv('trimmed_combined_table.csv')


inter.drop(['ncName','tarName','interDescription','experiment','reference','tag','class','level'], axis=1, inplace=True)

inter.drop(inter[inter.ncType!='lncRNA'].index, inplace=True)
inter.drop(inter[inter.tarType!='protein'].index, inplace=True)
inter.drop(inter[inter.organism != 'Homo sapiens'].index, inplace=True)

inter.drop(inter[inter.tarID=='-'].index, inplace=True)
inter.drop(inter[inter.ncID=='-'].index, inplace=True)

inter=inter[inter.datasource.str.contains("High-throughput", na=False)]
inter=inter[~inter.datasource.str.contains("prediction", na=False)]
inter.drop_duplicates()
print(inter['tissueOrCell'].unique())
print("Number of unique RNA: ")
print(len(inter['ncID'].unique()))
print("Number of unique protein: ")
print(len(inter['tarID'].unique()))
print("Number of unique interactions: ")
print(inter.shape[0])

inter.to_csv('HT_without_PRED.tsv',sep='\t') #output redirected to the file 'filtered_by_high_throughput_without_pred.txt'
print(inter.tail())


