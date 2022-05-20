import pandas as pd
import numpy as np

inter=pd.read_csv('CTAON_interactions.txt',sep='\t')
print(inter['RNA_NAME'].unique().shape)
arr=np.random.choice(inter['RNA_NAME'].unique(),size=990,replace=False)

inter=inter[inter['RNA_NAME'].isin(arr)]

print("Number of rna: ",end='')
print(inter['RNA_NAME'].unique().shape)

print("Number of prot: ",end='')
print(inter['PROT_NAME'].unique().shape)

print("Number of interactions: ",end='')
print(inter.shape)

f=open('sampled_rna.txt','w')
for itern_cnt,i in enumerate(inter['RNA_NAME'].unique()):
    f.write(str(itern_cnt)+"\t"+str(i)+"\n")
f.close()

f=open('sampled_prot.txt','w')
for itern_cnt,i in enumerate(inter['PROT_NAME'].unique()):
    f.write(str(itern_cnt)+"\t"+str(i)+"\n")
f.close()

inter.to_csv('sampled_inter.tsv',sep='\t')

# (18812,)
# Number of rna: (990,)
# Number of prot: (64,)
# Number of interactions: (3714, 4)

