import numpy as np
import pandas as pd

def random_walk(inter, rna_cnt, similarity_matrix, prot_id, Wq, Wu, rq, ru, output_file_name):
    
    associated_rna=inter[inter.PROT_ID==prot_id].RNA_ID.to_numpy() 
    # rna in associated_rna == True if rna is labelled node
    
    R=similarity_matrix.copy() # correlation matrix
    for i in range(rna_cnt):
        if i in associated_rna:
            R[i]=R[i]*Wq
        else:
            R[i]=R[i]*Wu
    for i in range(rna_cnt):
        sigma=np.sum(R[i])
        R[i]=R[i]/sigma
    
    Lq=R.copy()
    Lu=R.copy()
    
    for i in range(rna_cnt):
        if i in associated_rna:
            Lu[i]*=0
        else:
            Lq[i]*=0
    
    filt=np.zeros(rna_cnt)
    for i in associated_rna:
        filt[i]=1
    
    mod_Q=associated_rna.shape[0]
    X_init=np.zeros((rna_cnt,1))
    for i in associated_rna:
        X_init[i][0]=1/mod_Q
    
    X=X_init.copy()
    pq=np.matmul(filt,X)[0]
    pu=1-pq
    iter_cnt=0
    # print(pq)
    # print("i:"+str(iter_cnt)+" "+str(X[0][0]))
    
    while(True):
        prev=X.copy()
        X=rq*(np.matmul(Lq.transpose(),prev))+pq*(1-rq)*X_init+ru*(np.matmul(Lu.transpose(),prev))+pu*(1-ru)*X_init
        pq=np.matmul(filt,X)[0]
        pu=1-pq

        if(np.linalg.norm(X-prev,ord=1)<1e-10):
            break
        iter_cnt+=1
        # print("i:"+str(iter_cnt)+" "+str(X[0][0]))
    
    print(f"Iterations[{prot_id}]: {iter_cnt}")
    
    X=np.reshape(X,rna_cnt)
    ser=pd.Series(X)
    
    ser.sort_values(ascending=False, inplace=True)
    
    f=open(output_file_name,'a')
    f.write('\nPROTEIN ID:\t'+str(prot_id)+'\n')
    f.write('RNA\tCorr. score(Sr)\n')
    for index, value in ser.items():
        if index not in associated_rna:
            f.write(str(index)+'\t'+str(value)+'\n')
    f.close()
    return

rna_sim=pd.read_csv('rna_similarity.txt',sep='\t')
prot_sim=pd.read_csv('prot_similarity.txt',sep='\t')
inter=pd.read_csv('list_of_interactions.txt',sep='\t')

rna_cnt=inter['RNA_ID'].unique().shape[0]
prot_cnt=inter['PROT_ID'].unique().shape[0]

similarity_matrix=np.empty((rna_cnt,rna_cnt))
for i in rna_sim.index:
    similarity_matrix[rna_sim['RNA(i)'][i]][rna_sim['RNA(j)'][i]]=rna_sim['Sim(i,j)'][i]

output_file_name='random_walk_scores.txt'
Wq=0.8
Wu=0.4
rq=0.8
ru=0.4
open(output_file_name,'w').close()
for i in range(prot_cnt):
    random_walk(inter, rna_cnt,similarity_matrix,i,Wq,Wu,rq,ru,output_file_name)
    

