import numpy as np
import pandas as pd

def random_walk(Y, rna_cnt, Sl, prot_id, Wq, Wu, rq, ru):
    
    associated_rna=np.array([i for i in range(rna_cnt) if Y[i][prot_id]==1])
    # rna in associated_rna == True if rna is labelled node
    
    R=Sl.copy() # correlation matrix
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
    
    while(True):
        prev=X.copy()
        X=rq*(np.matmul(Lq.transpose(),prev))+pq*(1-rq)*X_init+ru*(np.matmul(Lu.transpose(),prev))+pu*(1-ru)*X_init
        pq=np.matmul(filt,X)[0]
        pu=1-pq

        if(np.linalg.norm(X-prev,ord=1)<1e-10):
            break
        iter_cnt+=1
    
    return X

def RW(rna_cnt,prot_cnt,Sl,Y,neg_prot_set=[],Wq=0.8,Wu=0.4,rq=0.8,ru=0.4):
    Sr=np.empty((rna_cnt,prot_cnt))
    for i in range(prot_cnt):
        if i in neg_prot_set:
            for j in range(rna_cnt):
                Sr[j][i]=1/rna_cnt
        else:
            X=random_walk(Y,rna_cnt,Sl,i,Wq,Wu,rq,ru)
            for j in range(rna_cnt):
                Sr[j][i]=X[j][0]
    return Sr

if __name__=='__main__':
    
    rna_sim=pd.read_csv('rna_similarity.txt',sep='\t')
    inter=pd.read_csv('list_of_interactions.txt',sep='\t')

    rna_cnt=inter['RNA_ID'].unique().shape[0] # m
    prot_cnt=inter['PROT_ID'].unique().shape[0] # n

    Sl=np.empty((rna_cnt,rna_cnt))

    for i in rna_sim.index:
        Sl[rna_sim['RNA(i)'][i]][rna_sim['RNA(j)'][i]]=rna_sim['Sim(i,j)'][i]

    Y=np.zeros((rna_cnt,prot_cnt))
    for i in inter.index:
        Y[inter['RNA_ID'][i]][inter['PROT_ID'][i]]=1

    Sr=RW(rna_cnt,prot_cnt,Sl,Y)

    f=open('RW_scores','w')
    f.close()
    
    for i in range(prot_cnt):
        ser=pd.Series(Sr[:,i])
        ser.sort_values(ascending=False, inplace=True)
        associated_rna=np.array([j for j in range(rna_cnt) if Y[j][i]==1])
        
        f=open('RW_scores','a')
        f.write('\nPROTEIN ID:\t'+str(i)+'\n')
        f.write('RNA\tCorr. score(Sr)\n')
        for index, value in ser.items():
            if index not in associated_rna:
                f.write(str(index)+'\t'+str(value)+'\n')
        f.close()