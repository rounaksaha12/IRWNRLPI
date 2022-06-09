#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt 


rna_sim=pd.read_csv('rna_similarity.txt',sep='\t')
prot_sim=pd.read_csv('prot_similarity.txt',sep='\t')
inter=pd.read_csv('list_of_interactions.txt',sep='\t')

rna_cnt=inter['RNA_ID'].unique().shape[0] # m
prot_cnt=inter['PROT_ID'].unique().shape[0] # n

Sl=np.empty((rna_cnt,rna_cnt))
Sp=np.empty((prot_cnt,prot_cnt))

Y=np.zeros((rna_cnt,prot_cnt))
for i in inter.index:
    Y[inter['RNA_ID'][i]][inter['PROT_ID'][i]]=1
# print(Y)

for i in rna_sim.index:
    Sl[rna_sim['RNA(i)'][i]][rna_sim['RNA(j)'][i]]=rna_sim['Sim(i,j)'][i]

for i in prot_sim.index:
    Sp[prot_sim['PROT(i)'][i]][prot_sim['PROT(j)'][i]]=prot_sim['Sim(i,j)'][i]

# print(Sl)
# print(Sl.shape)

# print(Sp)
# print(Sp.shape)

def construct_neighbourhood(S,M,K,cnt):
    for i in range(cnt):
        ser=pd.Series(S[i])
        ser.sort_values(ascending=False, inplace=True)
        for iter_cnt, j in enumerate(ser.index):
            if iter_cnt>K:
                break
            if j==i:
                continue
            M[i][j]=ser[j]

K1=5
K2=5

A=np.zeros((rna_cnt,rna_cnt))
B=np.zeros((prot_cnt,prot_cnt))

construct_neighbourhood(Sl,A,K1,rna_cnt)
construct_neighbourhood(Sp,B,K1,prot_cnt)

# construction of Ll
Dl=np.zeros((rna_cnt,rna_cnt))
Dl_=np.zeros((rna_cnt,rna_cnt))

for i in range(rna_cnt):
    Dl[i][i]=np.sum(A[i])
    Dl_[i][i]=np.sum(A[:,i])

Ll=Dl+Dl_-(A+A.transpose())

# construction of Lp
Dp=np.zeros((prot_cnt,prot_cnt))
Dp_=np.zeros((prot_cnt,prot_cnt))

for i in range(prot_cnt):
    Dp[i][i]=np.sum(B[i])
    Dp_[i][i]=np.sum(B[:,i])
    
Lp=Dp+Dp_-(B+B.transpose())


# print(Ll)
# print(Lp)

def prob(U,V):
    temp=np.exp(np.matmul(U,V.transpose()))
    return temp/(1+temp)

# P=prob(U,V)
# print(P)

def dL_dU(P,U,V,c,Y,Ll,Lambda_l,alpha,rna_cnt):
    return np.matmul(P,V)+np.matmul(np.multiply(Y,P),V)*(c-1)-np.matmul(Y,V)*c+np.matmul((Lambda_l*np.identity(rna_cnt)+alpha*Ll),U)

def dL_dV(P,U,V,c,Y,Lp,Lambda_p,beta,prot_cnt):
    return np.matmul(P.transpose(),U)+np.matmul(np.multiply(Y.transpose(),P.transpose()),U)*(c-1)-np.matmul(Y.transpose(),U)*c+np.matmul((Lambda_p*np.identity(prot_cnt)+beta*Lp),V)

def grad_des(U,V,Y,c,Ll,Lp,Lambda_l,Lambda_p,alpha,beta,gamma,rna_cnt,prot_cnt,max_iter):
    phi_l=np.zeros((rna_cnt,r))
    phi_p=np.zeros((prot_cnt,r))
    iter_cnt=1
    cont=True
    while(cont):
        P=prob(U,V)
        prev_U=U.copy()
        prev_V=V.copy()
        
        Gl=dL_dU(P,U,V,c,Y,Ll,Lambda_l,alpha,rna_cnt)
        phi_l=phi_l+np.multiply(Gl,Gl)
        U=U-gamma*(Gl/np.sqrt(phi_l))

        Gp=dL_dV(P,U,V,c,Y,Lp,Lambda_p,beta,prot_cnt)
        phi_p=phi_p+np.multiply(Gp,Gp)
        V=V-gamma*(Gp/np.sqrt(phi_p))
        
        # print(f'Iteration[{iter_cnt}]: ',end='')
        iter_cnt+=1
        # print(V[0])
        
        cont=False
        for i in range(rna_cnt):
            if(np.linalg.norm(U[i]-prev_U[i], ord=1)>1e-3):
                cont=True
                break
        
        for i in range(prot_cnt):
            if(np.linalg.norm(V[i]-prev_V[i], ord=1)>1e-3):
                cont=True
                break
    return U,V

def curve_smoothening(vec,sim,neg_set,K2,cnt,r):
    for i in range(cnt):
        if i in neg_set:
            ser=pd.Series(sim[i])
            ser.sort_values(ascending=False, inplace=True)
            num=np.zeros(r)
            den=0
            nbr_cnt=0
            for j in ser.index:
                if j in neg_set:
                    continue
                num+=sim[i][j]*vec[j]
                den+=sim[i][j]
                nbr_cnt+=1
                if(nbr_cnt>=K2):
                    break
            num/=den
            for k in range(r):
                vec[i][k]=num[k]
    return vec

max_iter=1000

r=20
phi_l=np.zeros((rna_cnt,r))
phi_p=np.zeros((prot_cnt,r))

Lambda_l=0.5
Lambda_p=0.5

sd=1/math.sqrt(r)

U=np.random.normal(0,sd,(rna_cnt,r))
V=np.random.normal(0,sd,(prot_cnt,r))
# P=prob(U,V)
# print(P)

c=5
alpha=1
beta=1
gamma=0.01

neg_rna_set=[]
neg_prot_set=[]
n=np.random.randint(0,989)
print(n)
print(Y[n])
true_positives=[i for i in range(prot_cnt) if Y[n][i]==1]
true_negatives=[i for i in range(prot_cnt) if Y[n][i]==0]
print(true_positives)
neg_rna_set.append(n)
modified_Y=Y.copy()
for i in range(prot_cnt):
    modified_Y[n][i]=0

U,V=grad_des(U,V,modified_Y,c,Ll,Lp,Lambda_l,Lambda_p,alpha,beta,gamma,rna_cnt,prot_cnt,max_iter)

U=curve_smoothening(U,Sl,neg_rna_set,K2,rna_cnt,r)
V=curve_smoothening(V,Sp,neg_prot_set,K2,prot_cnt,r)


P=prob(U,V)
ser=pd.Series(P[n])
ser.sort_values(ascending=False, inplace=True)
for i in ser.index:
    print(i," ",ser[i])
pred=ser.to_numpy()
print(pred)

# plotting roc
thresholds = np.linspace(0,1,100)
tpr=[]
fpr=[]
pos=len(true_positives)
fal=len(true_negatives)
for th in thresholds:
    tp=0
    fp=0
    for i in range(64):
        if pred[i]<th:
            break
        if i in true_positives:
            tp+=1
        else:
            fp+=1
    tpr.append(tp/pos)
    fpr.append(fp/fal)
plt.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r',label='Random', alpha=.8)
plt.plot(fpr,tpr, label="ROC Curve",color="blue")
plt.text(0.5, 0.5, "varying threshold scores (0-1)", rotation=0, size=12,ha="center", va="center",bbox=dict(boxstyle="rarrow"))
plt.xlabel("False Positve Rate")
plt.ylabel("True Positive Rate")
plt.legend()
plt.show()

