import pandas as pd
from itertools import combinations
from sklearn.utils import shuffle

df=pd.read_csv('C:/Users/blues/Desktop/PPI/All_CSV/physical.txt', delimiter=' ')
df_filtered = df[(df['experimental']>0) | (df['database']>0)]
df_filtered=df_filtered.reset_index(drop=True)

pos=df_filtered.loc[:, 'protein1':'protein2']
pos['protein1']=pos['protein1'].replace(r'\D', r'', regex=True)
pos['protein2']=pos['protein2'].replace(r'\D', r'', regex=True)
quer=0
for index, row in pos.iterrows():
    try:
        if row['protein1']>row['protein2']:
            quer=row['protein1']
            row['protein1']=row['protein2']
            row['protein2']=quer
    except:
        print(index)

pos=pos.sort_values(by=['protein1', 'protein2'])
pos=pos.drop_duplicates(keep='first').reset_index(drop=True)

nonrpt_prot=[]
cnt=0
for i in range(0, len(df_filtered['protein1'])):
    if df_filtered['protein1'][i] not in nonrpt_prot:
        nonrpt_prot.append(df_filtered['protein1'][i])
        cnt+=1

def combines(pro, r):
    return list(combinations(pro, r))

r=2
raw_list=combines(nonrpt_prot, r)
rd_smp=pd.DataFrame(raw_list, columns=['protein1', 'protein2'])

rd_smp['protein1']=rd_smp['protein1'].replace(r'\D', r'', regex=True)
rd_smp['protein2']=rd_smp['protein2'].replace(r'\D', r'', regex=True)

neg = pd.concat([pos, rd_smp])
neg = neg.drop_duplicates(keep=False)

pos=shuffle(pos).reset_index(drop=True)
for index, row in pos.iterrows():
    try:
        new1=row['protein1'].find('00000')
        row['protein1']=row['protein1'][:new1]+'.ENSP'+row['protein1'][new1:]
        new2=row['protein2'].find('00000')
        row['protein2']=row['protein2'][:new2]+'.ENSP'+row['protein2'][new2:]
    except:
        print(index)
        pos=pos.drop([index])
pos.insert(2, column="interact", value=1)

neg=shuffle(neg).reset_index(drop=True)
neg=neg.loc[0:937280, :]
for index, row in neg.iterrows():
    try:
        new1=row['protein1'].find('00000')
        row['protein1']=row['protein1'][:new1]+'.ENSP'+row['protein1'][new1:]
        new2=row['protein2'].find('00000')
        row['protein2']=row['protein2'][:new2]+'.ENSP'+row['protein2'][new2:]
    except:
        print(index)
        neg=neg.drop([index])
neg.insert(2, column="interact", value=0)

all_prot=pd.concat([pos, neg])
all_prot=shuffle(all_prot).reset_index(drop=True)
all_prot=all_prot.drop_duplicates(keep=False)
all_prot.to_csv('C:/Users/blues/Desktop/PPI/All_CSV/all_prot.csv')