from Bio import SeqIO
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Bio import Align
from Bio.Seq import Seq
import swalign
match=2
mismatch=-1
scoring = swalign.NucleotideScoringMatrix(match, mismatch)
sw = swalign.LocalAlignment(scoring)

for tem in SeqIO.parse('template.fasta','fasta'):
    template=tem.seq
    templatename=tem.id
filelist=os.listdir()
data=np.zeros((len(filelist),10000))
data=pd.DataFrame(data)

proteinid=np.zeros((len(filelist),10000))
proteinid=pd.DataFrame(proteinid)
proteinid.at[0,:]='a'

for i in range(len(filelist)):
    name=filelist[i]
    m=0
    for seq_record in SeqIO.parse(name,'fasta'):
        seq=str(seq_record.seq)
        seq=seq.replace('*','')
        seq=Seq(seq)
        alig=sw.align(template,seq)
        data.iat[i,m]=alig.score
        proteinid.iat[i,m]=seq_record.id
        m=m+1
    fig,ax=plt.subplots()
    ax.plot(range(10000),data.iloc[i,:],color='g',linewidth=4)
    ax.set_ylim(10,640)
    plt.savefig('C:\\Users\\Administrator\\Desktop\\kpc_protein\\'+str(name)+'.png',dpi=500)
data.to_csv('C:\\Users\\Administrator\\Desktop\\kpc_protein\\kpc_score.csv')
