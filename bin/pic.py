import numpy as np
import pandas as pd
import sys
import matplotlib.pyplot as plt
def smooth(y,box_pts):
    box=np.ones(box_pts)/box_pts
    y_smooth=np.convolve(y,box,mode = 'same')
    return y_smooth

file=sys.argv[1]
picout=sys.argv[2]
df=pd.read_csv(file,header=None,sep='\t')
df.columns=['chr','position','depth']
x=df['position']
y=df['depth']
plt.plot(x,smooth(y,100))
plt.title('smooth with box size 100')
plt.xlabel('genome position')
plt.ylabel('depth')

import time ,os
time1=time.strftime('%Y-%m-%d')
sv_path=picout+'/picture'
os.makedirs(sv_path,exist_ok=True)

plt.savefig(sv_path+'/depth.pdf')
plt.close()
