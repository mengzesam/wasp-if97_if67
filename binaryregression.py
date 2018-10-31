#!/usr/bin/python
# -*- coding: UTF-8 -*-
import numpy as np
from scipy.optimize import leastsq
import pandas as pd

file='fitdata.csv'
data=np.array(pd.read_csv(file,sep='\t',header=None))
x=data[:,0:6]
y=data[:,6]
func=lambda coeff,x:np.dot(x,coeff)
errFunc=lambda coeff,x,y:func(coeff,x)-y
coeff=[1.0,1.0,1.0,1.0,1.0,1.0]
ans,success=leastsq(errFunc,coeff,args=(x,y))
print(ans[0],ans[1],ans[2],ans[3],ans[4],ans[5])

