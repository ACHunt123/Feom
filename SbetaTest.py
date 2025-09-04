#!/usr/bin/env python3
'''
   +---------------------------------------+
   |   FEOM: Fortran heirarchical          |  
   |       Equations Of Motion             |
   |           By A. C. Hunt 2025          |
   +---------------------------------------+
'''
from Feom.setup import Setup
from Feom.parser import params
import matplotlib.pyplot as plt
import numpy as np
from Feom.baths.coth_decomp.cothAAA import get_coeffs

### Load all the parameters into the setup object 
sim = Setup(params)

bath=sim.bath
params=sim.params
params.mu=params.K

support=np.arange(-100,100,0.01)
fig,ax=plt.subplots()
ax.plot(support,bath.J(support))


get_coeffs(params, support, bath.J(support))
plt.show()






