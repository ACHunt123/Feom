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


### Load all the parameters into the setup object 
sim = Setup(params)

# bath=sim.bath

support=np.arange(-1000,1000,10000)


# plt.plot(bath.J(support),support)

# plt.show()






