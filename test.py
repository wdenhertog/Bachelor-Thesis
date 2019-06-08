import numpy as np
import matplotlib.pyplot as plt

x=[1,2]
indices = np.arange(2)
plt.bar(indices,x,color='r')
plt.xticks(indices,['FEM', 'EFR'])

plt.show()