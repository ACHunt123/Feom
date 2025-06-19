import matplotlib.pyplot as plt
import numpy as np

# Data for plotting
filename = "output.txt"
data = np.loadtxt(filename)
t=data[:,0]
ct=data[:,1]

fig, axes = plt.subplots(1, len(data[1,:])-2, figsize=(8, 6))
# get the headers from the first row
with open(filename, 'r') as f:
    header = f.readline().strip().split('\t')  # Split the header by tab
print(header)
for i,ax in enumerate(axes):
    ax.plot(t, data[:,i+1])
    ax.set_xlabel('Time (s)')
    ax.set_ylabel(f'${header[i+1]}$')

plt.tight_layout()
plt.show()  