#!/usr/bin/env python
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from scipy.stats import shapiro

# Analysing uniformly distributed numbers
numbers = []
with open('uniform.dat','r') as file:
    for line in file:
        numbers.append(float(line))
std = np.std(numbers)
average = np.average(numbers)
mu = (max(numbers)+min(numbers))/2
stdev = (max(numbers)-min(numbers))/np.sqrt(12)
print("""The standard deviation of a uniform distribution is given by
stdev = width_of_interval/sqrt(12)
""")
print("statistic", "from data", '"analytic"', "error")
print("mean", average, mu, (mu-average)/mu)
print("std.dev.", std, stdev, (stdev-std)/stdev)
n, bins, patches = plt.hist(numbers, 100, normed=1, facecolor='green', alpha=0.75)
plt.show()
plt.clf()

# Analysing normally distributed numbers
numbers = []
with open('normal.dat','r') as file:
    for line in file:
        numbers.append(float(line))

print('average is ' + str(np.average(numbers)))
print('st dev is ' + str(np.std(numbers)))

# normality test
stat, p = shapiro(numbers)
print('Statistics=%.9f, p=%.9f' % (stat, p))

n, bins, patches = plt.hist(numbers, 100, normed=1, facecolor='green', alpha=0.75)
plt.show()
plt.clf()

# Analysing normally distributed numbers 2
numbers = []
with open('normal2.dat','r') as file:
    for line in file:
        numbers.append(float(line))

print('average is ' + str(np.average(numbers)))
print('st dev is ' + str(np.std(numbers)))

# normality test
stat, p = shapiro(numbers)
print('Statistics=%.9f, p=%.9f' % (stat, p))

n, bins, patches = plt.hist(numbers, 100, normed=1, facecolor='green', alpha=0.75)
plt.show()
