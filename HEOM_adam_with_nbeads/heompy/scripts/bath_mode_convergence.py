#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate

eta = 1
gamma = 1
m = 1

x_min = 0
x_max = 20

dx = 0.01

bins = 200
n_b = 1000000
koef = 2.3
#******************************************************************888
n_max = 1000
nbs = np.arange(1, n_max+1, 1)

omgs = np.zeros([n_max])
omgs2 = np.zeros([n_max])
c = np.zeros([n_max])
c2 = np.zeros([n_max])
for i in range(n_max):
    omgs[i] = gamma*np.sqrt((nbs[i]/0.5 -1))
    c[i] = np.sqrt((omgs[i]**2+gamma**2)*m*omgs[i]*eta*gamma/(np.pi*nbs[i]))
    c[i] = (np.pi/2)*((c[i]**2)/(m*omgs[i]))
    omgs2[i] = -gamma*np.log(0.5/nbs[i])
    c2[i] = omgs2[i]*np.sqrt(2*eta*m*gamma/(np.pi*nbs[i]))
    c2[i] = (np.pi/2)*((c2[i]**2)/(m*omgs2[i]))

plt.title("Frequencies")
plt.plot(nbs, omgs, label="Debye")
plt.plot(nbs, omgs2, label="Ohm")
plt.legend()
plt.show()
plt.clf()

plt.title("Coefficients")
plt.plot(nbs, c, label="Debye")
plt.plot(nbs, c2, label="Ohm")
plt.legend()
plt.show()
plt.clf()
#*********************************************************************

omega = np.zeros([n_b])
cs = np.zeros([n_b])

omega2 = np.zeros([n_b])
cs2 = np.zeros([n_b])

def J(omega):
    J = eta*gamma*omega/(omega**2 + gamma**2)
    return J
def J2(omega):
    J2 = eta*gamma*np.exp(-omega/gamma)
    return J2

xs = np.arange(x_min, x_max, dx)
ys = J(xs)
ys2 = J2(xs)

for i in range(n_b):
    j = i+1
    omega[i] = gamma*np.sqrt((n_b/(j-0.5)) - 1)
    cs[i] = np.sqrt((omega[i]**2+gamma**2)*m*omega[i]*eta*gamma/(np.pi*n_b))
    cs[i] = (np.pi/2)*((cs[i]**2)/(m*omega[i]))
    omega2[i] = -gamma*np.log((j-0.5)/n_b)
    cs2[i] = omega2[i]*np.sqrt(2*eta*m*gamma/(np.pi*n_b))
    cs2[i] = (np.pi/2)*((cs2[i]**2)/(m*omega2[i]))

hist, bin_edges = np.histogram(omega, weights=cs, bins=100, density=True, range=(0, x_max))
hist2, bin_edges2 = np.histogram(omega2, weights=cs2, bins=100, density=True, range=(0, x_max))

bins = np.zeros([len(bin_edges)-1])

for i in range(len(bins)):
    bins[i] = 0.5*(bin_edges[i+1]+bin_edges[i])

bins2 = np.zeros([len(bin_edges2)-1])

for i in range(len(bins2)):
    bins2[i] = 0.5*(bin_edges2[i+1]+bin_edges2[i])

koef2 = dx*sum(ys2)

koef = dx*sum(ys)

#result = integrate.quad(lambda x: special.jv(2.5,x), 0, 4.5)

ax = plt.axes(xlim=(x_min,x_max))
plt.plot(xs, ys, label="Debye bath")
hist = koef*hist
plt.plot(bins,hist, label="histogram")
plt.legend()
plt.show()
plt.clf()

ax = plt.axes(xlim=(x_min,x_max))
plt.plot(xs, ys2, label="Ohmic bath")
hist2 = koef2*hist2
plt.plot(bins2,hist2, label="histogram")
plt.legend()
plt.show()
plt.clf()

print(np.amin(omega), "minimum Debye omega")
print(np.amax(omega), "maximum Debye omega")

print(np.amin(omega2), "minimum Ohm omega")
print(np.amax(omega2), "maximum Ohm omega")
plt.plot(omega, cs, "ro", label="Debye bath")
plt.plot(omega2, cs2, "bo", label="Ohmic bath")
plt.legend()
plt.show()


