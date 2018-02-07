import numpy as np
import matplotlib.pyplot as plt
import math
from pylab import *

class FDTD:

	def __init__(self):
		#Omega = ome
		#OmegaP = omep
		#Aj = A0,A1,...
		#vj = v0,v1,...
		self.omep = 9.03
		self.A = np.array([0.760, 0.024, 0.010, 0.071, 0.601, 4.384])
		self.omej = np.array([0.0, 0.415, 0.830, 2.969, 4.304, 13.32])
		self.v = np.array([0.053, 0.241, 0.345, 0.870, 2.494, 2.214])

	def showFig(self):
		#create x
		X = np.arange(300, 850, 50)
		E = self.calcuEpsilon(X)
		print E

		#plot y
		plt.plot(X, E.real)
		plt.plot(X, E.imag)
		plt.show()

	#Culculate Epsilon
	def calcuEpsilon(self, lam):
		ome = (4.135E-15*2.99E+8) / (lam*1.0E-9)
		#print ome
		epsilon = 1 - self.A[0]*self.omep*self.omep/(ome*(ome+self.v[0]*1.0j)) + self.calcuSigmaXi(5, ome)

		return epsilon

	#Calculate SigmaXi
	def calcuSigmaXi(self, j, ome):
		xi = 0
		for i in range(1, j+1):
			xi += self.A[i]*self.omep*self.omep / ((self.omej[i]*self.omej[i]-ome*ome)-ome*self.v[i]*1.0j)

		return xi



if __name__ == '__main__':
	print "==== Start ===="
	fdtd = FDTD()
	fdtd.showFig()
