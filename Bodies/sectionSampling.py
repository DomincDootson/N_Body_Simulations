import matplotlib.pyplot as plt

# we assume that vc =1 for this code

from math import *
VC = 1

def potential(r):
	return log(r)*VC

def energy(vr, angMom, radius):
	return 0.5*(vr**2 + (angMom**2)/(radius**2)) + potential(radius)

def jacobiIntegral(energy, angMom, patternSpeed):
	return energy - angMom*patternSpeed

def circularOrbitParams(radius):
	angMom = radius * VC
	return angMom, energy(0, angMom, radius) 



#First find the set of ICs for the circular orbit --> from this get the energy 
# take steps to decrease the angular momentum 


def vR(r, E, J):
	if (2*(E-potential(r)) - (J/r)**2 < 0): 
		print("Unphysical Value of vr")
	return sqrt(abs(2*(E-potential(r)) - (J/r)**2))

def radialMomentum(r, E, J):
	return (2 * (E - potential(r)) - (J/r)**2)

def rMin(E, J):
	x1, x2, nStep  = 0, exp(E-0.5), 20

	for i in range(nStep):
		holding = x1 + 0.5*(x2-x1)
		if radialMomentum(holding, E, J) <0:
			x1 = holding
		else:
			x2 = holding

	return x2 


def rMax(E, J):
	x1, x2, nStep  = exp(E-0.5), exp(E),  20

	for i in range(nStep):
		holding = x1 + 0.5*(x2-x1)
		if radialMomentum(holding, E, J)  > 0:
			x1 = holding
		else:
			x2 = holding

	return x1

def pow(x, y):
	x**y

def energyFuncRad(rApo, rPer):
	return (rApo*rApo*potential(rApo) -(rPer*rPer)*potential(rPer)) / (rApo*rApo - rPer*rPer)


def angMomFuncRad(rApo, rPer):
	return sqrt(2*(potential(rApo) - potential(rPer))/ (1/(rPer*rPer) - 1/(rApo*rApo)))


def lstEnergies(rCirc, patternSpeed, filename = "particleSamplesSections.out", jacobiFrac = 1.):
	jCirc, eCirc = circularOrbitParams(rCirc)
	jacobi = jacobiFrac*(eCirc - jCirc*patternSpeed)
	print("The jacobi Integral is: ", jacobi)


	step, nOrbits = 0.05, 24
	j = [jCirc-(i)*step for i in range(nOrbits)]

	f = open(filename, "w+")
	f.write(str(4*nOrbits))
	x, y = [], []
	for i in range(nOrbits):
		f.flush()
		vRadial = vR(rCirc, jacobi + j[i]*patternSpeed, j[i])

		f.write('\n' +"9.82208e-05" + " " + str(rCirc) + " " + str(0) + " " +  str(vRadial) + " " + str(j[i]/rCirc))
		f.write('\n' +"9.82208e-05" + " " + str(0) + " " + str(rCirc) + " " +   str(-j[i]/rCirc)+ " " + str(vRadial))
		f.write('\n' +"9.82208e-05" + " " + str(-rCirc) + " " + str(0) + " " +  str(-vRadial) + " " + str(-j[i]/rCirc))
		f.write('\n' +"9.82208e-05" + " " + str(0) + " " + str(-rCirc) + " " +   str(j[i]/rCirc)+ " " + str(-vRadial))

		e = energy(vRadial, j[i], rCirc)
		y.append(jacobiIntegral(e, j[i], patternSpeed))
		x.append(j[i])



	f.close()
	plt.plot(x, y) # Should be a vlat line 
	plt.xlabel("Ang Mom")
	plt.ylabel("Jacobi Integral")

	plt.show()


#lstEnergies(1/0.18, 0.18, "sormaniCR.out")
#lstEnergies((1/0.18), 0.18, "sormaniCR.out")
lstEnergies(0.3*(1/0.18), 0.18, "sormaniILR.out")

# lstEnergies((1/0.18) * (1-0.5*sqrt(2)), 0.18, "test.csv", jacobiFrac = 1  )
# lstEnergies((1/0.18), 0.18, "test.csv", jacobiFrac = 1  )