#!/usr/bin/python

# calculates Sigma and Pi from

# http://inspirehep.net/record/1703751

# using symbolic algebra

# Wed Feb 12 15:14:12 CST 2020 Jeff and Mark


from sympy import assoc_legendre
from sympy import cos,sin,Abs,factorial
from sympy import sqrt
from sympy.abc import x, m, n
from sympy import trigsimp
from sympy import Derivative

import sympy as sym

costheta=sym.Symbol('costheta')
z=sym.Symbol('z')
y=sym.Symbol('y')
x=sym.Symbol('x')

theta=sym.Symbol('theta',real=True,positive=True)
phi=sym.Symbol('phi',real=True,positive=True)
r=sym.Symbol('r',real=True,positive=True)

def c(ell,m):
    if(m>=0):
        return (factorial(ell-1)*(-2)**Abs(m))/factorial(ell+Abs(m))*cos(m*phi)
    else:
        return (factorial(ell-1)*(-2)**Abs(m))/factorial(ell+Abs(m))*sin(Abs(m)*phi)


ell=3
m=0

Sigma=c(ell,m)*r**ell*assoc_legendre(ell,Abs(m),cos(theta))
Sigma=sym.simplify(Sigma)
Sigma=Sigma.subs({Abs(sin(theta)):sin(theta)})
Sigma=Sigma.expand(trig=True)
print("Sigma in spherical coordinates is %s"%Sigma)
Sigma=Sigma.subs({r*cos(theta):z,
                  r*sin(theta)*sin(phi):y,
                  r*sin(theta)*cos(phi):x})
Sigma=Sigma.subs({r**2*sin(theta)**2:x**2+y**2})
Sigma=Sigma.subs({r**2:x**2+y**2+z**2})
Sigma=sym.simplify(Sigma.expand())
print("Sigma in cartesian coordinates is %s"%Sigma)

Pix=Derivative(Sigma,x)
Piy=Derivative(Sigma,y)
Piz=Derivative(Sigma,z)

Pix=Pix.doit()
Piy=Piy.doit()
Piz=Piz.doit()

Pix=sym.simplify(Pix.expand())
Piy=sym.simplify(Piy.expand())
Piz=sym.simplify(Piz.expand())

print(ell-1,m)
print("Pix is %s"%Pix)
print("Piy is %s"%Piy)
print("Piz is %s"%Piz)

