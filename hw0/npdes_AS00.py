#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 10 10:46:07 2017

@author: izabelaguiar

nPDES AS00
Homework 0, due 2017-09-13

Fork the class repository, clone, and create a directory hw0 inside the repository.
Add your source file(s) to that directory.

Write a function diffmat(x) that returns a matrix $D$ that computes first derivatives.

Write a function diff2mat(x) that returns a matrix $D_2$ that computes second derivatives.

Use test solutions to determine the order of accuracy of your methods for evenly and
 non-evenly spaced points. Which norm did you use?

Add README.md in the hw0 directory and summarize your results
 (one paragraph or a few bullet items is fine).

You may assume that the points x are monotonically increasing.

You'll need to think about what to do at the endpoints.

"""

import numpy as np
from matplotlib import pyplot as plt

def diffmat(x):
    """Compute first derivative matrix operator.
    
    Parameters
    ----------
    x : ndarray
        1-by-m array of x values for which the solution, u, is evaluated
    Returns
    -------
    D : ndarray
        m-by-m matrix that computes first derivatives for the system, u'=Du
    
    Notes
    -----
    Note that matrix used here to compute first derivatives is a second order accurate
    centered approximation. At the first and last points we implement a second order 
    accurate forward and backward approximation, respectively. 
    All approximations are derived in 
    
        [1] LeVeque, R. "Finite Difference Methods for Ordinary and Partial Differential
            Equations". SIAM (2007), Philadelphia, PA.
            
    The centered approximation is given by eq. (1.3) and the for/back(wards) by eq. (1.11)
    """
    m = len(x)
    
    h = np.ones(m-1)
    
    for i in range(m-1): 
        h[i] = x[i+1]-x[i]
        
    D = np.eye(m, k=1) - np.eye(m, k=-1)
    
    D[0, 0:3] = (1. / (h[0] + h[1])) * np.array([1., -4., 3. ]) #by 1.11
    
    #mth row
    D[-1, -3:] = (1. / (h[-1] + h[-2])) * np.array([1., -4. , 3. ]) #by 1.11
    
    #ith row
    for i in range(1, m-1):
        D[i, :] *= (1. / (h[i] + h[i-1]))
    #D *= (1. / (2*h) )
    return D




def diff2mat(x):
    """Compute first derivative matrix operator.
    
    Parameters
    ----------
    x : ndarray
        1-by-m array of x values for which the solution, u, is evaluated
    Returns
    -------
    D_2 : ndarray
        m-by-m matrix that computes first derivatives for the system, u'=Du
    
    Notes
    -----
    Note that matrix used here to compute second derivatives is a second order accurate
    centered approximation. At the first and last points we implement a second order 
    accurate forward and backward approximation, respectively. 
    All approximations are derived in 
    
        [1] LeVeque, R. "Finite Difference Methods for Ordinary and Partial Differential
            Equations". SIAM (2007), Philadelphia, PA.
            
    The centered approximation is given by eq. (1.14) and the for/back(wards) by eq. (1.11)
    """ 
    m = len(x)
    h = np.ones(m-1)

    for i in range(m-1): 
        h[i] = x[i+1]-x[i]
        
    D_2 = np.zeros((m,m))
    
    for i in range(1, m-1):
        c_1 = 2. / (h[i-1]*(h[i-1] + h[i]))
        c_2 = -2. / (h[i-1] * h[i])
        c_3 = 2. / (h[i]*(h[i-1] + h[i]))
        D_2[i, i-1] = c_1
        D_2[i, i] = c_2
        D_2[i, i+1] = c_3
    
    r = 1./(h[0]+h[1])
    a = 9./(h[0]+h[1])
    b = 12./(h[0]+h[1]) + 12./(h[1]+h[2])
    c = 3./(h[0]+h[1]) + 16./(h[1]+h[2]) + 3./(h[2]+h[3])
    d = 4./(h[1]+h[2]) + 4./(h[2]+h[3])
    e = 1./(h[2]+h[3])
    
    D_2[0, 0:5] = np.array([r*a, -r*b, r*c, -r*d, r*e])
    
    r = 1./(h[-1]+h[-2])
    a = 9./(h[-1]+h[-2])
    b = 12./(h[-1]+h[-2]) + 12./(h[-2]+h[-3])
    c = 3./(h[-1]+h[-2]) + 16./(h[-2]+h[-3]) + 3./(h[-3]+h[-4])
    d = 4./(h[-2]+h[-3]) + 4./(h[-3]+h[-4])
    e = 1./(h[-3]+h[-4])
    
    D_2[-1, -5:] = np.array([r*e, -r*d, r*c, -r*b, r*a])

    
    return D_2
    
""" Test function: u(x) = sin(x), u'(x) = cos(x), u''(x) = -sin(x) """


#Compute the errors 
ns = 2 ** np.arange(3, 9)
errors_1 = np.zeros(len(ns))
errors_2 = np.zeros(len(ns))   
h = np.zeros(len(ns)) 
i = 0
for n in ns:
    x = np.linspace(0, 2.*np.pi, n)
    h[i] = 1. / n
    du_analytic = np.cos(x)
    d2u_analytic = -np.sin(x)
    u = np.sin(x)
    D = diffmat(x)
    du_numerical = D.dot(u)
    D_2 = diff2mat(x)
    d2u_numerical = D_2.dot(u)
    errors_1[i] = np.linalg.norm(du_analytic - du_numerical, np.inf)
    errors_2[i] = np.linalg.norm(d2u_analytic - d2u_numerical, np.inf)
    i += 1
 
"""Test order of convergence for both methods"""
   
plt.loglog(h, errors_1, 'o', label='numerical')    
plt.loglog(h, h**2, label =r'$\mathcal{O}(h^2)$')  
plt.legend(loc='best')  
plt.xlabel('h')
plt.ylabel('Error')
plt.title('Cost vs. Accuracy, First Derivative')


plt.loglog(h, errors_2, 'o', label='numerical')    
plt.loglog(h, h**2, label =r'$\mathcal{O}(h^2)$')  
plt.legend(loc='best')  
plt.xlabel('h')
plt.ylabel('Error')
plt.title('Cost vs. Accuracy, Second Derivative')   


""" Plot the numerical vs. analytic solutions for even h """ 

n=2**4
x = np.linspace(0, 2.*np.pi, n)
du_analytic = np.cos(x)
d2u_analytic = -np.sin(x)
u = np.sin(x)
D = diffmat(x)
du_numerical = D.dot(u)
D_2 = diff2mat(x)
d2u_numerical = D_2.dot(u)
    
plt.plot(x, du_analytic, label= '$\cos(x)$')
plt.plot(x, du_numerical, 'o', label= '$Du$')
plt.plot(x, d2u_analytic, label='$-\sin(x)$')
plt.plot(x, d2u_numerical, 'o', label = r'$D_2u$')
plt.legend(loc='best')
plt.title('Numerical vs. Analytic Solution, even h')    
    
""" Plot the numerical vs. analytic solutions for uneven h"""
 
uneven_h = (1./5.)*np.ones(12)
uneven_h[0] *= 1/4
uneven_h[1] *= 1/4
uneven_h[2] *= 1/4
uneven_h[3] *= 1/4

uneven_h[5] *=1/2
uneven_h[4] *=3/4
uneven_h[6] *= 2/3

uneven_h[8] *= 1/5
uneven_h[9] *= 1/5
uneven_h[10] *= 1/5
uneven_h[11] *= 1/5

x = np.zeros(13)
for i in range(0, 12):
    x[i+1] = x[i] + uneven_h[i]

du_analytic = np.cos(x)
d2u_analytic = -np.sin(x)
u = np.sin(x)
D = diffmat(x)
du_numerical = D.dot(u)
D_2 = diff2mat(x)
d2u_numerical = D_2.dot(u)
    
plt.plot(x, du_analytic, label= '$\cos(x)$')
plt.plot(x, du_numerical, 'o', label= '$Du$')
plt.plot(x, d2u_analytic, label='$-\sin(x)$')
plt.plot(x, d2u_numerical, 'o', label = r'$D_2u$')
plt.legend(loc='best')
plt.title('Numerical vs. Analytic Solution, uneven h') 
  
    