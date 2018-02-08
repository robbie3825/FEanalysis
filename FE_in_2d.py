# -*- coding: utf-8 -*-
"""
Created on Fri Nov 11 14:57:15 2016
Author: Robert Cunningham



FINITE ELEMENT ANALYSIS, MODEL PROBLEM IN THE PLANE, POISSON'S EQ:

The following program determines the solution to Poisson's equation, 
-div(grad(u))=f on (0,1)x(0,1) with u=0 on the boundary thereof.
 
M is the number of blocks in each coordinate direction that make up the 
discretization of the unit square, and a triangulation is determined by the
addition of edges with positive slope inside each of the M^2 blocks.  The 
Barycentric quadrature rule is used to compute the load vector, b (cf. 
Galerkin's Method) and the element stiffness martrices of triangles type I and 
II are determined by a bottom left corner - ccw orientation and then used 
for the element-wise global stiffness matrix assemblage.  The triangles are 
enumerated from the bottom left and then the top right of the square. First
moving upward, with lower triangles being listed first, then the uppers.
The interior points of the grid are labeled similarly and all boundary nodes 
are given a "value" of -1 as they are not needed in the composition of the 
system; Au=b, which is a complete linear system of equations with (M-1)^2 
unkowns (values of the desired solution, u, at the intereior nodes of our 
domain). 

Running the program opens three figures, Figure 1 is the approximate
solution plotted against the square, Figure 2 is the exact solution which, 
assuming you have one, you must add into this program on line 354 in order for 
Figure 3 to make sense, which is a plot of the error between our approximate 
solution and the exact, analytical one.

All examples are given in green print under the assumption that M=3

"""


import numpy as np
import math as m 
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D



#establishes controls

M=20


h=1.0/M
pi=m.pi
def f(x):
    return 2*(pi**2*m.sin(pi*x[0])*m.sin(pi*x[1]))   
    


#equips vector arithmetic:

#addition

def add(a,b): 
    return [sum(x) for x in zip(a,b)]
    
"""
add([1,2,3],[4,5,6])
Out[78]: [5, 7, 9]

"""

#scalar multiplication
 
def scam(r,v):
    return [r*v[i] for i in range(len(v))]
    
"""
scam(3,[1,2,3])
Out[77]: [3, 6, 9]

"""



#creates a nodal mesh to extract indicial information about the grid

G=np.full((M+1,M+1),-1.)
for i in range(1,M):
    for j in range(1,M):
        G[j,i]=i+(j-1)*(M-1)
G=G.astype(int)

"""
G
Out[53]: 
array([[-1, -1, -1, -1],
       [-1,  1,  2, -1],
       [-1,  3,  4, -1],
       [-1, -1, -1, -1]])
       
"""



#builds a list of nodes from G

N=[]
for i in range(M):
    for j in range(M):
        N.append([G[i,j],G[i+1,j],G[i+1,j+1]])
for i in range(M):
    for j in range(M):
        N.append([G[i,j],G[i+1,j+1],G[i,j+1]])
             
"""
N
Out[59]: 
[[-1, -1, 1],
 [-1, 1, 2], ******* this is N[1] in the example below
 [-1, 2, -1],
 [-1, -1, 3],
 [1, 3, 4],
 [2, 4, -1],
 [-1, -1, -1],
 [3, -1, -1],
 [4, -1, -1],
 [-1, 1, -1],
 [-1, 2, -1],
 [-1, -1, -1],
 [-1, 3, 1],
 [1, 4, 2], ******** this is N[13]
 [2, -1, -1],
 [-1, -1, 3],
 [3, -1, 4],
 [4, -1, -1]]
 
"""



#determines which entries in some member of the above list are > -1

def l(x):
    L=[]
    for i in range(3):
        if x[i]>0:
            L.append(i)
    return L
    
"""
l(N[1])
Out[61]: [1, 2]

l(N[13])
Out[62]: [0, 1, 2]

"""



#creates a map between the positions and placement of the members of N above 
  
def map(k,x):
    a=N[k]
    b=l(N[k])
    if x in b: return a[x]-1
    else: return []
    
"""
[map(1,x) for x in l(N[1])]
Out[63]: [0, 1]

[map(13,x) for x in l(N[13])]
Out[64]: [0, 3, 1]

"""



#element stiffness matrices for lower and upper triangles, resp.

Ae=np.array([[.5,-.5,0],[-.5,1,-.5],[0,-.5,.5]])

Ae2=np.array([[.5,0,-.5],[0,.5,-.5],[-.5,-.5,1]])    

"""
Ae
Out[65]: 
array([[ 0.5, -0.5,  0. ],
       [-0.5,  1. , -0.5],
       [ 0. , -0.5,  0.5]])

Ae2
Out[66]: 
array([[ 0.5,  0. , -0.5],
       [ 0. ,  0.5, -0.5],
       [-0.5, -0.5,  1. ]])
       
"""



#builds the local element stifness matrices for each K in Th

def alpha(k):
    B=np.zeros(((M-1)**2,(M-1)**2))
    if k in range(M**2):
        for i in range(3):
            for j in range(3):
                
                B[map(k,i),map(k,j)]=Ae[i,j]
        return B
    elif k in range(M**2,2*M**2):
        for i in range(3):
            for j in range(3):
                
                B[map(k,i),map(k,j)]=Ae2[i,j]
        return B

"""
alpha(1)
Out[68]: 
array([[ 1. , -0.5,  0. ,  0. ],
       [-0.5,  0.5,  0. ,  0. ],
       [ 0. ,  0. ,  0. ,  0. ],
       [ 0. ,  0. ,  0. ,  0. ]])

alpha(13)
Out[69]: 
array([[ 0.5, -0.5,  0. ,  0. ],
       [-0.5,  1. ,  0. , -0.5],
       [ 0. ,  0. ,  0. ,  0. ],
       [ 0. , -0.5,  0. ,  0.5]])
       
"""        

#Assemblage (uneventful)

A=sum([alpha(i) for i in range(2*M**2)])



"""
A
Out[70]: 
array([[ 4., -1., -1.,  0.],
       [-1.,  4.,  0., -1.],
       [-1.,  0.,  4., -1.],
       [ 0., -1., -1.,  4.]])
       
"""



#builds a list of barycenters

a=round(h/3,5)
b=round(2*h/3,5)
c=[]
j=0
i=0
while j<M:
    for i in range(M):
        c.append([b+j*h,a+i*h])
        c.append([a+j*h,b+i*h])
    j+=1

C=[]
for i in range(M**2):
    C.append(c[2*i])
for j in range(M**2):
    C.append(c[2*j+1])

C=np.around(C,decimals=5)

"""
C
Out[71]: 
array([[ 0.22222,  0.11111],
       [ 0.22222,  0.44444],
       [ 0.22222,  0.77778],
       [ 0.55555,  0.11111],
       [ 0.55555,  0.44444],
       [ 0.55555,  0.77778],
       [ 0.88889,  0.11111],
       [ 0.88889,  0.44444],
       [ 0.88889,  0.77778],
       [ 0.11111,  0.22222],
       [ 0.11111,  0.55555],
       [ 0.11111,  0.88889],
       [ 0.44444,  0.22222],
       [ 0.44444,  0.55555],
       [ 0.44444,  0.88889],
       [ 0.77778,  0.22222],
       [ 0.77778,  0.55555],
       [ 0.77778,  0.88889]])
       
"""



#beta returns element load vectors

def beta(k):
    B=np.zeros((M-1)**2)
    for i in l(N[k]):
        B[N[k][i]-1]=f(C[k])
    return B

"""
beta(1)
Out[72]: array([ 0.85287674,  0.85287674,  0.        ,  0.        ])

beta(13)
Out[73]: array([ 0.63305595,  0.63305595,  0.        ,  0.63305595])

"""



#assembles b, the load vector

b=sum([beta(k) for k in range(2*M**2)])
b=scam((h**2)/6,b)

"""
b
Out[5]: 
[1.3325462291748869,
 1.2145798973376134,
 1.2145798973376134,
 1.3325487927477415]

"""



#Solves

Uint=np.linalg.solve(A,b)



#Reshapes

Uint=Uint.reshape((M-1,M-1))
U=np.zeros((M+1,M+1))
i=1
while i<M:
    U[i,1:M]=Uint[(i-1),:]
    i+=1

"""
Uint
Out[74]: 
array([[ 0.0104992 ,  0.00646052],
       [ 0.00801982,  0.00362009]])

U
Out[75]: 
array([[ 0.        ,  0.        ,  0.        ,  0.        ],
       [ 0.        ,  0.0104992 ,  0.00646052,  0.        ],
       [ 0.        ,  0.00801982,  0.00362009,  0.        ],
       [ 0.        ,  0.        ,  0.        ,  0.        ]])

"""


#plots the approximate solution

x=np.arange(0,1+h,h)
y=np.arange(0,1+h,h)
fig = plt.figure()
ax = fig.gca(projection='3d')
X,Y = np.meshgrid(x,y,indexing='ij')
surf = ax.plot_surface(X, Y, U, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
plt.xlabel('x')
plt.ylabel('y')
plt.show()


#gets a matrix of the exact solution on the (M+1)^2 grid points that make up 
#the closure of our domain, must be entered in order to utilize error analysis

x=np.arange(0,1+h,h)
y=np.arange(0,1+h,h)

u = np.zeros((M+1,M+1))
for i in range(M+1):
    for j in range(M+1):
        u[i,j] = m.sin(pi*x[i])*m.sin(pi*y[j])

      
#plots the exact solution
        
fig = plt.figure() 
ax = fig.gca(projection='3d')
X,Y = np.meshgrid(x,y,indexing='ij')   
surf = ax.plot_wireframe(Y,X, u.transpose(), rstride=1, cstride=1)
plt.xlabel('x')
plt.ylabel('y')
plt.show()



#gets the absolute value of the error

E=np.zeros((M+1,M+1))
for i in range(M+1):
    for j in range(M+1):
        E[i,j]=abs(U[i,j]-u[j,i])


#plots the absolute error over the mesh

fig = plt.figure()
ax = fig.gca(projection='3d')
X,Y = np.meshgrid(y,x,indexing='ij')
surf = ax.plot_surface(X, Y, E.transpose(), rstride=1, cmap=cm.bwr, cstride=1, linewidth=0, antialiased=False)
plt.xlabel('x')
plt.ylabel('y')
plt.show()


#computes and prints the maximum error

z=np.amax(E)
print('The maximum error along the mesh grid is')
print(z)

"""
runfile('C:/Users/Coffee/Desktop/FE in 2d.py', wdir='C:/Users/Coffee/Desktop')
The maximum error along the mesh grid is
0.123048782394

"""


















"""
Error Analysis

"""

#computes the approximate value of the L2 norm of E=|u-U|, e, by using a 2-
#dimensional trapezoidal rule for the integral over the interior points of
#the mesh grid.  An expression for this e is given by the square root of 
#quantity h^2/2 times quantity the sum over (M-1)^2 interior points P, of 
#|U-u|^2 evaluated at P. Since for each interior point P at (i*h,j*h) is 
#found to be the i,jth entry of E, we square E element wise and sum before
#taking this quatitiy's square root, times the square root of h^2/2.

E2=[]
for i in range(M+1):
    for j in range(M+1):
        E2.append(E[i,j]**2)
e=h*m.sqrt(.5*sum(E2))

#below is the values of this e for M's ewual to 4, 8, 16 and 32

"""
print(M)
print(e)

runfile('C:/Users/Coffee/Desktop/FE in 2d.py', wdir='C:/Users/Coffee/Desktop')
4
0.030342347652035372

runfile('C:/Users/Coffee/Desktop/FE in 2d.py', wdir='C:/Users/Coffee/Desktop')
8
0.0076384454369269475

runfile('C:/Users/Coffee/Desktop/FE in 2d.py', wdir='C:/Users/Coffee/Desktop')
16
0.001912828403938994

runfile('C:/Users/Coffee/Desktop/FE in 2d.py', wdir='C:/Users/Coffee/Desktop')
32
0.00047840217781461064

"""

#below we construct a table of log base 2's of the ratios of error estimate's
#e, the last row of which serves as an approximate value for the convergence
#rate of our solution method carried out in the above work.

a=[1./4,1./8,1./16,1./32]
b=[0.030342347652035372,0.0076384454369269475,0.001912828403938994,
   0.00047840217781461064]
c=['NA',b[0]/b[1],b[1]/b[2],b[2]/b[3]]
d=['NA']
d+=[m.log(c[i],2) for i in [1,2,3]]

table=np.array([a,b,c,d])
print(table)

"""
runfile('C:/Users/Coffee/Desktop/FE in 2d.py', wdir='C:/Users/Coffee/Desktop')
[['0.25' '0.125' '0.0625' '0.03125']
 ['0.030342347652035372' '0.0076384454369269475' '0.001912828403938994'
  '0.00047840217781461064']
 ['NA' '3.9723197478572967' '3.9932726956571067' '3.998368930252339']
 ['NA' '1.9899817558173631' '1.9975715951963318' '1.9994115959667695']]
 
 """
 
#the above table suggests that our method converges at the correct rate, 
#2, and we conclude that our code conatins no bugs, and that the Finite 
#Element method in this model problem agrees with the canonical Galerkin's
#method found in the literature.


 
