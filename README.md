# Optimality and Consencus


"""
Created on Fri Nov 11 14:57:15 2016\
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
\
\
\
#preliminaries\
\
import numpy as np\
import math as m \
import matplotlib.pyplot as plt\
from matplotlib import cm\
from mpl_toolkits.mplot3d import Axes3D\
\
#establishes controls\
\
#number of subdivisions in the dicretized unit plane, M.\
M=20\
#width of each subdivision\
h=1.0/M\
#lol, pi\
pi=m.pi\
\
\
\
"""
here we choose a particualr f, but you can define anyhting you'd like, \
just make sure f:R2 -> R, [x[0],x[1]] -> f\
"""
\
def f(x):\
    return 2*(pi**2*m.sin(pi*x[0])*m.sin(pi*x[1]))  \ 
 \
 """\
 a simpler example:\
 \
 def f(x): return x[0] + x[1]\
 """
