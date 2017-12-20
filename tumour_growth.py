
"""
tumour_growth.py
final project for numerical solutions for partial differential equations. 

@author: izabelaguiar
"""

from fenics import *
from dolfin import *
import numpy as np

#Initial Conditions
class InitialConditions(Expression):
    def eval(self, values, x):
        if x<=0.25 and x>=0:
            values[0] = np.exp(-x**2/.01)
        else:
            values[0] = 0
        values[1] = 1 - 0.5*values[0]
        values[2] = 0.5*values[0]
    def value_shape(self):
        return (3,)
    
class TumourGrowth(NonlinearProblem):
    def __init__(self, a, L):
        NonlinearProblem.__init__(self)
        self.L = L
        self.a = a
    def F(self, b, x):
        assemble(self.L, tensor=b)
    def J(self, A, x):
        assemble(self.a, tensor=A)
        
def boundary(x, on_boundary):
    return on_boundary

theta = 0.5
T = 10.0            # final time
num_steps = 10000    # number of time steps
dt = T / num_steps # time step size
d_n = 0.001         # tumor random motility coefficient
gamma = 0.005      # cancer cell proliferation rate
eta = 10           # degradation rate
d_m = 0.001         # MDE diffusion coefficient
alpha = 0.1        # MDE production rate
beta = 0.0         # natural MDE decay rate

# Form compiler options
parameters["form_compiler"]["optimize"]     = True
parameters["form_compiler"]["cpp_optimize"] = True

# Create mesh and define function spaces
mesh = IntervalMesh(96, 0, 1)

# Define function space for system of concentrations
P1 = FiniteElement('P', interval, 1)
element = MixedElement([P1, P1, P1])
V = FunctionSpace(mesh, element)

# Define test functions
v_1, v_2, v_3 = TestFunctions(V)
du = TrialFunction(V)
xi = FacetNormal(mesh)

# Define functions for state
u = Function(V)
u_n = Function(V)
u_0 = Function(V)

# Split system functions to access components
dn, df, dm = split(du)
n, f, m = split(u)
n_n, f_n, m_n = split(u_n)

# Enforce initial conditions
u.interpolate(InitialConditions(degree=5))
u_n.interpolate(InitialConditions(degree=5))
u_0.interpolate(InitialConditions(degree=5))

# Define Dirichlet Boundary Conditions
bc = DirichletBC(V, InitialConditions(degree=5), boundary)


dx = Measure('dx', domain=mesh)
ds = Measure('ds', domain=mesh)

# Define variational problem
F_1 = -d_n*grad(n) + n*gamma*grad(f)
F_2 = -d_m*grad(m)

L0 = ((n - n_n) / dt)*v_1*dx + dot(grad(v_1), F_1)*dx \
    
L1 = ((f - f_n) / dt)*v_2*dx + eta*m*f*v_2*dx 

L2 = ((m - m_n) / dt)*v_3*dx + dot(grad(v_3), F_2)*dx - alpha*n*v_3*dx + beta*m*v_3*dx

artdif = alpha * dx

L = L0 + L1 + L2
    
a = derivative(L, u, du)

#Define parameters for solving
problem = TumourGrowth(a, L)
solver = NewtonSolver()
solver.parameters["linear_solver"] = "lu"
solver.parameters["convergence_criterion"] = "incremental"
solver.parameters["relative_tolerance"] = 1e-6
solver.parameters["maximum_iterations"] = 50

t = 0.0
while (t<T):
    #print('over here!', t)
    t += dt
    u_n.vector()[:] = u.vector()
    A = assemble(a)
    b = assemble(L)
    bc.apply(A, b)
    solve(A, u.vector(), b)

