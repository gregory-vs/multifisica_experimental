import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
import math 

# u(x) = x^2*e^x+cos(3x)

# --- input data --- 
# boundary conditions
q = 1.728289
h = 0

# source f(x)
def f(x): 
    #f = 5.0 # Problema original da atividade avaliativa 03
    f = (x**2 + 4*x + 2)*math.e**x - 9*math.cos(3*x) # Um problema um pouco mais interessante 
    # teste aqui qualquer função f
    return f
# A solução analítica depende da integral da função f(x) no intervalo 0-y. 
# Essa integral poderia ser feita usando o pacote de integração simbólica do python.
# Entretanto, nesse exemplo ela é feita analiticamente por meio da função int_f
def int_f(z): # to avoid symbolic integation of f over z 
    #f = 5.0*z
    f = (z**2 + 2*z)*math.e**z - 3*math.sin(3*z)
    # para qualquer nova função f, é necessário criar a int_0^y f(z)dz
    return f

# degrees of freedom
dof = 7
samples = 201 # used to plot the figures

# --- Solution ---
# sample points
x = np.linspace(0, 1, samples)
u = np.zeros(np.size(x))

j = 0
for i in x:
    u[j] = q + (1-i)*h + quad(int_f,i,1)[0] # ,args=(f,i))[0]
    j = j+1
    
fig = plt.figure(figsize=(8,5))
plt.plot(x, u, 'k', linewidth=2);
plt.grid(True);
plt.xlabel('$x$');
plt.ylabel('$u$');
fig.suptitle('Fig. 1 - Exact solutions using quadratic integration');

# 1D mesh
xi = np.linspace(0, 1, dof+1) # do not necessarily have to a linear spaced vector
np.set_printoptions(precision=3)
print("\n")
print("Nós usados na discretização:")
print(xi)
print("\n")

# first order shape function
def N(xp,xi,i):
    last_i = len(xi)-1 
    if(i == 0 and xp <= xi[1]):
        y = (xi[1]-xp)/(xi[1]-xi[0])
    elif(i == last_i and xp >= xi[last_i-1]):
        y = (xp-xi[last_i-1])/(xi[last_i]-xi[last_i-1])
    elif(xp >= xi[i] and xp <= xi[i+1]):
        y = (xi[i+1]-xp)/(xi[i+1]-xi[i])
    elif(xp <= xi[i] and xp >= xi[i-1]):
        y = (xp-xi[i-1])/(xi[i]-xi[i-1])
    else:
        y = 0
    return y

# ploting shape functions
fig = plt.figure(figsize=(8,5))
for j in range(len(xi)):
    Ni = np.zeros(np.size(x))
    for i in range(len(x)):
        Ni[i] = N(x[i],xi,j)

    plt.plot(x, Ni, linewidth=2);
plt.grid(True);
plt.xlabel('$x$');
plt.ylabel('$Ni$');
fig.suptitle('Fig. 2 - First order shape functions');
plt.show()

# derivatives of the first order shape function
def dNdx(xp,xi,i):
    last_i = len(xi)-1 
    if(i == 0 and xp <= xi[1]):
        y = (-1)/(xi[1]-xi[0])
    elif(i == last_i and xp >= xi[last_i-1]):
        y = (1)/(xi[last_i]-xi[last_i-1])
    elif(xp >= xi[i] and xp <= xi[i+1]):
        y = (-1)/(xi[i+1]-xi[i])
    elif(xp <= xi[i] and xp >= xi[i-1]):
        y = (1)/(xi[i]-xi[i-1])
    else:
        y = 0
    return y

# ploting shape functions
fig = plt.figure(figsize=(8,5))
for j in range(len(xi)):
    dNidx = np.zeros(np.size(x))
    for i in range(len(x)):
        dNidx[i] = dNdx(x[i],xi,j)

    plt.plot(x, dNidx, linewidth=2);
plt.grid(True);
plt.xlabel('$x$');
plt.ylabel('$Ni$');
fig.suptitle('Fig. 3 - Derivatives of the first order shape functions');
plt.show()

# Initializing the variables
K_AB = np.zeros([dof,dof])
F_A = np.zeros(dof)

# Bilinear operetor a(w,v)
def innerprod_a(xp,xi,a,b):  
    y = dNdx(xp,xi,a)*dNdx(xp,xi,b)
    return y

# Bilinear operetor (w,f)
def innerprod(xp,xi,a):  
    y = N(xp,xi,a)*f(xp)
    return y

# loop over the N_A functions
for a in range(dof):
    # loop over the N_B functions
    for b in range(dof): # cold be considerably optimized due to the kronecker delta properties
        K_AB[a,b] = quad(innerprod_a,0,1,args=(xi,a,b))[0]
    
    F_A[a] = quad(innerprod,0,1,args=(xi,a))[0] + N(0,xi,a)*h - q*quad(innerprod_a,0,1,args=(xi,a,dof))[0]

print("\n")        
print("K_AB =")
print(K_AB)
print("\n")
print("F_A =")
print(F_A)

# solving the linea system
d = np.linalg.solve(K_AB, F_A)
print("\n System solution d = ")
print(d)

# Adding d_(n+1)
d = np.append(d,q)
print("\n Adding d_(n+1) = q")
print(d)

#ploting solution.
u_fem = np.zeros(np.size(x))
for i in range(len(x)): # loop over all sample points
    for j in range(len(xi)): # loop over all shape functios
    
        u_fem[i] = u_fem[i] + N(x[i],xi,j)*d[j]

# calculo do RMSE
RMSE = np.sqrt(np.mean((u_fem - u)**2))

print("\n")
print("Erro quadrático médio (RMSE) - FEM:")
print(RMSE)
        
fig = plt.figure(figsize=(8,5))
plt.plot(x, u_fem, 'r', linewidth=2);
plt.plot(x, u, 'k', linewidth=2);
plt.grid(True);
plt.xlabel('$x$');
plt.ylabel('$u$');
fig.suptitle('Fig. 3 - FEM solution');

plt.text(
    0.05, 0.90,
    f'RMSE = {RMSE:.6e}',
    transform=plt.gca().transAxes,
    bbox=dict(boxstyle='round', facecolor='white', edgecolor='black')
)

plt.show()