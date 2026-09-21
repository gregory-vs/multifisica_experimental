import numpy as np
import matplotlib.pyplot as plt

# input parameters
Omega = [0,1]
n = 6
h = 1/(n+1)
print(f'h = {h}')
xj = np.linspace(Omega[0]+h,Omega[1]-h,n)

# f(x)
def fx(x):
    return ((x**2 + 4*x + 2)*np.exp(x) - 9*np.cos(3*x))

fig = plt.figure(figsize=(8,3))
plt.plot(xj,fx(xj),'.b')
plt.plot(0,fx(0),'.r')
plt.plot(1,fx(1),'*r')
plt.legend(['f(xj)', 'f0', 'f1'])

plt.title('f(x)')
plt.xlabel('xj')
plt.ylabel('f(xj)')
plt.grid(True)
plt.show()

from scipy.sparse import diags

# exact solution
def ue(x):
    return (x**2)*np.exp(x) + np.cos(3*x)
xe = np.linspace(-1,1,100)

# Define the diagonals
main_diag = 2 * np.ones(n)
upper_diag = -1 * np.ones(n - 1)
lower_diag = -1 * np.ones(n - 1)

# Create the tridiagonal matrix A
A = diags(
    [lower_diag, main_diag, upper_diag], offsets=[-1, 0, 1], format='csr')

# condicoes de contorno
u0 = ue(0)
u1 = ue(1)
print(u0, " ", u1)

# Right-hand side
b = -h**2 * fx(xj)

# Incorporate boundary conditions
b[0] += u0
b[-1] += u1

print("Sparse Matrix A:\n", A)

plt.spy(A)
plt.show()

print("Vector b:\n", b)

from scipy.sparse.linalg import spsolve

# Solve the system Ax = b
vj = spsolve(A, b)

# calculo do RMSE

RMSE = np.sqrt(np.mean((vj - ue(xj))**2))

print("Erro quadrático médio (RMSE) - FDM:")
print(RMSE)

fig = plt.figure(figsize=(8,4))
plt.plot(xj,vj,'.b')
plt.plot(xe,ue(xe),'g')
plt.legend(['vj', 'ue'])

plt.title('Solution')
plt.xlabel('xj')
plt.ylabel('u(xj)')
plt.grid(True)

plt.text(
    0.05, 0.90,
    f'RMSE = {RMSE:.6e}',
    transform=plt.gca().transAxes,
    bbox=dict(boxstyle='round', facecolor='white', edgecolor='black')
)

plt.show()