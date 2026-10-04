import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve


# ============================================================
# Parâmetros gerais
# ============================================================

Omega = [0.0, 1.0]

# Valores exigidos pelo enunciado
n_values = [4, 8, 16, 32]

# Malha comum usada apenas para calcular o RMSE
# e comparar os métodos no domínio inteiro
x_eval = np.linspace(Omega[0], Omega[1], 2001)


# ============================================================
# Solução manufaturada
# ============================================================

def ue(x):
    """
    Solução analítica:
    u_e(x) = x² e^x + cos(3x)
    """
    return (x**2) * np.exp(x) + np.cos(3 * x)


def f(x):
    """
    Função fonte:
    f(x) = u_e''(x)
    """
    return (
        (x**2 + 4*x + 2) * np.exp(x)
        - 9 * np.cos(3 * x)
    )


# Condições de contorno
u0 = ue(Omega[0])
u1 = ue(Omega[1])


# ============================================================
# FDM
# ============================================================

def solve_fdm(n):
    """
    Resolve o problema pelo Método das Diferenças Finitas.

    n = número de pontos internos.

    Retorna:
        x_nodes  -> todos os nós, incluindo contornos
        u_nodes  -> solução FDM nesses nós
    """

    h = 1.0 / (n + 1)

    # pontos internos
    x_internal = np.linspace(
        Omega[0] + h,
        Omega[1] - h,
        n
    )

    # matriz tridiagonal
    main_diag = 2 * np.ones(n)
    upper_diag = -1 * np.ones(n - 1)
    lower_diag = -1 * np.ones(n - 1)

    A = diags(
        [lower_diag, main_diag, upper_diag],
        offsets=[-1, 0, 1],
        format="csr"
    )

    # 2u_i - u_(i-1) - u_(i+1) = -h² f(x_i)
    b = -h**2 * f(x_internal)

    # condições de contorno
    b[0] += u0
    b[-1] += u1

    # solução nos pontos internos
    u_internal = spsolve(A, b)

    # adiciona os pontos de contorno
    x_nodes = np.concatenate(
        ([Omega[0]], x_internal, [Omega[1]])
    )

    u_nodes = np.concatenate(
        ([u0], u_internal, [u1])
    )

    return x_nodes, u_nodes


# ============================================================
# FEM
# ============================================================

def N1(x, xa, xb):
    return (xb - x) / (xb - xa)


def N2(x, xa, xb):
    return (x - xa) / (xb - xa)


def solve_fem(n):
    """
    Resolve o problema pelo Método dos Elementos Finitos.

    n = número de pontos internos.

    Retorna:
        x_nodes  -> nós da malha
        u_nodes  -> solução FEM nos nós
    """

    n_nodes = n + 2
    n_elements = n_nodes - 1

    x_nodes = np.linspace(
        Omega[0],
        Omega[1],
        n_nodes
    )

    K = np.zeros((n_nodes, n_nodes))
    F = np.zeros(n_nodes)

    # --------------------------------------------------------
    # Montagem elemento por elemento
    # --------------------------------------------------------

    for e in range(n_elements):

        xa = x_nodes[e]
        xb = x_nodes[e + 1]

        he = xb - xa

        # matriz de rigidez local
        Ke = (1.0 / he) * np.array([
            [1.0, -1.0],
            [-1.0, 1.0]
        ])

        # vetor local
        # forma fraca:
        # integral w' u' dx = - integral w f dx

        Fe = np.zeros(2)

        Fe[0] = -quad(
            lambda x: N1(x, xa, xb) * f(x),
            xa,
            xb
        )[0]

        Fe[1] = -quad(
            lambda x: N2(x, xa, xb) * f(x),
            xa,
            xb
        )[0]

        nodes = [e, e + 1]

        for a in range(2):

            A = nodes[a]
            F[A] += Fe[a]

            for b in range(2):

                B = nodes[b]
                K[A, B] += Ke[a, b]

    # --------------------------------------------------------
    # Condições de contorno
    # --------------------------------------------------------

    F = F - K[:, 0] * u0
    F = F - K[:, -1] * u1

    K[0, :] = 0.0
    K[:, 0] = 0.0
    K[0, 0] = 1.0
    F[0] = u0

    K[-1, :] = 0.0
    K[:, -1] = 0.0
    K[-1, -1] = 1.0
    F[-1] = u1

    # --------------------------------------------------------
    # Solução
    # --------------------------------------------------------

    u_nodes = np.linalg.solve(K, F)

    return x_nodes, u_nodes


# ============================================================
# Cálculo dos erros
# ============================================================

errors_fdm = []
errors_fem = []
h_values = []

print("\nAnálise de convergência")
print("=" * 65)

print(
    f"{'n':>5} "
    f"{'h':>12} "
    f"{'RMSE FEM':>18} "
    f"{'RMSE FDM':>18}"
)

print("-" * 65)


for n in n_values:

    h = 1.0 / (n + 1)
    h_values.append(h)

    # --------------------------------------------------------
    # FDM
    # --------------------------------------------------------

    x_fdm, u_fdm_nodes = solve_fdm(n)

    # interpolação linear sobre a malha comum
    u_fdm_eval = np.interp(
        x_eval,
        x_fdm,
        u_fdm_nodes
    )

    rmse_fdm = np.sqrt(
        np.mean(
            (u_fdm_eval - ue(x_eval))**2
        )
    )

    errors_fdm.append(rmse_fdm)

    # --------------------------------------------------------
    # FEM
    # --------------------------------------------------------

    x_fem, u_fem_nodes = solve_fem(n)

    # Como o FEM é linear por partes,
    # np.interp reproduz exatamente a solução FEM.
    u_fem_eval = np.interp(
        x_eval,
        x_fem,
        u_fem_nodes
    )

    rmse_fem = np.sqrt(
        np.mean(
            (u_fem_eval - ue(x_eval))**2
        )
    )

    errors_fem.append(rmse_fem)

    print(
        f"{n:5d} "
        f"{h:12.6f} "
        f"{rmse_fem:18.8e} "
        f"{rmse_fdm:18.8e}"
    )


# Conversão para arrays
h_values = np.array(h_values)
errors_fem = np.array(errors_fem)
errors_fdm = np.array(errors_fdm)


# ============================================================
# Estimativa da ordem de convergência
# ============================================================

order_fem = np.polyfit(
    np.log(h_values),
    np.log(errors_fem),
    1
)[0]

order_fdm = np.polyfit(
    np.log(h_values),
    np.log(errors_fdm),
    1
)[0]


print("\nOrdens aproximadas de convergência:")
print(f"FEM: {order_fem:.4f}")
print(f"FDM: {order_fdm:.4f}")


# ============================================================
# Gráfico de convergência em função de n
# ============================================================

plt.figure(figsize=(8, 5))

plt.loglog(
    n_values,
    errors_fem,
    marker="o",
    label="FEM"
)

plt.loglog(
    n_values,
    errors_fdm,
    marker="s",
    label="FDM"
)

plt.xlabel("Número de pontos internos n")
plt.ylabel("RMSE")
plt.title("Convergência dos métodos FEM e FDM")

plt.grid(True, which="both")
plt.legend()

plt.savefig(
    "convergencia.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()


# ============================================================
# Gráfico de erro em função de h
# ============================================================

plt.figure(figsize=(8, 5))

plt.loglog(
    h_values,
    errors_fem,
    marker="o",
    label=f"FEM - ordem ≈ {order_fem:.2f}"
)

plt.loglog(
    h_values,
    errors_fdm,
    marker="s",
    label=f"FDM - ordem ≈ {order_fdm:.2f}"
)

plt.xlabel("Espaçamento da malha h")
plt.ylabel("RMSE")
plt.title("Erro em função do refinamento da malha")

plt.grid(True, which="both")
plt.legend()

# Inverte o eixo para a malha ficar mais refinada
# da esquerda para a direita.
plt.gca().invert_xaxis()

plt.savefig(
    "convergencia_h.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()


# ============================================================
# Salvar resultados em CSV
# ============================================================

results = np.column_stack(
    (
        n_values,
        h_values,
        errors_fem,
        errors_fdm
    )
)

np.savetxt(
    "convergencia.csv",
    results,
    delimiter=",",
    header="n,h,RMSE_FEM,RMSE_FDM",
    comments="",
    fmt=[
        "%d",
        "%.10e",
        "%.10e",
        "%.10e"
    ]
)


print("\nArquivos gerados:")
print("- convergencia.png")
print("- convergencia_h.png")
print("- convergencia.csv")