import math


# Parametros do problema
m = 0.073        # kg
g = 9.81         # m/s^2
k = 6.51e-5      # N*m^2/A^2
R = 11.0         # ohm

# Condicao inicial de referencia
x0 = 8.5e-3      # m
v0 = 0.0         # m/s


def L(x: float) -> float:
    """Indutancia equivalente."""
    return 2.0 * k / x


def fm(i: float, x: float) -> float:
    """Forca magnetica."""
    return k * i**2 / x**2


def Bl(i: float, x: float) -> float:
    """Fator eletromecanico."""
    return -2.0 * k * i / x**2


def equilibrio() -> tuple[float, float]:
    """Calcula corrente e tensao no ponto de equilibrio."""
    i0 = x0 * math.sqrt(m * g / k)
    u0 = R * i0
    return i0, u0


def main() -> None:
    i0, u0 = equilibrio()

    print("Ponto de equilibrio")
    print(f"x0 = {x0:.6e} m")
    print(f"v0 = {v0:.6e} m/s")
    print(f"i0 = {i0:.6e} A")
    print(f"u0 = {u0:.6e} V")
    print()
    print("Funcoes do modelo")
    print(f"L(x0) = {L(x0):.6e} H")
    print(f"fm(i0, x0) = {fm(i0, x0):.6e} N")
    print(f"Bl(i0, x0) = {Bl(i0, x0):.6e}")


if __name__ == "__main__":
    main()
