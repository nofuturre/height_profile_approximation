import numpy as np
from math import sqrt
from sympy import symbols,simplify
from matplotlib import pyplot as plt

points = 14   # ilosc punktów testowych
all_points = 513
span = int(all_points/points)


def f(x, datax, i):
    suma = 1
    for j in range(len(datax)):
        if i != j:
            suma *= (x - datax[j])/(datax[i] - datax[j])
    return suma


def F(x, datax, datay):
    suma = 0
    for i in range(len(datax)):
        suma += datay[i] * f(x, datax, i)
    return suma


def interpolacja_Lagrangea(datax, datay):
    L = []
    for i in range(all_points):
        L.append(F(i, datax, datay))
    return L


def S(x, x0, a, b, c, d):
    return a + b*(x - x0) + c*pow(x - x0, 2) + d*pow(x-x0, 3)


def interpolacja_splajtami(datax, datay):
    n = len(datax) - 1    # ilość wspolczynnikow a, b, c, d

    b = np.zeros(n)
    c = np.zeros(n+1)
    d = np.zeros(n)

    xdiff = np.diff(datax)
    ydiff = np.diff(datay)

    alfa = np.empty(n)
    L = np.ones(n)
    U = np.zeros(n)
    z = np.zeros(n)

    for i in range(1, n):
        alfa[i] = ((3/xdiff[i])*(ydiff[i]) - (3/xdiff[i-1])*(ydiff[i]))
        L[i] = (2*(datax[i+1] - datax[i-1]) - xdiff[i-1]*U[i-1])
        U[i] = xdiff[i]/L[i]
        z[i] = (alfa[i] - xdiff[i-1]*z[i-1])/L[i]

    for i in range(n-1, 0, -1):
        c[i] = z[i] - U[i]*c[i+1]
        b[i] = (ydiff[i])/xdiff[i] - xdiff[i]*(c[i+1]+ 2*c[i])/3
        d[i] = (c[i+1]-c[i])/3*xdiff[i]

    W = []
    for i in range(1, len(datax)):
        for j in range(datax[i-1], datax[i]):
            W.append(S(j, datax[i-1], datay[i-1], b[i-1], c[i-1], d[i-1]))
    return W


def show_plot(name):
    data = np.genfromtxt(f'{name}.csv', delimiter=",", names=["distance", "height"])    # wczytanie danych
    plt.plot(data["height"])

    test_data = [data["height"][i]for i in range(1, len(data["height"]), span)]     # wybór punktów testowych
    x = [i - 1 for i in range(1, len(data["height"]), span)]

    if data["height"][len(data["height"])-1] not in test_data:      # upewnienie sie ze ostani punkt zawira sie w punktach testowych
        test_data.append(data["height"][len(data["height"])-1])
        x.append(len(data["height"])-1)

    plt.plot(x, test_data, 'o', color='red')        # zaznaczenie punktów testowych na wykresie

    func1 = interpolacja_Lagrangea(x, test_data)
    plt.plot(func1, color='orange')

    func2 = interpolacja_splajtami(x, test_data)
    #plt.plot(func2, color='orange')
    plt.show()


if __name__ == '__main__':
    show_plot("many_peaks")
    show_plot("one_peak")
    show_plot("ascending")
