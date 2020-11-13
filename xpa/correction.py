import numpy as np
from math import sin, cos, pi, radians, tan, sqrt, log1p
from scipy import stats


def get_interval(x, y, start, end):
    init = x.index(start)
    end = x.index(end)

    return x[init:end], y[init:end]


def remove_kalpha2(x, y):
    novoy = []
    lambida2 = 1.541220
    lambida1 = 1.537400
    deltaL = lambida2 - lambida1
    deltaL = deltaL / lambida1
    diferenca = x[1] - x[0]

    for i in range(len(y)):
        deltasoma = x[1] - x[0]
        ase = np.tan(np.radians(x[i] / 2)) * 2 * deltaL / (diferenca)
        n = 1

        while (ase > deltasoma):
            deltasoma = deltasoma + diferenca
            n += 1
        try:
            yy = y[i] - 0.5 * y[i - n]

            if yy < 0: yy = (yy + y[i]) / 8

            novoy.append(yy)
        except:
            novoy.append(y[i])

    return novoy


def remove_background(x, y, m=5):
    min_value = np.mean(np.sort(y)[:10])
    for i in range(len(y)):
        y[i] = y[i] - min_value
    slope, intercept, r_value, p_value, std_err = stats.linregress(np.append(x[:m], x[-m:]), np.append(y[:m], y[-m:]))
    line_values = [slope * i + intercept for i in x]
    line_values = np.asarray(line_values)
    return y - line_values


def normalize(vector):
    list_temp = []
    max_value = max(vector)
    for value in vector:
        list_temp.append(value / max_value)
    return list_temp


def lorentz_and_polarization(x, y):
    for i in range(len(y)):
        y[i] /= (1 + pow(cos(radians(x[i])), 2)) / (cos(radians(x[i])) * pow(sin(radians(x[i])), 2))

    return y


def centralizer(y1):
    middle = list(y1).index(max(y1))
    left = list(y1).index(max(y1))
    right = abs(list(y1).index(max(y1)) - len(y1))
    dif = abs(left - right)
    if left < right:
        new_vector = dif * [y1[0]] + list(y1[:middle - 1]) + list([y1[middle]]) + list(y1[middle + 1:])
    else:
        new_vector = list(y1[:middle - 1]) + list([y1[middle]]) + list(y1[middle + 1:]) + dif * [y1[-1]]

    new_vector = 50 * [y1[0]] + new_vector + 50 * [y1[-1]]

    return new_vector
