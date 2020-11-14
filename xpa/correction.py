import numpy as np

from scipy import stats

from math import sin, asin, cos, pi, radians, tan, sqrt, log1p
from lmfit.models import VoigtModel,PseudoVoigtModel, LinearModel, GaussianModel
import matplotlib.pyplot as plt


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


def centralize_peak(x, y):
    y = centralizer(y)
    x, y = made_x_after_centralizer(x, y)
    return x, y


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


def made_x_after_centralizer(x, y1):
    center = list(y1).index(max(y1))
    precision = len(str(x[213]).split('.')[-1])
    delta = x[center + 1] - x[center]
    delta = float(f"{delta:.{precision}f}")

    B = y1[:len(y1) // 2]
    C = y1[len(y1) // 2:]

    x2 = []
    add_x = x[center]
    for i in C:
        x2.append(add_x)
        add_x = add_x + delta

    x1 = []
    add_x = x[center] - delta
    for i in reversed(B):
        x1.append(add_x)
        add_x = add_x - delta
    x1 = x1[::-1]

    x3 = x1 + x2

    return x3, y1


def calc_d_rachinguer(angle):
    # lambda1 = 1.78892  # Angstron
    # lambda2 = 1.79278  # Angstron
    # lambda0  = 1.79021

    lambda2 = 1.541220
    lambda1 = 1.537400
    lambda0 = 1.5406

    lambda1_1 = lambda1/lambda0
    lambda2_1 = lambda2/lambda0

    angle = radians(angle)

    s_angle = sin(angle)

    d = 2 * ( asin(lambda2_1 * s_angle) - asin(lambda1_1 * s_angle)  )

    return d


def kalpha_is_value(x, y, i, j):
    y1 = y[i] - 0.5*y[j]
    return y1


def kalpha_interpolation(x, y, j, k):
    pass


def get_position_x(d, angle, x):
    angle_i_2 = angle - d
    for i in x:
        if i == angle_i_2:
            return # metodo para valor que encontrado - raro
        if i > angle_i_2:
            return # metodo para valor intermediario - mais comum


def rachinger_correction(x,y):
    for angle in x:
        d = calc_d_rachinguer(angle)
        index_i2 = get_position_x(d, angle, x)


def ScherrerMethod(x, y, xs, ys):

    lambida=0.154056

    #mod = methodfunciont((comboBox.get()))

    mod = GaussianModel()


    pars = mod.guess(y, x=x)
    pars1 = mod.guess(ys, x=xs)
    out  = mod.fit(y, pars, x=x)
    out1  = mod.fit(ys, pars1, x=xs)


    print(out.best_values)
    print(out1.best_values)
    try:

        #if mod = 'GaussianModel'
        z1=out.best_values['sigma']
        z2=out1.best_values['sigma']
        #else:
            #z1=out.best_values['sigma']
            #z2=out1.best_values['sigma']

        z1=pow(z1,2)
        z2=pow(z1,2)

        wg=sqrt(z1-z2)
        print(wg, 'radians')
        wg=np.radians(wg)

    except:
        pass

    print(wg)

    cosseno = out.best_values['center']/2
    cosseno = cos(np.radians(cosseno))


    #pdb.set_trace()
    scherrervalue = (0.94*lambida)/(wg*cosseno)

    scherrervalue = int(scherrervalue)

    print(scherrervalue , ' nm')
    plt.plot(x,y,'k--',label ='data')
    #plt.plot(xs,ys)
    #plt.plot(xs, out1.best_fit, 'k--')
    plt.plot(x, out.best_fit, 'k-',label='best fit')#str(scherrervalue) + ' nm'
    #plt.plot(x, out.init_fit, 'k--',label = 'init fit' )

    plt.title('Amostra')
    plt.xlabel('$2\Theta$')
    plt.ylabel("Intensity")
    plt.grid()

    plt.legend()
    plt.show()