import matplotlib.pyplot as plt
from numpy import *

tau = 10  # [-]
D = 86.4 * 1e-3  # [m]
C = 2 * 45 * 1e-3  # [m]
L = 172.6 * 1e-3  # [m]
mpiston = 0.310  # [kg]
mbielle = 0.544  # [kg]
Q = 2800 * 1e3 # [J/kg_inlet gas]


def my_func(rpm, s, theta, thetaC, deltaThetaC):
    theta = radians(theta)
    thetaC = -radians(thetaC)
    deltaThetaC = radians(deltaThetaC)
    w = rpm / 60 * 2 * pi

    Rgaz = 287.1
    R = C / 2
    Vc = pi * (D ** 2) * R / 2
    beta = L / R
    gamma = 1.3
    h = theta[1] - theta[0]
    Tadmi = 303.15


    ### V_output ###

    V = lambda x: Vc / 2 * (1 - cos(x) + beta - sqrt(beta ** 2 - sin(x) ** 2)) + Vc / (tau - 1)
    V_output = V(theta)

    ### masse initial###

    m = s * 1e5 * V_output[0] / (Rgaz * Tadmi)

    ### Q_output ###
    Q_output = zeros_like(theta)

    start = min(where(abs(theta - thetaC) < 1e-2)[0])
    finish = max(where(abs(theta - (thetaC + deltaThetaC)) < 1e-2)[0])

    Q_output[start:finish + 1] = m * Q / 2 * (1 - cos(pi * (theta[start:finish + 1] - thetaC) / deltaThetaC))
    ### p_output ###

    dVdo = lambda x: Vc / 2 * sin(x) * (cos(x) / (sqrt(beta ** 2 - sin(x) ** 2)) + 1)


    f = lambda theta, p, index: -gamma * (p / V(theta)) * dVdo(theta) + (gamma - 1) * (1 / V(theta)) * dQdtheta(start,finish,theta, index, m, Q,thetaC, deltaThetaC)

    p_output = zeros_like(theta)
    p_output[0] = 1e5 * s

    for i in range(len(theta) - 1):
        K1 = f(theta[i], p_output[i],i)
        K2 = f(theta[i] + h, p_output[i] + h * K1, i)
        p_output[i + 1] = p_output[i] + h / 2 * (K1 + K2)



    ### F_pied_output ###

    F_pied_output = pi * D ** 2 / 4 * p_output - mpiston * R * w ** 2 * cos(theta)

    ### F_tete_output ###

    F_tete_output = - pi * D ** 2 / 4 * p_output + (mpiston + mbielle) * R * w ** 2 * cos(theta)

    ### F_compression ###

    F_compression = max(pi * D ** 2 / 4 * p_output + (mpiston + mbielle) * R * w ** 2 * cos(theta))

    ### t ###

    E = 200 * 1e9
    feuler = lambda I, K: pi ** 2 * I * E / ((K * L) ** 2)

    sigmaC = 450 * 1e6
    tt = arange(0.00005, 0.005, 0.00001)  # de 0.5 cm -> 10 cm, par pas de 1cm
    A = 11 * (tt ** 2)

    Fcrit_x = 1 / (1 / feuler(419 / 12 * (tt ** 4), 1) + 1 / (sigmaC * A))  ### selon x
    Fcrit_y = 1 / (1 / feuler(131 / 12 * (tt ** 4), 0.5) + 1 / (sigmaC * A))  ### selon y



    t = max(min(tt[abs((Fcrit_x - F_compression)) < 1000]), min(tt[abs((Fcrit_y - F_compression)) < 1000]))

    return (V_output, Q_output, F_pied_output, F_tete_output, p_output, t)


def dQdtheta(start,finish,theta, index, m, Q ,thetaC, deltaThetaC):
    if (index >= start and index <= finish):
        return m * pi * Q * sin(pi * (theta - thetaC) / deltaThetaC) / (2 * deltaThetaC)
    else: return 0

rpm = 4339
s = 2.9
thetaC = 16
deltaThetaC = 67

theta = arange(-180, 180 + 0.1, 0.1)

V, Q, Fp, Ft, p, t = my_func(rpm, s, theta, thetaC, deltaThetaC)

#plt.plot(radians(theta), Fp)
#plt.plot(radians(theta),Ft)
print(t)
plt.plot(radians(theta),Q)
plt.show()
