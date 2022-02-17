from numpy import *
import matplotlib.pyplot as plt


"""
Intégration numérique de la pression par la méthode de Rounge Kouta

"""

### Constante :

gamma = 1.3
L = 0.15
R = 0.05
D = 0.05
Vc = pi * D**2 * R/2
beta = L/R
tho = 10
Vmax = tho / (tho -1) * Vc
Vmin = 1/(tho-1) * Vc
o_d = -2*pi
do_conv = 2*pi
Qtot = 0.0050*2800 *1e3 # [kj/kg] 2800 pour essence / 1650 pour le diesel




def pression_RK(h):

    V = lambda x: Vc/2 * (1 - cos(x) + beta - sqrt(beta**2 - sin(x)**2)) + Vc / (tho -1)
    dVdo = lambda x: Vc/2 * sin(x) * (cos(x) / (sqrt(beta**2 - sin(x)**2)) + 1)
    dQdo = lambda o : Qtot/2 * (pi/do_conv * sin(pi*(o - o_d)/do_conv))
    f = lambda o,p: -gamma * p / V(o) * dVdo(o) + (gamma -1 ) * 1/ V(o) * dQdo(o)


    o = arange(-2*pi, 2*pi, h)
    p = zeros_like(o)

    for i in range(len(o)-1):
        K1 = f(o[i], p[i])
        K2 = f(o[i] + h, p[i] + h * K1)
        p[i + 1] = p[i] + h / 2 * (K1 + K2)


    plt.plot(o,p*1e-5,"or")
    plt.plot(o,p*1e-5)
    plt.xlabel("$\Theta$ [rad]")
    plt.ylabel("$p$ [bar]")
    plt.show()
    return p


pression_RK(0.1)

