from numpy import *
import matplotlib.pyplot as plt


"""
Intégration numérique de la pression par la méthode de Rounge Kouta

"""

### Constante :

gamma = 1.3
L = 152* 1e-3
R = 42.2 * 1e-3
D = 86.4 * 1e-3
Vc = pi * D**2 * R/2
beta = L/R
tho = 11.8
o_d = -2*pi
do_conv = 2*pi
Qtot = 2800 *1e3 # [kj/kg] 2800 pour essence / 1650 pour le diesel
Rgaz = 8.31415
Tadmi = 303.15
mmair = 114



def pression_RK(theta, h):
    theta = radians(theta)
    print(theta)
    print(arange(-pi,pi+0.01,0.01))
    V_output = Vc / 2 * (1 - cos(theta) + beta - sqrt(beta ** 2 - sin(theta) ** 2)) + Vc / (tho - 1)
    m = 1e5 * mmair* V_output[0] / ( Rgaz * Tadmi)  * 1e-3


    V = lambda x: Vc/2 * (1 - cos(x) + beta - sqrt(beta**2 - sin(x)**2)) + Vc / (tho -1)
    dVdo = lambda x: Vc/2 * sin(x) * (cos(x) / (sqrt(beta**2 - sin(x)**2)) + 1)
    dQdo = lambda o : m* Qtot/2 * (pi/do_conv * sin(pi*(o - o_d)/do_conv))
    f = lambda o,p: -gamma * p / V(o) * dVdo(o) + (gamma -1) * 1/ V(o) * dQdo(o)


    o = theta
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


pression_RK(arange(-180,180+0.01,0.01),radians(0.01))

