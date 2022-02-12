from numpy import *
import matplotlib.pyplot as plt


"""
Intégration numérique de la pression par la méthode de Rounge Kouta

"""

### Constante :

gamma = 1.3
L =
R =
D =
Vc = pi* D**2 * R/2
beta = L/R
Vmax =
Vmin =
o_d =
do_conv =
Qtot = 2800 # [kj/kg] 2800 pour essence / 1650 pour le diesel
tho  = Vmax/Vmin # à déterminer apres avoir obtenu (D,R,L,Vc)



def pression(h):

    V = lambda x: Vc/2 * (1 - cos(x) + beta - sqrt(beta**2 -  (sin(x)**2))) + Vc / (tho -1)
    dVdo = lambda x: Vc/2 * sin(x) * (cos(x) / (sqrt(beta**2 - (sin(x))**2)) + 1)
    dQdo = lambda o : Qtot/2 * (pi/do_conv * sin((o - o_d)/do_conv))
    f = lambda o,p: -gamma * p / V(o) * dVdo(o) + (gamma -1 ) * 1/ V(o) * dQdo(o)


    o = arange(-6, 6+h, h)
    p = zeros_like(o)

    for i in range(len(o)-1):
        K1 = f(o[i],p[i])
        K2 = f(o[i]+h,p[i]+h*K1)
        p[i+1] = p[i] + h/2 * (K1 + K2)


    plt.plot(o,p,"or")
    plt.plot(o,p)
    plt.show()
    return p





