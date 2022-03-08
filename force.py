import matplotlib.pyplot as plt
from numpy import *
from pression import *

mpiston = 0.310
mbielle = 0.544
w = 300
h = 0.1

def force() :


    theta = arange(-pi,pi +h,h)

    F_pied_output = pi * D ** 2 / 4 * pression_RK(theta,h) - mpiston * R * w ** 2 * cos(theta)
    F_tete_output = - pi * D ** 2 / 4 * pression_RK(theta,h) + (mpiston + mbielle) * R * w ** 2 * cos(theta)
    E = 200 *1e9
    feuler = lambda I, K: (pi ** 2) * I * E / ((K * L) ** 2)

    sigmaC = 450*1e6
    tt = arange(0.0001, 0.005, 0.001)  # de 0.5 cm -> 10 cm, par pas de 1cm
    A = 11 * (tt ** 2)

    Fcrit_x = array(1/ ((1 / feuler(419 * (tt ** 4)/12, 1)) + (1 / (sigmaC * A))) )### selon x
    Fcrit_y = array(1/ ((1 / feuler(131 * (tt ** 4)/12, 0.5)) + (1 / (sigmaC * A))) ) ### selon y

    F_compression = max(pi * D ** 2 / 4 * pression_RK(theta,h) + (mpiston + mbielle) * R * w** 2 * cos(theta))
    print(tt[abs(Fcrit_y-F_compression) > 1e3])
    plt.plot(tt*1e2,abs(Fcrit_y-F_compression))
    plt.plot(tt*1e2,abs(Fcrit_x-F_compression))
    #plt.plot([0,0.6],[F_compression,F_compression])
    #plt.plot(tt,1/ (sigmaC*A))
    #plt.plot(tt,Fcrit_y)
    #plt.plot(theta,F_tete_output)
    #plt.plot(theta,F_pied_output)

    plt.xlabel("$t$ [cm]")
    plt.ylabel("$F$ [N]")
    plt.show()


force()