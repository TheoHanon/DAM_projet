import matplotlib.pyplot as plt
from numpy import *

tau = 11.8#[-]
D = 86.4 * 1e-3 #[m]
C = 2* 42.2 * 1e-3 #[m]
L = 152* 1e-3 #[m]
mpiston = 0.310#[kg]
mbielle = 0.544 #[kg]
Q = 2800 * 1e3 #[J/kg_inlet gas]




def my_func(rpm, s, theta, thetaC, deltaThetaC):
    #theta = radians(theta)
    thetaC = -radians(thetaC)
    deltaThetaC = radians(deltaThetaC)
    w = rpm * 2 * pi / 60

    Rgaz = 8.31415
    R = C / 2
    Vc = pi * (D ** 2) * R / 2
    beta = L / R
    gamma = 1.3
    h = round(theta[1]-theta[0],4)
    print(h)
    Tadmi = 303.15
    mmair = 28.965

    ### V_output ###

    V_output = Vc / 2 * (1 - cos(theta) + beta - sqrt(beta ** 2 - sin(theta) ** 2)) + Vc / (tau - 1)

    ### masse initial###

    m = s* 1e5 *mmair *  V_output[0] / ( Rgaz * Tadmi) * 1e-3

    ### p_output ###

    V = lambda x: Vc / 2 * (1 - cos(x) + beta - sqrt(beta ** 2 - sin(x) ** 2)) + Vc / (tau - 1)
    dVdo = lambda x: Vc / 2 * sin(x) * (cos(x) / (sqrt(beta ** 2 - sin(x) ** 2)) + 1)
    dQdo = lambda theta: m*Q / 2 * (pi / deltaThetaC * sin(pi * (theta - thetaC) / deltaThetaC))
    f = lambda theta, p: -gamma * p / V(theta) * dVdo(theta) + (gamma - 1) * 1 / V(theta) * dQdo(theta)


    p_output = zeros_like(theta)

    for i in range(len(theta) - 1):
        K1 = f(theta[i], p_output[i])
        K2 = f(theta[i] + h, p_output[i] + h * K1)
        p_output[i + 1] = p_output[i] + h / 2 * (K1 + K2)



    ### Q_output ###

    Q_output = m*Q/2 * (1 - cos(pi* (theta - thetaC) / deltaThetaC ))

    ### F_pied_output ###

    F_pied_output = pi * D**2 /4 * p_output - mpiston * R*(w)**2 * cos(theta)

    ### F_tete_output ###

    F_tete_output = - pi * D**2 /4 * p_output + (mpiston +mbielle) * R*(w)**2 * cos(theta)

    ### F_compression ###

    F_compression = max(pi* D**2 /4 * p_output + (mpiston+mbielle)*R*(w)**2 * cos(theta))

    ### t ###

    E = 200 * 1e9
    feuler = lambda I,K: pi**2 *I * E / ((K *L)**2)

    sigmaC = 450 * 1e6
    tt = arange(0.005, 0.1, 0.01) #de 0.5 cm -> 10 cm, par pas de 1cm
    A = 11 * (tt ** 2)

    Fcrit_x = 1 / (1 / feuler(419/12 * (tt**4), 1)  + 1 / (sigmaC* A)) ### selon x
    Fcrit_y = 1 / (1 / feuler(131/12 * (tt**4), 0.5)  + 1 / (sigmaC* A)) ### selon y

    t = max(min(tt[abs((Fcrit_x - F_compression)) > 1e3]) , min(tt[abs((Fcrit_y - F_compression)) > 1e3 ]))


    return (V_output,Q_output,F_pied_output,F_tete_output,p_output, t)



rpm = 2412
s = 1.2
thetaC = 30.0
deltaThetaC = 60.0
theta = arange(-pi, pi + 0.001,0.001)

V, Q,Fp, Ft, p, t = my_func(rpm,s,theta,thetaC, deltaThetaC)


plt.plot(theta,p*1e-5)
plt.show()