import numpy as np
import scipy.special as sp
import matplotlib.pyplot as plt

def uniform_loading(r, x, a, rou_p):
    if x <= r:
        kerpa = sp.kerp(a)
        berx = sp.ber(x)
        keipa = sp.keip(a)
        beix = sp.beip(x)
        w = rou_basalt*uniform_h*(a*kerpa*berx-a*keipa*beix+1)/rou_p
    else:
        berpa = sp.berp(a)
        kerx = sp.ker(x)
        beipa = sp.beip(a)
        keix = sp.kei(x)
        w = rou_basalt*uniform_h*(a*berpa*kerx-a*beipa*keix)/rou_p

    return w


def calc():
    w = np.zeros((3, 3, 20), np.float64)
    for i in range(3):
        for j in range(3):
            for k in range(20):
                for step in steps:
                    r = A[k]*(10-step)/10
                    a = r/l[i]
                    w[i,j,k] += uniform_loading(r, x[j], a, rou_prime[i])

    return w

def plot_result():
    plt.figure()
    plt.plot(t, w[0,0,:])
    plt.plot(t, w[0,1,:])
    plt.plot(t, w[0,2,:])
    plt.show()

if __name__ == "__main__":
    ### parameters
    E = 7*1e16 # GPa, Young's modulus
    poisson = 0.25 # Poisson's ratios
    rou_basalt = 2700*9.8*1e9 # kg/m^3, basalt density
    rou_crust = 2900*9.8*1e9 # kg/m^3, crust density
    rou_mantle = 3370*9.8*1e9 # kg/m^3
    th_crust = np.array([30, 50, 70]) # km, crust thickness
    v = 300 # km/Ma, Basalts spread velocity
    h = 2 # km, basalts accumulated height
    T = 2 # Ma, Total time of volcanic eruption
    R = 6350 # km, 

    ### 
    D = E * pow(th_crust, 3) / (12 * (1-pow(poisson,2)))
    l = pow(D/((E*th_crust/pow(R, 2)) + rou_mantle), 1/4)
    rou_prime = rou_mantle + E*th_crust/pow(R, 2)
    x = np.array([150, 300, 450])
    t = np.linspace(0.1, 2, 20)
    steps = np.linspace(0.005, 2, 400)
    uniform_h = 0.005
    A = v*t
    
    w = calc()
    # print(w[0,0,0:5])
    plot_result()


