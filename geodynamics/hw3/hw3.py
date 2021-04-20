import numpy as np
import scipy.special as sp
import matplotlib.pyplot as plt

def uniform_loading(r, x, a, rou_p, l):
    if x <= r:
        kerpa = sp.kerp(a)
        berx = sp.ber(x/l)
        keipa = sp.keip(a)
        beix = sp.beip(x/l)
        w = rou_basalt*uniform_h*(a*kerpa*berx-a*keipa*beix+1)/rou_p
    else:
        berpa = sp.berp(a)
        kerx = sp.ker(x/l)
        beipa = sp.beip(a)
        keix = sp.kei(x/l)
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
                    w[i,j,k] += uniform_loading(r, x[j], a, rou_prime[i], l[i])

    return w

def plot_result():
    fontdict = {'family' : 'Times New Roman', 'weight' : 'normal', 'size' : 18}
    fontdict2 = {'family' : 'Times New Roman', 'weight' : 'normal', 'size' : 14}
    plt.figure(figsize=(13, 8))
    plt.plot(t, w[0,0,:], c='red', label='thickness=30km, r=150km')
    plt.plot(t, w[0,1,:], c='green', label='thickness=30km, r=300km')
    plt.plot(t, w[0,2,:], c='blue', label='thickness=30km, r=450km')
    plt.plot(t, w[1,0,:], c='red', linestyle='--', label='thickness=50km, r=150km')
    plt.plot(t, w[1,1,:], c='green', linestyle='--', label='thickness=50km, r=300km')
    plt.plot(t, w[1,2,:], c='blue', linestyle='--', label='thickness=50km, r=450km')
    plt.plot(t, w[2,0,:], c='red', linestyle=':', label='thickness=70km, r=150km')
    plt.plot(t, w[2,1,:], c='green', linestyle=':', label='thickness=70km, r=300km')
    plt.plot(t, w[2,2,:], c='blue', linestyle=':', label='thickness=70km, r=450km')
    plt.axis([0, 2.5, -1.8, 0.1])
    plt.xlabel("Time (Ma)", fontdict=fontdict)
    plt.ylabel("Combined topography at different location (km)", fontdict=fontdict)
    plt.xticks(fontproperties='Times New Roman', size=16)
    plt.yticks(fontproperties='Times New Roman', size=16)
    plt.legend(prop=fontdict2, loc=1)
    plt.savefig('result.pdf', bbox_inches='tight',dpi=600,pad_inches=0.0)
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
    steps = np.linspace(0.0005, 2, 4000)
    uniform_h = 0.0005
    A = v*t
    
    w = -calc()
    plot_result()


