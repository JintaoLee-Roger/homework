import numpy as np
import matplotlib.pyplot as plt
import os, re

def loadfromlines(lines):
    assert (len(lines[-1].split()) == 4)
    lines[-1] += '  0'*4
    d = np.loadtxt(lines).reshape(-1)
    return d[:-13]

def readsp(filename):
    f = open(filename)
    lines = f.readlines()

    loc = []
    for i, line in enumerate(lines):
        if line[:7] == 'Damping':
            loc.append(i)
    sd = np.zeros((3, N, 5), float)
    sv = np.zeros((3, N, 5), float)
    sa = np.zeros((3, N, 5), float)
    for k, i in enumerate(loc):
        sd[k//5, :, k%5] = loadfromlines(lines[i+1 : i+1+num])
        sv[k//5, :, k%5] = loadfromlines(lines[i+1+num : i+1+2*num])
        sa[k//5, :, k%5] = loadfromlines(lines[i+1+2*num : i+1+3*num])

    return sd, sv, sa

def plotfigure(sd, sv, sa):
    plt.figure(figsize=(9, 5))
    damp = [0, 0.02, 0.05, 0.10, 0.20]
    for i in range(5):
        plt.plot(t[:-9], sa[0,:, i], label='damp={}'.format(damp[i]))
    plt.yscale('log')
    # plt.xscale('log')
    plt.legend()
    plt.xlabel("period")
    plt.ylabel("Sa($cm/s^2$)")
    plt.title("Acceleration Response")
    plt.savefig(savep+"sa.pdf",bbox_inches='tight',dpi=300,pad_inches=0.0)

    plt.figure(figsize=(9, 5))
    for i in range(5):
        plt.plot(t[:-9], sv[0,:, i], label='damp={}'.format(damp[i]))
    plt.yscale('log')
    # plt.xscale('log')
    plt.legend()
    plt.xlabel("period")
    plt.ylabel("Sv(cm/s)")
    plt.title("Velocity Response")
    plt.savefig(savep+"sv.pdf",bbox_inches='tight',dpi=300,pad_inches=0.0)

    plt.figure(figsize=(9, 5))
    for i in range(5):
        plt.plot(t[:-9], sd[0,:, i], label='damp={}'.format(damp[i]))
    plt.yscale('log')
    # plt.xscale('log')
    plt.legend()
    plt.xlabel("period")
    plt.ylabel("Sd(cm)")
    plt.title("Displacement Response")
    plt.savefig(savep+"sd.pdf",bbox_inches='tight',dpi=300,pad_inches=0.0)


    plt.show()

if __name__ == "__main__":
    filename = 'D:\\Data\\eseis\\CICLC\\CICLC.v3'
    savep = 'D:\\homework\\engseis\\hw1\\figure\\'
    N = 8*12+4-9
    num = 13
    t = np.arange(0.04, 16, 0.16)
    sd,sv,sa = readsp(filename)
    plotfigure(sd, sv, sa)
    

