import numpy as np
from geopy.distance import geodesic
import os, re
import matplotlib.pyplot as plt

# 每个文件中有三个channels

def readpeak(filename):
    f = open(filename)
    lines = f.readlines()

    lonlat = lines[5].split()[3:5]
    # lonlat[1] = float(lonlat[1][:-1])
    # print(lonlat[1])
    # assert 0
    if lonlat[0][-2] == 'N':
        lonlat[0] = float(lonlat[0][:-2])
    elif lonlat[0][-2] == 'S':
        lonlat[0] = -float(lonlat[0][:-2])
    if lonlat[1][-1] == 'E':
        lonlat[1] = float(lonlat[1][:-1])
    elif lonlat[1][-1] == 'W':
        lonlat[1] = -float(lonlat[1][:-1])

    dis = geodesic((lonlat[0],lonlat[1]), (35.7695,-117.5993)).km

    pga = []
    pgv = []
    pgd = []
    loc = [0]
    for i, line in enumerate(lines):
        if line[:2] == '/&':
            loc.append(i+1)

    assert (len(loc) == 4)

    for i in loc[:3]:
        pga.append(float(lines[i + 17].split()[3]))
        pgv.append(float(lines[i + 18].split()[3]))
        pgd.append(float(lines[i + 19].split()[3]))

    return pga, pgv, pgd, dis

def readacc(lines):
    n = len(lines[-1].split())
    assert (n > 0)
    for i in range(len(lines)):
        lines[i] = " ".join(re.findall(".{10}", lines[i]))
    if n < 8:
        lines[-1] += ' 0'*(8-n)
    acc = np.loadtxt(lines).reshape(-1)

    if n == 8:
        return acc
    else:
        return acc[:n-8]

def duration(acc, T, N):
    t = np.linspace(0, T-T/N, N)
    a = np.trapz(acc**2, t)
    I = []
    for i in range(1, N):
        temp = np.trapz(acc[:i+1]**2, t[:i+1]) / a
        I.append(temp)
    I = np.array(I)
    d = N - (np.sum(I < 0.05) + np.sum(I > 0.95))
    return d*T/N

def calc_dur(filename):
    f = open(filename)
    lines = f.readlines()

    loc = [0]
    for i, line in enumerate(lines):
        if line[:2] == '/&':
            loc.append(i+1)
    assert (len(loc) == 4)
    loc = loc[:3]

    d = []
    for i in loc:
        assert ("points of accel data equally" in lines[i+45])
        if lines[i+11].split()[2] == "=":
            T = float(lines[i+11].split()[3])
        else:
            T = float(lines[i+11].split()[2][1:])
        N = int(lines[i+15].split()[0])
        num = N // 8 if N%8 == 0 else (N // 8) + 1
        acc = readacc(lines[i+46 : i+46+num])
        d.append(duration(acc, T, N))

    return d

def plotfigure(pga, pgv, pgd, dur, dis):
    # plot pga
    plt.figure(figsize=(9, 5))
    plt.scatter(dis, pga[:, 0], c='red', marker='o', label='channel: 360')
    plt.scatter(dis, pga[:, 1], c='blue', marker='v', label='channel: up')
    plt.scatter(dis, pga[:, 2], c='black', marker='*', label='channel: 90')
    plt.legend()
    plt.xlabel("Distance (km)")
    plt.ylabel("PGA ("+"$cm/s^2)$")
    plt.title("PGA with distance")
    plt.savefig("pga.pdf",bbox_inches='tight',dpi=300,pad_inches=0.0)
    # plt.show()

    # plot pgv
    plt.figure(figsize=(9, 5))
    plt.scatter(dis, pgv[:, 0], c='red', marker='o', label='channel: 360')
    plt.scatter(dis, pgv[:, 1], c='blue', marker='v', label='channel: up')
    plt.scatter(dis, pgv[:, 2], c='black', marker='*', label='channel: 90')
    plt.legend()
    plt.xlabel("Distance (km)")
    plt.ylabel("PGV ("+"$cm/s)$")
    plt.title("PGV with distance")
    plt.savefig("pgv.pdf",bbox_inches='tight',dpi=300,pad_inches=0.0)
    # plt.show()

    # plot PGD
    plt.figure(figsize=(9, 5))
    plt.scatter(dis, pgd[:, 0], c='red', marker='o', label='channel: 360')
    plt.scatter(dis, pgd[:, 1], c='blue', marker='v', label='channel: up')
    plt.scatter(dis, pgd[:, 2], c='black', marker='*', label='channel: 90')
    plt.legend()
    plt.xlabel("Distance (km)")
    plt.ylabel("PGD (cm)")
    plt.title("PGD with distance")
    plt.savefig("pgd.pdf",bbox_inches='tight',dpi=300,pad_inches=0.0)
    # plt.show()

    # plot duration
    plt.figure(figsize=(9, 5))
    plt.scatter(dis, dur[:, 0], c='red', marker='o', label='channel: 360')
    plt.scatter(dis, dur[:, 1], c='blue', marker='v', label='channel: up')
    plt.scatter(dis, dur[:, 2], c='black', marker='*', label='channel: 90')
    plt.legend()
    plt.xlabel("Distance (km)")
    plt.ylabel("Duration (s)")
    plt.title("Duration with distance")
    plt.savefig("duration.pdf",bbox_inches='tight',dpi=300,pad_inches=0.0)
    plt.show()


if __name__ == "__main__":
    path = 'D:\\Data\\eseis\\v2file\\'
    filenames = os.listdir(path)
    pga = np.zeros((len(filenames), 3), float)
    pgv = np.zeros((len(filenames), 3), float)
    pgd = np.zeros((len(filenames), 3), float)
    dur = np.zeros((len(filenames), 3), float)
    dis = np.array([])
    for i, filename in enumerate(filenames):
        print(filename)
        pa,pv,pd,ds = readpeak(path+filename)
        pga[i, :] = pa
        pgv[i, :] = pv
        pgd[i, :] = pd
        dis = np.append(dis, ds)
        dur[i, :] = calc_dur(path+filename)
    # np.save('pga.npy', pga)
    # np.save('pgv.npy', pgv)
    # np.save('pgd.npy', pgd)
    # np.save('dis.npy', dis)
    # np.save('dur.npy', dur)

    # pga = np.load('pga.npy')
    # pgv = np.load('pgv.npy')
    # pgd = np.load('pgd.npy')
    # dis = np.load('dis.npy')
    # dur = np.load('dur.npy')
    plotfigure(pga, pgv, pgd, dur, dis)
