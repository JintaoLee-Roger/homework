import numpy as np

if __name__ == "__main__":
    crustroot = '/home/roger/.gmt/GMTDB/crust1.0/'
    thkdir = crustroot + 'crust/out/'
    # 读取地壳厚度
    crsthk = np.loadtxt(thkdir+'crsthk.xyz')
    reg = crsthk[:,:2]
    crsthk = crsthk[:, 2]
    # 读取地壳密度
    ro = []
    thk = []
    ro_crust = np.zeros(crsthk.shape, np.float32)
    for i in range(9):
        ro.append(np.loadtxt(thkdir+'xyz-ro'+str(i+1))[:, 2])
        if i < 8:
            thk.append(np.loadtxt(thkdir+'xyz-th'+str(i+1))[:, 2])
        if i > 0 and i < 8:
            ro_crust += ro[i]*thk[i]
    # 地壳密度是第2层到第8层的平均
    density_crust = ro_crust / crsthk
    density_mantle = ro[8].mean()

    # 划分陆壳和洋壳 以地壳厚度小于12Km的区域为洋壳
    crsthk_continental = np.zeros(crsthk.shape, np.float32)
    crsthk_ocean = np.zeros(crsthk.shape, np.float32)
    crsthk_continental[crsthk > 12] = crsthk[crsthk > 12]
    crsthk_ocean[crsthk <= 12] = crsthk[crsthk <= 12]

    # 计算大陆和大洋的海拔高度
    elevation_continental = crsthk_continental*(1 - density_crust/density_mantle)
    elevation_ocean = (density_mantle - density_crust.mean())/(density_mantle - ro[0].mean())* \
                    crsthk_continental - \
                    (density_mantle - density_crust)/(density_mantle - ro[0])*crsthk_ocean

    # 合并
    elevation_isostasy = np.zeros(crsthk.shape, np.float32)
    elevation_isostasy[crsthk > 12] = elevation_continental[crsthk > 12]
    elevation_isostasy[crsthk <= 12] = elevation_ocean[crsthk <= 12]

    # etop1 海拔 (水层的下边界就是etop1的海拔)
    elevation_etop1 = np.loadtxt(thkdir+'xyz-bd2')[:, 2]

    # diffence
    diffence = elevation_isostasy - elevation_etop1
    print(diffence.max(), diffence.min())

    # save
    elevation_isostasy = elevation_isostasy[:, np.newaxis]
    diffence = diffence[:, np.newaxis]
    elevation_isostasy = np.concatenate((reg, elevation_isostasy), 1)
    diffence = np.concatenate((reg, diffence), 1)
    np.savetxt('elevation.xyz', elevation_isostasy)
    np.savetxt('diffence.xyz', diffence)
