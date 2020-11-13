import numpy as np
import matplotlib.pyplot as plt

x = np.array([6.0000, 10.1333, 14.2667, 18.4000, 22.5333, 26.6667]).reshape(6,1)
t = np.array([3.4935, 4.2853, 5.1374, 5.8181, 6.8632, 8.1841]).reshape(6,1)

t_m = np.array([3.3540,4.2645,5.1750,6.0855,6.9960,7.9065]).reshape(6,1)
r = t - t_m

plt.figure(figsize=(15,5))
plt.subplot(1, 2, 1)
plt.plot(x, t, 'ro', markerfacecolor='none')
plt.plot(x, t_m, 'b')
plt.xlabel("Distances/km", fontsize=16)
plt.ylabel("The first arrival times/s", fontsize=16)
plt.title("The data and fitted model", fontsize=16)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.legend(["data", "the fitted model"], fontsize=15)

plt.subplot(1, 2, 2)
plt.plot(x, r, '-bo', markerfacecolor='none')
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.xlabel("Distances/km", fontsize=16);
plt.ylabel("The differece of times/s", fontsize=16)
plt.title("The residual", fontsize=16)

plt.savefig('hw3e1a.pdf', bbox_inches='tight',dpi=600, pad_inches=0.0)
plt.show()


