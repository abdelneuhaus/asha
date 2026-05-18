import matplotlib.pyplot as plt
import numpy as np

data = [1817,1818,1819,1820,1821,1822,1828,1829,1830,1831,1832,1840,1841,1842,1843,1844,1845,1847,1848]

def plot_binary_on_off(df, frames_num):
    # df = df.to_numpy()
    tmp = []
    for i in range(1800, frames_num+1):
        if i in df:
            tmp.append(1)
        else:
            tmp.append(0)
    bbin = tmp
    bxaxis = np.arange(1800, 1861)
    byaxis = np.array(bbin)
    plt.step(bxaxis, byaxis, color="green")
    plt.yticks([0,1])
    plt.show()

plot_binary_on_off(data, 1, 1860)