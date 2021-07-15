import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import pickle


### u

def plot_log(filename, label="", BDIM=False):
    DNS = pickle.load(open(filename, "rb"))

    plt.subplots(1,1)

    ts= 25000
    sampling = 10
    range = 50000

    t = np.arange(0,DNS.shape[0]/10**4, 10**-4)
    # print(t)

    if not BDIM:
        xmin=0
        xmax=1
    else:
        xmin=-0.0078125
        xmax=+1.0078125

    X = np.linspace(xmin,xmax, DNS.shape[1])
    Y = t[ts:ts+range:sampling]
    Z = np.flip(DNS[ts:ts+range:sampling], 0)

    pcm = plt.pcolormesh(X, Y, Z,
                    norm=colors.SymLogNorm(linthresh=0.03, linscale=0.03,
                                              vmin=Z.min(), vmax=Z.max(), base=10),
                    cmap='coolwarm', shading='auto')

    cbar = plt.colorbar(pcm, extend='max')
    cbar.set_label(label, rotation=270)

    plt.ylabel("t")
    plt.xlabel("x")

    plt.show()
    plt.savefig("heatmap.png")

    del DNS

def plot_loglog(filename, label="", BDIM=False):
    DNS = pickle.load(open(filename, "rb"))

    plt.subplots(1,1)

    ts= 25000
    sampling = 1
    range = 50000

    t = np.arange(0,DNS.shape[0]/10**4, 10**-4)
    # print(t)

    if not BDIM:
        xmin=0
        xmax=1
    else:
        xmin=-0.0078125
        xmax=+1.0078125

    X = np.linspace(xmin,xmax, DNS.shape[1])
    Y = t[ts:ts+range:sampling]
    Z = np.flip(DNS[ts:ts+range:sampling], 0)

    pcm = plt.pcolormesh(X, Y, Z,
                    norm=colors.SymLogNorm(linthresh=0.03, linscale=0.03,
                                              vmin=Z.min(), vmax=Z.max(), base=10),
                    cmap='coolwarm', shading='auto')
    # plt.xticks()
    plt.xscale("log")
    plt.plot([1,1], [Y.min(), Y.max()], color='yellow', linewidth=1)

    cbar = plt.colorbar(pcm, extend='max')
    cbar.set_label(label, rotation=270)

    plt.ylabel("t")
    plt.xlabel("x")
    plt.xlim(xmax=xmax, xmin=0.5)
    plt.ylim(Y.min(), Y.max())

    plt.show()
    plt.savefig("heatmap.png")

    del DNS


def plot(filename, label="", BDIM=False):
    DNS = pickle.load(open(filename, "rb"))

    plt.subplots(1,1)

    ts= 25000
    sampling = 10
    range = 50000

    t = np.arange(0,DNS.shape[0]/10**4, 10**-4)
    # print(t)

    if not BDIM:
        xmin=0
        xmax=1
    else:
        xmin=-0.0078125
        xmax=+1.0078125
    
    plt.imshow(DNS[ts:ts+range:sampling], extent = [xmin,xmax,t[ts], min(t[ts+range-1], max(t))], aspect="auto", cmap="coolwarm")

    plt.ylabel("t")
    plt.xlabel("x")
    cbar = plt.colorbar()
    cbar.set_label(label, rotation=270)
    plt.show()
    plt.savefig("heatmap.png")

    del DNS

# plot("flowdata/channel_2048DNS.p", label="u")

plot("flowdata/perfect.p", label="Perfect closure", BDIM=True)

# plot_log("flowdata/tau.p", label="LES subgrid", BDIM=True)

# plot_loglog("flowdata/diff.p", label="Perfect - LES subgrid", BDIM=True)

# plot_loglog("flowdata/ImLES_u.p", label="ImLES u", BDIM=True)

# plot_loglog("flowdata/tau.p", label="LES subgrid", BDIM=True)

plot("flowdata/channel_ImLES_Dij_0.0078125.p", label="Dij", BDIM=True)


