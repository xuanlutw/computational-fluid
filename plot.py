import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator

def read_file(mu, rho, v, r):
    filename = "%.2lf-%.2lf-%.2lf-%d" % (mu, rho, v, r)
    f = open('result/' + filename)
    linelist = f.readlines()
    f.close()
    return linelist

def plot_force(mu, rho, v, r):
    linelist = read_file(mu, rho, v, r)
    res = np.zeros(50000)
    for i in range(0, 49999):
        res[i] = float(linelist[2 * i + 1].split(' ')[0])
    res[49999] = float(linelist[-1].split(' ')[0])

    fig, ax = plt.subplots()
    ax.plot(range(0, 50000), res)

    ax.set(xlabel = 'iter', ylabel = 'force', 
           title = 'mu = %.2lf, rho = %.2lf, v0 = %.2lf, r = %.2d' % (mu, rho, v, r))
    ax.grid()

    #fig.savefig("test.png")
    plt.show()

def plot_pressure(mu, rho, v, r, dim, frame):
    linelist = read_file(mu, rho, v, r)
    pressure = np.zeros((61, 61, 61))
    for i in range(0, 61):
        for j in range(0, 61):
            num = linelist[100001 + i * 61 + j].split(' ')
            for k in range(0, 61):
                pressure[i][j][k] = float(num[k])

    fig, ax = plt.subplots()
    if dim == 0:
        p = pressure[frame][:][:]
    if dim == 1:
        p = pressure[:][frame][:]
    if dim == 2:
        p = pressure[:][:][frame]

    levels = MaxNLocator(nbins = 30).tick_values(p.min(), p.max())
    cmap = plt.get_cmap('bwr')
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
    im = ax.pcolormesh(p, cmap=cmap, norm=norm)
    fig.colorbar(im, ax = ax)
    ax.set_title('Pressure along cut at axis %d frame %d' % (dim, frame))
    plt.show()

def plot_velocity(mu, rho, v, r, dim, frame):
    linelist = read_file(mu, rho, v, r)
    vx = np.zeros((61, 61, 61))
    vy = np.zeros((61, 61, 61))
    vz = np.zeros((61, 61, 61))
    for i in range(0, 61):
        for j in range(0, 61):
            num = linelist[100001 + 61 * 61  + i * 61 + j].split(' ')
            for k in range(0, 61):
                vx[i][j][k] = float(num[k])
    for i in range(0, 61):
        for j in range(0, 61):
            num = linelist[100001 + 2 * 61 * 61  + i * 61 + j].split(' ')
            for k in range(0, 61):
                vy[i][j][k] = float(num[k])
    for i in range(0, 61):
        for j in range(0, 61):
            num = linelist[100001 + 3 * 61 * 61  + i * 61 + j].split(' ')
            for k in range(0, 61):
                vz[i][j][k] = float(num[k])

    fig, ax = plt.subplots()
    if dim == 0:
        v1 = vy[frame][:][:]
        v2 = vz[frame][:][:]
    if dim == 1:
        v1 = vx[:][frame][:]
        v2 = vz[:][frame][:]
    if dim == 2:
        v1 = vx[:][:][frame][:][:]
        v2 = vy[:][:][frame]

    #levels = MaxNLocator(nbins = 30).tick_values(p.min(), p.max())
    #cmap = plt.get_cmap('bwr')
    #norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
    im = ax.streamplot(np.asarray(range(0, 61)), np.asarray(range(0, 61)), v1, v2)
    #fig.colorbar(im, ax = ax)
    ax.set_title('Pressure along cut at axis %d frame %d' % (dim, frame))
    plt.show()

#plot_force(5, 1, 1, 10)
#plot_pressure(5, 1, 1, 5, 0, 11)
plot_velocity(5, 2, 1, 5, 2, 31)


