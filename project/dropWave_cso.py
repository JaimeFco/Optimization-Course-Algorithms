# Cargando librerías a usar
from math import ceil, exp
import numpy as np
from random import choice, shuffle
import warnings

from matplotlib import pyplot as plt
import matplotlib.animation
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd

def animation(agents, function, lb, ub, sr=False):

    side = np.linspace(lb, ub, (ub - lb) * 5)
    X, Y = np.meshgrid(side, side)
    Z = np.array([np.array([function([X[i][j], Y[i][j]])
                            for j in range(len(X))])
                  for i in range(len(X[0]))])

    fig = plt.figure()
    plt.axes(xlim=(lb, ub), ylim=(lb, ub))
    plt.pcolormesh(X, Y, Z, shading='gouraud')
    plt.colorbar()

    x = np.array([j[0] for j in agents[0]])
    y = np.array([j[1] for j in agents[0]])
    sc = plt.scatter(x, y, color='black')

    plt.title(function.__name__, loc='left')

    def an(i):
        x = np.array([j[0] for j in agents[i]])
        y = np.array([j[1] for j in agents[i]])
        sc.set_offsets(list(zip(x, y)))
        plt.title('iteration: {}'.format(i), loc='right')

    ani = matplotlib.animation.FuncAnimation(fig, an, frames=len(agents) - 1)

    if sr:

        ani.save('result.mp4')

    plt.show()


def animation3D(agents, function, lb, ub, sr=False):

    side = np.linspace(lb, ub, 45)
    X, Y = np.meshgrid(side, side)
    zs = np.array([function([x, y]) for x, y in zip(np.ravel(X), np.ravel(Y))])
    Z = zs.reshape(X.shape)

    fig = plt.figure()

    ax = Axes3D(fig)
    surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap='jet',
                           linewidth=0, antialiased=False)
    ax.set_xlim(lb, ub)
    ax.set_ylim(lb, ub)

    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

    fig.colorbar(surf, shrink=0.5, aspect=5)

    iter = len(agents)
    n = len(agents[0])
    t = np.array([np.ones(n) * i for i in range(iter)]).flatten()
    b = []
    [[b.append(agent) for agent in epoch] for epoch in agents]
    c = [function(x) for x in b]
    a = np.asarray(b)
    df = pd.DataFrame({"time": t, "x": a[:, 0], "y": a[:, 1], "z": c})

    def update_graph(num):
        data = df[df['time'] == num]
        graph._offsets3d = (data.x, data.y, data.z)
        title.set_text(function.__name__ + " " * 45 + 'iteration: {}'.format(
            num))

    title = ax.set_title(function.__name__ + " " * 45 + 'iteration: 0')

    data = df[df['time'] == 0]
    graph = ax.scatter(data.x, data.y, data.z, color='black')

    ani = matplotlib.animation.FuncAnimation(fig, update_graph, iter,
                                             interval=50, blit=False)

    if sr:

        ani.save('result.mp4')

    plt.show()

class sw(object):

    def __init__(self):

        self.__Positions = []
        self.__Gbest = []

    def _set_Gbest(self, Gbest):
        self.__Gbest = Gbest

    def _points(self, agents):
        self.__Positions.append([list(i) for i in agents])

    def get_agents(self):
        """Returns a history of all agents of the algorithm (return type:
        list)"""

        return self.__Positions

    def get_Gbest(self):
        """Return the best position of algorithm (return type: list)"""

        return list(self.__Gbest)
    
warnings.filterwarnings("ignore")

# Chicken Swarm Optimization
class chso(sw):

    def __init__(self, n, function, lb, ub, dimension, iteration, G=5, FL=0.5):
        """
        n: number of agents
        function: test function
        lb: lower limits for plot axes
        ub: upper limits for plot axes
        dimension: space dimension
        iteration: number of iterations
        G: after what time relationship will be upgraded (default value is 5)
        FL: parameter, which means that the chick would follow its
            mother to forage for food (0 < FL < 2. Default value is 0.5)
        """

        super(chso, self).__init__()

        rn = ceil(0.15 * n)
        hn = ceil(0.7 * n)
        cn = n - rn - hn
        mn = ceil(0.2 * n)

        self.__agents = np.random.uniform(lb, ub, (n, dimension))
        pbest = self.__agents
        self._points(self.__agents)

        fitness = [function(x) for x in self.__agents]
        pfit = fitness

        Pbest = self.__agents[np.array(fitness).argmin()]
        Gbest = Pbest

        for t in range(iteration):

            if t % G == 0:

                chickens = self.__update_relationship(n, function, rn, hn,
                                                      cn, mn)
                roosters, hines, chicks = chickens

            for i in roosters:

                k = choice(roosters)
                while k == i:
                    k = choice(roosters)

                if pfit[i] <= pfit[k]:
                    sigma = 1
                else:
                    sigma = exp((pfit[k] - pfit[i]) / (abs(pfit[i]) + 0.01))

                self.__agents[i] = pbest[i] * (1 + np.random.normal(0, sigma,
                                                                    dimension))

            for i in hines:

                r1 = i[1]
                r2 = choice([choice(roosters), choice(hines)[0]])
                while r2 == r1:
                    r2 = choice([choice(roosters), choice(hines)[0]])

                s1 = exp((pfit[i[0]] - pfit[r1]) / (abs(pfit[i[0]]) + 0.01))

                try:
                    s2 = exp(pfit[r2] - pfit[i[0]])
                except OverflowError:
                    s2 = float('inf')

                rand1 = np.random.random((1, dimension))[0]
                rand2 = np.random.random((1, dimension))[0]

                self.__agents[i[0]] = pbest[i[0]] + s1 * rand1 * (
                    pbest[r1] - pbest[i[0]]) + s2 * rand2 * (
                    pbest[r2] - pbest[i[0]])

            for i in chicks:
                self.__agents[i[0]] = pbest[i[0]] * FL * (pbest[i[1]] -
                                                          pbest[i[0]])

            self.__kill(n, function, lb, ub, dimension)

            self.__agents = np.clip(self.__agents, lb, ub)
            self._points(self.__agents)

            fitness = [function(x) for x in self.__agents]

            for i in range(n):
                if fitness[i] < pfit[i]:
                    pfit[i] = fitness[i]
                    pbest[i] = self.__agents[i]

            Pbest = self.__agents[np.array(fitness).argmin()]
            if function(Pbest) < function(Gbest):
                Gbest = Pbest

        self._set_Gbest(Gbest)

    def __update_relationship(self, n, function, rn, hn, cn, mn):

        fitness = [(function(self.__agents[i]), i) for i in range(n)]
        fitness.sort()

        chickens = [i[1] for i in fitness]
        roosters = chickens[:rn]
        hines = chickens[rn:-cn]
        chicks = chickens[-cn:]

        shuffle(hines)
        mothers = hines[:mn]

        for i in range(cn):
            chicks[i] = chicks[i], choice(mothers)

        for i in range(hn):
            hines[i] = hines[i], choice(roosters)

        return roosters, hines, chicks

    def __kill(self, n, function, lb, ub, dimension):

        for i in range(n):

            fit = None

            try:
                fit = function(self.__agents[i])
            except OverflowError:
                for j in range(dimension):
                    self.__agents[i][j] = round(self.__agents[i][j])

            if str(fit) == 'nan':
                self.__agents[i] = np.random.uniform(lb, ub, (1, dimension))
def f(x):
    frac1 = 1. + np.cos(12.*np.sqrt(x[0]**2+x[1]**2))
    frac2 = 0.5*(x[0]**2+x[1]**2) + 2.0
    return (-frac1/frac2)


opt = chso(32, f, -5., 5., 2, 200, G=5, FL=0.4)
print(opt.get_Gbest())
animation(opt.get_agents(), f, -5, 5)
