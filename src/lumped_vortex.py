import numpy as np
import matplotlib.pyplot as plt

def vor2D(gamma, x, z, x0, z0):
    u = (gamma/(2.0*np.pi))*(z-z0)/((x-x0)**2 + (z - z0)**2)
    w = -(gamma/(2.0*np.pi))*(x-x0)/((x-x0)**2 + (z - z0)**2)

    return u, w



if __name__ == "__main__":

    c = 1.0
    v_inf = 1.0
    alpha = np.radians(5.0)

    plt.figure()

    num = [2,8,24,96,384]
    for nPanels in num:

        #nPanels = 2

        cp = np.zeros((nPanels,2))
        r_0 = np.zeros((nPanels,2))

        for i in range(nPanels):
            cp[i,0] = c/nPanels * (i + 0.75)
            r_0[i,0] = c/nPanels * (i + 0.25)

        print(cp)
        print(r_0)
    
        #x0 = np.array([c/8.0, 5*c/8.0])
        #xc = np.array([3*c/8.0, 7*c/8.0])
        #z0 = np.zeros(2)
        #zc = np.zeros(2)

        n = np.array([0.0,1.0])

        a = np.zeros((nPanels, nPanels))
        for i in range(nPanels):
            for j in range(nPanels):
                u, w = vor2D(1.0, cp[i,0], cp[i,1], r_0[j,0], r_0[j, 1])
                a[i,j] = np.dot([u,w], n)

        print(a)

        V_inf = v_inf*np.array([np.cos(alpha), np.sin(alpha)])
        b = np.zeros(nPanels)
        b[:] = np.dot(-V_inf, n)

        print(b)

        x = np.linalg.solve(a, b)

        print(x)

        plt.plot(r_0[:,0], x)

    plt.xlim(0,1)
    plt.show()
