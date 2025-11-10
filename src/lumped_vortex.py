import numpy as np

def vor2D(gamma, x, z, x0, z0):
    u = (gamma/(2.0*np.pi))*(z-z0)/((x-x0)**2 + (z - z0)**2)
    w = -(gamma/(2.0*np.pi))*(x-x0)/((x-x0)**2 + (z - z0)**2)

    return u, w



if __name__ == "__main__":

    c = 1.0
    v_inf = 1.0
    alpha = np.radians(5.0)
    
    x0 = np.array([c/8.0, 5*c/8.0])
    xc = np.array([3*c/8.0, 7*c/8.0])
    z0 = np.zeros(2)
    zc = np.zeros(2)

    n = np.array([0.0,1.0])

    a = np.zeros((2,2))
    for i in range(2):
        for j in range(2):
            u, w = vor2D(1.0, xc[i], zc[i], x0[j], z0[j])
            a[i,j] = np.dot([u,w], n)

    print(a)

    V_inf = v_inf*np.array([np.cos(alpha), np.sin(alpha)])
    b = np.zeros(2)
    b[:] = np.dot(-V_inf, n)

    print(b)

    x = np.linalg.solve(a, b)

    print(x)