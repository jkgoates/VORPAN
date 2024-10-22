import numpy as np
import matplotlib.pyplot as plt


def get_NACA_4_camber(x,c,series):
    x_mc = 0.1*series[1]*c
    y_mc = 0.01*series[0]*c

    y = np.zeros(len(x))

    for i, x_i in enumerate(x):
        if (x_i<x_mc):
            y[i] = y_mc*(2*(x_i/x_mc) - (x_i/x_mc)**2)
        else:
            y[i] = y_mc*(2*((c-x_i)/(c-x_mc)) - ((c-x_i)/(c-x_mc))**2)

    return y

def get_NACA_4_camber_der(x,c,series):
    x_mc = 0.1*series[1]*c
    y_mc = 0.01*series[0]*c

    dydx = np.zeros(len(x))

    for i, x_i in enumerate(x):
        if (x_i<x_mc):
            dydx[i] = y_mc*(2*(1/x_mc) - 2*(x_i/x_mc)*(1/x_mc))
        else:
            dydx[i] = y_mc*(2*((-1)/(c-x_mc)) - 2*((c - x_i)/(c-x_mc))*((-1.0)/(c-x_mc)))

    return dydx


def get_NACA_4_thickness(x,c,series):
    t_m = 0.1*series[2]*c + 0.01*series[3]*c

    t = np.zeros_like(x)

    t = t_m*(2.969*np.sqrt(x/c) - 1.260*(x/c) - 3.516*(x/c)**2 + 2.843*(x/c)**3 - 1.015*(x/c)**4)

    return t

def get_NACA_4_surface(x,c,series):

    y_c = get_NACA_4_camber(x,c,series)
    dy_cdx = get_NACA_4_camber_der(x,c,series)
    t = get_NACA_4_thickness(x,c,series)

    x_u = np.zeros_like(x)
    x_u = x - (t/(2*np.sqrt(1 + (dy_cdx**2))))*dy_cdx

    x_l = np.zeros_like(x)
    x_l = x + (t/(2*np.sqrt(1 + (dy_cdx**2))))*dy_cdx

    y_u = np.zeros_like(x)
    y_u = y_c + (t/(2*np.sqrt(1 + (dy_cdx**2))))

    y_l = np.zeros_like(x)
    y_l = y_c - (t/(2*np.sqrt(1 + (dy_cdx**2))))

    return x_u, x_l, y_u, y_l

def get_NACA_4_nodes(n,c,series):
    if (n%2 == 0): #even    
        dtheta = 2*np.pi/(n-1)
        x = np.zeros(n//2)
        for i in range(n//2):
            x[i] = 0.5*(1.0 - np.cos((i+0.5)*dtheta))
    else: # odd
        dtheta = 2*np.pi/(n-1)
        x = np.zeros(n//2)
        for i in range(n//2):
            x[i] = 0.5*(1.0 - np.cos((i+1)*dtheta))


    x_u, x_l, y_u, y_l = get_NACA_4_surface(x,c,series)


    nodes = np.zeros((n,2))

    if (n%2 == 0): #even    
        nodes[:n//2,0] = x_l[::-1]
        nodes[:n//2,1] = y_l[::-1]
        nodes[n//2:,0] = x_u
        nodes[n//2:,1] = y_u
    else: # odd
        nodes[:n//2,0] = x_l[::-1]
        nodes[:n//2,1] = y_l[::-1]
        nodes[(n//2)+1:,0] = x_u
        nodes[(n//2)+1:,1] = y_u
        nodes[n//2,:] = 0.0

    return nodes



if __name__ == "__main__":
    
    series = [0,0,1,0]
    c = 1.0
    n = 50 # Number of nodes


    nodes = get_NACA_4_nodes(n=n, c=c, series=series)

    plt.figure()

    #plt.plot(x,y_c, label='Camber')
    #plt.plot(x,t, label='Thickness')
    #plt.plot(x_u,y_u, label='Upper')
    #plt.plot(x_l,y_l, label='Lower')
    plt.plot(nodes[:,0],nodes[:,1],'k--')
    plt.plot(nodes[:,0],nodes[:,1],'ko')

    plt.ylim(-c/2,c/2)
    plt.xlim(0,c)
    #plt.legend()

    plt.show()