import numpy as np
import json
import sys
import datetime


class Geometry:

    """Airfoil Geometry class
    
    Parameters
    ----------
    file : string
        Airfoil Geometry file

    """


    def __init__(self, **kwargs):
        self.file = kwargs["file"]
        #self.__series = kwargs["series"]

        self._read_airfoil_file()

        self._calc_geom()


    def _read_airfoil_file(self):

        print('Reading airfoil...', end='')
        self.nodes = np.loadtxt(self.file, dtype=float)
        self.N = self.nodes.shape[0]
        print(' Done. Found ', self.N, ' nodes.')


    def _calc_geom(self):
        print('Calculating Airfoil Geometry...', end='')
        self.cp = 0.5*(self.nodes[1:] + self.nodes[:-1])
        self.dx = (self.nodes[1:,0] - self.nodes[:-1,0])
        self.dy = (self.nodes[1:,1] - self.nodes[:-1,1])
        self.l = np.sqrt(self.dx**2 + self.dy**2)
        print(' Done.')



class Solver:
    """Vortex Solver Class
    
    Parameters
    ----------
    airfoil : Geometry
    
    """

    def __init__(self, **kwargs):

        self._airfoil = kwargs["airfoil"]
        self._alpha = kwargs["alpha"]
        self._V_inf = kwargs["V_inf"]
        self._N = self._airfoil.N

        self.A = np.zeros((self._N,self._N))
        self.b = np.zeros(self._N)


    def calc_influence(self):

        print('Calculating influences...', end='')

        for i in range(self._N-1):
            for j in range(self._N-1):
                xi = (self._airfoil.dx[j] * (self._airfoil.cp[i,0] - self._airfoil.nodes[j,0]) + self._airfoil.dy[j] * (self._airfoil.cp[i,1] - self._airfoil.nodes[j,1]))/self._airfoil.l[j]
                eta = (-self._airfoil.dy[j] * (self._airfoil.cp[i,0] - self._airfoil.nodes[j,0]) + self._airfoil.dx[j] * (self._airfoil.cp[i,1] - self._airfoil.nodes[j,1]))/self._airfoil.l[j]

                phi = np.arctan2(eta*self._airfoil.l[j], eta**2 + xi**2 - xi*self._airfoil.l[j])
                psi = 0.5*np.log((xi**2 + eta**2)/((xi-self._airfoil.l[j])**2 + eta**2))

                P1 = np.array([[self._airfoil.dx[j], -self._airfoil.dy[j]], [self._airfoil.dy[j], self._airfoil.dx[j]]])
                P2 = np.array([[(self._airfoil.l[j]-xi)*phi + eta*psi, xi*phi - eta*psi],
                            [eta*phi - (self._airfoil.l[j]-xi)*psi - self._airfoil.l[j], -eta*phi - xi*psi + self._airfoil.l[j]]])
            
                P = (1/(2*np.pi*self._airfoil.l[j]**2))*np.matmul(P1,P2)

                self.A[i,j] += self._airfoil.dx[i]*P[1,0]/self._airfoil.l[i] - self._airfoil.dy[i]*P[0,0]/self._airfoil.l[i]
                self.A[i,j+1] += self._airfoil.dx[i]*P[1,1]/self._airfoil.l[i] - self._airfoil.dy[i]*P[0,1]/self._airfoil.l[i]

        
        self.A[-1,0] = 1.0
        self.A[-1,-1] = 1.0

        self.b[:-1] = self._V_inf*(self._airfoil.dy*np.cos(np.deg2rad(self._alpha)) - self._airfoil.dx*np.sin(np.deg2rad(self._alpha)))/self._airfoil.l
        self.b[-1] = 0.0

        print('Done.')


    def solve(self):
        
        print('Solving...', end='')
        self.gamma = np.linalg.solve(self.A,self.b)
        print('Done.')
        

    def calc_forces(self):
        print('Caculating forces...', end='')
        self.C_L = 0.0
        self.C_m_le = 0.0
        for i in range(self._N-1):
            self.C_L += self._airfoil.l[i]*(self.gamma[i] + self.gamma[i+1])/self._V_inf
            self.C_m_le += self._airfoil.l[i]*(((2*self._airfoil.nodes[i,0]*self.gamma[i] + self._airfoil.nodes[i,0]*self.gamma[i+1] + self._airfoil.nodes[i+1,0]*self.gamma[i] + 2*self._airfoil.nodes[i+1,0]*self.gamma[i+1])/self._V_inf)*np.cos(np.deg2rad(self._alpha)) 
                                            + ((2*self._airfoil.nodes[i,1]*self.gamma[i] + self._airfoil.nodes[i,1]*self.gamma[i+1] + self._airfoil.nodes[i+1,1]*self.gamma[i] + 2*self._airfoil.nodes[i+1,1]*self.gamma[i+1])/self._V_inf)*np.sin(np.deg2rad(self._alpha)))
        self.C_m_le = -1/3 * self.C_m_le

        self.C_m_qtr = self.C_m_le + 0.25*self.C_L*np.cos(np.deg2rad(self._alpha))
        print('Done.')

        print("       C_L: ",self.C_L)
        print("    C_m_le: ",self.C_m_le)
        print("   C_m_qtr: ",self.C_m_qtr)

        return self.C_L, self.C_m_le, self.C_m_qtr

    def write_report(self, report_filename):

        report_dict = {
            "info": {
                "generated_by": "VORPAN 2024 Jeremiah Goates"
            },
            "input": {
                "geometry": self._airfoil.file,
                "alpha[deg]" : self._alpha,
                "freestream_velocity": self._V_inf
            },
            "forces": {
                "C_L": self.C_L,
                "C_m_le": self.C_m_le,
                "C_m_qtr": self.C_m_qtr,
            }

        }

        with open(report_filename, 'w') as report_handle:
            json.dump(report_dict, report_handle, indent=4)


def _parse_json(input_filename):
    
    with open(input_filename) as json_file_handle:
        input_dict = json.load(json_file_handle)

    geometry = input_dict["geometry"]
    V_inf = input_dict["freestream_velocity"]
    alpha_deg = input_dict["alpha[deg]"]

    return geometry, V_inf, alpha_deg



if __name__ == "__main__":



    # Print logo

    print(' _______________________________________________')
    print()
    print(' | | \ \     / _ \   _ \   _ \   \     \  |  | |  ')
    print(' | |  \ \   / |   | |   | |   | _ \     \ |  | |  ')
    print(' | |   \ \ /  |   | __ <  ___/ ___ \  |\  |  | |  ')
    print(' | |    \_/  \___/ _| \_\_|  _/    _\_| \_|  | |  ')
    print(' _______________________________________________')
    print('                 Jeremiah Goates                  ')
    print('                    (c) 2024                 ')

    # Check for interactive mode
    if "-i" in sys.argv:
        pass
    else:
        # Get input filename from command line arguments
        input_filename = sys.argv[-1]

        # Check for valid input
        if ".json" not in input_filename:
            raise IOError("Please specify a .json input file (got {0}).".format(input_filename))

    geometry, V_inf, alpha_deg = _parse_json(input_filename)


    myAirfoil = Geometry(file=geometry)
    mySolver = Solver(airfoil=myAirfoil, alpha=alpha_deg, V_inf=V_inf)

    mySolver.calc_influence()

    mySolver.solve()

    mySolver.calc_forces()

    mySolver.write_report("dev/reports/report.json")




