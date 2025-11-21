import numpy as np
import json
import matplotlib.pyplot as plt
import sys
import time

# Geometry class that holds all information relevant to the specific geometry including velocity functions
class Geometry:
    def __init__(self, file_name):
        with open(file_name, 'r') as json_geometry:
            input_vals = json.load(json_geometry)

        self.__cylinder_radius = input_vals["geometry"]["cylinder_radius"]

        self.__x_start = input_vals["plot_options"]["x_start"]
        self.__x_lower_limit = input_vals["plot_options"]["x_lower_limit"]
        self.__x_upper_limit = input_vals["plot_options"]["x_upper_limit"]
        self.__delta_s = input_vals["plot_options"]["delta_s"]
        self.__n_lines = input_vals["plot_options"]["n_lines"]
        self.__delta_y = input_vals["plot_options"]["delta_y"]

        self.__freestream_velocity = input_vals["operating"]["freestream_velocity"]
        self.__alpha_deg = input_vals["operating"]["alpha[deg]"]
        self.__vortex_strength = input_vals["operating"]["vortex_strength"]

    def cylinder_radius(self):
        return self.__cylinder_radius
    
    def x_start(self):
        return self.__x_start
    
    def x_lower_limit(self):
        return self.__x_lower_limit
    
    def x_upper_limit(self):
        return self.__x_upper_limit
    
    def delta_s(self):
        return self.__delta_s
    
    def n_lines(self):
        return self.__n_lines
    
    def delta_y(self):
        return self.__delta_y
    
    def freestream_velocity(self):
        return self.__freestream_velocity
    
    def alpha_deg(self):
        return self.__alpha_deg
    
    def vortex_strength(self):
        return self.__vortex_strength
    
    def camber(self, x):
        self.__camber = [None,None]

        if x <= self.__cylinder_radius and x >= -self.__cylinder_radius:
            self.__camber = [x, 0]
    
        return self.__camber
    
    def upper_surface(self, x):
        self.__upper_surface = [None,None]

        if x <= self.__cylinder_radius and x >= -self.__cylinder_radius:
            self.__upper_surface = np.array([x, (self.__cylinder_radius * np.sin(np.arccos(x/self.__cylinder_radius)))])

        return self.__upper_surface
    
    def lower_surface(self, x):
        self.__lower_surface = [None,None]

        if x <= self.__cylinder_radius and x >= -self.__cylinder_radius:
            self.__lower_surface = np.array([x, -(self.__cylinder_radius * np.sin(np.arccos(x/self.__cylinder_radius)))])

        return self.__lower_surface

    def velocity(self, point):
        s = [np.sqrt(point[0]**2 + point[1]**2), np.arctan2(point[1],point[0])]

        v_p = [None, None]

        radius = self.__cylinder_radius
        V_inf = self.__freestream_velocity
        alpha_deg = self.__alpha_deg
        vortex_strength = self.__vortex_strength
        alpha = alpha_deg*np.pi/180

        v_p = [ V_inf * ( 1 - ((radius**2)/(s[0]**2)) ) * np.cos(s[1] - alpha), -(V_inf * ( 1 + ((radius**2)/(s[0]**2)) ) * np.sin(s[1] - alpha) + vortex_strength/(2*np.pi*s[0]))]

        v_c = [ v_p[0]*(np.cos(s[1])) - (np.sin(s[1]))*v_p[1] , v_p[0]*(np.sin(s[1])) + (np.cos(s[1]))*v_p[1] ]

        return np.array(v_c)

    def surface_tangential_velocity(self,x):
    
        radius = self.__cylinder_radius
        upper_surface = self.upper_surface(x)
        lower_surface = self.lower_surface(x)

        if x >= -radius and x<= radius:
            upper_tangential_velocity = self.velocity(upper_surface)
            lower_tangential_velocity = self.velocity(lower_surface)

        else:
            upper_tangential_velocity = [None, None]
            lower_tangential_velocity = [None, None]

        return upper_tangential_velocity, lower_tangential_velocity
    
    def geo_plot(self):

        print("Plotting geometry... Done")

        x_array = []
        step = 0.001
        x = self.__x_lower_limit
        while x <= self.__x_upper_limit:
            x_array.append(x)
            x += step
        
        c_plot_x = []
        c_plot_y = []

        u_plot_x = []
        u_plot_y = []

        l_plot_x = []
        l_plot_y = []

        for i in x_array:
            camber = self.camber(i)
            upper_surface = self.upper_surface(i)
            lower_surface = self.lower_surface(i)

            c_plot_x.append(camber[0])
            c_plot_y.append(camber[1])

            u_plot_x.append(upper_surface[0])
            u_plot_y.append(upper_surface[1])

            l_plot_x.append(lower_surface[0])
            l_plot_y.append(lower_surface[1])

        plt.ylabel("y")
        plt.xlabel("x")
        plt.title("Streamlines")
        plt.axis([self.__x_lower_limit,self.__x_upper_limit, self.__x_lower_limit, self.__x_upper_limit])
        plt.plot(c_plot_x, c_plot_y, color = "red")
        plt.plot(u_plot_x, u_plot_y, color = "blue")
        plt.plot(l_plot_x, l_plot_y, color = "blue")
        plt.gca().set_aspect("equal")

# Class with all necessary functions for calculating the streamlines over the Geometry
class Streamlines(Geometry):
    def __init__(self, file_name):
        super().__init__(file_name)
        # Lower computation cost
        self.__one_sixth = 1/6

    # Give the magnitude of a vector
    def magnitude(self, vector):
    
        magnitude = np.sqrt(vector[0]**2 + vector[1]**2)

        return magnitude
    
    # RK4 method for integrating over the velocity field
    def runge_kutta(self, point, forward=True):
        if forward:
            delta_s = super().delta_s()
        else:
            delta_s = -super().delta_s()
        k_1 = super().velocity(point)/self.magnitude(super().velocity(point))
        k_2 = super().velocity(point + .5*delta_s*k_1)/self.magnitude(super().velocity(point + .5*delta_s*k_1))
        k_3 = super().velocity(point + .5*delta_s*k_2)/self.magnitude(super().velocity(point + .5*delta_s*k_2))
        k_4 = super().velocity(point + delta_s*k_3)/self.magnitude(super().velocity(point + delta_s*k_3))
        point = point + delta_s*self.__one_sixth*(k_1 + 2*k_2 + 2*k_3 + k_4)
        return point

    # Calculate the stagnation points on the geometry
    def stagnation(self):

        print("Calculating stagnation points...", end=" ")

        forward_stagnation_point = [None, None]
        aft_stagnation_point = [None, None]

        radius = super().cylinder_radius()
        x = 0
        omega = np.pi/radius

        # Calculate Forward stagnation point
        while x >= -radius:
            # Variable step size
            step = 0.001 * (np.cos(omega * x) + 1.1)

            upper_tangential_velocity, lower_tangential_velocity = super().surface_tangential_velocity(x)
            upper_surface = super().upper_surface(x)
            lower_surface = super().lower_surface(x)
        
            if self.magnitude(upper_tangential_velocity) < 0.01:
                forward_stagnation_point = upper_surface
                break
            elif self.magnitude(lower_tangential_velocity) < 0.01:
                forward_stagnation_point = lower_surface  
                break     
            else:
                x -= step

        x = 0

        # Calculate Aft stagnation point
        while x <= radius:
            step = 0.001 * (np.cos(omega * x) + 1.1)
            upper_tangential_velocity, lower_tangential_velocity = self.surface_tangential_velocity(x)
            upper_surface = super().upper_surface(x)
            lower_surface = super().lower_surface(x)
            if self.magnitude(upper_tangential_velocity) <  0.01:
                aft_stagnation_point = upper_surface
                break
            elif self.magnitude(lower_tangential_velocity) < 0.01:
                aft_stagnation_point = lower_surface
                break
            else:
                x += step
        
        # Check if points were found
        if forward_stagnation_point[0] == None:
            print("ERROR: Unable to find stagnation points....Quitting")
            sys.exit(0)
        else:
            print("Done")
        return np.array(forward_stagnation_point), np.array(aft_stagnation_point)
    
    # Calculates the stagnation streamlines
    def stagnation_streamlines(self):
        
        forward_stagnation_point, aft_stagnation_point = self.stagnation()
        forward_streamline_array = [forward_stagnation_point]
        aft_streamline_array = [aft_stagnation_point]

        delta_s = super().delta_s()

        start = forward_stagnation_point
        x_lower_limit = super().x_lower_limit()

        point = np.array([start[0] + (delta_s*start[0]/self.magnitude(start)), start[1] + (delta_s*start[1]/self.magnitude(start))])
        while point[0] >= x_lower_limit:
            point = self.runge_kutta(point, False)
            forward_streamline_array.insert(0,point)

        start = aft_stagnation_point
        x_upper_limit = super().x_upper_limit()

        point = np.array([start[0] + (delta_s*start[0]/self.magnitude(start)), start[1] + (delta_s*start[1]/self.magnitude(start))])
        while point[0] <= x_upper_limit:
            point = self.runge_kutta(point, True)
            aft_streamline_array.append(point)


        return np.array(forward_streamline_array), np.array(aft_streamline_array)
    
    # Calculate streamline up and downstream of a specified starting location
    def streamline(self, start):

        streamline_array = [start]

        x_lower_limit = super().x_lower_limit()
        x_upper_limit = super().x_upper_limit()

        point = start
        while point[0] <= x_upper_limit:
                point = self.runge_kutta(point, True)
                streamline_array.append(point)
        
        point = start
        while point[0] >= x_lower_limit:
                point = self.runge_kutta(point,False)
                streamline_array.insert(0,point)
        
        return np.array(streamline_array)

    # Plots a collection of streamlines
    def plot(self):
        super().geo_plot()
        forward_stagnation_line , aft_stagnation_line= self.stagnation_streamlines()

        
        plt.plot(forward_stagnation_line[:,0], forward_stagnation_line[:,1], color = "black")
        plt.plot(aft_stagnation_line[:,0], aft_stagnation_line[:,1], color = "black")

        start = np.array([super().x_start(),forward_stagnation_line[0,1]])

        upper_line_count = 0
        while upper_line_count < super().n_lines():

            start[1] += super().delta_y()
            streamline_array = self.streamline(start)
            plt.plot(streamline_array[:,0], streamline_array[:,1], color = "black")
            upper_line_count += 1

        start = np.array([super().x_start(),forward_stagnation_line[0,1]])
        lower_line_count = 0
        while lower_line_count < super().n_lines():
            start[1] -= super().delta_y()
            streamline_array = self.streamline(start)
            plt.plot(streamline_array[:,0], streamline_array[:,1], color = "black")
            lower_line_count += 1

        print()
        
def main():

    print()
    start_time = time.time()
    stream = Streamlines("geometry.json")
    stream.plot()    
    print()
    print("Program executed in %6.3f seconds" % (time.time() - start_time))
    plt.show()


if __name__ == '__main__':
    main()