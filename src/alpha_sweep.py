import matplotlib.pyplot as plt
import numpy as np
import os
import subprocess as sp
import json
import sys

def write_input_file(input_dict, input_filename):
    with open(input_filename, 'w') as input_handle:
        json.dump(input_dict, input_handle, indent=4)

def run_vorpan(input_filename, delete_input=True, run=True):
    """Runs VORPAN with the given input and returns the report if VORPAN generated one."""

    # Run
    if run:
        sp.call(["python","src/main.py", input_filename])
    
    with open("dev/reports/report.json") as report_handle:
        report = json.load(report_handle)

    if delete_input:
        os.remove(input_filename)

    return report

def get_data_from_csv(csv_file, remove_csv=True):
    """Pulls in data from a csv file.
    
    Parameters
    ----------
    csv_file : str
        File to pull the data in from.

    remove_csv : bool, optional
        Whether to delete the file the data are pulled from. Defaults to True.

    Returns
    -------
    column_headers : ndarray
        Array of column headers.

    data : ndarray
        Data array. Columns correspond to column headers.
    """

    # Read into arrays
    with open(csv_file, 'r') as data_file:
        column_headers = [x.strip().replace('"', '') for x in data_file.readline().split(',')]
    cell_data = np.genfromtxt(csv_file, delimiter=',', skip_header=1)

    # Remove csv
    if remove_csv:
        os.remove(csv_file)

    return column_headers, cell_data

def get_data_column_from_array(headers, data, col_des):
    """Returns the data vector from data with the header col_des.
    
    Parameters
    ----------
    headers : list of str
        Lsit of column headers.

    data : ndarray
        Data array.

    col_des : str
        Header of column to pull data from.

    Returns
    -------
    vector
        Data in desired column.

    Raises
    ------
    ValueError
        If col_des is not found in headers.
    """

    # Get column index
    try:
        ind = headers.index(col_des)
    except ValueError:
        raise ValueError("Desired column header not found. Headers available: {0}".format(headers))

    return data[:,ind].flatten()

if __name__ == "__main__":

    
    alpha_deg = np.linspace(-10,10,21)
    geometries = ["2412", "2421", "0015"]
    alpha_L_0 = [-2.077, -2.077, 0.0]
    #experimental = 
    #geometries = ["dev/input/2412_200.txt","dev/input/2421_200.txt", "dev/input/0015_200.txt"]

    for j, geometry in enumerate(geometries):

        geometry_file = "dev/input/"+geometry+"_200.txt"
        experimental_file = "dev/data/"+geometry+"_exp.csv"

        headers, data = get_data_from_csv(experimental_file, remove_csv=False)

        alpha_exp = get_data_column_from_array(headers, data, "alpha")
        C_L_exp = get_data_column_from_array(headers, data, "C_L")

        C_L_thin = 2*np.pi*(np.deg2rad(alpha_deg - alpha_L_0[j]))

        C_L = np.zeros_like(alpha_deg)
        C_m_le = np.zeros_like(alpha_deg)
        C_m_qtr = np.zeros_like(alpha_deg)
        
        for i, alpha in enumerate(alpha_deg):
            input_dict = {
                "geometry": geometry_file,
                "alpha[deg]": alpha,
                "freestream_velocity": 10.0
            }

            input_filename = "input.json"
            write_input_file(input_dict, input_filename)
            report = run_vorpan(input_filename)

            C_L[i] = report["forces"]["C_L"]
            C_m_le[i] = report["forces"]["C_m_le"]
            C_m_qtr[i] = report["forces"]["C_m_qtr"]


        plt.figure()
        if j == 2:
            plt.plot(alpha_exp, C_L_exp, 'ko', label='XFOIL')
        else:
            plt.plot(alpha_exp, C_L_exp, 'ko', label='Experimental')
        plt.plot(alpha_deg, C_L, 'k-', label='VORPAN')
        plt.plot(alpha_deg, C_L_thin, 'k-.', label='Thin Airfoil Theory')
        plt.legend(fontsize=6)
        plt.xlabel("$\\alpha$")
        plt.ylabel("$C_L$")
        plt.savefig(geometry+".png")
        plt.savefig(geometry+".pdf")
        #plt.show()
        plt.close()

        print("C_L: ")
        for i, alpha in enumerate(alpha_deg):
            print(C_L[i])
        print()
        print("C_m_le")
        for i, alpha in enumerate(alpha_deg):
            print(C_m_le[i])
        print()
        print("C_m_qtr")
        for i, alpha in enumerate(alpha_deg):
            print(C_m_qtr[i])
            


