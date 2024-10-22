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


if __name__ == "__main__":

    
    alpha_deg = np.linspace(-10,10,21)
    geometries = ["2412", "2421", "0015"]
    experimental = 
    #geometries = ["dev/input/2412_200.txt","dev/input/2421_200.txt", "dev/input/0015_200.txt"]

    for geometry in geometries:

        geometry_file = "dev/input/"+geometry+"_200.txt"

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

        print(C_L)
        print(C_m_le)
        print(C_m_qtr)
            


