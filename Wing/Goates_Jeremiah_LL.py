import json
import numpy as np
import matplotlib.pyplot as plt


if __name__ == "__main__":

    
    # Parse input
    with open("input.json", "r") as input_handle:
        input_settings = json.load(input_handle)

    R_A = input_settings["wing"]["planform"]["aspect_ratio"]
    R_T = input_settings["wing"]["planform"]["taper_ratio"]
    N = 1+2*input_settings["wing"]["nodes_per_semispan"]
    c_la = input_settings["wing"]["airfoil_lift_slope"]



    # Calculate planform
    theta = np.linspace(0,np.pi,N)
    z_b = -0.5*np.cos(theta)
    c_b = 2*(1-(1-R_T)*np.abs(np.cos(theta)))/(R_A*(1+R_T))

    # Calculate washout
    washout_distribution = input_settings['wing']['washout']['distribution']
    
    if washout_distribution == 'optimum':
        omega = 1-np.sin(theta)/(c_b/c_b[N//2])
    elif washout_distribution == 'linear':
        omega = abs(np.cos(theta))
    elif washout_distribution == 'none':
        omega = 0.0
    else:
        pass

    # Calculate ailerons
    hinge_efficiency = input_settings['wing']['aileron']['hinge_efficiency']
    z_b_begin = input_settings['wing']['aileron']['begin[z/b]']
    z_b_end = input_settings['wing']['aileron']['end[z/b]']
    cf_c_begin = input_settings['wing']['aileron']['begin[cf/c]']
    cf_c_end = input_settings['wing']['aileron']['end[cf/c]']

    cf_c = np.zeros(N)
    x_b_ail = np.zeros(N)
    
    c_b_begin = 2*(1-(1-R_T)*np.abs(2*z_b_begin))/(R_A*(1+R_T))
    c_b_end = 2*(1-(1-R_T)*np.abs(2*z_b_end))/(R_A*(1+R_T))
    x_b_begin = cf_c_begin*c_b_begin
    x_b_end = cf_c_end*c_b_end
    for i in range(N):
        if abs(z_b[i]) > z_b_end:
            cf_c[i] = 0.0
        elif abs(z_b[i]) < z_b_begin:
            cf_c[i] = 0.0
        else:
            x_b_ail[i] = (1- (z_b_end - abs(z_b[i])) / (z_b_end - z_b_begin))*x_b_end + ((z_b_end - abs(z_b[i])) / (z_b_end - z_b_begin)) * x_b_begin
            cf_c[i] = ( x_b_ail[i])/c_b[i]

    theta_f = np.arccos(2*cf_c - 1)
    eps_fi = 1 - (theta_f - np.sin(theta_f))/np.pi
    eps_f = eps_fi*hinge_efficiency

    chi = -1*eps_f*np.sign(z_b)


    C = np.zeros((N,N))
    for i in range(N):
        for j in range(N):
            n = j+1
            if i == 0:
                C[i,j] = n**2
            elif i==N-1:
                C[i,j] = ((-1)**(n+1))*n**2
            else:
                C[i,j] = ( (4/(c_la*c_b[i])) + (n/np.sin(theta[i])) )*np.sin(n*theta[i])
    
    C_inv = np.linalg.inv(C)

    a_n = np.matmul(C_inv, np.ones(N))
    b_n = np.matmul(C_inv, omega)
    c_n = np.matmul(C_inv, chi)
    d_n = np.matmul(C_inv, np.cos(theta))

    # Write text files
    np.savetxt('C.txt', C)
    np.savetxt('C_inv.txt', C_inv)
    np.savetxt('a_vec.txt', a_n)
    np.savetxt('b_vec.txt', b_n)
    np.savetxt('c_vec.txt', c_n)
    np.savetxt('d_vec.txt', d_n)
    
    # Calculate constants
    kappa_L = (1 - (1 + np.pi*R_A/c_la)*a_n[0])/( (1 + np.pi*R_A/c_la)*a_n[0])
    C_La = np.pi*R_A*a_n[0]
    eps_omega = b_n[0]/a_n[0]
    kappa_D = 0.0
    kappa_DL = 0.0
    kappa_D_omega = 0.0
    for i in range(1,N):
        j = i+1
        kappa_D += j*a_n[i]**2/a_n[0]**2
        kappa_DL += j*(a_n[i]/a_n[0])*((b_n[i]/b_n[0]) - (a_n[i]/a_n[0]))
        kappa_D_omega += j*( b_n[i]/b_n[0] - a_n[i]/a_n[0] )**2
    
    kappa_DL = 2*b_n[0]/a_n[0] * kappa_DL
    kappa_D_omega = (b_n[0]/a_n[0])**2 * kappa_D_omega

    C_l_delta_a = - np.pi*R_A/4 * c_n[1]
    C_l_p = - np.pi*R_A/4 * d_n[1]

    print()
    print("CONSTANTS")
    print()
    print('kappa_L')
    print(kappa_L)
    print('C_La')
    print(C_La)
    print('eps_omega')
    print(eps_omega)
    print('kappa_D')
    print(kappa_D)
    print('kappa_DL')
    print(kappa_DL)
    print('kappa_D_omega')
    print(kappa_D_omega)
    print('C_l_delta_a')
    print(C_l_delta_a)
    print('C_l_pbar')
    print(C_l_p)



    if input_settings['wing']['washout']['amount[deg]'] == 'optimum':
        CL_design = input_settings['wing']['washout']['CL_design']
        washout_deg = kappa_DL*CL_design/(2*kappa_D_omega*C_La)
    else:
        washout_deg = np.deg2rad(input_settings['wing']['washout']['amount[deg]'])

    delta_a = np.deg2rad(input_settings['condition']['aileron_deflection[deg]'])

    if input_settings['condition']['pbar'] == 'steady':
        pbar = - C_l_delta_a/C_l_p * delta_a
    else:
        pbar = input_settings['condition']['pbar']

    if input_settings['condition']['alpha_root[deg]'] == 'CL':
        CL = input_settings['condition']['CL']
        alpha_root = CL/C_La + eps_omega*washout_deg
    else:
        alpha_root = np.deg2rad(input_settings['condition']['alpha_root[deg]'])

    # Compute A series
    A = np.zeros(N)
    A = a_n*alpha_root - b_n*washout_deg + c_n*delta_a + d_n*pbar
    
    # Compute forces and moments
    CL = np.pi*R_A*A[0]
    CDi = 0.0
    for i in range(N):
        j = i+1
        CDi += j*A[i]**2
    CDi = np.pi*R_A*CDi - np.pi*R_A*pbar/2 * A[1]
    Cl = -np.pi*R_A/4 * A[1]
    Cn = 0.0
    for i in range(1,N):
        j = i+1
        Cn += (2*j-1)*A[i-1]*A[i]
    Cn = np.pi*R_A/4 * Cn - np.pi*R_A*pbar/8 * (A[0] + A[2])

    print()
    print('FORCES AND MOMENTS')
    print()
    print('CL')
    print(CL)
    print('CDi')
    print(CDi)
    print('Cl')
    print(Cl)
    print('Cn')
    print(Cn)
    print('pbar_steady')
    print(pbar)



    
    # Compute lift distributions
    C_hat_L_plan = np.zeros(N)
    C_hat_L_wash = np.zeros(N)
    C_hat_L_ail =  np.zeros(N)
    C_hat_L_roll = np.zeros(N)
    for i in range(N):
        j = i+1
        C_hat_L_plan += a_n[i]*np.sin(j*theta)
        C_hat_L_wash += b_n[i]*np.sin(j*theta)
        C_hat_L_ail += c_n[i]*np.sin(j*theta)
        C_hat_L_roll += d_n[i]*np.sin(j*theta)
        
    C_hat_L_plan = 4*alpha_root*C_hat_L_plan 
    C_hat_L_wash = -4*washout_deg*C_hat_L_wash 
    C_hat_L_ail  = 4*delta_a*C_hat_L_ail  
    C_hat_L_roll = 4*pbar*C_hat_L_roll 
    C_hat_L = C_hat_L_ail+C_hat_L_plan+C_hat_L_roll+C_hat_L_wash

    c_L_plan = C_hat_L_plan/c_b
    c_L_wash = C_hat_L_wash/c_b
    c_L_ail  = C_hat_L_ail/c_b
    c_L_roll = C_hat_L_roll/c_b
    c_L = C_hat_L/c_b

    # Plot geometry
    plt.figure()
    plt.plot(z_b,c_b/4, 'k-', label='Planform')
    plt.plot(z_b,-3*c_b/4, 'k-')
    for i, z in enumerate(z_b):
        if i == 1:
            plt.plot([z, z],[c_b[i]/4,-3*c_b[i]/4], 'k--', label='Wing sections')
        else:
            plt.plot([z, z],[c_b[i]/4,-3*c_b[i]/4], 'k--')
    plt.plot(z_b, -3*c_b/4+x_b_ail, color='red', label='Aileron')
    plt.legend(fontsize=6)
    plt.xlabel('$z/b$')
    plt.ylabel('$c/b$')
    plt.axis('equal')
    plt.show()

    # Plot lift distribution
    plt.figure()
    plt.plot(z_b, C_hat_L, label='$\hat{ C}_{L_{total}}$')
    plt.plot(z_b, C_hat_L_plan, label='$\hat{ C}_{L_{planform}}$')
    plt.plot(z_b, C_hat_L_wash, label='$\hat{ C}_{L_{washout}}$')
    plt.plot(z_b, C_hat_L_ail, label='$\hat{ C}_{L_{aileron}}$')
    plt.plot(z_b, C_hat_L_roll, label='$\hat{ C}_{L_{roll}}$')
    plt.xlabel('$z/b$')
    plt.ylabel('$\hat{ C } _L$')
    plt.legend(fontsize=6)
    plt.show()

    # Plot section lift distribution
    plt.figure()
    plt.plot(z_b, c_L, label='$\\tilde{ C}_{L_{total}}$')
    plt.plot(z_b, c_L_plan, label='$\\tilde{ C}_{L_{planform}}$')
    plt.plot(z_b, c_L_wash, label='$\\tilde{ C}_{L_{washout}}$')
    plt.plot(z_b, c_L_ail, label='$\\tilde{ C}_{L_{aileron}}$')
    plt.plot(z_b, c_L_roll, label='$\\tilde{ C}_{L_{roll}}$')
    plt.xlabel('$z/b$')
    plt.ylabel('$\\tilde{ C} _L$')
    plt.legend(fontsize=6)
    plt.show()


    plt.figure()
    plt.plot(z_b, chi, 'k-')
    plt.xlabel('$z/b$')
    plt.ylabel('$\chi (z)$')
    plt.show()

    plt.figure()
    plt.plot(z_b, omega, 'k-')
    plt.xlabel('$z/b$')
    plt.ylabel('$\omega (z)$')
    plt.show()



