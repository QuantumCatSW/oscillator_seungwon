# my_package/module.py

import matplotlib as mpl
import matplotlib.pyplot as plt
from qutip import *
import numpy as np

def plot_wignerfunction(state, xvec, pvec, ax=None):
    # input parameter example:
        # xvec = np.linspace(-3, 3, 100)
        # pvec = np.linspace(-3, 3, 100)

    W = wigner(state, xvec, pvec)
    X, P = np.meshgrid(xvec, pvec)

    if ax is None:
        fig, ax = plt.subplots(1, 1)
    
    ax.axhline(0, color='grey', lw=0.5)
    ax.axvline(0, color='grey', lw=0.5)

    # Plot the Wigner function
    c=ax.contourf(X, P, W, 100, norm=mpl.colors.Normalize(-0.25, 0.25), cmap=plt.get_cmap('RdBu_r'))
    ax.set_xlim([xvec[0], xvec[-1]])
    ax.set_ylim([pvec[0], pvec[-1]])
    ax.set_xlabel(r'$x: real$', fontsize=16)
    ax.set_ylabel(r'$p: img$', fontsize=16)
    # ax.set_title('Wigner function', fontsize=16)

    # Add a colorbar to the plot
    cbar = plt.colorbar(c, ax=ax)
    cbar.set_label('Wigner Function Intensity')    

    # plt.show()
    return ax


def plot_fock_dist(cavity_dm, dim_c, max_fock_state):
    # Calculate populations for each Fock state |n⟩ from |0⟩ to |max_fock_state⟩
    fock_populations = []
    for n in range(max_fock_state + 1):
        # Create the Fock state |n⟩
        fock_state_n = basis(dim_c, n)

        # Calculate the population (overlap with the Fock state)
        population_n = expect(fock_state_n * fock_state_n.dag(), cavity_dm)

        # Store the result
        fock_populations.append(population_n)

    # Plot the Fock state population distribution as a histogram
    fock_states = np.arange(max_fock_state + 1)
    plt.figure(figsize=(10, 6))
    plt.bar(fock_states, fock_populations, color='b', alpha=0.7)

    # Calculate the mean photon number
    mean_photon_number = expect(create(dim_c).dag() * destroy(dim_c), cavity_dm)
    # print(f"Mean photon number: {mean_photon_number}")

    # Plot settings
    plt.xlabel('Fock State |n⟩')
    plt.ylabel('Population')
    plt.title(f'Fock State Population Distribution\nMean Photon Number = {mean_photon_number:.3f}')
    plt.grid(True)
    plt.show()

def Hamiltonian_dispersive(dim_c, dim_q, w_c, w_q, K_c, K_q, chi_cq, kappa_c, kappa_q, plotter):
    # Hamiltonian
    H_cavity = w_c*create(dim_c)*destroy(dim_c) - 0.5*K_c*create(dim_c)*create(dim_c)*destroy(dim_c)*destroy(dim_c)
    H_qubit = w_q*create(dim_q)*destroy(dim_q) - 0.5*K_q*create(dim_q)*create(dim_q)*destroy(dim_q)*destroy(dim_q)

    H0 = (
        tensor(H_cavity,identity(dim_q)) + tensor(identity(dim_c),H_qubit)
        + chi_cq*tensor(create(dim_c)*destroy(dim_c),create(dim_q)*destroy(dim_q))
    )

    c_ops = [
    np.sqrt(kappa_c)*tensor(destroy(dim_c),identity(dim_q))
    ,np.sqrt(kappa_q)*tensor(identity(dim_c),destroy(dim_q))
    ]

    # Charge drive Hamiltonian
    H_drive_cavity = tensor((create(dim_c)+destroy(dim_c)),identity(dim_q))
    H_drive_qubit  = tensor(identity(dim_c),(create(dim_q)+destroy(dim_q)))

    if plotter==True:
        # Assume these eigenenergy arrays (replace with actual calculations)
        eigenenergies_cavity = np.array(H_cavity.eigenenergies()/(2*np.pi*1e9))   # Replace with H_cavity.eigenenergies()/(2*np.pi*1e9)
        eigenenergies_qubit = np.array(H_qubit.eigenenergies()/(2*np.pi*1e9))    # Replace with H_qubit.eigenenergies()/(2*np.pi*1e9)
    
        # Define a cutoff frequency
        cutoff_frequency = 100  # GHz
        bottom_frequency = 0  # GHz

        # Apply cutoff filter for each component
        eigenenergies_cavity = eigenenergies_cavity[eigenenergies_cavity <= cutoff_frequency]
        eigenenergies_cavity = eigenenergies_cavity[bottom_frequency <= eigenenergies_cavity]
        eigenenergies_qubit = eigenenergies_qubit[eigenenergies_qubit <= cutoff_frequency]
        eigenenergies_qubit = eigenenergies_qubit[bottom_frequency <= eigenenergies_qubit]

        # Calculate the eigenenergies of the full Hamiltonian H0
        eigenenergies_H0 = np.array(H0.eigenenergies()/(2*np.pi*1e9))  # Replace with actual eigenenergies calculation

        # Filter the eigenenergies of H0 based on the cutoff frequency
        eigenenergies_H0 = eigenenergies_H0[eigenenergies_H0 <= cutoff_frequency]
        eigenenergies_H0 = eigenenergies_H0[bottom_frequency <= eigenenergies_H0]

        # Offsets for each component for better visualization
        offsets = [0, 1]  # Just horizontal offsets to separate them visually on the plot
        labels = ['Cavity', 'Qubit']

        # Create a plot
        plt.figure(figsize=(6, 10))

        # Plot H0 eigenenergies as background bars (spanning the whole plot)
        for e in eigenenergies_H0:
            plt.hlines(y=e, xmin=-0.5, xmax=2, color='gray', alpha=0.3, label='H0' if e == eigenenergies_H0[0] else "")
            # if (eigenenergies_H0[e] < cutoff_frequency):
            #     plt.text(offsets[1] + 0.5, e, f'{e:.2f} GHz', color='r', va='center')  # Add text label next to each level

        # Plot Cavity Energy Levels
        for e in eigenenergies_cavity:
            plt.hlines(y=e, xmin=offsets[0]-0.2, xmax=offsets[0]+0.2, color='r', linewidth=2, label='Cavity' if e == eigenenergies_cavity[0] else "")

        # Plot Qubit Energy Levels
        for e in eigenenergies_qubit:
            plt.hlines(y=e, xmin=offsets[1]-0.2, xmax=offsets[1]+0.2, color='b', linewidth=2, label='Qubit' if e == eigenenergies_qubit[0] else "")

        # Customize the plot
        plt.ylim(-1, 40) # Set y-axis limits to focus on specific region (e.g., between 0 and 30 GHz)
        plt.xticks(ticks=offsets, labels=labels)
        plt.ylabel('Energy Levels (GHz)')
        plt.title('Energy Levels of Bare Components')
        plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
        plt.grid(True)

        # Show the plot
        plt.show()

    return H0, c_ops, H_drive_cavity, H_drive_qubit

def Hamiltonian_lamor(dim_c, dim_q, w_c, w_q, chi_cq):
    # Hamiltonian
    H_cavity = w_c*create(dim_c)*destroy(dim_c)
    H_qubit = w_q*create(dim_q)*destroy(dim_q)

    H0 = (
        tensor(H_cavity,identity(dim_q)) + tensor(identity(dim_c),H_qubit)
        + chi_cq*tensor(create(dim_c)*destroy(dim_c),create(dim_q)*destroy(dim_q))
    )

    return H0

def Observables(H0):
    # # eigenstates
    # es_0g = H0.eigenstates()[1][0]
    # es_0e = H0.eigenstates()[1][1]
    # es_1g = H0.eigenstates()[1][2]
    # es_0f = H0.eigenstates()[1][3]
    # es_1e = H0.eigenstates()[1][4]
    # es_2g = H0.eigenstates()[1][5]
    # es_0h = H0.eigenstates()[1][6]
    # es_1f = H0.eigenstates()[1][7]
    # es_2e = H0.eigenstates()[1][8]
    # es_3g = H0.eigenstates()[1][9]
    # es_0i = H0.eigenstates()[1][10]

    # # observables
    # observables = [
    #     es_0g*es_0g.dag(), es_0e*es_0e.dag(), es_1g*es_1g.dag(), 
    #     es_0f*es_0f.dag(), es_1e*es_1e.dag(), es_2g*es_2g.dag(), 
    #     es_0h*es_0h.dag(), es_1f*es_1f.dag(), es_2e*es_2e.dag(), 
    #     es_3g*es_3g.dag(), es_0i*es_0i.dag()
    #     ]
    
    observables = []
    for i in range(0, 20):
        observables.append(H0.eigenstates()[1][i]*H0.eigenstates()[1][i].dag())

    return observables