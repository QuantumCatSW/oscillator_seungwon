import numpy as np
import matplotlib.pyplot as plt

# Define the envelopes for qubit and cavity using flat-top pulses with carrier frequencies
def envelope_qubit(t, args):
    A = args['A_qubit']
    carrier_frequency = args['w_qubit']  # Carrier frequency for the qubit
    phi=args['phi_qubit']

    flat_time = args['flat_time_qubit']  # Duration of the flat top
    slope_time = args['slope_time_qubit']  # Duration of the rising/falling edges
    total_time = slope_time * 2 + flat_time  # Total duration of the pulse

    if flat_time == 0 and slope_time == 0:
        return 0
    elif t < slope_time:
        return A * np.cos(np.pi * t / slope_time / 2-np.pi/2) * np.cos(carrier_frequency * t + phi)
    elif t < slope_time + flat_time:
        return A * np.cos(carrier_frequency * t + phi)
    elif t < slope_time + flat_time + slope_time:
        return A * np.sin(np.pi * (total_time - t) / slope_time / 2) * np.cos(carrier_frequency * t + phi)
    else:
        return 0

def envelope_qubit_2nd(t, args):
    A = args['A_qubit_2nd']
    carrier_frequency = args['w_qubit_2nd']  # Carrier frequency for the qubit
    phi=args['phi_qubit_2nd']

    flat_time = args['flat_time_qubit_2nd']  # Duration of the flat top
    slope_time = args['slope_time_qubit_2nd']  # Duration of the rising/falling edges
    total_time = slope_time * 2 + flat_time  # Total duration of the pulse

    if flat_time == 0 and slope_time == 0:
        return 0
    elif t < slope_time:
        return A * np.cos(np.pi * t / slope_time / 2-np.pi/2) * np.cos(carrier_frequency * t + phi)
    elif t < slope_time + flat_time:
        return A * np.cos(carrier_frequency * t + phi)
    elif t < slope_time + flat_time + slope_time:
        return A * np.sin(np.pi * (total_time - t) / slope_time / 2) * np.cos(carrier_frequency * t + phi)
    else:
        return 0

def envelope_qubit_3rd(t, args):
    A = args['A_qubit_3rd']
    carrier_frequency = args['w_qubit_3rd']  # Carrier frequency for the qubit
    phi=args['phi_qubit_3rd']

    flat_time = args['flat_time_qubit_3rd']  # Duration of the flat top
    slope_time = args['slope_time_qubit_3rd']  # Duration of the rising/falling edges
    total_time = slope_time * 2 + flat_time  # Total duration of the pulse

    if flat_time == 0 and slope_time == 0:
        return 0
    elif t < slope_time:
        return A * np.cos(np.pi * t / slope_time / 2-np.pi/2) * np.cos(carrier_frequency * t + phi)
    elif t < slope_time + flat_time:
        return A * np.cos(carrier_frequency * t + phi)
    elif t < slope_time + flat_time + slope_time:
        return A * np.sin(np.pi * (total_time - t) / slope_time / 2) * np.cos(carrier_frequency * t + phi)
    else:
        return 0

def envelope_qubit_4th(t, args):
    A = args['A_qubit_4th']
    carrier_frequency = args['w_qubit_4th']  # Carrier frequency for the qubit
    phi=args['phi_qubit_4th']

    flat_time = args['flat_time_qubit_4th']  # Duration of the flat top
    slope_time = args['slope_time_qubit_4th']  # Duration of the rising/falling edges
    total_time = slope_time * 2 + flat_time  # Total duration of the pulse

    if flat_time == 0 and slope_time == 0:
        return 0
    elif t < slope_time:
        return A * np.cos(np.pi * t / slope_time / 2-np.pi/2) * np.cos(carrier_frequency * t + phi)
    elif t < slope_time + flat_time:
        return A * np.cos(carrier_frequency * t + phi)
    elif t < slope_time + flat_time + slope_time:
        return A * np.sin(np.pi * (total_time - t) / slope_time / 2) * np.cos(carrier_frequency * t + phi)
    else:
        return 0
    
def envelope_qubit_5th(t, args):
    A = args['A_qubit_5th']
    carrier_frequency = args['w_qubit_5th']  # Carrier frequency for the qubit
    phi=args['phi_qubit_5th']

    flat_time = args['flat_time_qubit_5th']  # Duration of the flat top
    slope_time = args['slope_time_qubit_5th']  # Duration of the rising/falling edges
    total_time = slope_time * 2 + flat_time  # Total duration of the pulse

    if flat_time == 0 and slope_time == 0:
        return 0
    elif t < slope_time:
        return A * np.cos(np.pi * t / slope_time / 2-np.pi/2) * np.cos(carrier_frequency * t + phi)
    elif t < slope_time + flat_time:
        return A * np.cos(carrier_frequency * t + phi)
    elif t < slope_time + flat_time + slope_time:
        return A * np.sin(np.pi * (total_time - t) / slope_time / 2) * np.cos(carrier_frequency * t + phi)
    else:
        return 0
    
def envelope_cavity(t, args):
    A = args['A_cavity']
    carrier_frequency = args['w_cavity']  # Carrier frequency for the cavity
    phi=args['phi_cavity']

    flat_time = args['flat_time_cavity']  # Duration of the flat top
    slope_time = args['slope_time_cavity']  # Duration of the rising/falling edges
    total_time = slope_time * 2 + flat_time  # Total duration of the pulse

    if flat_time == 0 and slope_time == 0:
        return 0
    elif t < slope_time:
        return A * np.cos(np.pi * t / slope_time / 2-np.pi/2) * np.cos(carrier_frequency * t + phi)
    elif t < slope_time + flat_time:
        return A * np.cos(carrier_frequency * t + phi)
    elif t < slope_time + flat_time + slope_time:
        return A * np.sin(np.pi * (total_time - t) / slope_time / 2) * np.cos(carrier_frequency * t + phi)
    else:
        return 0







# Deprecated
# Plot the pulse sequence
def pulse_plotter(tlist,pulse_sequence):
    plt.figure(figsize=(10, 5))
    plt.plot(tlist*1e9,pulse_sequence*1e-6,label='output pulse')
    plt.legend()
    plt.xlabel('Time (ns)')
    plt.ylabel('Amplitude (MHz)')

def envelope_qubit_plotter(t, args):
    A = args['A_qubit']
    w = args['w_qubit']  # Carrier frequency for the qubit
    phi=args['phi_qubit']

    flat_time = args['flat_time_qubit']  # Duration of the flat top
    slope_time = args['slope_time_qubit']  # Duration of the rising/falling edges
    total_time = slope_time * 2 + flat_time  # Total duration of the pulse

    # Define boolean masks for regions
    rising_edge = (t >= 0) & (t < slope_time)
    flat_edge = (t >= slope_time) & (t < slope_time + flat_time)
    falling_edge = (t >= slope_time + flat_time) & (t <= total_time)
    # Initialize the pulse array with zeros
    pulse = np.zeros_like(t)
    # Rising edge: Sinusoidal increase
    pulse[rising_edge] = A * np.cos(np.pi * t[rising_edge] / slope_time /2-np.pi/2) 
    # Flat region: constant amplitude
    pulse[flat_edge] = A 
    # Falling edge: Sinusoidal decrease
    pulse[falling_edge] = A * np.sin(np.pi * (total_time - t[falling_edge]) / slope_time/2)
    # return pulse 
    return pulse * np.cos(w * t + phi)

def envelope_cavity_plotter(t, args):
    A = args['A_cavity']
    w = args['w_cavity']  # Carrier frequency for the qubit
    phi=args['phi_cavity']

    flat_time = args['flat_time_cavity']  # Duration of the flat top
    slope_time = args['slope_time_cavity']  # Duration of the rising/falling edges
    total_time = slope_time * 2 + flat_time  # Total duration of the pulse

    # Define boolean masks for regions
    rising_edge = (t >= 0) & (t < slope_time)
    flat_edge = (t >= slope_time) & (t < slope_time + flat_time)
    falling_edge = (t >= slope_time + flat_time) & (t <= total_time)
    # Initialize the pulse array with zeros
    pulse = np.zeros_like(t)
    # Rising edge: Sinusoidal increase
    pulse[rising_edge] = A * np.cos(np.pi * t[rising_edge] / slope_time /2-np.pi/2) 
    # Flat region: constant amplitude
    pulse[flat_edge] = A 
    # Falling edge: Sinusoidal decrease
    pulse[falling_edge] = A * np.sin(np.pi * (total_time - t[falling_edge]) / slope_time/2)
    # return pulse 
    return pulse * np.cos(w * t + phi)
