import scipy 
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import math
import numpy as np
import pandas as pd
import csv
import re
from uncertainties import ufloat

def read_data():

    #not_indices = [2,3,22,24,28,31,34,38]

    stop_voltage_file = open('stopping_voltages.csv')
    csv_reader = csv.reader(stop_voltage_file, delimiter = ',')

    stop_voltages = []
    trial_numbers = []

    start = True
    for row in csv_reader:
        if start:
            start = False
            continue
        trial_number = int(row[0])
        stop_voltages.append(int(row[1]))
        trial_numbers.append(trial_number)

    data = []
    count = 0
    for i in trial_numbers:
        file = open('Milk - Copy (' + str(i) + ').csv')
        csv_reader = csv.reader(file, delimiter = ',')
        values = []
        start = True
        try:
            for row in csv_reader:
                if start:
                    start = False
                    continue
                try:
                    value = ufloat(int(row[0]), 0.1)
                    values.append(value)
                except ValueError:
                    continue
            try:
                data_pair = (stop_voltages[count], values)
            except IndexError:
                print("INDEX error from file " + str(i))
            data.append(data_pair)
        except Exception as e:
            print("ERROR from file " + str(i))
        
        count += 1
    return data, trial_numbers

def plot_data(data, critical_points):
    print("DATA", data)
    pos = []
    for i in range(len(data[1])):
        pos.append(data[1][i].nominal_value)

    
    time = []
    for i in range(len(pos)):
        time.append(i*0.2)
    plt.plot(time, pos)
    for i in range(len(critical_points)):
        plt.axvline(x = critical_points[i]*0.2, color = 'b')
    plt.title('Time [s] v. Position [pixels]')
    plt.xlabel('Time [s]')
    plt.ylabel('Position [pixels]')
    plt.show()

def find_critical_points(pos, tolerance):
    prev_velocity = 0
    curr_velocity = 0
    critical_points = []

    for i in range(len(pos) - 2):
        prev_velocity = (pos[i + 1] - pos[i]).nominal_value 
        curr_velocity = (pos[i + 2] - pos[i + 1]).nominal_value

        #print(prev_velocity, curr_velocity, i + 1)

        if prev_velocity <= 0 and curr_velocity >= 0:
            critical_points.append(i + 1)
        elif prev_velocity >= 0 and curr_velocity <= 0:
            critical_points.append(i + 1)

        # if start_found == False:
        #     if prev_velocity < 0 & curr_velocity > 0 :
        #         segment_starts.append(i)
        #         start_found = True
        # if start_found == True:
        #     if curr_velocity < prev_velocity:
        #         if i < segment_starts[-1] + tolerance:
        #             continue
        #         segment_ends.append(i)
        #         start_found = False

    return critical_points

def cluster_points(pos, critical_points, tolerance):
    curr = critical_points[0]

    cluster = [curr]
    for i in range(len(critical_points) - 1):
        if critical_points[i] - curr > tolerance:
            cluster.append(critical_points[i])
        curr = critical_points[i]
    
    cluster.append(critical_points[-1])
    cluster.append(len(pos))

    return cluster


def linear(x, m, b):
    return m*x + b
def find_terminal_velocity(pos, critical_points):

    #print(pos)
    real_values = []
    uncertainties = []
    slope_uncertainty = []
    intercept_uncertainty =[]
    chi_array = []
    for i in range(len(pos)):
        real_values.append(pos[i].nominal_value)
        uncertainties.append(pos[i].std_dev)
    
    #print(uncertainties)
    #print(real_values)

    intervals = []
    for i in range(len(critical_points) - 1):
        intervals.append((critical_points[i], critical_points[i + 1]))

    for interval in intervals:
        if interval[0] == interval[1] -1:
            intervals.remove(interval)
    
    velocities = []
    for i in range(len(intervals)):
        start = intervals[i][0]
        end = intervals[i][1]
        time = []
        for j in range(int(end - start)):
            time.append((start + j)*0.2)
        popt, pcov = scipy.optimize.curve_fit(linear, time, real_values[start:int(end)], sigma = uncertainties[start:int(end)])
        sum = 0
        square_sum = 0
        diff_sum = 0
        R_total = 0
        chi = 0
        for j in range(len(time)):
            sum += time[j]
            square_sum += time[j]**2
        
        delta = (int(end) - start) * square_sum - sum**2

        syx = 0
        for j in range(len(time)):
            syx += (real_values[start + j] - linear(time[j], popt[0], popt[1]))**2
        

        if (int(end) - start -2 != 0):
            syx = syx / (int(end) - start - 2)     
        else:
            continue

        for j in range(len(time)):
            diff_sum += (real_values[start + j] - sum/len(time))**2
        for j in range(len(time)):
            chi += (real_values[start + j] - linear(time[j], popt[0], popt[1]))**2/1.83534/len(time)
        if chi < 3 and  chi > 0.5:
            chi_array.append(chi)
        R_squared = 1- (int(end) - start - 2)*syx/diff_sum
        #print('R=Squared', R_squared)
        #R_total += R_squared
        s_m = (int(end) - start) * syx/delta
        s_m = s_m **0.5

        s_b = syx*square_sum/delta
        s_b = s_b**0.5

        slope_uncertainty.append(s_m)
        intercept_uncertainty.append(s_b)

        velocities.append(popt[0]*0.001/543)
    
    method1 = []

    method2 = []


    for i in range(len(velocities)):
        if velocities[i] < 0:
            method2.append(-velocities[i])
        else:
            method1.append(velocities[i])


    chi_sum = 0
    for j in range(len(chi_array)):
        chi_sum += chi_array[j]
    print("ChI array", chi_array)
    if len(chi_array) != 0:
        print("Avergae chi", chi_sum/len(chi_array))
    
    #plot_data(real_values, critical_points)
    return max(method1), max(method2), max(slope_uncertainty), max(intercept_uncertainty)
    #popt, pcov = scipy.curve_fit(linear, time , pos)


        
data, trial_numbers = read_data()

def method(data):
    method1_charges = []
    method2_charges = []
    slope_unc = 0
    intercept_unc = 0
    method1_uncertainty = 0
    method2_uncertainty = 0
    const1 = 2.02 *10**(-10)
    radius = 0

    radius_array = []
    
    nu = 1.827*10**(-5)
    rho_oil = 875.3
    rho_air = 1.204
    g = 9.8
    
    radius_uncertainty = 0
    for i in range(len(data)):

        critical_points = find_critical_points(data[i][1], tolerance = 5)
        cluster = cluster_points(data[i][1], critical_points, tolerance=3)
        terminal_velocity_1, terminal_velocity_2, slope_uncertainty, intercept_uncertainty = find_terminal_velocity(data[i][1], cluster)
        # Uncertainty Propogation
        term1 = (-const1*terminal_velocity_1**1.5/(data[i][0])**2 * 5)**2
        term2 = (const1**0.5/(data[i][0])*3/2 *3.38*10**(-6))**2

        term3 = (-const1*(terminal_velocity_1 + terminal_velocity_2)*terminal_velocity_1**0.5/data[i][0]**2*5)**2
        term4 = ((const1*terminal_velocity_1**0.5/data[i][0]*3/2 + const1*terminal_velocity_2*terminal_velocity_1**(-0.5)/data[i][0]*0.5)*3.38*10**(-6))**2
        term5 = (const1*terminal_velocity_1**0.5/data[i][0] *3.38*10**(-6))**2

        r = ((4.5 *nu* terminal_velocity_1)/(g *(rho_oil - rho_air)))**0.5


        r_unc = (((4.5 *nu)/(g *(rho_oil - rho_air)))**0.5 * terminal_velocity_1**(-0.5) * 0.5 * 3.38*10**(-6))**2

        radius_uncertainty += r_unc

        method1_uncertainty += term1
        method1_uncertainty += term2

        method2_uncertainty += term3
        method2_uncertainty += term4
        method2_uncertainty += term5

        radius += r

        radius_array.append(r)

        print(terminal_velocity_1)
        slope_unc += slope_uncertainty
        intercept_unc += intercept_uncertainty
        method1_charges.append(method_1(terminal_velocity_1, data[i][0]))
        method2_charges.append(method_2(terminal_velocity_1, terminal_velocity_2, data[i][0]))
    
    print("Method 1 Uncertainty: ", method1_uncertainty**0.5)
    print("Method 2 Uncertainty: ", method2_uncertainty**0.5)
    print("Radius Value: ", radius/len(data))
    print("Radius Uncertainty: ", radius_uncertainty**0.5)

    slope_unc = slope_unc/len(data)
    intercept_unc = intercept_unc/len(data)
    return method1_charges, method2_charges, slope_unc, intercept_unc, radius_array


def method_1(terminal_velocity, stop_voltage):
    charge = terminal_velocity**1.5/(float(stop_voltage))
    const = 2.023*10.0**(-10)
    charge = charge*const
    return charge

def method_2(terminal_velocity, velocity_2, stop_voltage):
    charge = (terminal_velocity + velocity_2) * (terminal_velocity)**0.5 /stop_voltage
    const = 2.02 * 10.0**(-10)
    charge = charge*const
    return charge


# for i in range(len(trial_numbers)):
#     critical_points = find_critical_points(data[i][1], tolerance=5)
#     cluster = cluster_points(data[i][1], critical_points, tolerance=3)
#     #plot_data(data[i], cluster)

# critical_points = find_critical_points(data[1][1], tolerance=5)
# cluster = cluster_points(critical_points, tolerance=2)
# plot_data(data[1], cluster)
#terminal_velocity = find_terminal_velocity(data[1][1], segment_start, segment_end)
#print(terminal_velocity)

charges1, charges2, slope_uncertainty, intercept_uncertainty, radius_array = method(data)

print(slope_uncertainty*0.001/543)
print(intercept_uncertainty)
charges1.sort()
charges2.sort()

print(charges1)
print(charges2)

charges1_new = [x*10**19 for x in charges1]
charges2_new = [x*10**19 for x in charges2]

radius_new = [x*10**7 for x in radius_array]

print(charges2_new)
x1 = np.linspace(0, len(charges1), num = len(charges1))
x2 = np.linspace(0, len(charges2), num = len(charges2))

# plt.scatter(x1, charges1_new)
# plt.xlabel('Data Point Number')
# plt.ylabel('Charge [C 1e-19]')
# plt.title('Calculated Charges for Method 1')
# plt.show()

plt.scatter(x2, charges2_new)
plt.xlabel('Data Point Number')
plt.ylabel('Charge [C 1e-19]')
plt.title('Calculated Charges for Method 2')
plt.show()


# plt.scatter(radius_new, charges1_new, c = 'b', label = 'Method 1')
# plt.scatter(radius_new, charges2_new, c = 'r', label = 'Method 2')
# plt.xlabel('radius of the oil drops [m 1e-7]')
# plt.ylabel('Charge [C 1e-19] ')
# plt.legend(loc = 'upper left')
# plt.title('Calculated Charges of Methods 1(blue) and Method 2(red) [C] v. Radius of Droplet [m] ')
# plt.show()


def gcd(data, init):
    gcd = init
    a = 0.01
    for i in range(500):
        gcd = gcd - a* gradientloss(gcd, data)
    return gcd

def gradientloss(gcd, data):
    gradient = 0
    for i in data:
        estimate = 0
        index = 0
        for j in range(100):
            if (gcd*j -i) **2 < (estimate*j - i) **2:
                estimate = gcd
                index = j
        if index != 0:
            gradient += (estimate*index - i)/index
        else:
            gradient += estimate*index - i
        gradient = gradient/len(data)
    return gradient

charge = gcd(charges2_new, 2.5)

print(charge)


# plt.scatter(x1, charges1_new)
# plt.show()
# plt.scatter(x2, charges2_new)
# plt.show()
#print(charges)





