import numpy as np 
import matplotlib.pyplot as plt
import csv

def read_output_file(file_path):
    '''This function reads the output textfile and converts it into numpy arrays in dictionaries'''
    with open(file_path, 'r') as file:
        reader = csv.reader(file)
        headers = next(reader)
        data = np.array([list(map(float, row)) for row in reader])

    data_dict = {header: data[:, i] for i, header in enumerate(headers)}
    return data_dict

def plot_data(x, y, title = "Plot", xlabel = "x-axis", ylabel = "y-axis"):
    '''This function takes in two sets of scalar arrays and spits out a plot'''
    plt.figure(figsize = (8, 6))
    plt.plot(x, y, color = 'b', label = 'Data')
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.show()

def proper_time_calculator(time, lapse):
    '''Calculates the proper time from the coordinate time given the lapse'''
    proper_time_array = np.array([0])
    tau = 0
    time_diff = time[1]-time[0]
    for i in np.arange(0, len(time) - 1, 1):
        tau += time_diff * lapse[i]
        proper_time_array = np.append(proper_time_array, tau)

    return proper_time_array

# Read file
file_path = "/home/mosny/Desktop/GRChombo/Examples/Analysis/output_data.txt"
data_dict = read_output_file(file_path)

# Print to see if it worked
for name, array in data_dict.items():
    print(f"{name}: {array}")

# Create the data arrays as numpy arrays
time = np.array(data_dict["time"])
e_fold = np.array(data_dict["e-fold"])
Hubble = np.array(data_dict["Hubble"])
lapse = np.array(data_dict["lapse"])
Average_Phi = np.array(data_dict["Average_Phi"])
Diff_Phi = np.array(data_dict["Diff_Phi"])
Average_Pi = np.array(data_dict["Average_Pi"])
Diff_Pi = np.array(data_dict["Diff_Pi"])
Rho = np.array(data_dict["Rho"])
Delta = np.array(data_dict["Delta"])
Ham = np.array(data_dict["Ham"])

# Calculate the proper time
proper_time = proper_time_calculator(time, lapse)
proper_time_Hubble = 1/(3*Hubble)
tau0 = 1/(3*Hubble[0])
N0 = e_fold[0]

# Check how close the Pi is to Kination evolution
# plot_data(proper_time, (tau0+proper_time)*Average_Pi/tau0)

# Check the scale factor power p, which for kination should be 1/3
# plot_data(proper_time+tau0, (e_fold - N0)/(np.log((tau0+proper_time)/tau0)))
# plot_data(proper_time+tau0, proper_time + tau0 - proper_time_Hubble)

print(proper_time + tau0)
print(proper_time_Hubble)
plot_data(proper_time, Delta / np.exp(e_fold))
# plot_data(proper_time, Diff_Pi)
# plot_data(e_fold, np.sqrt(8*np.pi)*(Average_Phi - Average_Phi[0])/(np.sqrt(6)*(e_fold - e_fold[0])))