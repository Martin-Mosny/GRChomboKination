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

# Read file
file_path = "/home/mosny/Desktop/GRChombo/Examples/Analysis/output_data.txt"
data_dict = read_output_file(file_path)

# Print to see if it worked
for name, array in data_dict.items():
    print(f"{name}: {array}")