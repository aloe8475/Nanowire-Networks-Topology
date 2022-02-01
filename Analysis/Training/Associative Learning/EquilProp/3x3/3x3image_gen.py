# -*- coding: utf-8 -*-
"""
Created on Thu Jul  1 13:16:32 2021

@author: gnate
"""

import numpy as np

# Define some noise functions
def zero_bit_noise(N):
    noise_list = []
    noise_vec = np.zeros((1,N))
    noise_list.append(noise_vec)
    return noise_list

def one_bit_noise(N):
    noise_list = []
    for i in range(0,N):
        noise_vec = np.zeros((1,N))
        noise_vec[0,i] = 1
        noise_list.append(noise_vec)
    return noise_list

def xor_array(x,y):
    a = x.reshape((1,9))
    b = y.reshape((1,9))
    output_vec = np.zeros((1,a.shape[1]), dtype = 'uint')
    for i in range(0, a.shape[1]):
        output_vec[0, i] = a[0, i] ^ b[0, i]
    return output_vec     

def hamming_dist(x,y):
    dist = 0
    output_vec = xor_array(x,y)
    dist = sum(output_vec.T)
    return dist

## Main program
  
#Generate 3x3 noiseless data
N = 9
orig_inputs = ['111010111', '101010010', '101111101']
orig_inputs_matrix = np.zeros((len(orig_inputs), 1, N), dtype = 'uint')
for i in range(0,len(orig_inputs)):
    input_tuple = tuple(orig_inputs[i])
    for j in range(0, N):
        if input_tuple[j] == '1': orig_inputs_matrix[i,0,j] = 1
#print('Input Matrix -', orig_inputs_matrix)
               
# Generate 5x5 noisy data
noise_vectors = []
noise_vectors.append(zero_bit_noise(N))
noise_vectors.append(one_bit_noise(N))
#print('Noise vectors - ', noise_vectors)

count_vec = np.zeros(2, dtype = 'uint')
for i in range(0,2):
    count_vec[i] = np.math.factorial(N)/(np.math.factorial(i)*np.math.factorial(N-i))

noisy_data_matrix = np.zeros((len(orig_inputs)*sum(count_vec), 1, N), dtype = 'uint')
row_count = 0
for i in range(0,len(orig_inputs)):
    for j in range(0,2):
        bit_noise = np.array(noise_vectors[j], dtype = 'uint')
        for k in range(0, int(count_vec[j])):
            noisy_data_matrix[row_count] = xor_array(orig_inputs_matrix[i], bit_noise[k])
            row_count += 1
            
final_data_matrix = np.zeros((len(orig_inputs)*sum(count_vec), 1, N+1), dtype = 'uint')
final_data_matrix[0:len(orig_inputs)*sum(count_vec), 0, 0:N] = noisy_data_matrix.reshape(len(orig_inputs)*sum(count_vec), N)

for i in range(0,len(orig_inputs)):
    final_data_matrix[sum(count_vec)*i: sum(count_vec)*(i+1), 0, N] = i
    
 
# Save data as csv file
np.savetxt("3x3image_gen.csv", final_data_matrix.reshape(len(orig_inputs)*sum(count_vec), N+1), delimiter=",")
      


            
