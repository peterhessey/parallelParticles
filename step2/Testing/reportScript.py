import subprocess
import os
import math
import numpy as npdi
from matplotlib import pyplot as plt

#  0.01 100.0 1e-8 0 0 0 0 0 0 4 3 0 0 0 0 0 5 3 4 0 0 0 0 3
#  7e-6 is largest time step that yields collisions

# go from high time step to high time step for testing purposes

args = ['reportOutput.exe', '-0.01', '100.0', '1e-6', '0', '0', '0', '0', '0',
        '0', '4', '3', '0', '0', '0', '0', '0', '5', '3', '4', '0', '0', '0',
        '0', '3']


distances = []
times = []

min_distance = float('inf')
max_distance = float('-inf')

min_time = float('inf')
max_time = float('-inf')

with open('./Results/output.txt', 'w') as results_file:
    for exponent in range(7,10,1):
        for mantissa in range(99, 9, -1):
                
            time = (mantissa/10) * (10 ** (-1 * exponent))
            args[3] = str(time)    
            print('Testing with time step %s' % time)

            new_time = (math.log10(-1 * math.log10(time)))

            if new_time < min_time:
                min_time = new_time

            if new_time > max_time:
                max_time = new_time


            process = subprocess.Popen(args, stdout=subprocess.PIPE)
            test_output = (process.communicate()[0].decode("utf-8")).split(',')[:-1]

            distance = math.sqrt(
                (float(test_output[0]) - float(test_output[3]))**2 +
                (float(test_output[1]) - float(test_output[4]))**2 +
                (float(test_output[2]) - float(test_output[5]))**2
            )
            print('Distance between collision points: %s' %distance)

            new_distance = -1*math.log(-1 * math.log(distance))

            if new_distance < min_distance:
                min_distance = new_distance
        
            if new_distance > max_distance:
                max_distance = new_distance

            results_file.write(str(time) + ',' + str(distance) + '\n')

            times.append(new_time)
            distances.append(new_distance)


result_graph = plt.figure()
axes = result_graph.add_axes([0.1,0.1,0.8,0.8])
axes.scatter(times, distances)


# axes.set_xlim(min_time,max_time)
axes.set_ylim(min_distance - (max_distance - min_distance) * 0.1, max_distance  + (max_distance - min_distance) * 0.1)

plt.xlabel('Log-Log Time step size')
plt.ylabel('Log-Log Distace between collision points')

plt.show()