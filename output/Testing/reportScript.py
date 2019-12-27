import subprocess
import os
import math
from matplotlib import pyplot as plt

#  0.01 100.0 1e-8 0 0 0 0 0 0 4 3 0 0 0 0 0 5 3 4 0 0 0 0 3
#  7e-6 is largest time step that yields collisions

# go from high time step to high time step for testing purposes

args = ['reportOutput.exe', '-0.01', '100.0', '1e-6', '0', '0', '0', '0', '0',
        '0', '4', '3', '0', '0', '0', '0', '0', '5', '3', '4', '0', '0', '0',
        '0', '3']


distances = []
times = []
# will increase to 9 once testing complete!
for exponent in range(5,7,1):
    for mantissa in range(10, 99, 1):

        time = (mantissa/10) * (10 ** (-1 * exponent))
        args[4] = str(time)    
        print('Testing with time step %s' % time)

        times.append(math.log10(-1 * math.log10(time)))


        process = subprocess.Popen(args, stdout=subprocess.PIPE)
        test_output = (process.communicate()[0].decode("utf-8")).split(',')[:-1]

        distance = math.sqrt(
            (float(test_output[0]) - float(test_output[3]))**2 +
            (float(test_output[1]) - float(test_output[4]))**2 +
            (float(test_output[2]) - float(test_output[5]))**2
        )
        print('Distance between collision points: %s' %distance)

        distances.append(math.log(-1 * math.log(distance)))


plt.scatter(distances, times)
plt.title('Time step size against Distance between collision points')
plt.show()