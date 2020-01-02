import subprocess
import os
import math
import time
import numpy as np
from pandas import Series, DataFrame
from matplotlib import pyplot as plt

numberOfBodiesToTest = [50000, 100000, 200000]
numberOfThreadsToTest = [1,2,4]

results = []

for i in range(len(numberOfThreadsToTest)):
    results_entry = []
    threadsArgument = numberOfThreadsToTest[i]

    for j in range(len(numberOfBodiesToTest)):
        args = []
        bodiesArgument = numberOfBodiesToTest[j]

        if threadsArgument == 1:
            args.append('assignment4NoParallel.exe')
            args.append(str(bodiesArgument))
        else:
            args.append('assignment4.exe')
            args.append(str(bodiesArgument))
            args.append(str(threadsArgument))

        print('Args: ' , args)
        total_time = 0
        number_of_tests = 1
        # loop 5 times and take average
        for k in range(number_of_tests):
            
            start_time = time.time()

            process = subprocess.Popen(args, stdout=subprocess.PIPE)

            test_output = (process.communicate()[0].decode("utf-8"))

            end_time = time.time()

            time_taken = end_time - start_time

            print('Test %s on %s bodies using %s threads took %ss' % (k, bodiesArgument, threadsArgument, time_taken))

            total_time += (time_taken)
        
        average_runtime = total_time / number_of_tests
        
        print('Average run time for %s bodies on %s threads: %s' %(bodiesArgument, threadsArgument, average_runtime))

        results_entry.append(average_runtime)

    results.append(results_entry)

width = 0.3

for i in range(len(results)):

    bar_chart_data = results[i]
    print(bar_chart_data)

    plt.bar(np.arange(len(bar_chart_data))+width*i, bar_chart_data, width=width)

plt.xticks(range(len(numberOfBodiesToTest)), numberOfBodiesToTest)
plt.xlabel('Number of bodies')
plt.ylabel('Average runtime (seconds)')
plt.legend(numberOfThreadsToTest)
plt.show()