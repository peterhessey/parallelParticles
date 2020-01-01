import subprocess
import os
import math
import time
import numpy as np
import random
from matplotlib import pyplot as plt

args = []

for i in range(2):
    if i == 0:
        args.append('assignment1.exe')
    else:
        args.append('assignment4.exe')

    time_step = 0.001

    args.append('-1.0')
    args.append('1.0')
    args.append(str(time_step))

    for j in range(400):
        for k in range(40):
            args.append(str(i*5))
            args.append(str(j*5))
            args.append('0')

            for _ in range(3):
                args.append('0')

            args.append(str(random.randint(1,5)))

    start_time = time.time()

    process = subprocess.Popen(args, stdout=subprocess.PIPE)
    test_output = process.communicate()[0].decode('utf-8')

    end_time = time.time()
    total_time = end_time - start_time

    print('Test on %s ran in %s seconds.' % (args[0], str(total_time)))