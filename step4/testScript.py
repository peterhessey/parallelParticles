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
        args.append('assignment4NoParallel.exe')
    else:
        args.append('assignment4.exe')

    start_time = time.time()

    process = subprocess.Popen(args, stdout=subprocess.PIPE)
    test_output = process.communicate()[0].decode('utf-8')

    end_time = time.time()
    total_time = end_time - start_time

    print('Test on %s ran in %s seconds.' % (args[0], str(total_time)))