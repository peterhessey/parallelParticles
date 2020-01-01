import random
with open('particles.txt', 'w') as particleOutput:
    for i in range(500):
        for j in range(40):
            particleOutput.write(str(i*5) + ' ')
            particleOutput.write(str(j*5) + ' ')
            particleOutput.write('0 ')
            for _ in range(3):
                particleOutput.write('0 ')
            particleOutput.write(str(random.randint(1,5)) + ' ')