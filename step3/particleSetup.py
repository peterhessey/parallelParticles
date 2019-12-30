import random
with open('particles.txt', 'w') as particleOutput:
    for i in range(10):
        for j in range(20):
            particleOutput.write(str(i) + ' ')
            particleOutput.write(str(j) + ' ')
            particleOutput.write('0 ')
            for _ in range(3):
                particleOutput.write('0 ')
            particleOutput.write(str(random.randint(1,5)) + ' ')