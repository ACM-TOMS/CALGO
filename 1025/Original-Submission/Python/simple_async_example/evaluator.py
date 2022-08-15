import sys
import numpy as np
import time

xs = [float(x) for x in sys.argv[1:]]
val = 1
for x in xs:
    val *= np.exp(-(x - 2.) ** 2) + np.exp(-(x - 6.) ** 2 / 10) + 1 / (x ** 2 + 1)

val *= -1.0

# time.sleep(np.random.random() * 10.)

time.sleep(5)

with open('result.txt', 'w') as f:
    f.write(str(val) + '\n')
