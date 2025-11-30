from matplotlib import pyplot as plt
import numpy as np

data = []
with open("points.dat", 'r') as f:
    for line in f:
        data.append(float(line))

logbins = np.logspace(-7, 0, 25)
fig, ax = plt.subplots()
ax.hist(data, bins=logbins)
plt.xscale('log')
plt.show()
