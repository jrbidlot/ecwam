import numpy as np
import matplotlib.pyplot as plt

filename = "input"
data = np.loadtxt(filename)

x = data[:, 0]
y = data[:, 1]

plt.plot(x, y, color="red", linestyle="-", label=filename)
plt.legend()
plt.show()
