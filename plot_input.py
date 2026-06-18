import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

filename = "input"
data = np.loadtxt(filename)

x = data[:, 0]
y = data[:, 1]

fig, ax = plt.subplots()
ax.plot(x, y, color="red", linestyle="-", label=filename)
ax.legend()

# Save as PostScript
fig.savefig(filename + ".ps")

# Save as PNG
fig.savefig(filename + ".png")
