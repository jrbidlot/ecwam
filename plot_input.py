import io
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

filename = "input"
data = np.loadtxt(filename)

x = data[:, 0]
y = data[:, 1]

fig, ax = plt.subplots()
ax.plot(x, y, color="red", linestyle="-", label=filename)
ax.legend()

# Save as PostScript
fig.savefig(filename + ".ps")

# Save as GIF via PIL
buf = io.BytesIO()
fig.savefig(buf, format="png")
buf.seek(0)
Image.open(buf).save(filename + ".gif")
