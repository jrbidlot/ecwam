import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime

filename = "input"

dates = []
values = []
with open(filename) as f:
    for line in f:
        cols = line.split()
        if len(cols) < 2:
            continue
        dates.append(datetime.strptime(cols[0], "%Y%m%d%H%M"))
        values.append(float(cols[1]))

x = dates
y = values

fig, ax = plt.subplots()
ax.plot(x, y, color="red", linestyle="-", label=filename)
ax.legend()
fig.autofmt_xdate()
ax.xaxis.set_major_locator(mdates.HourLocator(byhour=[0, 12]))
ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d %Hz"))

# Save as PostScript
fig.savefig(filename + ".ps")

# Save as PNG
fig.savefig(filename + ".png")
