import pandas as pd
import numpy as np
import gzip 
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go

filepath = r"C:\Users\salma\Downloads\2027.dat.gz"

read_orbit = pd.read_csv(filepath, sep = r"\s+", comment = "#")

pd.set_option("display.max_rows", None)
pd.set_option("display.max_columns", None)
pd.set_option("display.width", None)

JDtime = read_orbit.iloc[:,0].to_numpy()
xGSE = read_orbit.iloc[:,1].to_numpy()
yGSE = read_orbit.iloc[:,2].to_numpy()
zGSE = read_orbit.iloc[:,3].to_numpy()

firstcolumn = read_orbit.index.to_numpy()

fig, ax  = plt.subplots(3)

ax[0].plot(JDtime, xGSE, label = "x-direction", color = "blue")
ax[0].set_xlabel("JD")
ax[0].set_ylabel("GSE (km)")
ax[0].set_title("Moon x-position 2027")
ax[0].legend()
ax[0].grid()

ax[1].plot(JDtime, xGSE, label = "y-direction", color = "green")
ax[1].set_xlabel("JD")
ax[1].set_ylabel("GSE (km)")
ax[1].set_title("Moon y-position 2027")
ax[1].legend()
ax[1].grid()

ax[2].plot(JDtime, xGSE, label = "z-direction",color = "cyan")
ax[2].set_xlabel("JD")
ax[2].set_ylabel("GSE (km)")
ax[2].set_title("Moon z-position 2027")
ax[2].legend()
ax[2].grid()

plt.tight_layout()
plt.show()

plt.plot(xGSE,zGSE, label = "xz-direction", color = "black")
plt.xlabel("JD")
plt.ylabel("GSE (km)")
plt.grid()
plt.title("Moon xz-plane 2027")
plt.legend()

plt.show()


fig = px.scatter(x = xGSE, y = zGSE, animation_frame = JDtime)
fig.update_traces(marker = dict(size=30, color  = "red" ))
fig.add_scatter(x = xGSE, y = zGSE, mode = "lines", line= dict(color = "blue", width = 2), name = "Orbit")
fig.update_layout(xaxis_range  = [xGSE.min(),xGSE.max()], yaxis_range  = [zGSE.min()*1.1,zGSE.max()*1.1],xaxis_title = "X (km)", yaxis_title = "Z (km)")
fig.show()
