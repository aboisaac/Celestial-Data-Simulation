import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import gzip

filepath = r"\Users\salma\Downloads\2026.dat.gz"

print(gzip.open(filepath,"rt").read(1200))

df = pd.read_csv(filepath, sep = r"\s+", comment = "#")


pd.set_option("display.max_rows",None)
pd.set_option("display.max_columns",None)
pd.set_option("display.width", None)

timelist = df.iloc[:,0].to_numpy()
xlist = df.iloc[:,1].to_numpy()
zlist = df.iloc[:,2].to_numpy()
firstcolumn = df.index.to_numpy()

plt.plot(timelist, xlist)
plt.xlabel("unit")
plt.ylabel("km")
plt.title("Moon Orbit Earth xy-plane")
plt.grid()
plt.show()



