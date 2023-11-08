import sys
import matplotlib

import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('psi_tutta.dat', sep="  ")
df.set_index("x")
df["psi"].plot()
plt.show()