import sys
import matplotlib

import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('output.txt').abs()

df.plot()
plt.yscale("log")
plt.show()