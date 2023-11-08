import sys
import matplotlib

import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('grafico.txt', delimiter=';', index_col='x')

df.plot()
#plt.yscale("log")
plt.show()