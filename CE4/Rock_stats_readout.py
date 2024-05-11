import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.optimize import curve_fit
from scipy.special import factorial
from scipy.stats import poisson

df = pd.read_csv('./RockNums_8a8b10d_0.03D_0.064142Dm.csv')
print(len(df.iloc[:]))
test = np.array(df.iloc[121][1:])
print(test)
plt.hist(test, density=True)
plt.show()
std = []
mean = []
for ii in np.arange(len(df.iloc[:])):
      std.append(np.std(df.iloc[ii][1:]))
      mean.append(np.mean(df.iloc[ii][1:]))

plt.plot(df['Diameter'], std)
plt.show()
#print(df.columns)
x = 0
for i in df.columns[1:]:
    if np.max(np.nonzero(np.array(df[i]))) > 100:
            x += 1

std_phi_D = []
for ii in range(len(std)):
      std_phi_D.append(std[ii]*np.pi*(df['Diameter'][ii])**3/(6*8*8*10))
plt.plot(df['Diameter'], std_phi_D)
plt.show()

df_phi_D = pd.DataFrame(std_phi_D, columns = ['STD(phi_D)'])
df_phi_D.to_csv('./std_phi_D_8a8b10d.csv')
print(x/10000)