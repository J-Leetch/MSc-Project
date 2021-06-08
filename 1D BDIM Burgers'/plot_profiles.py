import sys
import pandas as pd
import matplotlib.pyplot as plt

plt.figure()

nameslist = ["DNS_Re500_255.csv", "test_Re500_255.csv"]

for arg in nameslist:
    data = pd.read_csv("1D BDIM Burgers'/"+arg, header=None)
    print(data)
    plt.plot(data.iloc[:,0], data.iloc[:,1], label=arg[:-4])

plt.legend()
plt.show()

