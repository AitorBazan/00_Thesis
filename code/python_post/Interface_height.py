import os
import re
from tkinter import Tk
from tkinter.filedialog import askdirectory
import csv
import matplotlib.pyplot as plt
import numpy as np
from glob import iglob
from cmath import pi, sqrt, exp

path = askdirectory(title='Select Folder') 
for root, dirs, files in os.walk(path):
    numL2 = []
    denL2 = []
    L2 = []
    numL2_sol2 = []
    denL2_sol2 = []
    L2_sol2 = []
    dir = []
    n = 0
    for i in files:
        with os.scandir(root) as it:
            alphaField = []
            posField = []
            t = []

            j = 0
            for entry in it:
                if entry.name.endswith(".csv") and entry.is_file():

                    t = np.append(t,int(re.search(r'\d+', entry.name).group()))
                    with open(entry.path, 'r') as f:
                        fields = []
                        alpha = []
                        x = []
                        csvreader = csv.reader(f)
                        fields = next(csvreader)
                        for row in csvreader:
                            alpha = np.append(alpha, float(row[0]))
                            x = np.append(x,float(row[1]))
                        idx = (np.abs(alpha - 0.5)).argmin()
                        alp = alpha[idx]
                        loc = x[idx]
                        # if alp == 1:
                        #     loc = 0.0
                        #     alp = 1

                posField = np.append(posField, loc)
                alphaField = np.append(alphaField, alp)

                j = j + 1
            m = t.argsort()
            alphaField = alphaField[m]
            posField = posField[m]
            # Analytical Interface position 
            cps = 2050
            k1 = 2.22
            rho1 = 916.2
            a1 = k1/(rho1*cps)
            lambd = 0.00032622525325939834
            lambd1 = 0.2299545377262345
            psi = []
            psi2 = []
            t = np.sort(t)
            for x in range(0,len(t)):
                psi = np.append(psi, lambd*sqrt(t[x]))
                psi2 = np.append(psi2, 2*lambd1*sqrt(a1*t[x]))

        numL2 = []
        denL2 = []
        L2 = []
        for k in range(0,len(t)):
            numL2 = np.append(numL2, np.sum(np.power((posField[k]-psi2[k]),2)))
            L2 = numL2

    f1 = plt.figure()
    f2 = plt.figure()

    ax1 = f1.add_subplot(111)
    ax1.plot(t, psi2.real, 'r--', label='Neumann solution')
    ax1.plot(t, posField, 'g--', label='Numerical solution')
    ax1.set(xlabel='t [s]', ylabel= 'Interface position [m]')
    ax1.grid(True)
    ax1.legend()
    L2 = L2.real
    ax2 = f2.add_subplot(111)
    ax2.plot(t, L2, 'r--', label='Relative error')
    ax2.set(xlabel='t [s]', ylabel= 'Relative error')
    ax2.grid(True)

    ax2.legend()
    plt.show()
 
 