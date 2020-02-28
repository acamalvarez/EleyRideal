import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import optimize

from functionsEleyRideal import (HC_effect, kineticData, r_SCR, r_SCR2, save_fig,
                                 select_data, set_labels, set_style)

set_style()
kineticData = kineticData()


# HC_effect(kineticData, ['NO_2', 'NO', 'NO'], [[0, 1, 2, 3, 4, 5, 6], [0, 1, 2, 3, 4, 5, 6], [
#           1, 2, 3, 4, 5, 6]]).plot('Cu12AEZ', 'NO_in', ['ln(NO)', 'ln(TOR)'])

# HC_effect(kineticData, ['NO', 'NO', 'NO'], [[0, 1, 2, 3, 4, 5, 6], [0, 1, 2, 3, 4, 5, 6], [
#           0, 1, 2, 3, 4, 5, 6]]).plot('Cu14SSZ1325', 'NO_in', ['ln(NO)', 'ln(TOR)'])

# HC_effect(kineticData, ['O2', 'O2', 'O2'], [[0, 1, 2, 3, 4, 5, 6], [0, 1, 2, 3, 4, 5, 6], [
#           0, 1, 2, 3, 4, 5]]).plot('Cu12AEZ', 'nominal', ['ln(O$_2$)', 'ln(TOR)'])

# HC_effect(kineticData, ['O2', 'O2', 'O2'], [[1, 2, 3, 4, 5, 6], [0, 1, 2, 3, 4, 5, 6], [
#           0, 1, 2, 3, 4, 5]]).plot('Cu14SSZ1325', 'nominal', ['ln(O$_2$)', 'ln(TOR)'])


# HC_effect(kineticData, ['NH3', 'NH3', 'NH3'], [[0, 1, 2, 3, 4], [0, 1, 2, 3, 4, 5, 6], [
#           0, 1, 2, 3, 4, 5]]).plot('Cu12AEZ', 'nominal', ['ln(O$_2$)', 'ln(TOR)'])

# HC_effect(kineticData, ['NH3', 'NH3', 'NH3'], [[0, 1, 2, 3, 4, 5, 6], [0, 1, 2, 3, 4, 5, 6], [
#           0, 1, 2, 3, 4, 5]]).plot('Cu14SSZ1325', 'nominal', ['ln(O$_2$)', 'ln(TOR)'])


# x = kineticData[(kineticData['catalyst'] == 'Cu12AEZ') & (kineticData['compound'] == 'NH3') & (kineticData['HC'] == 'None')]['NH3_in']
# y = kineticData[(kineticData['catalyst'] == 'Cu12AEZ') & (kineticData['compound'] == 'NH3') & (kineticData['HC'] == 'None')]['rateMeasured']

# x = select_data(np.array(x), [1,2,3,4])
# y = select_data(np.array(y), [1,2,3,4])

# plt.plot(x, y, 's', fillstyle='none')

# plt.show()


Cu12AEZ = kineticData[kineticData.catalyst.isin(['Cu12AEZ']) & kineticData.HC.isin(
    ['None']) & kineticData.compound.isin(['NO_2', 'NH3', 'O2'])]

Cu12AEZ_O2_in = np.array(Cu12AEZ[['O2_in', 'NO_in', 'NH3_in']])[:, 0]
Cu12AEZ_NO_in = np.array(Cu12AEZ[['O2_in', 'NO_in', 'NH3_in']])[:, 1]
Cu12AEZ_NH3_in = np.array(Cu12AEZ[['O2_in', 'NO_in', 'NH3_in']])[:, 2]


popt1, pcov = optimize.curve_fit(r_SCR, [Cu12AEZ_O2_in, Cu12AEZ_NO_in, Cu12AEZ_NH3_in], np.array(Cu12AEZ['rateMeasured']), p0=[1e-5, 1e-5, 1e-5, 1e-5])

print(popt1)


new_rate = r_SCR([Cu12AEZ_O2_in, Cu12AEZ_NO_in, Cu12AEZ_NH3_in], *popt1)

line = [0.7e-8, 2.7e-8]

plt.plot(line, line, '-')
plt.plot(Cu12AEZ['rateMeasured'], new_rate, 's', fillstyle='none')

plt.xlim(line)
plt.ylim(line)

set_labels('r measured', 'r fitted')
plt.show()


sstot = np.sum((Cu12AEZ['rateMeasured'] - np.mean(Cu12AEZ['rateMeasured']))**2)
ssres = np.sum((Cu12AEZ['rateMeasured'] - new_rate)**2)
R2 = 1 - ssres / sstot

print(R2)

