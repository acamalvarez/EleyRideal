# popt2, pcov = optimize.curve_fit(r_SCR2, [Cu12AEZ_O2_in, Cu12AEZ_NO_in, Cu12AEZ_NH3_in], np.array(
#     Cu12AEZ['rateMeasured']), p0=[1e-5, 1e-5, 1e-5])

# print(popt2)


# new_rate2 = r_SCR2([Cu12AEZ_O2_in, Cu12AEZ_NO_in, Cu12AEZ_NH3_in], *popt2)

# line = [0.7e-8, 2.7e-8]

# plt.plot(line, line, '-')
# plt.plot(Cu12AEZ['rateMeasured'], new_rate2, 's', fillstyle='none')

# plt.xlim(line)
# plt.ylim(line)

# set_labels('r measured', 'r fitted')
# plt.show()


# sstot = np.sum((Cu12AEZ['rateMeasured'] - np.mean(Cu12AEZ['rateMeasured']))**2)
# ssres = np.sum((Cu12AEZ['rateMeasured'] - new_rate2)**2)
# R2 = 1 - ssres / sstot

# print(R2)
