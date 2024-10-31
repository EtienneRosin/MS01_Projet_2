import numpy as np
import matplotlib.pyplot as plt

file_name = "validation/jacobi_asymptotic_acuracy_measurement.csv"

file_props = dict(dtype = float, delimiter = ',', names=True)

data = np.genfromtxt(fname = file_name, **file_props)

print(data['h'])

fig = plt.figure()
ax = fig.add_subplot()


ax.plot(np.log(data['h']), np.log(data['error']))
ax.plot(np.log(data['h']), np.log(data['h']**2))
plt.show()