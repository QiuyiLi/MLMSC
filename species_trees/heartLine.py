import numpy as np

import matplotlib.pyplot as plt

from math import sin



x = np.arange(-2, 2, 0.01)

y = abs(x)**(2/3)+.8*abs(4-x**2)**.5*np.sin(31.41593*x)

fig = plt.figure('heart line')

axes = fig.add_subplot(111)

axes.plot(x, y,'-r')

axes.axis('equal')

plt.title('heart line')

plt.show()