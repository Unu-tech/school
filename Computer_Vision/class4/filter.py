import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(0.0, 10.0, num=101)
y = np.random.rand(101)

y = 25 * y
y[20:40] = 200 + y[20:40]
y[60:80] = 120 + y[60:80]
y = y.astype(int)

low_pass = [1 / 3, 1 / 3, 1 / 3]
high_pass = [1/2, 0, -1/2]

y_low = np.convolve(y, low_pass, 'same')
y_high = np.convolve(y, high_pass, 'same')

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))
fig.suptitle(f"Visualization of low/high pass filtering.")

ax1.set_title("A slice of a 2d image")
ax2.set_title("Low pass filter (Box/Blur)")
ax3.set_title("High pass filter (1st order derivative)")

ax1.plot(x,y)
ax2.plot(x,y_low)
ax3.plot(x,y_high)

plt.subplots_adjust(wspace=0, hspace=0)
plt.savefig("Filtering.png")
plt.show()
