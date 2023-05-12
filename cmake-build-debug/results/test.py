import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Generate some data
x = np.arange(0, 2*np.pi, 0.01)
y = np.sin(x)

# Create the figure and axis
fig, ax = plt.subplots()

# Create a line object with the initial data
line, = ax.plot(x, y)

# Define the animation function
def animate(frame):
    # Update the y-data of the line
    line.set_ydata(np.sin(x + frame/10.0))
    return [line]

# Create the animation object
ani = animation.FuncAnimation(fig, animate, frames=200, interval=20, blit=True)

# Show the plot
plt.show()
