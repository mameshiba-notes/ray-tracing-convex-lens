import numpy as np
import matplotlib.pyplot as plt

d = 3
l = 8

alpha = l**2
beta = l**2-d**2
theta = np.linspace(-np.pi/2,np.pi/2,400)
x = np.sqrt(alpha)*np.cos(theta)
y = np.sqrt(beta)*np.sin(theta)
plt.plot(x,y,'k-')
plt.axis([-10,12,-6,6])
plt.xticks(np.arange(-10, 14, 2))
plt.yticks(np.arange(-8, 10, 2))

a0 = -1*d
b0 = 0

theta0_list = np.linspace(-2*np.pi/6, 2*np.pi/6, 5)
p0 = np.sqrt(alpha) * np.cos(theta0_list)
q0 = np.sqrt(beta) * np.sin(theta0_list)

for px, qy in zip(p0, q0):
    plt.plot([a0, px], [b0, qy],'b-')
    nx = px/alpha
    ny = qy/beta
    norm = np.sqrt(nx**2+ny**2)
    nx_hat = nx/norm
    ny_hat = ny/norm
    dx = px-a0
    dy = qy-b0
    dot = dx*nx_hat+dy*ny_hat
    rx = dx-2*dot*nx_hat
    ry = dy-2*dot*ny_hat
    t = np.linspace(0, 5, 100)
    rx_line = px + t * rx
    ry_line = qy + t * ry
    plt.plot(rx_line, ry_line,'c-')

plt.scatter(a0, b0,color='black')
plt.text(d+0.5,0.5,f"({d},0)")
plt.scatter(d, 0,color='black')
plt.text(a0-2,b0+0.5,f"({a0},{b0})")

plt.text(
    0.05, 0.95,
    f"d = {d}\nl = {l}",
    transform=plt.gca().transAxes,
    verticalalignment='top',
    bbox=dict(facecolor='white', alpha=0.8)
)
plt.plot([], [], color='blue', label='Incident ray')
plt.plot([], [], color='cyan', label='Reflected ray')
plt.legend()

plt.gca().set_aspect('equal',adjustable='box')
plt.show()
