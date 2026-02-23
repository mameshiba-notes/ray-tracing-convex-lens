import numpy as np
import matplotlib.pyplot as plt
#凸レンズを用意
a = 11.0
r = 13.0

theta = np.arccos(a/r)
theta_list=np.linspace(-1*theta,theta,400)
x0 = r*np.cos(theta_list)-a
y0 = r*np.sin(theta_list)
plt.plot(x0,y0,'k-')
x1 = r*np.cos(np.pi-theta_list)+a
y1 = r*np.sin(np.pi+theta_list)
plt.plot(x1,y1,'k-')
plt.axis([-40,40,-10,10])
plt.plot([-40.0,40.0],[0.0,0.0],'k-')
plt.gca().set_aspect('equal',adjustable='box')
#点光源の位置
p0=-27.0
q0=5.12
#r*np.sin(np.pi-theta*2/3)

plt.plot(p0,q0,".")
plt.text(p0-5,q0+1,f"({p0:.1f},{q0:.2f})")

phi_list=np.linspace(np.pi-theta*2/3,np.pi+theta*2/3,3)
p1 = r*np.cos(phi_list)+a
q1 = r*np.sin(phi_list)

n0 = 1.0
n1 = 1.5
#光路の交点を求める用に、それぞれの点を保管する箱を用意
lines = []
#凸レンズから出る光の始点の座標を保管する用の箱を用意
p2_list=[]
q2_list=[]

for p1x,q1y in zip(p1,q1):
    plt.plot([p0,p1x],[q0,q1y],color='blue')
    d0x = p0-p1x
    d0y = q0-q1y
    d0norm = np.sqrt(d0x**2+d0y**2)
    d0x_hat = d0x/d0norm
    d0y_hat = d0y/d0norm
    nx = p1x-a
    ny = q1y
    nnorm = np.sqrt(nx**2+ny**2)
    nx_hat = nx/nnorm
    ny_hat = ny/nnorm
    d0dot = d0x_hat*nx_hat+d0y_hat*ny_hat
    d1x = (-1*n0/n1*(d0x_hat-d0dot*nx_hat)-
           nx_hat*np.sqrt(1-((n0/n1)**2)
                          *(1-d0dot**2)))  
    d1y = (-1*n0/n1*(d0y_hat-d0dot*ny_hat)-
           ny_hat*np.sqrt(1-((n0/n1)**2)
                          *(1-d0dot**2)))      
    t1_line = np.linspace(0,20,100)
    d1x_line = p1x+t1_line*d1x
    d1y_line = q1y+t1_line*d1y
    mask1 = (d1x_line+a)**2+d1y_line**2<r**2
    plt.plot(d1x_line[mask1],d1y_line[mask1],
             color='blue')
    A = d1x**2+d1y**2
    B = 2*(d1x*(p1x+a)+q1y*d1y)
    C = (p1x+a)**2+q1y**2-r**2
    #判別式
    D = B**2-4*A*C
    if(D>=0):
        t1_1 = (-1*B+np.sqrt(B**2-4*A*C))/(2*A)
        t1_2 = (-1*B-np.sqrt(B**2-4*A*C))/(2*A)
        t1 = min(t for t in [t1_1,t1_2] if t>0)
    p2 = p1x+t1*d1x
    q2 = q1y+t1*d1y
    
    p2_list.append(p2)
    q2_list.append(q2)
    
    d1x_prime = p1x-p2
    d1y_prime = q1y-q2
    d1norm = np.sqrt(d1x_prime**2+d1y_prime**2)
    d1x_hat = d1x_prime/d1norm
    d1y_hat = d1y_prime/d1norm
    nx_prime = -(p2+a)
    ny_prime = -q2
    nprimenorm = np.sqrt(nx_prime**2+ny_prime**2)
    nxprime_hat = nx_prime/nprimenorm
    nyprime_hat = ny_prime/nprimenorm
    d1dot = d1x_hat*nxprime_hat+d1y_hat*nyprime_hat
    d2x = (-1*n1/n0*(d1x_hat-d1dot*nxprime_hat)-
           nxprime_hat*np.sqrt(1-((n1/n0)**2)
                          *(1-d1dot**2))) 
    d2y = (-1*n1/n0*(d1y_hat-d1dot*nyprime_hat)-
           nyprime_hat*np.sqrt(1-((n1/n0)**2)
                          *(1-d1dot**2))) 
    t2_line = np.linspace(0,50,100)
    d2x_line = p2+t2_line*d2x
    d2y_line = q2+t2_line*d2y
    plt.plot(d2x_line,d2y_line,color='blue')
    lines.append((p2,q2,d2x,d2y))
#光路の交点を求める
(x1,y1,v1x,v1y) = lines[0]
(x2,y2,v2x,v2y) = lines[1]
E = np.array([[v1x,-v2x],
             [v1y,-v2y]])             
F = np.array([x2-x1,y2-y1])
t, s = np.linalg.solve(E, F)
x = x1 + t*v1x
y = y1 + t*v1y
plt.plot(x,y,'k-')
plt.text(x-5, y-3, f"({x:.1f}, {y:.1f})")

#座標を表示
#for px_1, py_1 in zip(p1, q1):
#    plt.text(px_1-15, py_1+2, 
#             f"({px_1:.3f}, {py_1:.3f})")

#for px_2, py_2 in zip(p2_list, q2_list):
#    plt.text(px_2+1, py_2+2, 
#             f"({px_2:.3f}, {py_2:.3f})")

plt.text(
    0.03,0.55,
    f"n0={n0}\nn1={n1}\nr={r}\na={a}",
    transform = plt.gca().transAxes,
    verticalalignment='top',
    bbox = dict(facecolor = 'white',alpha=0.8)
)
 
plt.show()

