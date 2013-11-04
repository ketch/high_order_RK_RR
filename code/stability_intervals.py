from nodepy import rk, loadmethod, stability_function
import matplotlib.pyplot as plt
import numpy as np


s_start=4
s_end=12
order = range(s_start,s_end+1)

I_R = np.zeros((4,len(order)))
I_I = np.zeros((4,len(order)))
I_R_scaled = np.zeros((4,len(order)))
I_I_scaled = np.zeros((4,len(order)))

for i,p in enumerate(order):

    print p

    if p == 4:
        ork = rk.loadRKM('Merson43')
    if p == 6:
        ork = rk.loadRKM('CMR6')
    elif p == 8:
        ork = rk.loadRKM('DP8')
    elif p == 10:
        ork = loadmethod.load_rkpair_from_file('rk108curtis.txt')
    elif p == 12:
        ork = loadmethod.load_rkpair_from_file('rk129hiroshi.txt')

    pp,q = ork.stability_function()
    I_R[0,i] = stability_function.real_stability_interval(pp,q)
    I_I[0,i] = stability_function.imaginary_stability_interval(pp,q)
    I_R_scaled[0,i] = I_R[0,i]/len(ork)
    I_I_scaled[0,i] = I_I[0,i]/len(ork)

    exe = rk.extrap_pair(p)
    pp,q = exe.stability_function()
    I_R[1,i] = stability_function.real_stability_interval(pp,q)
    I_I[1,i] = stability_function.imaginary_stability_interval(pp,q)
    I_R_scaled[1,i] = I_R[1,i]/len(exe)
    I_I_scaled[1,i] = I_I[1,i]/len(exe)

    if p%2==0:
        exm = rk.extrap_pair(p/2,'midpoint')
        pp,q = exm.stability_function()
        I_R[2,i] = stability_function.real_stability_interval(pp,q)
        I_I[2,i] = stability_function.imaginary_stability_interval(pp,q)
        I_R_scaled[2,i] = I_R[2,i]/len(exm)
        I_I_scaled[2,i] = I_I[2,i]/len(exm)

    dc  = rk.DC_pair(p-1)
    pp,q = dc.stability_function(formula='pow',use_butcher=True)
    I_R[3,i] = stability_function.real_stability_interval(pp,q)
    I_I[3,i] = stability_function.imaginary_stability_interval(pp,q)
    I_R_scaled[3,i] = I_R[3,i]/len(dc)
    I_I_scaled[3,i] = I_I[3,i]/len(dc)

styles = ['k-D','b-p','g-s','r-o']
leg=[]
leg.append('Reference RK')
leg.append('Ex-Euler ')
leg.append('Ex-Midpoint ')
leg.append('IDC-Euler ')

plt.figure()
plt.hold(True)
for i in range(4):
    if i == 2:
        plt.plot(order[::2],I_R[i,::2],styles[i],linewidth=3,markersize=14)
    else:
        plt.plot(order,I_R[i,:],styles[i],linewidth=3,markersize=14)
plt.hold(False)
plt.grid()
plt.xlabel('Order (p)')
plt.ylabel('Real stability interval')
plt.legend(leg,loc='best')
plt.ylim((0.0,6.5))
plt.savefig('real_stability_intervals.pdf')


plt.figure()
plt.hold(True)
for i in range(4):
    if i == 2:
        plt.plot(order[::2],I_I[i,::2],styles[i],linewidth=3,markersize=14)
    else:
        plt.plot(order,I_I[i,:],styles[i],linewidth=3,markersize=14)
plt.hold(False)
plt.grid()
plt.xlabel('Order (p)')
plt.ylabel('Imaginary stability interval')
plt.legend(leg,loc='best')
plt.savefig('imaginary_stability_intervals.pdf')

# Scaled figures
plt.figure()
plt.hold(True)
for i in range(4):
    if i == 2:
        plt.plot(order[::2],I_R_scaled[i,::2],styles[i],linewidth=3,markersize=14)
    else:
        plt.plot(order,I_R_scaled[i,:],styles[i],linewidth=3,markersize=14)
plt.hold(False)
plt.grid()
plt.xlabel('Order (p)')
plt.ylabel('$I_R/s$')
plt.legend(leg,loc='best')
plt.ylim((0.0,1.0))
plt.savefig('real_stability_intervals_scaled.pdf')

plt.figure()
plt.hold(True)
for i in range(4):
    if i == 2:
        plt.plot(order[::2],I_I_scaled[i,::2],styles[i],linewidth=3,markersize=14)
    else:
        plt.plot(order,I_I_scaled[i,:],styles[i],linewidth=3,markersize=14)
plt.hold(False)
plt.grid()
plt.xlabel('Order (p)')
plt.ylabel('$I_I/s$')
plt.legend(leg,loc='best')
plt.savefig('imaginary_stability_intervals_scaled.pdf')
