#This code shows crazy behavior of DC

from nodepy import runge_kutta_method as rk
from nodepy import ivp
from nodepy import conv
import pylab as pl
import numpy as np

print 'initializing...'
methods = []
for p in (6,8,10,12):
    methods.append(rk.extrap_pair(p))
print 'done'

myivp = ivp.detest('SB1');
tols = 10.**(-np.arange(2,14))
[work,EX]=conv.ptest(methods,myivp,tols,verbosity=1);

pl.clf()
pl.hold(True)
pl.loglog(EX[0,:],work[0,:],'b-o',linewidth=2.5)
pl.loglog(EX[1,:],work[1,:],'y-d',linewidth=2.5)
pl.loglog(EX[2,:],work[2,:],'c-p',linewidth=2.5)
pl.loglog(EX[3,:],work[3,:],'m-s',linewidth=2.5)
pl.xticks([1.e0,1.e-2,1.e-4,1.e-6,1.e-8,1.e-10,1.e-12], fontsize=15)
pl.yticks(fontsize=15)
pl.grid()
pl.xlabel('Error at $\mathbf{t}_\mathrm{\mathbf{final}}$', fontsize=15,fontweight='bold')
pl.ylabel('Cost', fontsize=18,fontweight='bold')
#pl.ylim([1.e2,1.e5])
pl.hold(False)
pl.legend(('Ex-Euler 6(5)','Ex-Euler 8(7)','Ex-Euler 10(9)','Ex-Euler 12(11)'),loc="best",prop={'size':15,'weight':'bold'})
#pl.title('ptest() on Orbit Equations SB1')
pl.savefig('extrap_internalstability.pdf')





