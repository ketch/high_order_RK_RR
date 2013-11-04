import matplotlib.pyplot as plt

def plot_efficiency(methods,err,work):
    styles = ['k-D','b-p','g-s','r-o']
    plt.figure()
    plt.hold(True)
    for i in range(len(methods)):
        plt.loglog(err[i],work[i],styles[i],linewidth=3,markersize=14)

    plt.grid()
    plt.xlabel(r'Error at $t_{final}$')
    plt.ylabel('Derivative evaluations')
    plt.hold(False)
    leg=[]
    leg.append(ork.name)
    leg.append('Ex-Euler '+str(p)+'('+str(p-1)+')')
    leg.append('Ex-Midpoint '+str(p)+'('+str(p-2)+')')
    leg.append('DC-Euler '+str(p)+'('+str(p-1)+')')
    plt.legend(leg,loc='best')

    plt.xticks(np.array([1.e-10,1.e-8,1.e-6,1.e-4,1.e-2,1.e0]))
    


def runtests(p,problem,parallel=False,tol=None):
    from nodepy import rk, ivp, conv
    import numpy as np
    from nodepy import loadmethod

    if tol is None:
        tol = [10.**(-m) for m in range(2,11)]

    # Load methods
    print 'constructing Ex-Euler'
    ex  = rk.extrap_pair(p);
    print 'constructing Ex-midpoint'
    exm = rk.extrap_pair(p/2,'midpoint');
    print 'constructing DC-Euler'
    dc  = rk.DC_pair(p-1);

    if p == 6:
        ork = rk.loadRKM('CMR6')
    elif p == 8:
        ork = rk.loadRKM('DP8')
    elif p == 10:
        ork = loadmethod.load_rkpair_from_file('rk108curtis.txt')
    elif p == 12:
        ork = loadmethod.load_rkpair_from_file('rk129hiroshi.txt')

    if problem == 'suite':
        myivp = ivp.detest_suite_minus()
    else:
        myivp = ivp.detest(problem);

    methods = [ork,ex,exm,dc]

    [work,err]   = conv.ptest(methods,myivp,tol,verbosity=1,parallel=parallel);

    plot_efficiency(methods, err, work)
    
    fname = 'eff_'+problem+'_'+str(p)
    if parallel:
        fname = fname + '_par'
        plt.ylabel('Sequential derivative evaluations')
    plt.savefig(fname+'.pdf')


