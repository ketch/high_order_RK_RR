#!/usr/bin/env python
# encoding: utf-8
import numpy as np

def setaux(x,rhoB=4,KB=4,rhoA=1,KA=1,alpha=0.5,xlower=0.,xupper=600.,bc=2):
    aux = np.empty([3,len(x)],order='F')
    xfrac = x-np.floor(x)
    #Density:
    aux[0,:] = rhoA*(xfrac<alpha)+rhoB*(xfrac>=alpha)
    #Bulk modulus:
    aux[1,:] = KA  *(xfrac<alpha)+KB  *(xfrac>=alpha)
    aux[2,:] = 0. # not used
    return aux

    
def stegoton(rkm,tol,iplot=0,htmlplot=0,outdir='./_output'):
    """
    Stegoton problem.
    Nonlinear elasticity in periodic medium.
    See LeVeque & Yong (2003).

    $$\\epsilon_t - u_x = 0$$
    $$\\rho(x) u_t - \\sigma(\\epsilon,x)_x = 0$$
    """

    from clawpack import pyclaw

    solver = pyclaw.SharpClawSolver1D()
    solver.time_integrator = 'RK'
    solver.A = rkm.A
    solver.b = rkm.b
    solver.b_hat = rkm.bhat
    solver.c = rkm.c
    solver.error_tolerance = tol
    solver.dt_variable = True
    solver.cfl_max = 2.5
    solver.cfl_desired = 1.5
    solver.weno_order=5

    from clawpack import riemann
    solver.rp = riemann.rp1_nonlinear_elasticity_fwave

    solver.bc_lower[0] = pyclaw.BC.periodic
    solver.bc_upper[0] = pyclaw.BC.periodic

    #Use the same BCs for the aux array
    solver.aux_bc_lower = solver.bc_lower
    solver.aux_bc_upper = solver.bc_lower

    xlower=0.0; xupper=300.0
    cellsperlayer=6; mx=int(round(xupper-xlower))*cellsperlayer
    x = pyclaw.Dimension('x',xlower,xupper,mx)
    domain = pyclaw.Domain(x)
    num_eqn = 2
    state = pyclaw.State(domain,num_eqn)

    #Set global parameters
    alpha = 0.5
    KA    = 1.0
    KB    = 4.0
    rhoA  = 1.0
    rhoB  = 4.0
    state.problem_data = {}
    state.problem_data['t1']    = 10.0
    state.problem_data['tw1']   = 10.0
    state.problem_data['a1']    = 0.0
    state.problem_data['alpha'] = alpha
    state.problem_data['KA'] = KA
    state.problem_data['KB'] = KB
    state.problem_data['rhoA'] = rhoA
    state.problem_data['rhoB'] = rhoB
    state.problem_data['trtime'] = 25000.0
    state.problem_data['trdone'] = False

    #Initialize q and aux
    xc=state.grid.x.centers
    state.aux=setaux(xc,rhoB,KB,rhoA,KA,alpha,xlower=xlower,xupper=xupper)

    a2=1.0
    sigma = a2*np.exp(-((xc-xupper/2.)/5.)**2.)
    state.q[0,:] = np.log(sigma+1.)/state.aux[1,:]
    state.q[1,:] = 0.

    tfinal=100.; num_output_times = 10;

    solver.max_steps = 5000000
    solver.fwave = True 
    solver.num_waves=2

    solver.lim_type = 2
    solver.char_decomp=0

    claw = pyclaw.Controller()
    claw.keep_copy = True
    claw.output_style = 1
    claw.num_output_times = num_output_times
    claw.tfinal = tfinal
    claw.solution = pyclaw.Solution(state,domain)
    claw.solver = solver

    # Solve
    status = claw.run()
    print claw.solver.status['totalsteps']

    if htmlplot:  pyclaw.plot.html_plot(outdir=outdir)
    if iplot:     pyclaw.plot.interactive_plot(outdir=outdir)

    epsilon = claw.frames[-1].q[0,:]
    aux = claw.frames[1].aux
    from clawpack.riemann.rp_nonlinear_elasticity import sigma 
    stress = sigma(epsilon,aux[1,:])
    dx = state.grid.delta[0]
    work = solver.status['totalsteps']

    return stress, dx, work


if __name__=="__main__":
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(stegoton)

def plot_efficiency(err,work,ork):
    import matplotlib.pyplot as plt
    p = ork.order()
    styles = ['k-s','b-o','g-d','r-p']
    plt.figure()
    plt.hold(True)
    for i in range(err.shape[0]):
        plt.loglog(err[i],work[i],styles[i],linewidth=3,markersize=11)

    plt.grid()
    plt.xlabel(r'Error',fontsize=15)
    plt.ylabel('Derivative evaluations',fontsize=15)
    plt.hold(False)
    leg=[]
    leg.append(ork.name)
    leg.append('Ex-Euler '+str(p)+'('+str(p-1)+')')
    leg.append('Ex-Midpoint '+str(p)+'('+str(p-2)+')')
    leg.append('DC-Euler '+str(p)+'('+str(p-1)+')')
    plt.legend(leg,loc='best',prop={'size':15,'weight':'bold'})
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    

    plt.xticks(np.array([1.e-10,1.e-8,1.e-6,1.e-4]))

def run_tolerances(tols,p=8):
    from nodepy import rk
    from nodepy import loadmethod
    
    err = np.zeros((4,len(tols)))
    work = np.zeros((4,len(tols)))

    rkm = rk.loadRKM('BS5').__num__()

    q_ref, dx, iwork = stegoton(rkm,tol=1.e-13)

    rkm = rk.extrap_pair(p/2,'midpoint').__num__()
    rke = rk.extrap_pair(p).__num__()
    rkd = rk.DC_pair(p-1,grid='cheb').__num__()

    if p == 6:
        ork = rk.loadRKM('CMR6')
    elif p == 8:
        ork = rk.loadRKM('DP8')
    elif p == 10:
        ork = loadmethod.load_rkpair_from_file('rk108curtis.txt')
    elif p == 12:
        ork = loadmethod.load_rkpair_from_file('rk129hiroshi.txt')

    for i, tol in enumerate(tols):
        q, dx, iwork = stegoton(ork,tol)
        err[ 0,i] = np.linalg.norm(q-q_ref,1)*dx
        work[0,i] = iwork*len(rkm)

        q, dx, iwork = stegoton(rke,tol)
        err[ 2,i] = np.linalg.norm(q-q_ref,1)*dx
        work[2,i] = iwork*len(rke)

        q, dx, iwork = stegoton(rkm,tol)
        err[ 1,i] = np.linalg.norm(q-q_ref,1)*dx
        work[1,i] = iwork*len(rkm)

        q, dx, iwork = stegoton(rkd,tol)
        err[ 3,i] = np.linalg.norm(q-q_ref,1)*dx
        work[3,i] = iwork*len(rkd)

    return err, work, ork
    
