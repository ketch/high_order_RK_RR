This needs to be run against the PyClaw branch error_control.
To produce the plot in the paper, do

    err, work, ork = stegoton.run_tolerances([1.e-7,1.e-8,1.e-9,1.e-10])
    stegoton.plot_efficiency(err,work,ork)
