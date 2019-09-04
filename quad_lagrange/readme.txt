This script calculates stress and displacment of scaled boundary finite element for square plate with square hole problem. 
Run quad_lagrange_script.m to run the script. There is no INPUT, the OUTPUT (displacement and stress) is stored in variable named 'local'.
After finish, it will plot the stress (the last two lines).
To plot the stress, run plot_stress(local, local) with local is the result from the solver.
To plot the displacement, run plot_disp(local).