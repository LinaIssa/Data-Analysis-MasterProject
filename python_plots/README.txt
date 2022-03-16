This directory contains several python routines to make easy first plots of the shock results
and test conservations laws, energetics, and chemistry

To use these scripts, go in the directory containing the results of a run and just type
   python ../python_plots/script_plot_XXXX.py

Any script can be used by default or with optional arguments
   -h      shows an help message and exit
   -type   output format: file (pdf) or interacticve plot (itv)
   -axis   x axis choice: time (time) or distance (dist)
   -xmin   min x value
   -xmax   max x value

Example 
   python ../python_plots/script_plot_main.py -type pdf -axis dist -xmin 1e11 -xmax 1e19

List of available scripts
   script_plot_carr.py        =>   plot probable main Carbon, Oxygen, and Silicon carriers
   script_plot_cons.py        =>   test chemical conservation, conservation of H2 populations, and gas neutrality
   script_plot_ener.py        =>   display energetics output, test energy conservation
   script_plot_main.py        =>   plot temperature, density, H and H2, and ions and electron abundance profiles
   script_plot_profiles.py    =>   plot the abundances of a few chemical species
   script_plot_thchm.py       =>   plot thermal balance and focus on chemical   contributions
   script_plot_thmec.py       =>   plot thermal balance and focus on mechanical contributions
   script_plot_thrad.py       =>   plot thermal balance and focus on radiative  contributions
