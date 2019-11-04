This repo contains simple codes that are used for practicing and testing various things
The codes available are 
1. Python GUI to make a cfg file for SU2 (outdated now)
This is a very crude attempt at making a GUI to write config files for SU2.

Need python's tkinter package.

At the moment, the GUI can only write cfg files for some limited scenarios.

To use:
python muli_option_gui.py confile_file_name.cfg

A config file with the name specified in the command will be genrated. If no such file exists, a ne one is created and if a file by hat name already exists, it will be overwritten.

2. A 1-D scalar equation solver used to test various convection schemes. Convection schemes currently added include first and second order upwind, QUICK, central difference scheme and gamma difference (a type of bounded central difference)
