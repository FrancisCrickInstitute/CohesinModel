This is 3D MMC Model of Cohesin


Main files:

main_run		- runs a simulation
main_continue		- if a simulations was stopped, this file will resume it from the same place
end_run			- calculates and plots output



Classes:

CohesinDNA		- Main class. Describes cohesin, DNA and their interaction
Cohesin5		- Class that describes cohesin 
DNA			- Class that describes DNA
RandomRotation		- A class of random rotations used in MMC process


Output functions:

scroll_cdsave(cdsave,0)	- plots results of simulations in a figure in laboratory frame of reference. Use right and left arrows to scroll through iterations.
scroll_cdsave(cdsave,2)	- plots results of simulations in a figure in local reference frame of cohesin. Use right and left arrows to scroll through iterations.
save_movie		- saves results as .avi file






