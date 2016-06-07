# S16work
Scripts from working with KMOS data during spring 2016 work term.

These are the scripts designed for use on my laptop.


./manual\_sky\_subtraction\_scripts/ contains python scripts that:

	1) Calculate the fraction of T_lambda/S_lambda for sky arms and target arms 

	2) Plot the spectrum of manual sky subtraction carried out in the scripts


./sky\_subtraction\_look\_scripts/ contains the python scripts that look into sky and object cubes and plot their residual fluxes on the same plot to compare the two


./total\_flux\_vs\_wavelength\_scripts/ contains the python scripts that will plot the spectrum of the data cubes. Read the comments at the beginning of each script to see exactly what each script does.

./diagnostic\_of\_objects/ contains the python script that will create the plots used by diagnostic\_report\_OB(1|2).tex to generate diagnostic reports of the pipeline data comparing the calculated Z magnitude values to the expected values and place the plots of the HST image of the target, the pipeline's output of the target before sky substraction, and the pipeline's output of the target after sky subtraction.


./figures/ contain all the figues that are saved by each of the above scripts.


./headers/ contain the headers of some files for easy viewing of arm names, which arms are locked, and other observation information.


./txt\_tables/ contains the txt tables provided by Dr. Balogh which designate how bright the target in each arms is, their redshift, and other user data. It also contains a script on how to extract info from the txt tables if needed.


Fwd_ Re_ kmos data.zip is the compressed .zip file from which the txt tables came from.


write\_obs\_mode\_to\_header.py is a script used to change the header information of all raw fits files from nod\_to\_sky mode to stare mode.
