# S16work
Python scripts from working with KMOS and SAMI data during spring 2016 work term on my laptop workstation. Assume all scripts work with KMOS data unless it indicates otherwise. I apologize in advance if my scripts are disgusting, I tend to write code on a need-to-use/do basis without concern of future readability, to the disadvantage of my future self and others.You will most likely not need to use these scripts as they only produce products as Balogh requested, so each are very specific to the request, but they can be used as a guideline for any future requests to see how to work with the data, produce plots, write fits files after some data manipulation, etc.


A good starting point for understanding everything I've done would be to read through my formal documentation in ./formal\_documentation/. There I describe in a kind of formal journal how I got started with the data and what I did on a day-by-day basis.


Related documents/articles to read:

Proposal: in this directory, called proposal.pdf

KMOS instrument explanation (overview, description, etc.): https://www.eso.org/sci/facilities/paranal/instruments/kmos.html

esoreflex installation (uses KEPLER to visually portray esorex pipeline workflows. Used to reduce KMOS data.): http://www.eso.org/sci/software/pipelines/reflex_workflows/

KMOS esoreflex tutorial: ftp://ftp.eso.org/pub/dfs/pipelines/kmos/kmos-reflex-tutorial-1.6.pdf

KMOS pipeline manual: ftp://ftp.eso.org/pub/dfs/pipelines/kmos/kmos-pipeline-manual-2.18.pdf

KMOS pipeline cookbook: ftp://ftp.eso.org/pub/dfs/pipelines/kmos/kmos-pipeline-cookbook-1.5.pdf

SAMI instrument front-page: http://sami-survey.org/

SAMI Early Data Release (this is the data used): http://sami-survey.org/edr

ALFALFA survey data: http://egg.astro.cornell.edu/alfalfa/data/

K98 (SFR from Halpha relation) article: http://www.annualreviews.org/doi/pdf/10.1146/annurev.astro.36.1.189


See bookmarks.html for more documents/articles. I would recommend importing it to firefox.


These are the scripts designed for use on my laptop. Will port them to the desktop later for others to use. Some scripts may not work on desktop, most likely due to me not changing the directory names in the scripts, which is mostly the same between my laptop and the desktop except for the following: on laptop, my working directory is /home/rburnet/ whereas on the desktop my working directory is /home/rcburnet/. All you'd most likely need to change is that directory naming (rburnet to rcburnet) to use the scripts on the desktop.


./manual\_sky\_subtraction\_scripts/ contains python scripts that:

1) Calculate the fraction of T\_lambda/S\_lambda for sky arms and target arms (target flux from target arms over sky flux from corresponding sky arms)

2) Plot the spectrum of manual sky subtraction carried out in the scripts


./sky\_subtraction\_look\_scripts/ contains the python scripts that look into sky and object cubes and plot their residual fluxes on the same plot to compare the two


./total\_flux\_vs\_wavelength\_scripts/ contains the python scripts that will plot the spectrum of the data cubes. Read the comments at the beginning of each script to see exactly what each script does.

./diagnostic\_of\_objects/ contains the python script that will create the plots used by diagnostic\_report\_OB(1|2).tex to generate diagnostic reports of the pipeline data comparing the calculated Z magnitude values to the expected values and place the plots of the HST image of the target, the pipeline's output of the target before sky substraction, and the pipeline's output of the target after sky subtraction. They also contain scripts for local sky subtraction (local\_sky\_subtraction) where the scripts and figures of my attempt at using the local sky around an elliptical aperture for each target and subtracting the flux of that sky from the each spaxel to ensure the aperture flux is positive and not negative (like it is right now for CL0034-IZ target 11) are. Started June 22. See Balogh's email on June 15 for information. It also contains the Magnitude plots (M\_IZ vs M\_Z) scripts and figures and a script to calculate expected Halpha flux from SFR(OII) as detailed in Sean's txt table files (SFR\_to\_Halpha.py) as well as the new report with Halpha flux values reported and spectrum plots integrated (new\_report).


./figures/ contain all the figues that are saved by each of the above scripts except when indicated that the figures are in the scripts' respective directories.


./headers/ contain the headers of some files for easy viewing of arm names, which arms are locked, and other observation information.


./txt\_tables/ contains the txt tables provided by Dr. Balogh which designate how bright the target in each arms is, their redshift, and other user data. It also contains a script on how to extract info from the txt tables if needed.

./SAMI/ contains the scripts for my second project with SAMI survey data.

./formal\_documentation/ contains the documentation of everything that I have done to the best of my memory. Basically a rewriting of my journal for the term.


Fwd_ Re_ kmos data.zip is the compressed .zip file from which the txt tables came from.


write\_obs\_mode\_to\_header.py is a script used to change the header information of all raw fits files from nod\_to\_sky mode to stare mode.
