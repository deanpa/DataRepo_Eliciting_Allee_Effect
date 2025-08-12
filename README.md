# DataRepo_Eliciting_Allee_Effect
Supporting material for manuscript submitted to "Ecology and Evolution", entitled, "Can a mate-finding Allee effect be elicited to help eradicate invasive mammals?"

DATA
    (1) ressy.img:  Erdas Imagine raster file of Resolution Island at 10-metre resolution projected in epsg:27200.
    (2) ressyalldatatraploc5.csv:  Eastings and northing of traps used in study.

COMPUTER CODE
Python 3.7.3 was used on linux with the following packages:
	(1) Packages
 		(A) Numpy		1.16.2
		(B) Scipy		1.2.1
		(C) Numba		0.43.1
		(D) RIOS		1.4.8
		(E) gdal		2.4.0	
		(F) matplotlib	3.0.3
		(G) prettytable	0.7
	(2) Scripts
	    (A) startSimulation.py: script to set data and results directories, number of iteration, and intiate the simulation.
    	(B) pheromone/params.py: script to set parameters for simulation.
    	(C) pheromone/calculation.py: script that runs the simulation and writes results to directory. 
    	(D) postSimulation.py: script to initiate the processing of results from simulation.
    	(E) pheromone/calcresults.py: script with functions to process results.

MOVIES 
    (1) movie.wmv: Video file that demonstrates movement, survival, trapping, reproduction and decoy deployment. The following symbols are used in the video:
        (A) Green squares are traps
        (B) Small blue immobile squares that fade with time are the pheromone decoys
        (C) Small black filled circles are males exhibiting home range movement
        (D) Large black filled circles are males exhibiting mate-searching movement
        (E) Small blue filled circles are non-pregnant females exhibiting home range movement
        (F) Large blue filled circles are non-pregnant females exhibiting mate-searching movement
        (G) Small red filled circles are pregnant females exhibiting home range movement
        (H) Large stationary black or blue filled circles are nests with non-fertilised female kits
        (I) Stationary red squares are nests with non-fertilised female kits
    (2) movieNoDecoy.wmv:  Video files as above, but without deployment of pheromone decoys.
