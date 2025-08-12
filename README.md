# DataRepo_Eliciting_Allee_Effect

Supporting material for the manuscript submitted to *Ecology and Evolution*, entitled:  
**"Can a mate-finding Allee effect be elicited to help eradicate invasive mammals?"**

## Data

1. **ressy.img**  
   - Erdas Imagine raster file of Resolution Island at 10-metre resolution.  
   - Projected in EPSG:27200.

2. **ressyalldatatraploc5.csv**  
   - Eastings and northings of traps used in the study.

## Computer Code

Python **3.7.3** was used on Linux with the following packages and scripts.

### Packages

1. **Numpy** – 1.16.2  
2. **Scipy** – 1.2.1  
3. **Numba** – 0.43.1  
4. **RIOS** – 1.4.8  
5. **GDAL** – 2.4.0  
6. **Matplotlib** – 3.0.3  
7. **PrettyTable** – 0.7

### Scripts

1. **startSimulation.py**  
   - Sets data and results directories, number of iterations, and initiates the simulation.

2. **pheromone/params.py**  
   - Sets parameters for simulation.

3. **pheromone/calculation.py**  
   - Runs the simulation and writes results to the directory.

4. **postSimulation.py**  
   - Initiates processing of results from the simulation.

5. **pheromone/calcresults.py**  
   - Functions to process results.

## Movies

1. **movie.wmv**  
   - Demonstrates movement, survival, trapping, reproduction, and decoy deployment.  
   - Symbols used in the video:  
     - Green squares – Traps  
     - Small blue immobile squares (fade with time) – Pheromone decoys  
     - Small black filled circles – Males with home range movement  
     - Large black filled circles – Males with mate-searching movement  
     - Small blue filled circles – Non-pregnant females with home range movement  
     - Large blue filled circles – Non-pregnant females with mate-searching movement  
     - Small red filled circles – Pregnant females with home range movement  
     - Large stationary black or blue filled circles – Nests with non-fertilised female kits  
     - Stationary red squares – Nests with non-fertilised female kits  

2. **movieNoDecoy.wmv**  
   - As above, but without deployment of pheromone decoys.

