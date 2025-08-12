#!/usr/bin/env python

import os
import pickle
from pheromone import calcresults
#from pheromone import params
from scipy.stats.mstats import mquantiles
import pylab
import numpy as np


def processResults():
    ### Output data path to write graphs, imagees and movies
    outputDataPath = os.path.join(os.getenv('PROJDIR', default='.'), 
            'pheromoneWork', 'Results', 'd1_multi')

    resultsDataPath = os.path.join(outputDataPath, 'results.pkl')

    simResultsFile = os.path.join(outputDataPath, 'simulationResults.csv')    



    results = calcresults.PheromoneResults.unpickleFromFile(resultsDataPath)

    calcresults.PheromoneResults.writeToFileFX(results, simResultsFile)


if __name__ == '__main__':
    processResults()
