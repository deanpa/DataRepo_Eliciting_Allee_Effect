#!/usr/bin/env python

import os
import multiprocessing
import pickle
from pheromone import calculation
from pheromone import calcresults
from pheromone import params
from rios.parallel import jobmanager
import resource

NITERATIONS = 400

# Use the same environment variable as RIOS to define the type of
# parallel processing.
# Default to the multiprocessing type.
JOBMGR_TYPE = os.getenv('RIOS_DFLT_JOBMGRTYPE', default='multiprocessing')

TMP_DIR = 'XXX'

def parallelRunModel(pars, results):
    """
    A slight variation on pheromone.runModel
    which makes the results a parameter so it 
    can be used with rios.parallel.
    """

    (eradicated, nAdd, decoySpacing, nDecoyDeplyment, alphaK, COA_decay_spatial,
        COA_decay_temporal, habituationDays, pDaySurv) = calculation.runModel(pars, save=False)
    results.eradicated = eradicated
    results.nAdd = nAdd
    results.decoySpacing = decoySpacing
    results.nDecoyDeplyment = nDecoyDeplyment
    results.alphaK = alphaK
    results.COA_decay_spatial = COA_decay_spatial
    results.COA_decay_temporal = COA_decay_temporal
    results.habituationDays = habituationDays
    results.pDaySurv = pDaySurv 
    results.iter = NITERATIONS

class PheromoneJobInfo(jobmanager.JobInfo):
    """
    Contains an implementation of RIOS's jobmanager.JobInfo
    for the pheromone model.
    """
    def __init__(self, pars):
        self.pars = pars

    def getFunctionParams(self):
        "make input suitable for parallelRunModel"
        results = calcresults.PheromoneResults()
        return self.pars, results

    def getFunctionResult(self, params):
        "output was the last parameter"
        return params[-1]
        
def runMultipleJobs(pars, resultsDataPath):
    # if using multiprocessing, run a job per cpu
    # otherwise (assume SLURM) run a job per iteration
    # not sure if this is correct
    if JOBMGR_TYPE == 'multiprocessing':
        nThreads = multiprocessing.cpu_count()
    else:
        nThreads = NITERATIONS
    jobmgrClass = jobmanager.getJobManagerClassByType(JOBMGR_TYPE)
    jobmgr = jobmgrClass(nThreads)

    # Home dir runs out of quota
    # I couldn't find a cluster-wide temp var created on Pan
    # so simply use this dir if it exists (assume we are on Pan)
    # otherwise leave as default.
    if os.path.isdir(TMP_DIR):
        jobmgr.setTempdir(TMP_DIR)

    jobInputs = []
    for i in range(NITERATIONS):
        jobInfo = PheromoneJobInfo(pars)
        jobInputs.append(jobInfo)
    # run all in parallel and collect results
    results = jobmgr.runSubJobs(parallelRunModel, jobInputs)

    # pickle results
    fileobj = open(resultsDataPath, 'wb')
    pickle.dump(results, fileobj, protocol=4) # so we get large file support
    fileobj.close()
    
if __name__ == '__main__':

    # DATA PATHS
    inputDataPath = os.path.join(os.getenv('PROJDIR', default='.'), 'pheromoneWork','Data')
    outputDataPath = os.path.join(os.getenv('PROJDIR', default='.'), 'pheromoneWork', 
            'Results', 'd1_multi')
    pars = params.PheromoneParams()

#    pars.decoySpacing = [250, 250, 250, 250]

    pars.setExtentMask(os.path.join(inputDataPath, 'ressy.img'))
    print('path', os.path.join(inputDataPath, 'ressy.img'))
    pars.setTrapsFile(os.path.join(inputDataPath, 'ressyalldatatraploc5.csv'))

    # TODO: should this be an environment variable?
    resultsDataPath = os.path.join(outputDataPath, 'results.pkl')

    runMultipleJobs(pars, resultsDataPath)

    maxMem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    print('Max Mem Usage', maxMem)