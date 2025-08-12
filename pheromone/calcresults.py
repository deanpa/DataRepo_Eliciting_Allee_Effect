
import pickle
import numpy as np


class PheromoneResults(object):
    """
    Dummy class to take the parameters for the rios
    parallel processing. Needed in a separate module
    so the pickling still works.

    Contains the info for a single iteration
    """
    eradicated = None
    nAdd = None
    decoySpacing = None
    nDecoyDeplyment = None
    alphaK = None
    COA_decay_spatial = None
    COA_decay_temporal = None
    habituationDays = None
    pDaySurv = None
    iter = None

    def pickleSelf(self, fname):
        fileobj = open(fname, 'wb')
        pickle.dump(self, fileobj, protocol=4) # so we get large file support
        fileobj.close()

    @staticmethod
    def unpickleFromFile(fname):
        fileobj = open(fname, 'rb')
        data = pickle.load(fileobj)
        fileobj.close()
        return data


    def writeToFileFX(results, simResultsFile):
        n = len(results)
        # create new structured array with columns of different types
        structured = np.empty((n,), dtype=[('Eradicated', 'bool'), 
            ('nAdd', np.integer), ('Decoy spacing', np.float), ('N deployments', np.integer),
            ('Alpha K', np.float), ('COA spatial decay', np.float),
            ('COA temporal decay', np.float), ('Habituation days', np.float),
            ('Daily P(survive)', np.float)])
        for i in range(n):
            # copy data over
            structured['Eradicated'][i] = results[i].eradicated
            structured['nAdd'][i] = results[i].nAdd
            structured['Decoy spacing'][i] = results[i].decoySpacing
            structured['N deployments'][i] = results[i].nDecoyDeplyment
            structured['Alpha K'][i] = results[i].alphaK
            structured['COA spatial decay'][i] = results[i].COA_decay_spatial
            structured['COA temporal decay'][i] = results[i].COA_decay_temporal
            structured['Habituation days'][i] = results[i].habituationDays
            structured['Daily P(survive)'][i] = results[i].pDaySurv
            np.savetxt(simResultsFile, structured, fmt=['%s', '%.0f', '%.0f', '%.0f', '%.4f', '%.6f', '%.4f', '%.1f', '%.4f'],
                comments = '', delimiter = ',', 
                header='eradicated, nAdd, decoySpacing, nDeploy, alphaK, coaSpatialDecay, coaTempDecay, habituationDays, pSurvive')

        pEradication = np.sum(structured['Eradicated']) / n

        print('########### PROBABILITY OF ERADICATION: ', pEradication)