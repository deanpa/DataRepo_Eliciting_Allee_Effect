
"""
Parameters for the pheromone model
"""

import datetime



class PheromoneParams(object):
    def __init__(self):
        self.extentMask = None
        self.trapsFile = None

        self.startDate = datetime.date(2021, 8, 1)
#        self.endDate = datetime.date(2022, 8, 15)
        self.endDate = datetime.date(2024, 2, 15)

        self.pheromoneReleaseDayMonths = [(15, 9), (15, 11)]
#        self.pheromoneReleaseDayMonths = []

        self.decoySpacing = [300, 301] # [300, 1000]
        self.estrousStartDayMonth = (15, 9)
        self.estrousEndDayMonth = (15, 1)
        self.dispersalDateDayMonth = (16, 1)

        self.hoursPerDay = 16           # 30 min steps - change to 16
        
        ## Range of initial stoats to add
        self.meanNAdd = [6,7] #[6, 12]

        self.stepScale = 50.0               # Weibull scale
        self.stepShape = 0.9                # Weibull shape

        # directional or biased random walk to find mates
        self.directionalVM = 3.5 #1.5

        self.alphaK = [0.02, 0.0201] # [0.01, 0.1]    # 0.075

        self.COA_radius = 2500      # 1000 # metres

        self.COA_decay_spatial =  [0.01, 0.0101]  # [0.001, 0.01]    # [0.001, 0.02]      # delta parameter
        self.COA_decay_temporal = [0.005, 0.006] # [0.005, 0.05]      # gamma parameter
        
        self.minK = 0.005
        
        self.meanRecruits = 9
        
        self.habituationDays =  [19, 20] # [8, 20] 

        self.encounterDistance = 25  # 25 metres

        ## probability of pregnacy given an encounter
        self.probPregnacy = 0.9

        self.birthDayMonth = (30, 10)                   # when multi-years need to populate
        self.nDaysPregnantBeforeBirth = 180 # 6 months

#        self.trappingDayMonths = None
        self.trappingDayMonths = [(20, 1), (20, 7), (20, 11)]
        self.nTrapDays = 14
        self.trapEncDist = 15.0
        self.trapProbRemoval = 0.20
        self.PAnnualSurv = [0.5, 0.501] # [0.45, 0.61]        # [0.45, 0.55]

    def setDecoySpacing(self, minRes, maxRes):
        self.decoySpacing = [minRes, maxRes]

    def setExtentMask(self, mask):
        self.extentMask = mask
        
    def setTrapsFile(self, filename):
        self.trapsFile = filename

    def setResultsFile(self, filename):
        self.resultsFile = filename