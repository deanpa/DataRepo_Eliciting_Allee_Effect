#!/usr/bin/env python

import os
import datetime
import tempfile
import subprocess
import pickle
import shutil
from numba import njit, objmode
import numpy as np
from osgeo import gdal


# the dtype we use for the stoat array
STOAT_DTYPE = [('deleted', np.bool), ('male', np.bool), ('pregnant', np.bool), 
                    ('x', np.float64), ('y', np.float64), ('home_x', np.float64),
                    ('home_y', np.float64), 
                    ('id', np.int32), # uniquely generated id for this stoat
                    ('pregnant_day', np.int32), # if pregnant, day she became inpregnated
                    ('parentid', np.int32), # -1 when individual has dispersed
                    ('bearingT_1', np.float64), # bearing of previous time step for directed search
                    ('homerange', np.bool) # True if homerange state, or False for searching. For movie.
                    ]
INITIAL_STOAT_ARRAY_SIZE = 1000

PHEROMONE_DTYPE = [('x', np.float64), ('y', np.float64)]
# keep a track of which pheromones each stoat has intereacted with
PHEROMONE_INTERACTION_DTYPE = [('pheromoneid', np.int32), # index of pheromone in pheromoneArray
                    ('stoatid', np.int32), # 'id' from stoatArray
                    ('ndays', np.int32) # number of days until interaction no longer applies
                    ]

# keep a track of which stoats have mated
MATING_DTYPE = [('maleid', np.int32), # 'id' from stoatArray
                ('femaleid', np.int32), # 'id' from stoatArray
                ('ndays', np.int32)  # number of days until interaction no longer applies
                ]

def readTrapsFile(filename):
    """
    Read the traps file and return an array with the data
    """
    data = np.loadtxt(filename, skiprows=1, usecols=(3,4), delimiter=',')
#    data = np.loadtxt(filename, skiprows=1, usecols=(0, 1), delimiter=',')
    return data

def makePheromoneArray(mask, tlx, tly, brx, bry, pixsize, spacing, transform):
    tinverse = gdal.InvGeoTransform(transform)

    xList = []
    yList = []
    for y in range(int(bry + spacing), int(tly), int(spacing)):
        for x in range(int(tlx), int(brx), int(spacing)):
            xPix, yPix = gdal.ApplyGeoTransform(tinverse, x, y)
            xPix = int(xPix)
            yPix = int(yPix)
            if mask[yPix, xPix] > 0:
                xList.append(x)
                yList.append(y)

    data = np.empty(len(xList), dtype=PHEROMONE_DTYPE)
    data['x'] = xList
    data['y'] = yList
    
    return data

@njit
def checkLocationIsOnIsland(mask, tlx, tly, brx, bry, pixsize, x, y):
    # inside the masked file?
    isInsideIsland = (x >= tlx and x <= brx and y <= tly and y >= bry)
    if isInsideIsland:
        # now check the actual mask
        xPix = int(np.round((x - tlx) / pixsize))
        yPix = int(np.round((tly - y) / pixsize))
        isInsideIsland = (mask[yPix, xPix] > 0)
    return isInsideIsland

@njit
def createRandomLocationOnIsland(mask, tlx, tly, brx, bry, pixsize):
    xsize = brx - tlx
    ysize = tly - bry

    x = 0
    y = 0
    isInsideIsland = False
    while not isInsideIsland:
        # choose a spot at random.
        x = tlx + np.random.random() * xsize
        y = bry + np.random.random() * ysize
        # inside the masked file?
        isInsideIsland = checkLocationIsOnIsland(mask, tlx, tly, 
                                brx, bry, pixsize, x, y)
    return x, y

@njit
def createInitialStoats(stoatArray, nAdd, mask, tlx, tly, brx, bry, pixsize, 
                            nDaysPregnantBeforeBirth):
    """
    Put down some initial stoats within the masked area.
    """
    Current_Id = 0

    for nStoats in range(nAdd):
        x, y = createRandomLocationOnIsland(mask, tlx, tly, brx, bry, pixsize)
        male = np.random.random() < 0.5
        ## make first two a female and male
        if nStoats == 0:
            stoatArray[nStoats]['male'] = False     # female
        elif nStoats == 1:
            stoatArray[nStoats]['male'] = True     # male
        else:
            stoatArray[nStoats]['male'] = male

        stoatArray[nStoats]['deleted'] = False
        if not stoatArray[nStoats]['male']:
            stoatArray[nStoats]['pregnant'] = True # all females start pregnant 
            # make them pregnant ~6 months ago so they give birth straight away
            stoatArray[nStoats]['pregnant_day'] = -nDaysPregnantBeforeBirth
        else:
            stoatArray[nStoats]['pregnant'] = False
            stoatArray[nStoats]['pregnant_day'] = -1
        stoatArray[nStoats]['x'] = x
        stoatArray[nStoats]['y'] = y
        stoatArray[nStoats]['home_x'] = x
        stoatArray[nStoats]['home_y'] = y
        stoatArray[nStoats]['homerange'] = False
        stoatArray[nStoats]['id'] = Current_Id
        Current_Id += 1
        stoatArray[nStoats]['parentid'] = -1  # these 'existing' stoats have parentid = -1

        ## INITIAL NUMBER OF STOATS
        nStoats = Current_Id

    return nStoats, Current_Id

@njit
def checkForEstrousMatesAndDecoysInRadius(lookingForMale, x, y, stoatArray, nStoats, stoatid, COA_radius, 
            COA_decay_spatial, COA_decay_temporal, minK, daysSincePheromoneRelease, pheromoneArray,
            pheromoneInteractionArray, matingArray):

    arraySize = INITIAL_STOAT_ARRAY_SIZE + pheromoneArray.shape[0]
    Kappac = np.zeros(arraySize, dtype=np.float64)
    # temp arrays
    xCoords = np.zeros(arraySize, dtype=np.float64)
    yCoords = np.zeros(arraySize, dtype=np.float64)

    result = None
    nStoatsInTmp = 0
    for i in range(nStoats):
        if (not stoatArray[i]['deleted'] and stoatArray[i]['parentid'] == -1 
                    and stoatArray[i]['male'] == lookingForMale):
            # estrous other gender
            xdist = x - stoatArray[i]['x']
            ydist = y - stoatArray[i]['y']
            dist = np.sqrt(xdist * xdist + ydist * ydist)
            if dist < COA_radius:
                # check if already mated
                alreadyMated = False
                for m in range(matingArray.shape[0]):
                    if (lookingForMale and matingArray[m]['femaleid'] == stoatid
                            and matingArray[m]['maleid'] == stoatArray[i]['id']):
                        alreadyMated = True
                        break
                    elif (not lookingForMale and matingArray[m]['maleid'] == stoatid
                            and matingArray[m]['femaleid'] == stoatArray[i]['id']):
                        alreadyMated = True
                        break

                if not alreadyMated:

                    xCoords[nStoatsInTmp] = stoatArray[i]['x']
                    yCoords[nStoatsInTmp] = stoatArray[i]['y']
                    # daysSincePheromoneRelease is 0 for an actual stoat so we can ignore
                    Kappac[nStoatsInTmp] = ((1.0 / minK) * np.exp(-COA_decay_spatial * dist))

                    nStoatsInTmp += 1
                    if nStoatsInTmp > xCoords.shape[0]:
                        raise ValueError('Too many items for tmp array 1')

    if daysSincePheromoneRelease != -1:
        # now search the pheromones
        for i in range(pheromoneArray.shape[0]):
            xdist = x - pheromoneArray[i]['x']
            ydist = y - pheromoneArray[i]['y']
            dist = np.sqrt(xdist * xdist + ydist * ydist)
            if dist < COA_radius:

                # not recently interacted with this pheromone.
                alreadyInteracted = False
                for p in range(pheromoneInteractionArray.shape[0]):
                    if (pheromoneInteractionArray[p]['stoatid'] == stoatid and 
                            pheromoneInteractionArray[p]['pheromoneid'] == i):
                        # already interacted
                        alreadyInteracted = True
                        break

                if not alreadyInteracted:

                    xCoords[nStoatsInTmp] = pheromoneArray[i]['x']
                    yCoords[nStoatsInTmp] = pheromoneArray[i]['y']
                    Kappac[nStoatsInTmp] = ((1.0 / minK) * np.exp(-COA_decay_spatial * dist) * 
                            np.exp(-COA_decay_temporal * daysSincePheromoneRelease))

                    nStoatsInTmp += 1
                    if nStoatsInTmp > xCoords.shape[0]:
                        raise ValueError('Too many items for tmp array 2')

    if nStoatsInTmp > 0:
        # we found some
        Kappac = Kappac[:nStoatsInTmp]
        Psel = Kappac / Kappac.sum()

        COAi = np.random.multinomial(1, Psel)
        for i in range(nStoatsInTmp):
            if COAi[i] > 0:
                result = (xCoords[i], yCoords[i], Kappac[i])
                break
        
    return result

@njit
def doMaleMating(x, y, stoatid, stoatArray, nStoats, encounterDistance, matingArray,
                habituationDays, day, probPregnacy):
    mated = False
    for i in range(nStoats):
        # males still look for pregnant females
        if (not stoatArray[i]['deleted'] and not stoatArray[i]['male'] 
                    and stoatArray[i]['parentid'] == -1): # must be mature
            # check not already mated with this one
            alreadyMated = False
            for m in range(matingArray.shape[0]):
                if (matingArray[m]['maleid'] == stoatid and 
                        matingArray[m]['femaleid'] == stoatArray[i]['id']):
                    alreadyMated = True
                    break

            if alreadyMated:
                # try next female stoat
                continue

            xdist = x - stoatArray[i]['x']
            ydist = y - stoatArray[i]['y']
            dist = np.sqrt(xdist * xdist + ydist * ydist)
            if dist < encounterDistance:
                # females can mate multiple times - but only update the flag if not already pregnant
                if not stoatArray[i]['pregnant']:
                    if np.random.binomial(1, probPregnacy) == 1:
                        stoatArray[i]['pregnant'] = True
                        stoatArray[i]['pregnant_day'] = day

                        # if pregnant, make all her non-adult female offspring also pregnant with probabil
                        for n in range(nStoats):
                            if (stoatArray[n]['parentid'] == stoatArray[i]['id'] and 
                                    not stoatArray[n]['male'] and 
                                    np.random.binomial(1, probPregnacy) == 1):
                                stoatArray[n]['pregnant'] = True
                                stoatArray[n]['pregnant_day'] = day

                addedMating = False
                for m in range(matingArray.shape[0]):
                    if matingArray[m]['maleid'] == -1:
                        # unused slot
                        matingArray[m]['maleid'] = stoatid
                        matingArray[m]['femaleid'] = stoatArray[i]['id']
                        matingArray[m]['ndays'] = habituationDays 
                        addedMating = True
                        break
                if not addedMating:
                    raise ValueError('Unable to add mating')

                mated = True
                break

    return mated

@njit
def doPheromoneInteraction(x, y, stoatid, stoatArray, nStoats, pheromoneArray, encounterDistance,
            pheromoneInteractionArray, habituationDays):
    """
    Interact with a pheromone (if one within encounterDistance). Update
    the pheromoneInteractionArray.
    """
    # for both males and females
    for p in range(pheromoneArray.shape[0]):
        xdist = pheromoneArray[p]['x'] - x
        ydist = pheromoneArray[p]['y'] - y
        dist = np.sqrt(xdist * xdist + ydist * ydist)
        if dist < encounterDistance:
            # not recently interacted with this pheromone.
            for i in range(pheromoneInteractionArray.shape[0]):
                if (pheromoneInteractionArray[i]['stoatid'] == stoatid and 
                        pheromoneInteractionArray[i]['pheromoneid'] == p):
                    # already interacted
                    continue

            # interact, but won't be attracted to this pheromone for another habituationDays
            # find slot
            added = False
            for i in range(pheromoneInteractionArray.shape[0]):
                if pheromoneInteractionArray[i]['stoatid'] == -1:
                    pheromoneInteractionArray[i]['stoatid'] = stoatid
                    pheromoneInteractionArray[i]['pheromoneid']  = p
                    pheromoneInteractionArray[i]['ndays'] = habituationDays
                    added = True
                    break
            if not added:
                raise ValueError('Unable to add to pheromoneInteractionArray')

            break

@njit
def doBirth(stoatArray, nStoats, stoat, meanRecruits, Current_Id):
    """
    Have the stoat give birth. Tries to re-use 'deleted' slots in stoatArray first,
    otherwise adds onto the end.
    """

    x = stoatArray[stoat]['x']
    y = stoatArray[stoat]['y']
    nKits = np.random.poisson(meanRecruits)
#    print('birth', nKits)
    i = 0
    while nKits > 0 and i < nStoats:
        # go through and find deleted slots for these baby stoats
        if stoatArray[i]['deleted']:
            male = np.random.random() < 0.5
            stoatArray[i]['deleted'] = False
            stoatArray[i]['x'] = x
            stoatArray[i]['y'] = y
            stoatArray[i]['home_x'] = x
            stoatArray[i]['home_y'] = y
            stoatArray[i]['male'] = male
            stoatArray[i]['pregnant'] = False
            stoatArray[i]['pregnant_day'] = -1
            stoatArray[i]['id'] = Current_Id
            stoatArray[i]['homerange'] = False
            Current_Id += 1
            stoatArray[i]['parentid'] = stoatArray[stoat]['id']
            nKits -= 1
        i += 1

    # ok add the rest onto the end
    while nKits > 0:
        male = np.random.random() < 0.5
        stoatArray[nStoats]['deleted'] = False
        stoatArray[nStoats]['x'] = x
        stoatArray[nStoats]['y'] = y
        stoatArray[nStoats]['home_x'] = x
        stoatArray[nStoats]['home_y'] = y
        stoatArray[nStoats]['male'] = male
        stoatArray[nStoats]['pregnant'] = False
        stoatArray[nStoats]['pregnant_day'] = -1
        stoatArray[nStoats]['id'] = Current_Id
        stoatArray[nStoats]['homerange'] = False
        Current_Id += 1
        stoatArray[nStoats]['parentid'] = stoatArray[stoat]['id']
        nKits -= 1
        nStoats += 1
        if nStoats > INITIAL_STOAT_ARRAY_SIZE:
            print('Too many stoats')
            raise ValueError('Too many stoats')
    
    # now not pregnant
    stoatArray[stoat]['pregnant'] = False
    stoatArray[stoat]['pregnant_day'] = -1

    return nStoats, Current_Id
    
@njit
def checkWithinDistanceOfTraps(x, y, trapsArray, trapEncDist):
    """
    Returns True if the given X, y is withing trapEncDist of a trap
    """
    foundTrap = False
    for i in range(trapsArray.shape[0]):
        xdist = trapsArray[i, 0] - x
        ydist = trapsArray[i, 1] - y
        dist = np.sqrt(xdist * xdist + ydist * ydist)
        #print('in trap', x, y, trapsArray[i, 0], trapsArray[i, 1], dist)
        if dist < trapEncDist:
            foundTrap = True
            break
    return foundTrap
    
@njit
def checkStoatHasKitsInNest(stoatArray, nStoats, stoatid):
    """
    Returns True if given stoat still has juvenile offspring
    """
    hasKits = False
    for i in range(nStoats):
        if (not stoatArray[i]['deleted'] and 
                stoatArray[i]['parentid'] == stoatid):
            hasKits = True
            break
        
    return hasKits

@njit
def inArray(value, array):
    """
    Helper function
    """
    found = False
    ## CHECK IF WE DO TRAPPING; IF NONE, RETURN FALSE
    if array is not None:
        for i in range(array.shape[0]):
            if value == array[i]:
                found = True
                break
    return found

@njit
def runRealisation(nDays, hoursPerDay, stoatArray, nStoats, stepScale, stepShape, 
            alphaK, minK, pheromoneReleaseDays, COA_radius, COA_decay_spatial, COA_decay_temporal,
            estrousStartDays, estrousEndDays, pheromoneArray, pheromoneInteractionArray, habituationDays,
            encounterDistance, birthDays, meanRecruits, trappingDays, trapsArray,
            trapEncDist, trapProbRemoval, pDaySurv, dispersalDays,
            mask, tlx, tly, brx, bry, pixSize, matingArray, Current_Id, stoatDebugEstrous,
            stoatDebugDaysSincePheromone, stoatDebugFrame, stoatDebugTrapping, directionalVM,
            nDaysPregnantBeforeBirth, probPregnacy):
    """
    Main function - iterates through all the days, hours etc
    """
    debugIndex = 0

    daysSincePheromoneRelease = -1
    inEstrous = False
    Kappac = 0.0
    for day in range(nDays):

        ############################
        ##
        ## Break out with failure if a lot of stoats
        ##
        if nStoats > 175:
            eradication = False
            print('nStoats > 75 and failed eradication: ', nStoats)
            break
        ##
        ############################


        if inArray(day, pheromoneReleaseDays):
            daysSincePheromoneRelease = 0

        if inArray(day, estrousStartDays):
            inEstrous = True
        elif inArray(day, estrousEndDays):
            inEstrous = False

        if inArray(day, dispersalDays):
            # disperse males and juvenile stoats to shake thing up a bit
            for i in range(nStoats):
                if not stoatArray[i]['deleted'] and (stoatArray[i]['parentid'] != -1 or
                                stoatArray[i]['male']):
                    x, y = createRandomLocationOnIsland(mask, tlx, tly, brx, bry, pixSize)
                    stoatArray[i]['x'] = x
                    stoatArray[i]['y'] = y
                    stoatArray[i]['home_x'] = x
                    stoatArray[i]['home_y'] = y
                    stoatArray[i]['bearingT_1'] = np.random.uniform(-np.pi, np.pi)
                # if a child, set so now an adult
                stoatArray[i]['parentid'] = -1

        for hour in range(hoursPerDay):
            count = 0
            #print('sdsds', nStoats, day, hour, stoatArray)
            eradication = True
            for stoat in range(nStoats):
                #if stoat in stoatDict:
                #    stoatDict[stoat].append((stoatArray[stoat]['x'], stoatArray[stoat]['y']))
                #else:
                #    stoatDict[stoat] = [(stoatArray[stoat]['x'], stoatArray[stoat]['y'])]
                if not stoatArray[stoat]['deleted']:
                    eradication = False  # at least one individual exists
                    count += 1

                    COA_x = stoatArray[stoat]['home_x']
                    COA_y = stoatArray[stoat]['home_y']
                    x = stoatArray[stoat]['x']
                    y = stoatArray[stoat]['y']

                    homerangeBehaviour = True # the default
                    if inEstrous:
                        if stoatArray[stoat]['male']:
                            # searching behaviour in this period
                            homerangeBehaviour = False
                        elif not stoatArray[stoat]['pregnant']:
                            # for females, searching behaviour unless pregnant
                            # or has no kits
                            if not checkStoatHasKitsInNest(stoatArray, nStoats, 
                                        stoatArray[stoat]['id']):
                                homerangeBehaviour = False
                        
                    # store, for movie
                    stoatArray[stoat]['homerange'] = homerangeBehaviour
                            
                    # if male search for non pregnant females
                    mated = False
                    if stoatArray[stoat]['male']:
                        # check inEstrous and mature male
                        if inEstrous and stoatArray[stoat]['parentid'] == -1:
                            mated = doMaleMating(x, y, stoatArray[stoat]['id'], 
                                stoatArray, nStoats, encounterDistance, matingArray, 
                                habituationDays, day, probPregnacy)
                    elif (inArray(day, birthDays) and hour == 0 and stoatArray[stoat]['pregnant'] and
                            (day - stoatArray[stoat]['pregnant_day']) > nDaysPregnantBeforeBirth):
                        # note: only one hour on this day results in giving birth
                        # female will give birth
                        #print('adding stoats', meanRecruits)
                        nStoats, Current_Id = doBirth(stoatArray, nStoats, stoat, meanRecruits, Current_Id)

                    # interact with pheromones
                    # TODO: just mated females won't do pheromones?
                    if not mated:
                        doPheromoneInteraction(x, y,  stoatArray[stoat]['id'],
                            stoatArray, nStoats, 
                            pheromoneArray, encounterDistance, pheromoneInteractionArray, habituationDays)

                    # trapping and mortality
                    # trapping done at each hour
                    killed = False
                    if inArray(day, trappingDays):
                        if checkWithinDistanceOfTraps(x, y, trapsArray, trapEncDist):
###                            print('within trap distance')
                            killed = np.random.binomial(1, trapProbRemoval) == 1
                            if killed:
###                                print('killed by traps')
                                if stoatDebugTrapping is not None:
                                    stoatDebugTrapping[day] += 1

                    # mortality done on the last hour per day
                    if not killed and hour == hoursPerDay-1:
                        killed = np.random.binomial(1, 1.0 - pDaySurv) == 1
#                        if killed:
#                            print('Natural Mortality')

                    if killed:
#                        print('killed')
                        stoatArray[stoat]['deleted'] = True
                        # now kill immature all children
                        for i in range(nStoats):
                            if stoatArray[i]['parentid'] == stoatArray[stoat]['id']:
                                stoatArray[i]['deleted'] = True

                        # and remove from matingArray
                        for i in range(matingArray.shape[0]):
                            if (matingArray[i]['maleid'] == stoatArray[stoat]['id'] or 
                                    matingArray[i]['femaleid'] == stoatArray[stoat]['id']):
                                matingArray[i]['maleid'] = -1

                        # phermone array
                        for i in range(pheromoneInteractionArray.shape[0]):
                            if pheromoneInteractionArray[i]['stoatid'] == stoatArray[stoat]['id']:
                                pheromoneInteractionArray[i]['stoatid'] = -1

                        # don't bother with movement now
                        continue
                            

                    # Now do movement
                    # skip this bit if it is offspring that hasn't dispered yet
                    if stoatArray[stoat]['parentid'] != -1:
                        continue

                    mateResult = None # stays None if homerange



                    if not homerangeBehaviour:
                        # look for other COA
                        lookingForMale = not stoatArray[stoat]['male']
                        mateResult = checkForEstrousMatesAndDecoysInRadius(lookingForMale, x, y, 
                                    stoatArray, nStoats, stoatArray[stoat]['id'], COA_radius, 
                                    COA_decay_spatial, COA_decay_temporal, minK, 
                                    daysSincePheromoneRelease, pheromoneArray, pheromoneInteractionArray,
                                    matingArray)
#                        if mateResult is None:
#                            print('Was unable to find new COA')

                    if mateResult is not None:
                        COA_x, COA_y, Kappac = mateResult






                    newPosOK = False # keep looping until new location on land
                    while not newPosOK:

                        stepLength = stepScale * np.random.weibull(stepShape)

                        xdistToCOA = COA_x - x
                        ydistToCOA = COA_y - y
                        distToCOA = np.sqrt(xdistToCOA * xdistToCOA + ydistToCOA * ydistToCOA)
                        if distToCOA == 0:
                            bearing = 0.0
                        else:
                            if mateResult is not None and distToCOA < (stepScale * 1.0):
                                stepLength = (distToCOA * 0.95) * np.random.weibull(stepShape)
                        
                            # CALC THE ARCSIN
                            asinTH = np.arcsin(xdistToCOA / distToCOA)
                            # CONDITIONS WITH ZERO CHANGE IN X DIRECTION
                            if xdistToCOA == 0:
                                bearing = np.arccos(ydistToCOA)
                            # CONDITION WITH CHANGES IN X DIRECTION
                            if ydistToCOA > 0.0:
                                bearing = asinTH
                            if xdistToCOA > 0.0 and ydistToCOA <= 0.0:
                                bearing = np.pi - asinTH
                            if xdistToCOA < 0.0 and ydistToCOA <= 0.0:
                                bearing = -np.pi - asinTH

                        # draw the actual bearing that the stoat moves in from vonMises
                        # note +1 on distance to ensure we don't end up with negative from the log
                        #print(bearing, distToCOA, stepLength)
                        if mateResult is None:
                            # IF NOT HOMERANGE BEHAVIOUR - RANDOM WALK MOVEMENT SEARCH MATE
                            if not homerangeBehaviour:
                                bearing = stoatArray[stoat]['bearingT_1']
                                bearing = np.random.vonmises(bearing, directionalVM)
#                                print('no mates and not HomeRangeBehaviour, bearing=', bearing)

                            # IF HOMERANGE BEHAVIOUR - NOT ESTROUS PERIOD
                            else:
                                bearing = np.random.vonmises(bearing, 
                                    np.log(np.power(distToCOA + 1, alphaK)))
                        else:
                            # new COA calculated - Kappac should be set
                            bearing = np.random.vonmises(bearing, Kappac)

                        # UPDATE 'bearingT_1' for directed search for next step
                        stoatArray[stoat]['bearingT_1'] = bearing

                        # GET DELTA X AND Y, AND NEW X AND Y
                        xDistToMove = np.sin(bearing) * stepLength
                        yDistToMove = np.cos(bearing) * stepLength
                        newx = stoatArray[stoat]['x'] + xDistToMove
                        newy = stoatArray[stoat]['y'] + yDistToMove

                        # is this new location on the island?
                        # otherwise start from the original pos and try again
                        newPosOK = checkLocationIsOnIsland(mask, tlx, tly, 
                                brx, bry, pixSize, newx, newy)
                        if newPosOK:
                            stoatArray[stoat]['x'] = newx
                            stoatArray[stoat]['y'] = newy

#            print(count, nStoats, day, hour)
            if stoatDebugEstrous is not None:
                stoatDebugEstrous[debugIndex] = inEstrous
                stoatDebugDaysSincePheromone[debugIndex] = daysSincePheromoneRelease
                for i in range(stoatArray.shape[0]):
                    stoatDebugFrame[debugIndex, i] = stoatArray[i]

            debugIndex += 1
            if eradication:
                print('eradicated')
                return True


        if daysSincePheromoneRelease != -1:
            daysSincePheromoneRelease += 1

        # decrement pheromoneInteractionArray
        for i in range(pheromoneInteractionArray.shape[0]):
            if pheromoneInteractionArray[i]['stoatid'] != -1:
                pheromoneInteractionArray[i]['ndays'] -= 1
                if pheromoneInteractionArray[i]['ndays'] == 0:
                    pheromoneInteractionArray[i]['stoatid'] = -1

        # decrement matingArray
        for i in range(matingArray.shape[0]):
            if matingArray[i]['maleid'] != -1:
                matingArray[i]['ndays'] -= 1
                if matingArray[i]['ndays'] == 0:
                    matingArray[i]['maleid'] = -1
                

    return eradication

def dayMonthToDays(startDate, endDate, dayMonths, ndays=1):
    result = []
    for day, month in dayMonths:
        for year in range(startDate.year, endDate.year+1):
            date = datetime.date(year, month, day)
            if date < startDate or date > endDate:
                # before the startDate
                continue
            startday = (date - startDate).days
            for n in range(ndays):
                result.append(startday + n)
    return np.array(result)

def getRandomVariates(params):
    """
    Get random variates for this realisation of model
    """
    ## get number of stoats to add to initial population 2 + nAdd
    nAdd = int(np.random.uniform(params.meanNAdd[0], params.meanNAdd[1]))
    ## get dates of pheromone release
    twoPheroPerYear = (np.random.binomial(1, .5) == 1)

    twoPheroPerYear = True

    if twoPheroPerYear:
        pheromoneReleaseDayMonths = params.pheromoneReleaseDayMonths
    else:
        if len(params.pheromoneReleaseDayMonths) > 0:
            pheromoneReleaseDayMonths = [params.pheromoneReleaseDayMonths[0]]
        else:
            pheromoneReleaseDayMonths = []
#    ## decide whether to use 250m, 500m, 750m or 1000m grid
#    result = int(np.random.uniform(0, len(params.decoySpacing)))
#    spacing = params.decoySpacing[result]

    ## GET DECOY SPACING 
    spacing = np.random.uniform(params.decoySpacing[0], params.decoySpacing[1])
    ## ROUND TO NEAREST 10 METRES
    spacing = np.round(spacing, -1)

    ## get spatial decay parameter for pheromone attractions (k)
    COA_decay_spatial = np.random.uniform(params.COA_decay_spatial[0],
                            params.COA_decay_spatial[1])
    ## get temporal decay parameter for pheromone attractions (k)
    COA_decay_temporal = np.random.uniform(params.COA_decay_temporal[0],
                            params.COA_decay_temporal[1])
    ## get random number of habituation days
    habituationDays = np.random.randint(params.habituationDays[0],
                            params.habituationDays[1])

    ## get random daily survivorship
    pAnnSurv = (np.random.uniform(params.PAnnualSurv[0], params.PAnnualSurv[1])) 
    pDaySurv = np.power(pAnnSurv, 1.0 / 365.0)

    alphaK = np.random.uniform(params.alphaK[0], params.alphaK[1])

    ## Print for assessing maps
    print('nAdd', nAdd, 'decoy spacing=', spacing, 'alphaK=', alphaK, 'COA_decay_spatial=', COA_decay_spatial, 
        'COA_decay_temporal=', COA_decay_temporal, 'habituationdays=',
        habituationDays, 'Daily surv prob', pDaySurv)


    return(nAdd, pheromoneReleaseDayMonths, spacing, alphaK, COA_decay_spatial,
            COA_decay_temporal, habituationDays, pDaySurv)


def runModel(params, save=True, savePath='.'):
    """
    Main function

    """

    # Open the mask and read it
    ds = gdal.Open(params.extentMask)
    transform = ds.GetGeoTransform()
    tlx, tly = gdal.ApplyGeoTransform(transform, 0, 0)
    brx, bry = gdal.ApplyGeoTransform(transform, ds.RasterXSize, ds.RasterYSize)
    pixSize = transform[1]
    mask = ds.GetRasterBand(1).ReadAsArray()
    del ds
    
    # create an array to handle the stoats
    stoatArray = np.zeros((INITIAL_STOAT_ARRAY_SIZE,), dtype=STOAT_DTYPE)
    stoatArray['deleted'] = True # empty

    # for keeping a track of which stoats have mated with which 
    matingArray = np.empty(INITIAL_STOAT_ARRAY_SIZE * 2, dtype=MATING_DTYPE)
    matingArray['maleid'] = -1  # unused flag

    # Draw random variates of parameters for this realisation
    (nAdd, pheromoneReleaseDayMonths, spacing, alphaK, COA_decay_spatial, 
        COA_decay_temporal, habituationDays, pDaySurv) = getRandomVariates(params)
  

    # initial stoats
    nStoats, Current_Id = createInitialStoats(stoatArray, nAdd, mask, tlx, tly, 
                    brx, bry, pixSize, params.nDaysPregnantBeforeBirth)

    # read in the traps
    trapsArray = readTrapsFile(params.trapsFile)

    # convert datetime object to number of days for numba code
    nDays = (params.endDate - params.startDate).days

    pheromoneReleaseDays = dayMonthToDays(params.startDate, 
                params.endDate, pheromoneReleaseDayMonths)


#    print('top of runModel fx, pheromoneReleaseDays', pheromoneReleaseDays)


    estrousStartDays = dayMonthToDays(params.startDate, 
                params.endDate, [params.estrousStartDayMonth])
    estrousEndDays = dayMonthToDays(params.startDate, 
                params.endDate, [params.estrousEndDayMonth])

    birthDays = dayMonthToDays(params.startDate, 
                params.endDate, [params.birthDayMonth])

    pheromoneArray = makePheromoneArray(mask, tlx, tly, brx, bry, pixSize,
                                spacing, transform)
    # clobber it
    #pheromoneArray = np.empty(0, dtype=PHEROMONE_DTYPE)

    pheromoneInteractionArray = np.empty(INITIAL_STOAT_ARRAY_SIZE * pheromoneArray.shape[0], 
                dtype=PHEROMONE_INTERACTION_DTYPE)
    pheromoneInteractionArray['stoatid'] = -1

    trappingDays = None
    if params.trappingDayMonths is not None: 
        trappingDays = dayMonthToDays(params.startDate, params.endDate, 
            params.trappingDayMonths, params.nTrapDays)



    dispersalDays = dayMonthToDays(params.startDate, 
                params.endDate, [params.dispersalDateDayMonth])

#    print('Trapping days', trappingDays, 'Dispersaldays', dispersalDays,
#        'pheromone Days', pheromoneReleaseDays)



#    pDaySurv = np.power(params.PAnnualSurv, 1.0 / 365.0)

    nhours = nDays * params.hoursPerDay
    if save:
        stoatDebugInEstrous = np.zeros(nhours, dtype=np.bool)
        stoatDebugDaysSincePheromone = np.zeros(nhours, dtype=np.int32)
        stoatDebugFrame = np.zeros((nhours, INITIAL_STOAT_ARRAY_SIZE), dtype=STOAT_DTYPE)
        stoatDebugTrapping = np.zeros(nDays, dtype=np.int32)
    else:
        stoatDebugInEstrous = None
        stoatDebugDaysSincePheromone = None
        stoatDebugFrame = None
        stoatDebugTrapping = None

    eradicated = runRealisation(nDays, params.hoursPerDay, stoatArray, nStoats, 
            params.stepScale, params.stepShape, alphaK, params.minK,
            pheromoneReleaseDays, params.COA_radius, COA_decay_spatial, COA_decay_temporal,
            estrousStartDays, estrousEndDays, pheromoneArray, pheromoneInteractionArray, habituationDays,
            params.encounterDistance, birthDays, params.meanRecruits, trappingDays, trapsArray,
            params.trapEncDist, params.trapProbRemoval, pDaySurv, dispersalDays,
            mask, tlx, tly, brx, bry, pixSize, matingArray, Current_Id, stoatDebugInEstrous,
            stoatDebugDaysSincePheromone, stoatDebugFrame, stoatDebugTrapping, params.directionalVM,
            params.nDaysPregnantBeforeBirth, params.probPregnacy)

    if save:
        outname = os.path.join(savePath, 'stoats.npz')
        np.savez_compressed(outname, inEstrous=stoatDebugInEstrous,
                daysSincePheromone=stoatDebugDaysSincePheromone, debugInfo=stoatDebugFrame,
                trappingDays=trappingDays, pheromones=pheromoneArray,
                traps=trapsArray, trappingCount=stoatDebugTrapping)
                
        outname = os.path.join(savePath, 'stoatsparams.pkl')
        paramsFile = open(outname, 'wb')
        # cache some other things 
        params.COA_decay_spatial = COA_decay_spatial
        params.COA_decay_temporal = COA_decay_temporal
        pickle.dump(params, paramsFile)
        paramsFile.close()

    return (eradicated, nAdd, spacing, len(pheromoneReleaseDayMonths), alphaK, 
            COA_decay_spatial, COA_decay_temporal, habituationDays, pDaySurv)
