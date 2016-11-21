from __future__ import division
import matplotlib.pyplot as plt
import math
from numpy import random

def scalarMultiple(lst, scalar):
    newList = list(lst)
    length = len(newList)
    for i in range(0, length):
        newList[i] = scalar * newList[i]
    return newList

def addLists(lst1, lst2):
    newList = list(lst1)
    length = len(lst2)
    for i in range(0, length):
        newList[i] += lst2[i]
    return newList

def dotMap(functions, arguments):
    vals = []
    length = len(functions)
    for i in range(0, length):
        vals.append(functions[i](arguments))
    return vals

def verticalAppend(listOfLists, listOfValues):
    N = len(listOfLists)
    M = len(listOfValues)
    if N > M:
        length = M
    else:
        length = N
    for i in range(0, length):
        listOfLists[i].append(listOfValues[i])

def horizontalExtract(listOfLists, index):
    length = len(listOfLists)
    cut = []
    for i in range(0, length):
        cut.append(listOfLists[i][index])
    return cut

def normalize(waveValues):
    length = len(waveValues)
    sum = 0
    for i in range(0, length):
        sum += abs(waveValues[i]) ** 2
    norm = 1/sum**0.5
    for i in range(0, length):
        waveValues[i] *= norm

def RungeKuttaStep(derivatives, baseArguments, additionalArguments, stepSize=0.01, methodSettings={}):
    currentX = baseArguments[0]
    currentYs = baseArguments[1:]
    k1s = dotMap(derivatives, baseArguments + additionalArguments)
    k1s = scalarMultiple(k1s, stepSize)
    k2s = dotMap(derivatives, [currentX + stepSize/2] + addLists(currentYs, scalarMultiple(k1s, 0.5)) + additionalArguments)
    k2s = scalarMultiple(k2s, stepSize)
    currentYs = addLists(currentYs, k2s)
    currentX += stepSize
    return [currentX]+currentYs

def getSpaceDerivatives(spacePoints, spaceValues, N, periodicDerivatives=True):
        
    currentSpaceDerivatives = [spaceValues]
    
    for derivativeOrder in range(1, N+1):
        
        currentSpaceDerivatives.append([])
        
        for i in range(0, len(spacePoints)):
            
            pointsAvailable = len(currentSpaceDerivatives[derivativeOrder-1])
            spacePointsAvailable = len(spacePoints)
            
            if i + 1 >= pointsAvailable or i - 1 < 0:
                
                if periodicDerivatives:
                    
                    YDiff = ( currentSpaceDerivatives[derivativeOrder-1][(i+1) % pointsAvailable] - \
                        currentSpaceDerivatives[derivativeOrder - 1][i - 1 % pointsAvailable])
                    
                    if i+2 >= spacePointsAvailable:
                        
                        derVal = YDiff / (spacePoints[-1] - spacePoints[-3])
                        
                    else:
                        
                        derVal = YDiff / (spacePoints[i+1] - spacePoints[i - 1])
                        
                else:
                    
                    derVal = currentSpaceDerivatives[derivativeOrder][-1]
                    
            else:
                
                derVal = ( currentSpaceDerivatives[derivativeOrder-1][i+1] - \
                currentSpaceDerivatives[derivativeOrder-1][i-1] ) / \
                (spacePoints[i+1] - spacePoints[i-1])
            
            currentSpaceDerivatives[derivativeOrder].append(derVal)
    
    return currentSpaceDerivatives

def timeSpaceDiffEq(timeDerivative, spaceDerivativesN, initialSpaceValues, timeLimits, timeStep=0.01, additionalArguments=[], spaceLimits=False, periodicDerivatives=True):
    '''Initial space values is a two list list, i. e. [[x1, ...], [y1, ...]]'''
    
    if not spaceLimits:
        spaceLimits = [0, len(initialSpaceValues[0])]
    elif spaceLimits[1] < 0:
        spaceLimits[1] = len(initialSpaceValues[0]) + spaceLimits[1]
    
    spacePoints = initialSpaceValues[0]
    currentSpaceDerivatives = getSpaceDerivatives(spacePoints, initialSpaceValues[1], spaceDerivativesN, periodicDerivatives)
    
    currentTime = timeLimits[0]
    while currentTime < timeLimits[1]:
        if spaceLimits[0] > 0:
            newValues = currentSpaceDerivatives[0][0:spaceLimits[0]]
        else:
            newValues = []
        
        for i in range(spaceLimits[0], spaceLimits[1]):
            baseArguments = [currentTime, currentSpaceDerivatives[0][i]]
            adArgs = [spacePoints[i]] + horizontalExtract(currentSpaceDerivatives[1:], i) + additionalArguments
            newValuePair = RungeKuttaStep([timeDerivative], baseArguments, adArgs, timeStep)
            
            newValues.append(newValuePair[1])
        currentTime += timeStep
        if spaceLimits[1] < len(currentSpaceDerivatives[0]):
            newValues = newValues + currentSpaceDerivatives[0][spaceLimits[1]:]
        # normalization
        normalize(newValues)
        currentSpaceDerivatives = getSpaceDerivatives(spacePoints, newValues, spaceDerivativesN, periodicDerivatives)
        print currentTime
    return [spacePoints] + currentSpaceDerivatives, currentTime

# hbar/2m
C1 = 0.6*10**(-2)
# V0/hbar
C2 = 30

def voltage(x):
    return voltageDict[x]

def timeDer(args):
    # t = args[0]
    # y = args[1]
    # x = args[2]
    # derYbyX = args[3]
    # 2derYbyX2 = args[4]
    return 1j*C1*args[4] - 1j*C2*voltage(args[2])*args[1]

X = []
initY = []
L = 10**-10
for i in range(0, 501):
    X.append(L * i/500)
    #~ initY.append(math.exp(-((i - 243) ** 2) / (2*100)) + 0j)
    initY.append(0+0j)

initY[241] = 0+0j
initY[242] = 0+0j
initY[243] = 1+0j
initY[244] = 0+0j
initY[245] = 0+0j
normalize(initY)
displayInitY = []
for i in range(0, len(initY)):
    displayInitY.append(abs(initY[i]) ** 2)
# generate random voltages
voltageDict = {}
voltageList = []
maxVoltage = 0.025
random.seed(231542)
for i in range(0, len(X)):
    #~ voltageDict[X[i]] = -random.random()*maxVoltage
    #~ voltageDict[X[i]] = (math.sin(2*math.pi*20*X[i]/L)-1)*maxVoltage
    voltageDict[X[i]] = -0.025
    voltageList.append(voltageDict[X[i]])

result, endTime = timeSpaceDiffEq(timeDer, 2, [X, initY], [0, 100], spaceLimits=[1, -1], periodicDerivatives=True, timeStep=0.1)

for i in range(0, len(result[1])):
    result[1][i] = abs(result[1][i])**2

g1, = plt.plot(X, displayInitY, label="Original")
g2, = plt.plot(X, voltageList, label="Potential energy")
g3, = plt.plot(result[0], result[1], label="Final")

plt.legend([g1, g3, g2], ["Original", "Final", "Potential energy"])
plt.xlabel("Distance [Angstrom]")
plt.ylabel("Probability density")
plt.ylim([-0.05, 0.05])

plt.savefig("constant_100s_30C2_0.025V.png", transparent=True)