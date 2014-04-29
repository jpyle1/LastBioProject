import optparse as optParse
import random,math,array
import matplotlib.pyplot as plt

def main(jOne,jTwo,offset,rOne,rTwo):
	"""
		Four runs should be performed, and the results should be averaged.
	"""
	#Simple caching mechanism...
	radCache = dict()
	distRadCache = dict()
	results  = []
	for runNum in xrange(4):
		results.append(computeSteps(jOne,jTwo,offset,rOne,rTwo,radCache,distRadCache,runNum))
	#Average the results...
	avgWl = sum(map(lambda x: x['wL'],results))/4.0
	avgEn = sum(map(lambda x: x['en'],results))/4.0	
	avgMI = map(lambda x: sum(x)/4.0,zip(*map(lambda x: x['mI'],results)))
	avgJE = map(lambda x: sum(x)/4.0,zip(*map(lambda x: x['jE'],results)))
	avgCo = map(lambda x: sum(x)/4.0,zip(*map(lambda x: x['cor'],results)))	
	fig1 = plt.figure(1)
	plt.subplot(211)
	plt.imshow(results[0]['disp'])
	plt.subplot(212)
	plt.ylim(0,2)
	plt.plot(xrange(14),avgMI)
	plt.plot(xrange(14),avgJE)
	plt.plot(xrange(14),avgCo)
	plt.figtext(0,0,"CorrWaveLength="+str(avgWl)+","+"En="+str(avgEn))
	plt.show()

def getRadius(index,rOne,rTwo):
	"""
		Returns the indices of the values for the summation.
	"""
	lessThan = []
	between = []
	other = []
	for i in xrange(30):
		for j in xrange(30):
			xDist = math.fabs(index[0] - i)
			xDist = xDist if xDist<=15 else 30-xDist
			yDist = math.fabs(index[1]-j)
			yDist = yDist if yDist<=15 else 30-yDist
			dist = xDist+yDist
			if dist < rOne:
				lessThan.append((i,j,dist))
			elif dist >= rOne and dist< rTwo:
				between.append((i,j,dist))
			else:
				other.append((i,j,dist))	
	return (lessThan,between,other)

def addToCache(radCache,distRadCache,index,item):
	"""
	Adds an item to the cache.
	"""		
	radCache[index] = item	
	#Add items to the distance cache...
	for val in item[0]:
		if val[2] in distRadCache:
			distRadCache[val[2]].append((index,(val[0],val[1])))
		else:
			distRadCache[val[2]] = []
			distRadCache[val[2]].append((index,(val[0],val[1])))
	for val in item[1]:
		if val[2] in distRadCache:
			distRadCache[val[2]].append((index,(val[0],val[1])))
		else:
			distRadCache[val[2]] = []
			distRadCache[val[2]].append((index,(val[0],val[1])))
	for val in item[2]:
		if val[2] in distRadCache:
			distRadCache[val[2]].append((index,(val[0],val[1])))
		else: 
			distRadCache[val[2]] = []
			distRadCache[val[2]].append((index,(val[0],val[1])))

def checkDistRadCache(distRadCache,length):
	"""
	Checks to see if an item is in the dist rad cache.
	"""	
	return length in distRadCache

def getItemFromDistRadCache(distRadCache,length):
	"""
		Returns an item from the cache.
	"""
	if(checkDistRadCache(distRadCache,length)):
		return distRadCache[length]
	else:
		return None

def checkCache(radCache,index):
	"""
	Checks to see if an item is in the cache.
	"""	
	return index in radCache

def getItemFromCache(radCache,index):
	"""
	Returns an item from the cache.
	"""
	if(checkCache(radCache,index)):
		return radCache[index] 
	return None

def computeSteps(jOne,jTwo,offset,rOne,rTwo,radCache,distRadCache,runNumber):
	#Initialize the cells to a random value.
	cells=[[float(random.randint(0,1)) for i in xrange(30)] for j in xrange(30)]
	#Attempt 5 runs at first.
	for iterNum in xrange(30):
		#Indices of the unupdated cells.
		unUpdated = range(900)
		#Now, update all the sells at a different point at time, one at at time.	
		for timeStep in xrange(900):
			#Get a random cell...
			toBeUpdated = unUpdated[random.randint(0,900-timeStep-1)]
			unUpdated.remove(toBeUpdated)		
			cellRowVal = toBeUpdated%30	
			cellColVal = toBeUpdated/30
			cell = (cellRowVal,cellColVal)
			radius = None
			if checkCache(radCache,cell):
				radius = getItemFromCache(radCache,cell)
			else:	
				radius = getRadius(cell,rOne,rTwo)
				addToCache(radCache,distRadCache,cell,radius)
			lessThanStates = map(lambda x: cells[x[0]][x[1]],radius[0])
			betweenStates =  map(lambda x: cells[x[0]][x[1]],radius[1])
			val = offset + jOne*sum(lessThanStates) + jTwo*sum(betweenStates)
			cells[cellRowVal][cellColVal] = 1 if val>=1 else -1 

	spatCorr = []
	#Compute the first spatial correlation value.
	spatCorr.append(math.fabs(1-math.pow((1/900.0)*sum([sum(cells[i]) for i in xrange(30)]),2)))
	#Now, firstly, compute the spatial correlation.
	for radiusLength in xrange(1,15):
		mult = (1.0/(3600*radiusLength))
		act = mult*sum(map(lambda x:cells[x[0][0]][x[0][1]]*cells[x[1][0]][x[1][1]]
			,getItemFromDistRadCache(distRadCache,radiusLength)))
		multInh = 1.0/900
		inh = multInh*sum([sum(cells[i]) for i in xrange(30)])
		inh = math.pow(inh,2)
		spatCorr.append(math.fabs(act-inh))
	firstWave = spatCorr[0]/math.exp(1)	
	closestIdx = 0
	closestVal = 1.5
	#Now, compute the wavelength.
	for radiusLength in xrange(1,15):
		val = math.fabs(firstWave-spatCorr[radiusLength])
		if(val < closestVal):
			closestIdx = radiusLength
			closestVal = val
	charWaveLength = closestIdx
	#Convert the cells.
	converted = [[ (cells[i][j]+1)/2.0 for j in xrange(30)] for i in xrange(30)]
	convertedNeg = [[ (cells[i][j]-1)/2.0 for j in xrange(30)] for i in xrange(30)]
	prPlusOne = 1/900.0*sum([sum(converted[i]) for i in xrange(30)])
	prMinusOne = 1 - prPlusOne
	val1 = prPlusOne*math.log(prPlusOne,2) if prPlusOne !=0 else 0
	val2 = prMinusOne*math.log(prMinusOne,2) if prMinusOne !=0 else 0
	#Calc the entropy.
	entropy = -1*(val1+val2)
	#Calc the joint entropy.
	posSums = [ sum(map(lambda x:converted[x[0][0]][x[0][1]]*converted[x[1][0]][x[1][1]], 
		getItemFromDistRadCache(distRadCache,radiusLength) )) for radiusLength in xrange(0,15)]
	posSums = [1 if posSums[radiusLength] == 0.0 else posSums[radiusLength] for radiusLength in
		xrange(0,15)] 
	negSums = [ sum(map(lambda x:convertedNeg[x[0][0]][x[0][1]]*convertedNeg[x[1][0]][x[1][1]],
		getItemFromDistRadCache(distRadCache,radiusLength) )) for radiusLength in xrange(0,15)]
	negSums = [1 if negSums[rl] == 0.0 else negSums[rl] for rl in xrange(0,15)]
	pPlusPlus = [ 1/(900.0*4*radiusLength)*posSums[radiusLength] for radiusLength in xrange(1,15)]
	pMinusMinus = [1/(900.0*4*radiusLength)*negSums[radiusLength] for radiusLength in xrange(1,15)]
	pPlusMinus = [1-pPlusPlus[rl]-pMinusMinus[rl] for rl in xrange(0,14)]
	pPlusMinus = [1 if pPlusMinus[rl] <= 0.0 else pPlusMinus[rl] for rl in xrange(0,14)]
	jointEntropy = [-1*(pPlusPlus[radiusLength]*math.log(pPlusPlus[radiusLength],2)+
		pMinusMinus[radiusLength]*math.log(pMinusMinus[radiusLength],2)+
		pPlusMinus[radiusLength]*math.log(pPlusMinus[radiusLength],2)) for radiusLength in xrange(0,14)]
	mutualInformation = [2*entropy-jointEntropy[radiusLength] for radiusLength in xrange(0,14)]
	spatCorr.pop(0)
	cells = [[ 1 if cells[i][j] == 1 else 0 for j in xrange(30)] for i in xrange(30)]
	return {'mI': mutualInformation, 'jE': jointEntropy, 'en': entropy, 'wL': charWaveLength,
		'cor': spatCorr, 'disp': cells} 

if __name__ == "__main__":
	#Parse the Options from the user.
	parser = optParse.OptionParser()
	parser.add_option("-j","--jOne")
	parser.add_option("-k","--jTwo")
	parser.add_option("-o","--offset")
	parser.add_option("-r","--rOne")
	parser.add_option("-s","--rTwo")
	(options,args) = parser.parse_args()
	#Default vals...
	jOne = 1
	jTwo = -.1
	offset = 0
	rOne = 1
	rTwo = 2
	if(options.jOne is not None):
		jOne = float(options.jOne)
	if(options.jTwo is not None):
		jTwo = float(options.jTwo)
	if(options.rOne is not None):
		rOne = float(options.rOne)
	if(options.rTwo is not None):
		rTwo = float(options.rTwo)
	if(options.offset is not None):
		offset = float(options.offset)			
	main(jOne,jTwo,offset,rOne,rTwo)
