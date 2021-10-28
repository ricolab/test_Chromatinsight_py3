###########################
### Chromatinsight v2.4 ###
###########################
#
# a set of methods
# used to use a random forest algorithm
# to detect differential patterns between two sets of samples
# analysed with ChIP-seq histone modifications
# and pre-binarised by ChromHMM
#
# Author: Marco Trevisan-Herraz, PhD
# Computational Epigenomics Laboratory
# Newcastle University, UK
# 2019-2021
#
# requirements:
# Python 2.7
# pandas installed, see https://pandas.pydata.org/pandas-docs/stable/getting_started/install.html
# sklearn installed, see https://scikit-learn.org/stable/install.html
# Quick way to install both:
# pip install pandas
# pip install sklearn
#

import glob
import os
import pandas
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.ensemble import RandomForestClassifier
import pdb
import random

#######################################################

# load any _tab separated values_ file
def load2stringList(fileName, removeCommas = False, splitChar = "\t"):
	
	reader = open(fileName, "r")
	fullList = []
	
	for myRow in reader:
	
		myRowStrip = myRow.strip()
		
		if len(myRowStrip) > 0:
			thisRow = myRowStrip.split(splitChar)
			
			for i in xrange(len(thisRow)):
				thisRow[i] = thisRow[i].strip()
				
			if removeCommas:
				for i in xrange(len(thisRow)):
					if thisRow[i].endswith('"') and thisRow[i].startswith('"'):
						thisRow[i] = thisRow[i][1:len(thisRow[i]) - 1]
			
			fullList.append(thisRow)
	
	return fullList

#------------------------------------------------------

def stringList2inputDataFile(input, format = ['s', 'f', 'f'], fillEmptyPositions = False, emptyFiller = ""):
	
	result = []
	
	counter = 0
	for myRow in input:
		# if it is an inputRawDataFile, it should be str - float - float
		if fillEmptyPositions or len(myRow) >= len(format):
			resultRow = []
			
			if len(format) > 0:
				for i in xrange(len(format)):
					if i > len(myRow) - 1:
						resultRow.append(emptyFiller)
					else:
						stringy = myRow[i].strip()
						if format[i] == 's':
							resultRow.append(stringy)
						elif format[i] == 'f' or format[i] == 'i':
							if len(stringy) > 0 and \
								((stringy[0] >= '0' and stringy[0] <= '9') or \
								(stringy[0] == '-' and (stringy[1] >= '0' and stringy[1] <= '9'))):
								if format[i] == 'f': resultRow.append(float(stringy))
								if format[i] == 'i': resultRow.append(int(float(stringy)))
							else:
								# if a row that is supposed to contain a float or an int
								# is empty, or there is an error while reading it,
								# the row is not read but the program keeps going
								resultRow = []
								break
			
			if len(resultRow) > 0:
				result.append(resultRow)
	
	return result

#------------------------------------------------------

def removeHeader(myList):
	
	if len(myList) > 0:
		myList.remove(myList[0])
	
	return myList

#------------------------------------------------------

def saveFile(fileName, list, header = ""):
	
	writer = open(fileName, "w")
	if len(header) > 0:
		writer.write(header + "\n")
	
	for row in list: saveRow(writer, row)
	
	writer.close()
	
	return
	
#------------------------------------------------------

def saveRow(writer, rowList):
	
	line = ""
	for element in rowList:
		line += str(element) + "\t"
		
	line = line[:-1] # to remove last \t
	line += "\n"
	writer.write(line)
	
	return
	
#------------------------------------------------------

def mergeRegionFiles(regionFileFolder = "", minDistance = 1000, regionFileId = "", outputFile = ""):
	
	# the default output does not have a *.bed extension,
	# in order to avoid reprocessing the result if the program is run again
	if len(outputFile.strip()) == 0: outputFile = os.path.join(regionFileFolder, "mergeRegionFiles_output.bed.txt")
	
	chromList = [str(x) for x in range(23)[1:]] + ["X", "Y"]
	
	# load files
	regionFiles = glob.glob(os.path.join(regionFileFolder, "*.bed"))
	allFileData = []
	
	for regionFile in regionFiles:
		regionList = stringList2inputDataFile(removeHeader(load2stringList(regionFile)), format = ["s", "i", "i", "f"])
		allFileData.append(regionList)
		
	# make list for each chromosome
	allChromDictionary = {}
	for chrom in chromList: allChromDictionary[chrom] = []
	for file in allFileData:
		for fileRow in file:
			allChromDictionary[fileRow[0]].extend(fileRow[1:3])
	for chrom in chromList:
		allChromDictionary[chrom].sort()
		n = 0
		while not n > len(allChromDictionary[chrom]) - 2:
			if abs(allChromDictionary[chrom][n] - allChromDictionary[chrom][n + 1]) < minDistance:
				# the int is not needeed, but I use it to make it more readable
				allChromDictionary[chrom][n] = int((allChromDictionary[chrom][n] + allChromDictionary[chrom][n + 1]) / 2)
				del(allChromDictionary[chrom][n + 1])
			else: n += 1
			
	finalList = []
	for chrom in chromList:
		n = 0
		while not n > len(allChromDictionary[chrom]) - 2:
			# print chrom, n, len(allChromDictionary[chrom])
			finalList.append([chrom, allChromDictionary[chrom][n], allChromDictionary[chrom][n + 1]])
			n += 1
	
	saveFile(outputFile, finalList, "chr\tstart\tend")
	
	message = "\n\nProcessed %i files.\nGenerated %i subregions.\nNew file saved at %s" % (len(regionFiles), len(finalList), outputFile)
	print message
	
	return

#------------------------------------------------------

def joinData(pattern = "", histmod = "ac", direc = "", verbose = False, removePAR = False, removeXist = False):

	PAR1chrXend = 13497 # this included, is in the PAR#1, one-based
	PAR2chrXstart = 774653
	XistStart = 365300 # 73,060,000 b
	XistEnd = 365400 # 73,080,000 b
	
	if histmod == "ac" or histmod == "H3K27ac":
		histmodPos = 0
		histmod = "H3K27ac"
	if histmod == "me1" or histmod == "H3K4me1":
		histmodPos = 1
		histmod = "H3K4me1"

	myDirec = direc
	myFilePath = os.path.join(myDirec, "*" + pattern + "*binary.txt")
	myFiles = glob.glob(myFilePath)
	myList = []
	if verbose: print "Loading files..."
	for myFile in myFiles:
		badFile = False
		if verbose: print myFile
		reader = open(myFile, "r")
		myFileContents = [os.path.basename(myFile)]
		skipThis = True
		counter = 0
		while True:
			counter += 1
			myLine = reader.readline().strip()
			if not myLine: break
			if not skipThis:
				myLine = myLine.split("\t")
				if counter == 2:
					if len(myLine) >= histmodPos+1:
						if not myLine[histmodPos] == histmod:
							badFile = True
							break
					else:
						badFile = True
						break
				if (counter > 2) \
					and ((not removePAR) or (counter > PAR1chrXend + 2 and counter < PAR2chrXstart + 2)) \
					and ((not removeXist) or (counter < XistStart + 2) or (counter > XistEnd + 2)):
					myValue = int(myLine[histmodPos])
					if myValue == 2:
						badFile = True
						break
					myFileContents.append(myValue)
			else: skipThis = False
		reader.close()
		if not badFile: myList.append(myFileContents)
		else:
			if verbose: print "Warning: %s is a bad file, skipping." % myFile
	
	if verbose: print "Generating main data frame..."
	myListP = pandas.DataFrame(myList)
	if verbose: print "Main data frame generated..."
	myListPi = myListP.set_index(0)
	
	return myListPi
	
#------------------------------------------------------

def testPrediction(prefix = "",
					testSize = 0.3,
					totRandomStates = 11,
					chrom = "",
					histmod = "ac",
					verbose = False,
					regionFile = "",
					interRegionTested = True,
					binSize = 200,
					inputFolder = "",
					outputFolder = "",
					output = "output.txt",
					randomize = False,
					randomizeMethod = "scramble",
					label_seed = None,
					RF_seed = None):

# randomizeMethod can be
# coin -> 50% chance of getting "fem" or "mal"
# scramble -> just scramble the existing labels (preserving their ratios)

	outputFile = os.path.join(outputFolder, output)
	
	if len(chrom) == 0:
		chromList = ["chr" + str(chrom) for chrom in range(23)[1:]]
		chromList.append("chrX")
		chromList.append("chrXnoPAR")
		chromList.append("chrXnoPARnoXist")
	else: chromList = [chrom]
	medianPos = totRandomStates // 2 # only for odd values
	
	regionList = []
	if(len(regionFile) > 0):
		# regionList = stringList2inputDataFile(removeHeader(load2stringList(regionFile)), format = ["s", "i", "i", "f"])
		regionList = stringList2inputDataFile(removeHeader(load2stringList(regionFile)), format = ["s", "i", "i"])
	else:
		# regionList = [["0", 0, 0, 0.0]]
		regionList = [["0", 0, 0]]
	
	myScoreList = []
	for chrom in chromList:
		thisChrom = chrom
		removePAR = False
		removeXist = False
		
		if chrom == "chrXnoPAR":
			thisChrom = "chrX"
			removePAR = True
		
		if chrom == "chrXnoPARnoXist":
			thisChrom = "chrX"
			removePAR = True
			removeXist = True
			
		myData = joinData(pattern = prefix + "*_" + thisChrom + "_", histmod = histmod, verbose = verbose, removePAR = removePAR, removeXist = removeXist, direc = inputFolder)
		if verbose: print "Data joined."
		
		sampleLabelList = [myData.index.str.split("_")[i][2] for i in range(len(myData))]
		if randomize:
			if verbose: print "Randomising labels, as requested..."
			random.seed(label_seed)
			if randomizeMethod == "coin": sampleLabelList = [["fem", "mal"][random.randint(0,1)] for i in range(len(myData))]
			if randomizeMethod == "scramble": random.shuffle(sampleLabelList)
						
		myScoreChrom = []
		previousRegionEnd = 0
		interTADLabel = "Starting"
		for chromRegion in regionList:
			thisChrom = "chr" + chromRegion[0]
			if thisChrom == chrom:
				
				# we check the interTAD region
				regionID = "%s_%i-%i_%s" % (chrom, previousRegionEnd, chromRegion[1], interTADLabel)
				regionStart = previousRegionEnd // binSize
				regionEnd = chromRegion[1] // binSize
				
				if regionStart < regionEnd - 1 and interRegionTested:
					if verbose: print "Getting patterns in inter-region %s" % regionID
					if chromRegion[2] == 0: regionEnd = len(myData.iloc[0,:]) # the last bin
					regionCoordinates = "%s:%s-%s" % (chrom, format(regionStart * binSize, ","), format(regionEnd * binSize, ","))
					
					thisData = myData.iloc[:,regionStart:regionEnd - 1]
					thisData.loc[:, "sex"] = sampleLabelList
					
					
					myScores = []
					for randomState in range(totRandomStates):
						myScore = getScore(thisData, testSize, randomState, RF_seed)
						myScores.append(myScore)
					
					myScoreChrom.append([regionID] + myScores)
					print ("%s %s: %s, median = %f" % (interTADLabel, regionCoordinates, myScores, sorted(myScores)[medianPos]))
					interTADLabel = "interTAD"
				
				# we check the TAD region
				regionID = "%s_%i-%i_TAD" % (chrom, chromRegion[1], chromRegion[2])
				regionStart = chromRegion[1] // binSize
				regionEnd = chromRegion[2] // binSize
				
				if regionStart < regionEnd - 1:
					if verbose: print "Getting patterns in region %s" % regionID
					if chromRegion[2] == 0: regionEnd = len(myData.iloc[0,:]) # the last bin
					regionCoordinates = "%s:%s-%s" % (chrom, format(regionStart * binSize, ","), format(regionEnd * binSize, ","))
					
					thisData = myData.iloc[:,regionStart:regionEnd - 1]
					thisData.loc[:, "sex"] = sampleLabelList
					
					myScores = []
					for randomState in range(totRandomStates):
						myScore = getScore(thisData, testSize, randomState, RF_seed)
						myScores.append(myScore)
					
					myScoreChrom.append([regionID] + myScores)
					
					print ("TAD %s: %s, median = %f" % (regionCoordinates, myScores, sorted(myScores)[medianPos]))
				previousRegionEnd = chromRegion[2]
				
		# check the last part of the chromosome
		regionStart = previousRegionEnd // binSize
		regionEnd = len(myData.iloc[0,:]) # the last bin
		if regionStart < regionEnd - 1 and interRegionTested:
			regionID = "%s_%i-%i_Ending" % (chrom, previousRegionEnd, regionEnd * binSize)
			regionCoordinates = "%s:%s-%s" % (chrom, format(regionStart * binSize, ","), format(regionEnd * binSize, ","))
			
			thisData = myData.iloc[:,regionStart:regionEnd - 1]
			thisData.loc[:, "sex"] = sampleLabelList
			
			myScores = []
			for randomState in range(totRandomStates):
				myScore = getScore(thisData, testSize, randomState, RF_seed)
				myScores.append(myScore)
			
			myScoreChrom.append([regionID] + myScores)
			print ("Ending %s: %s, median = %f" % (regionCoordinates, myScores, sorted(myScores)[medianPos]))
		
		myScoreList.append(myScoreChrom)
		
	if verbose: print "Saving results in file %s" % outputFile
	saveFile(outputFile, myScoreList[0], header = "chrom_init-end_region")
		
	return myScoreList

#------------------------------------------------------

def getScore(myData, testSize, randomState = None, RF_seed = None):
	train_index, test_index = next(StratifiedShuffleSplit(test_size = testSize, random_state=randomState).split(myData, myData.sex))
	# myData.sex.value_counts().plot.bar(x="sex")
	myData_train = myData.iloc[train_index,:]
	myData_test = myData.iloc[test_index,:]
	myData_train
	
	rf = RandomForestClassifier(n_estimators = 200, random_state = RF_seed)
	myDataLength = myData.shape[1] - 1
	myFit = rf.fit(myData_train.iloc[:,0:myDataLength], myData_train.loc[:,"sex"])
	myScore = rf.score(myData_test.iloc[:,0:myDataLength], myData_test.loc[:,"sex"])
	
	return myScore

#######################################################