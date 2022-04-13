import chromatinsight as ci
import os

myRegionFile = "./utest/input/regionFile_chrX_selection.bed"
myGroupingFile = "./utest/input/grouping.txt"
myOutputFolder = "./utest/output"
if not os.path.isdir(myOutputFolder): os.mkdir(myOutputFolder)

# RF_seed = 0 and label_seed = 0 are set to generate a deterministic random behaviour
resultObserved = ci.testPrediction(groupingFile = myGroupingFile,
						regionFile = myRegionFile,
						histmod = "ac",
						chrom = "chrX",
						outputFolder = myOutputFolder,
						output = "mono_ac_chrX_real.txt",
						totRandomStates = 11,
						randomize = False,
						verbose = True,
						RF_seed = 0)
resultRnd = ci.testPrediction(groupingFile = myGroupingFile,
						regionFile = myRegionFile,
						histmod = "ac",
						chrom = "chrX",
						outputFolder = myOutputFolder,
						output = "mono_ac_chrX_rnd.txt",
						totRandomStates = 11,
						randomize = True,
						verbose = True,
						label_seed = 0,
						RF_seed = 0)