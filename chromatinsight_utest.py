import chromatinsight as ci
import os

myRegionFile = "./utest/input/regionFile_chrX_selection.bed"
myInputFolder = "./utest/input"
myOutputFolder = "./utest/output"
if not os.path.isdir(myOutputFolder): os.mkdir(myOutputFolder)

# RF_seed = 0 and label_seed = 0 are set to generate a deterministic random behaviour
resultObserved = ci.testPrediction(prefix = "mono",
						histmod = "ac",
						regionFile = myRegionFile,
						chrom = "chrX",
						inputFolder = myInputFolder,
						outputFolder = myOutputFolder,
						output = "mono_ac_chrX_real.txt",
						totRandomStates = 11,
						randomize = False,
						verbose = True,
						RF_seed = 0)
resultRnd = ci.testPrediction(prefix = "mono",
						histmod = "ac",
						regionFile = myRegionFile,
						chrom = "chrX",
						inputFolder = myInputFolder,
						outputFolder = myOutputFolder,
						output = "mono_ac_chrX_rnd.txt",
						totRandomStates = 11,
						randomize = True,
						verbose = True,
						label_seed = 0,
						RF_seed = 0)