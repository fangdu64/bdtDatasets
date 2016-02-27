import sys, traceback, array, csv, random,re

#load previously determined DNase bin idxs,
#only consider those high signal bins
inFile="E:/Biostatistics/EncodeCluster/results/DNase_HighSignalBins/highSignalBins.txt"
featureIdxs=[]
for line in open(inFile):
    fields=line.rstrip('\n').split('\t')
    featureIdxs.append(int(fields[0]))
sampleCnt=250000
finalDNaseFeatureIdxs=random.sample(featureIdxs,sampleCnt)
finalDNaseFeatureIdxs.sort()
#==========================
# Export DNase
#==========================
outdir="E:/iBS/trunk/analysis/iBS.Projects/BDVD/DukeUWDNase/100bp/comman/FeatureIdxs/SignalRandom"
outFile="{0}/DNase_FeatureIdxs.txt".format(outdir)
with open(outFile, 'w') as f:
    for idx in finalDNaseFeatureIdxs:
        rt=f.write(str(idx)+'\n')
