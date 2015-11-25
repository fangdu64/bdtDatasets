import sys, traceback, array, csv, random,re
import iBSFCDClient as fcdc
import iBS
import ExonArray.Module_TTS_TCIDs as ExonTTS


def sample_wr(population, k): 
    n = len(population)
    sps=[]
    for i in range(k):
        j=int(n*random.random())
        sps.append(population[j])
    return sps

(TTSBinIdxs,TCIDIdxs)=ExonTTS.GetTTSwithUniqueTCIDs100bp()

radiusBinCnt=20 #+-2k bp

DNaseFeatureIdxFrom=0
DNaseFeatureIdxTo=30363049 #whole genome

DNaseFeatureCnt=DNaseFeatureIdxTo - DNaseFeatureIdxFrom
DNaseFeatureIdx2NearestTCIDIdx=[-1]*DNaseFeatureCnt
DNaseFeatureIdx2NearestDist=[radiusBinCnt+100]*DNaseFeatureCnt

for i in range(len(TTSBinIdxs)):
    ttsBinIdx=TTSBinIdxs[i]
    tcidIdx=TCIDIdxs[i]

    # valid DNase bin idxs that are near to this TSS
    minBinIdx=ttsBinIdx-radiusBinCnt
    maxBinIdx=ttsBinIdx+radiusBinCnt
    if maxBinIdx < DNaseFeatureIdxFrom or minBinIdx > DNaseFeatureIdxTo:
        continue
    # convert to relative bin idx ( with regard to DNaseFeatureIdxFrom)
    idxFrom=max(minBinIdx,DNaseFeatureIdxFrom) - DNaseFeatureIdxFrom
    idxTo=min(maxBinIdx,DNaseFeatureIdxTo) - DNaseFeatureIdxFrom
    midIdx=ttsBinIdx-DNaseFeatureIdxFrom
    # a DNase bin can be near to multiple TSS, determine the nearest TTS a DNase is to
    for j in range(idxFrom,idxTo):
        d=abs(j-midIdx)
        if DNaseFeatureIdx2NearestDist[j]>d:
            DNaseFeatureIdx2NearestDist[j]=d
            DNaseFeatureIdx2NearestTCIDIdx[j]=tcidIdx

#load previously determined DNase bin idxs,
#only consider those high signal bins
inFile="E:/Biostatistics/EncodeCluster/results/DNase_HighSignalBins/highSignalBins.txt"
featureIdxs=[]
for line in open(inFile):
    fields=line.rstrip('\n').split('\t')
    featureIdxs.append(int(fields[0]))

validDNaseFeatureIdxs=[]
validDNaseFeatureIdxs_inHighSignal=[]
validNearestTCIDIdxs=[]
idx=-1
for featureIdx in featureIdxs:
    idx=idx+1
    i=featureIdx-DNaseFeatureIdxFrom # convert ot local idx
    tcidIdx=DNaseFeatureIdx2NearestTCIDIdx[i]
    if tcidIdx!=-1:
        validDNaseFeatureIdxs.append(featureIdx)
        validNearestTCIDIdxs.append(tcidIdx)
        validDNaseFeatureIdxs_inHighSignal.append(idx)

rowIdxs=list(range(len(validDNaseFeatureIdxs))) # len =254371
finalDNaseFeatureIdxs=[]
finalNearestTCIDIdxs=[]
finalHighDNaseFeatureIdxs=[]
for i in rowIdxs:
    finalDNaseFeatureIdxs.append(validDNaseFeatureIdxs[i])
    finalNearestTCIDIdxs.append(validNearestTCIDIdxs[i])
    finalHighDNaseFeatureIdxs.append(validDNaseFeatureIdxs_inHighSignal[i])

finalNearestTCIDIdxs=sample_wr(sorted(set(TCIDIdxs)),len(finalNearestTCIDIdxs))
#==========================
# Export Transcript Cluster Idxs
#==========================
outdir="E:/iBS/trunk/analysis/iBS.Projects/BDVD/DNaseExonCorrelation/100bp/s02-Random-PairIdxs"
outFile="{0}/Exon_FeatureIdxs.txt".format(outdir)
with open(outFile, 'w') as f:
    for idx in finalNearestTCIDIdxs:
        rt=f.write(str(idx)+'\n')

uniqueIdxs=sorted(set(finalNearestTCIDIdxs))
outFile="{0}/Exon_UniqueFeatureIdxs.txt".format(outdir)
with open(outFile, 'w') as f:
    for idx in uniqueIdxs:
        rt=f.write(str(idx)+'\n')

rowmap={}
for idx in uniqueIdxs:
        rowmap[idx]=len(rowmap)+1 # 1-based

outFile="{0}/Exon_RowIDs.txt".format(outdir) #row indx in exported matrix
with open(outFile, 'w') as f:
    for idx in finalNearestTCIDIdxs:
        rt=f.write(str(rowmap[idx])+'\n')

