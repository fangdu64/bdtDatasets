import sys, traceback, array, csv, random,re
requiredPaths=["E:\\iBS\\trunk\\core\\src\\iBSPyTest",
    "D:\\CS\\Ice\\Ice-3.5.1\\Install\\Ice-3.5.1\\python"]
for path in requiredPaths:
    if path not in sys.path:
        sys.path.append(path)

import iBSFCDClient as fcdc
import iBS
import ExonArray.Module_TTS_TCIDs as ExonTTS


(TTSBinIdxs,TCIDIdxs,TTSRefs)=ExonTTS.GetTTSwithUniqueTCIDs100bp_FullInfo()
sortedIdxs = [i[0] for i in sorted(enumerate(TTSBinIdxs), key=lambda x:x[1])]

#
# TTS wiht distinct Locations and corresponding to unqiue TCIDs
#
sortedIdx_uniqueTTSBinIdxs=[]
preTTBBinIdx=-1
for i in sortedIdxs:
    if TTSBinIdxs[i] !=preTTBBinIdx:
        sortedIdx_uniqueTTSBinIdxs.append(i)
        preTTBBinIdx=TTSBinIdxs[i]

TTSBinIdxs_temp=[TTSBinIdxs[i] for i in sortedIdx_uniqueTTSBinIdxs]
TTSBinIdxs = TTSBinIdxs_temp

TCIDIdxs_temp=[TCIDIdxs[i] for i in sortedIdx_uniqueTTSBinIdxs]
TCIDIdxs = TCIDIdxs_temp

TTSRefs_temp=[TTSRefs[i] for i in sortedIdx_uniqueTTSBinIdxs]
TTSRefs = TTSRefs_temp
del TTSBinIdxs_temp,TCIDIdxs_temp,TTSRefs_temp




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
validNearestTCIDIdxs=[]
idx=-1
for featureIdx in featureIdxs:
    idx=idx+1
    i=featureIdx-DNaseFeatureIdxFrom # convert ot local idx
    tcidIdx=DNaseFeatureIdx2NearestTCIDIdx[i]
    if tcidIdx!=-1:
        validDNaseFeatureIdxs.append(featureIdx)
        validNearestTCIDIdxs.append(tcidIdx)

rowIdxs=list(range(len(validDNaseFeatureIdxs))) # len =254371
finalDNaseFeatureIdxs=[]
finalNearestTCIDIdxs=[]
finalHighDNaseFeatureIdxs=[]
for i in rowIdxs:
    finalDNaseFeatureIdxs.append(validDNaseFeatureIdxs[i])
    finalNearestTCIDIdxs.append(validNearestTCIDIdxs[i])

#==========================
# Export Ordered TTSs
#==========================
outdir="E:/iBS/trunk/analysis/iBS.Projects/BDVD/DNaseExonCorrelation/100bp/s04-NearbyTSS"



#==========================
# Export Transcript Cluster Idxs
#==========================
outFile="{0}/Exon_FeatureIdxs.txt".format(outdir)
with open(outFile, 'w') as f:
    for idx in TCIDIdxs:
        rt=f.write(str(idx)+'\n')

uniqueIdxs=sorted(set(TCIDIdxs))
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

TCIDIdx2LocationTblRowID={}
outFile="{0}/OrderedTTS_Tbl.txt".format(outdir)
outf = open(outFile, "w")
outf.write("{0}\t{1}\t{2}\t{3}\n".format("TTSBinIdx","DataRowID","TCIDIdx","TTSRef"))
for i in range(len(TTSBinIdxs)):
    rt=outf.write("{0}\t{1}\t{2}\t{3}\n".format(TTSBinIdxs[i],rowmap[TCIDIdxs[i]],TCIDIdxs[i],TTSRefs[i]))
    if TCIDIdxs[i] not in TCIDIdx2LocationTblRowID:
        TCIDIdx2LocationTblRowID[TCIDIdxs[i]]=i+1
outf.close()

outFile="{0}/Exon_LocationRowIDs.txt".format(outdir) #rowID in location table
with open(outFile, 'w') as f:
    for idx in finalNearestTCIDIdxs:
        rt=f.write(str(TCIDIdx2LocationTblRowID[idx])+'\n')

#==========================
# Export DNase
#==========================

outFile="{0}/DNase_FeatureIdxs.txt".format(outdir)
with open(outFile, 'w') as f:
    for idx in finalDNaseFeatureIdxs:
        rt=f.write(str(idx)+'\n')

uniqueIdxs=sorted(set(finalDNaseFeatureIdxs))
outFile="{0}/DNase_UniqueFeatureIdxs.txt".format(outdir)
with open(outFile, 'w') as f:
    for idx in uniqueIdxs:
        rt=f.write(str(idx)+'\n')

rowmap={}
for idx in uniqueIdxs:
        rowmap[idx]=len(rowmap)+1 # 1-based

outFile="{0}/DNase_RowIDs.txt".format(outdir) #row indx in exported matrix
with open(outFile, 'w') as f:
    for idx in finalDNaseFeatureIdxs:
        rt=f.write(str(rowmap[idx])+'\n')
