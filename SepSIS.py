from recycle_utils import *
from utils import *
import os, sys, collections, pysam, copy, itertools, argparse


############# HANDLE INPUTS START #############
parser = argparse.ArgumentParser()
parser.add_argument("--RUNMODE", "-RM", help="Sets the manner in which a subsequence is evaluated as strain-unique. Options: 'ORGANIC_Z', 'ORGANIC_P', or 'SYNTH'")
parser.add_argument("--SUBMODE", "-SM", help="Sets the graph subset to be evaluated. CYCLIC is recommended in most cases. Options: 'CYCLIC', 'ISOLATED', 'BOTH'")
parser.add_argument("--fastgFileIn", "-F", help="The path to the SPAdes output .fastg file.")
parser.add_argument("--MaxScoreValue", "-MAXS", help="The max thresholding value used to evaluate if a subsequence is strain-specific.")
parser.add_argument("--MinScoreValue", "-MINS", help="The min thresholding value used to evaluate if a subsequence is strain-specific.")
#--need min score value and max score value
parser.add_argument("--outDirectory", "-O", help="The path to the output directory for the output .fasta sequence files.")
parser.add_argument("--kmerLength", "-K", help="The length of the kmer setting used by SPAdes.")
parser.add_argument("--outSuffix", "-OS", help="An optional string to add to the output files.")
parser.add_argument("--bamFileIn", "-B", help="The path to the .BAM file used in the SYNTH RUNMODE.")
parser.add_argument("--maxPathNodeLength", "-MP", help="The maximum number of nodes making up a potential strain-specifc sequence.")
#parser.add_argument("--decidingCoverageSection", "-DC", help="The section from which the thresholding Z-Score or percentile is calculated.")
args = parser.parse_args()

if args.RUNMODE:
    if args.RUNMODE in ['ORGANIC_Z', 'ORGANIC_P', 'SYNTH']:
        RUNMODE = args.RUNMODE
    else:
        raise Exception('Enter a valid RUNMODE. You have entered %s' % args.RUNMODE)
else:
    raise Exception('RUNMODE not entered.')

#if args.decidingCoverageSection:
#    if(args.RUNMODE == 'SYNTH'):
#        raise Exception('Coverage section selection not needed for SYNTH RUNMODE.')
#    if args.decidingCoverageSection == 'all' or args.decidingCoverageSection == 'branch':
#        decidingCoverageSection = args.decidingCoverageSection
#    else:
#        raise Exception('Please enter either \'all\' or \'branch\' for the value of --decidingCoverageSection.')

decidingCoverageSection = 'all'

if args.SUBMODE:
    if args.SUBMODE in ['CYCLIC', 'ISOLATED', 'BOTH']:
        SUBMODE = args.SUBMODE
    else:
        raise Exception('Enter a valid SUBMODE. You have entered %s' % args.SUBMODE)
else:
    raise Exception('SUBMODE not entered.')

if args.fastgFileIn:
    if not os.path.exists(args.fastgFileIn):
        raise Exception('The fastg file does not exist at the designated location.')
    else:
        fastgFileIn = args.fastgFileIn
else:
    raise Exception('fastgFileIn not entered.')

#Zscore cutoffs are not checked for validity. Please make sure to unter a valid one
#Examples:
# 90% = 1.282
# 50% = 0
# 10% = -1.282
if args.MaxScoreValue:
    if RUNMODE == 'ORGANIC_Z':
        MaxScoreValue = float(args.MaxScoreValue)
    elif RUNMODE == 'ORGANIC_P':
        if float(args.MaxScoreValue) < 0 or float(args.MaxScoreValue) > 100:
            raise Exception('Given RUNMODE ORGANIC_P, MaxScoreValue must  within the range [0,100]. You have entered %s' % args.MaxScoreValue)
        else:
            MaxScoreValue = float(args.MaxScoreValue)
    elif RUNMODE == 'SYNTH':
        if float(args.MaxScoreValue) < 0 or float(args.MaxScoreValue) > 1:
            raise Exception('Given RUNMODE SYNTH, MaxScoreValue must  within the range [0,1]. You have entered %s' % args.MaxScoreValue)
        else:
            MaxScoreValue = float(args.MaxScoreValue)
    else:
        raise Exception('ERROR: RUNMODE Unrecognized Late')
else:
    raise Exception('MaxScoreValue not entered.')

if args.MinScoreValue:
    if RUNMODE == 'ORGANIC_Z':
        MinScoreValue = float(args.MinScoreValue)
    elif RUNMODE == 'ORGANIC_P':
        if float(args.MinScoreValue) < 0 or float(args.MinScoreValue) > 100:
            raise Exception('Given RUNMODE ORGANIC_P, MinScoreValue must  within the range [0,100]. You have entered %s' % args.MinScoreValue)
        else:
            MinScoreValue = float(args.MinScoreValue)
    elif RUNMODE == 'SYNTH':
        if float(args.MinScoreValue) < 0 or float(args.MinScoreValue) > 1:
            raise Exception('Given RUNMODE SYNTH, MinScoreValue must  within the range [0,1]. You have entered %s' % args.MinScoreValue)
        else:
            MinScoreValue = float(args.MinScoreValue)
    else:
        raise Exception('ERROR: RUNMODE Unrecognized Late')
else:
    raise Exception('MinScoreValue not entered.')

if args.outDirectory:
    if not os.path.exists(args.outDirectory):
        raise Exception('The output directory does not exist.')
    else:
        outDirectory = args.outDirectory
else:
    raise Exception('The output directory was not entered.')

#outSuffix is optional
if args.outSuffix:
    outSuffix = args.outSuffix
else:
    outSuffix = ''

if args.bamFileIn:
    if RUNMODE != 'SYNTH':
        raise Exception('BAM file argument not empty when RUNMODE is not SYNTH')
    if not os.path.exists(args.bamFileIn):
        raise Exception('The BAM file does not exist at the designated location.')
    bamFileIn = args.bamFileIn
    bamFile = pysam.AlignmentFile(bamFileIn)
elif RUNMODE == 'SYNTH':
    raise Exception('BAM file argument empty when RUNMODE is SYNTH.')

if args.kmerLength:
    if int(args.kmerLength) > 1:
        kmerLength = int(args.kmerLength)
    else:
        raise Exception('Invalid kmer length enterd')
else:
    raise Exception('The kmer size used by the SPAdes assembly was not entered.')

if args.maxPathNodeLength:
    if int(args.maxPathNodeLength) > 1:
        maxPathLength = int(args.maxPathNodeLength)
    else:
        raise Exception('Please enter a valid path node length.')
else:
    raise Exception('The maximum number of nodes in strain specific sequence path was not entered.')
############# HANDLE INPUTS END #############

#The wholegraph must exist across all submodes
digraphG = get_fastg_digraph(fastgFileIn)

WholeGraph = copy.deepcopy(digraphG)
WholeGraph_Components = nx.strongly_connected_component_subgraphs(WholeGraph)
listOfWholeGraph_Components = list(WholeGraph_Components)

sequenceDictionary = get_fastg_seqs_dict(fastgFileIn, digraphG)

if (SUBMODE == 'CYCLIC'):
    ####
    #create the second copy of the digraph for the multiple node strongly connected components
    GraphOfCSCCs = copy.deepcopy(digraphG)
    SCC_Group_Components = nx.strongly_connected_component_subgraphs(GraphOfCSCCs)
    #Create a subgraph of all of the single nodes in the graph that are a part of a SCC of 2 or more nodes.
    for component in list(SCC_Group_Components):
        if len(component.nodes()) == 1:
            #component.nodes()
            GraphOfCSCCs.remove_nodes_from(component.nodes())
    #GraphOfCSCCs is now pruned
    #Create a list of each SCC for later use
    GraphOfCSCCs_Components = nx.strongly_connected_component_subgraphs(GraphOfCSCCs)
    #ABSOLUTELY NECESSARY TO CONVERT TO COMPONENTS TO LIST ONLY ONCE, AS FOR SOME REASON THE COMPONENTS
    ## WILL DISAPPEAR AFTER 1 ITERATION
    listOfGraphOfCSCCs_Components = list(GraphOfCSCCs_Components)


if (SUBMODE == 'ISOLATED'):
    ####
    #create the first copy of the digraph for the singles in the list strongly connected components
    GraphOfISCCs = copy.deepcopy(digraphG)
    #Create a subgraph of all of the single nodes in the graph that are not a part of a SCC of 2 or more nodes.
    SCC_Singles_Components = nx.strongly_connected_component_subgraphs(GraphOfISCCs)
    for component in list(SCC_Singles_Components):
        if len(component.nodes()) > 1:
            #component.nodes()
            GraphOfISCCs.remove_nodes_from(component.nodes())
    #GraphOfISCCs is now pruned
    #Create a list of each SCC for later use
    #the list is necessary due to the components only able to be iterated over once in subgraph
    GraphOfISCCs_Components = nx.strongly_connected_component_subgraphs(GraphOfISCCs)
    listOfGraphOfISCCs_Components = list(GraphOfISCCs_Components)



listOfPredPaths = list()
listOfSuccPaths = list()
listOfStartPaths = list()
listOfEndPaths = list()

listOfSuccessorNodes = list()
listOfPredecessorNodes = list()
listOfStartNodes = list()
listOfEndNodes = list()

uniqueListOfPaths = list()

print('GRAPHS done')

if (SUBMODE == 'CYCLIC'):
    for component in listOfGraphOfCSCCs_Components:
        listOfSuccessorNodes += getAllSuccessorNodes(component, GraphOfCSCCs)
        listOfPredecessorNodes += getAllPredecessorsNodes(component, GraphOfCSCCs)
    for node in listOfSuccessorNodes:
        listOfSuccPaths += traverseFromCSCCNode(node, True, 'successor', GraphOfCSCCs, listOfSuccessorNodes, listOfPredecessorNodes)
    for node in listOfPredecessorNodes:
        listOfPredPaths += traverseFromCSCCNode(node, True, 'predecessor', GraphOfCSCCs, listOfSuccessorNodes, listOfPredecessorNodes)
    listOfCSCCPaths = listOfSuccPaths + listOfPredPaths
    listOfPrimaryNodes = list(set(listOfSuccessorNodes + listOfPredecessorNodes))
    #removes duplicates
    uniqueListOfPaths = list(listOfCSCCPaths for listOfCSCCPaths,_ in itertools.groupby(listOfCSCCPaths))

if (SUBMODE == 'ISOLATED'):
    for component in listOfGraphOfISCCs_Components:
        listOfStartNodes += getAllNodesWithoutPredecessors(component, GraphOfISCCs)
        listOfSuccessorNodes += getAllSuccessorNodes(component, GraphOfISCCs)
        listOfPredecessorNodes += getAllPredecessorsNodes(component, GraphOfISCCs)
        listOfEndNodes += getAllNodesWithoutSuccessors(component, GraphOfISCCs)
    for node in listOfStartNodes:
        listOfStartPaths += traverseFromISCCNode(node, True, 'start', GraphOfISCCs, listOfSuccessorNodes, listOfPredecessorNodes, listOfStartNodes, listOfEndNodes)
    for node in listOfSuccessorNodes:
        listOfSuccPaths += traverseFromISCCNode(node, True, 'successor', GraphOfISCCs, listOfSuccessorNodes, listOfPredecessorNodes, listOfStartNodes, listOfEndNodes)
    for node in listOfPredecessorNodes:
        listOfPredPaths += traverseFromISCCNode(node, True, 'predecessor', GraphOfISCCs, listOfSuccessorNodes, listOfPredecessorNodes, listOfStartNodes, listOfEndNodes)
    for node in listOfEndNodes:
        listOfEndPaths += traverseFromISCCNode(node, True, 'end', GraphOfISCCs, listOfSuccessorNodes, listOfPredecessorNodes, listOfStartNodes, listOfEndNodes)
    listOfISCCPaths = listOfSuccPaths + listOfPredPaths + listOfStartPaths + listOfEndPaths
    listOfPrimaryNodes = list(set(listOfStartNodes + listOfSuccessorNodes + listOfPredecessorNodes + listOfEndNodes))
    #removes duplicates
    uniqueListOfPaths = list(listOfISCCPaths for listOfISCCPaths,_ in itertools.groupby(listOfISCCPaths))

if (SUBMODE == 'BOTH'):
    for component in listOfWholeGraph_Components:
        listOfStartNodes += getAllNodesWithoutPredecessors(component, WholeGraph)
        listOfSuccessorNodes += getAllSuccessorNodes(component, WholeGraph)
        listOfPredecessorNodes += getAllPredecessorsNodes(component, WholeGraph)
        listOfEndNodes += getAllNodesWithoutSuccessors(component, WholeGraph)
    for node in listOfStartNodes:
        listOfStartPaths += traverseFromISCCNode(node, True, 'start', WholeGraph, listOfSuccessorNodes, listOfPredecessorNodes, listOfStartNodes, listOfEndNodes)
    for node in listOfSuccessorNodes:
        listOfSuccPaths += traverseFromISCCNode(node, True, 'successor', WholeGraph, listOfSuccessorNodes, listOfPredecessorNodes, listOfStartNodes, listOfEndNodes)
    for node in listOfPredecessorNodes:
        listOfPredPaths += traverseFromISCCNode(node, True, 'predecessor', WholeGraph, listOfSuccessorNodes, listOfPredecessorNodes, listOfStartNodes, listOfEndNodes)
    for node in listOfEndNodes:
        listOfEndPaths += traverseFromISCCNode(node, True, 'end', WholeGraph, listOfSuccessorNodes, listOfPredecessorNodes, listOfStartNodes, listOfEndNodes)
    listOfBothPaths = listOfSuccPaths + listOfPredPaths + listOfStartPaths + listOfEndPaths
    listOfPrimaryNodes = list(set(listOfStartNodes + listOfSuccessorNodes + listOfPredecessorNodes + listOfEndNodes))
    # removes duplicates
    uniqueListOfPaths = list(listOfBothPaths for listOfBothPaths, _ in itertools.groupby(listOfBothPaths))


print('Paths done')


if (RUNMODE == 'ORGANIC_Z'):
    if decidingCoverageSection == 'branch':
        cov_vals = [get_cov_from_spades_name(n) for n in listOfPrimaryNodes]
    elif (SUBMODE == 'CYCLIC'):
        cov_vals = [get_cov_from_spades_name(n) for n in GraphOfCSCCs.nodes()]
    elif (SUBMODE == 'ISOLATED'):
        cov_vals = [get_cov_from_spades_name(n) for n in GraphOfISCCs.nodes()]
    elif (SUBMODE == 'BOTH'):
        cov_vals = [get_cov_from_spades_name(n) for n in WholeGraph.nodes()]
    else:
        raise Exception("Specify either CYCLIC, ISOLATED, or BOTH for argument 2.")
    mean_cov = np.mean(cov_vals)
    std_cov = np.std(cov_vals)
    MergedSequences = mergeZTestSequences(uniqueListOfPaths, mean_cov, std_cov, MinScoreValue, MaxScoreValue, maxPathLength)
    print('MergeDone')
    ThresholdedSequences_BothEnds = thresholdSubsequencesAfterSplit_BothEnds_ZTest(MergedSequences, mean_cov,
                                                                                   std_cov, MinScoreValue, MaxScoreValue)
    print('bothDone')
    ThresholdedSequences_FrontEnds = thresholdSubsequencesAfterSplit_FrontEnd_ZTest(MergedSequences, mean_cov,
                                                                                    std_cov, MinScoreValue, MaxScoreValue)
    print('frontDone')
    ThresholdedSequences_BackEnds = thresholdSubsequencesAfterSplit_BackEnd_ZTest(MergedSequences, mean_cov,
                                                                                  std_cov, MinScoreValue, MaxScoreValue)
    print('backDone')
elif (RUNMODE) == 'ORGANIC_P':
    if decidingCoverageSection == 'branch':
        cov_vals = [get_cov_from_spades_name(n) for n in listOfPrimaryNodes]
    if (SUBMODE == 'CYCLIC'):
        cov_vals = [get_cov_from_spades_name(n) for n in GraphOfCSCCs.nodes()]
    elif (SUBMODE == 'ISOLATED'):
        cov_vals = [get_cov_from_spades_name(n) for n in GraphOfISCCs.nodes()]
    elif (SUBMODE == 'BOTH'):
        cov_vals = [get_cov_from_spades_name(n) for n in WholeGraph.nodes()]
    else:
        raise Exception("Specify either CYCLIC, ISOLATED, or BOTH for argument 2.")
    MinQualityTolerance = np.percentile(cov_vals, MinScoreValue)
    MaxQualityTolerance = np.percentile(cov_vals, MaxScoreValue)
    #print(qualityTolerance)
    MergedSequences = mergePercentileSequences(uniqueListOfPaths, MinQualityTolerance, MaxQualityTolerance, maxPathLength)
    ThresholdedSequences_BothEnds = thresholdSubsequencesAfterSplit_BothEnds_Percentile(MergedSequences,
                                                                                        MinQualityTolerance,
                                                                                        MaxQualityTolerance)
    ThresholdedSequences_FrontEnds = thresholdSubsequencesAfterSplit_FrontEnd_Percentile(MergedSequences,
                                                                                         MinQualityTolerance,
                                                                                         MaxQualityTolerance)
    ThresholdedSequences_BackEnds = thresholdSubsequencesAfterSplit_BackEnd_Percentile(MergedSequences,
                                                                                       MinQualityTolerance,
                                                                                       MaxQualityTolerance)
elif (RUNMODE) == 'SYNTH':
    if (SUBMODE == 'CYCLIC'):
        Dictionary = {}
        for component in listOfGraphOfCSCCs_Components:
            for node in component:
                nodeRP = getReadPercentage(node, bamFile)
                nodeDR = getDominantRead(node, bamFile)
                Dictionary[node] = [nodeDR,nodeRP]
    elif (SUBMODE == 'ISOLATED'):
        Dictionary = {}
        for component in listOfGraphOfISCCs_Components:
            for node in component:
                nodeRP = getReadPercentage(node, bamFile)
                nodeDR = getDominantRead(node, bamFile)
                Dictionary[node] = [nodeDR,nodeRP]
    elif (SUBMODE == 'BOTH'):
        Dictionary = {}
        for component in listOfWholeGraph_Components:
            for node in component:
                nodeRP = getReadPercentage(node, bamFile)
                nodeDR = getDominantRead(node, bamFile)
                Dictionary[node] = [nodeDR, nodeRP]
    else:
        raise Exception("Specify either CYCLIC, ISOLATED, or BOTH for argument 2.")
    MergedSequences = mergeReadPercentageSequences(uniqueListOfPaths, MinScoreValue, MaxScoreValue, Dictionary,maxPathLength)
    ThresholdedSequences_BothEnds = thresholdSubsequencesAfterSplit_BothEnds_ReadPercentage(MergedSequences,
                                                                                            MinScoreValue,
                                                                                            MaxScoreValue,
                                                                                            Dictionary)
    ThresholdedSequences_FrontEnds = thresholdSubsequencesAfterSplit_FrontEnd_ReadPercentage(MergedSequences,
                                                                                             MinScoreValue,
                                                                                             MaxScoreValue,
                                                                                             Dictionary)
    ThresholdedSequences_BackEnds = thresholdSubsequencesAfterSplit_BackEnd_ReadPercentage(MergedSequences,
                                                                                           MinScoreValue,
                                                                                           MaxScoreValue,
                                                                                           Dictionary)
    FinalDictionary = Dictionary
else:
    raise Exception("Specify either ORGANIC_Z, ORGANIC_P, or SYNTH for argument 1.")

#ThresholdedSequences_BothEnds_U = list(ThresholdedSequences_BothEnds for ThresholdedSequences_BothEnds, _ in itertools.groupby(ThresholdedSequences_BothEnds))
#ThresholdedSequences_FrontEnds_U = list(ThresholdedSequences_FrontEnds for ThresholdedSequences_FrontEnds, _ in itertools.groupby(ThresholdedSequences_FrontEnds))
#ThresholdedSequences_BackEnds_U = list(ThresholdedSequences_BackEnds for ThresholdedSequences_BackEnds, _ in itertools.groupby(ThresholdedSequences_BackEnds))

ThresholdedSequences_BothEnds = sorted(ThresholdedSequences_BothEnds)
ThresholdedSequences_BothEnds_U = [ThresholdedSequences_BothEnds[i] for i in range(len(ThresholdedSequences_BothEnds)) if i == 0 or ThresholdedSequences_BothEnds[i] != ThresholdedSequences_BothEnds[i-1]]

ThresholdedSequences_FrontEnds = sorted(ThresholdedSequences_FrontEnds)
ThresholdedSequences_FrontEnds_U = [ThresholdedSequences_FrontEnds[i] for i in range(len(ThresholdedSequences_FrontEnds)) if i == 0 or ThresholdedSequences_FrontEnds[i] != ThresholdedSequences_FrontEnds[i-1]]

ThresholdedSequences_BackEnds = sorted(ThresholdedSequences_BackEnds)
ThresholdedSequences_BackEnds_U = [ThresholdedSequences_BackEnds[i] for i in range(len(ThresholdedSequences_BackEnds)) if i == 0 or ThresholdedSequences_BackEnds[i] != ThresholdedSequences_BackEnds[i-1]]




print('Merges and Thresholds done')
#raise Exception("OOF")

BackEndsListBoolean = [False for i in range(len(ThresholdedSequences_BackEnds_U))] #need to remove backlist from list if merged int new list
MergedList = []
PostEndsMerge_FrontEndsList = []
PostEndsMerge_BackEndsList = []
innerIndex = 0
for FrontSeq in ThresholdedSequences_FrontEnds_U:
    #print(FrontSeq)
    innerIndex += 1
    #print(str(innerIndex) + ' of ' + str(len(ThresholdedSequences_FrontEnds_U)))
    FrontSeqMergedBoolean = False
    backCounter = 0
    for BackSeq in ThresholdedSequences_BackEnds_U:
        bothList = testIfMatchedEndsSubset(FrontSeq, BackSeq)
        if (bothList != []):
            FrontSeqMergedBoolean = True
            MergedList = MergedList + [bothList]
            BackEndsListBoolean[backCounter] = True
        backCounter = backCounter + 1
    if not FrontSeqMergedBoolean:
        PostEndsMerge_FrontEndsList = PostEndsMerge_FrontEndsList + [FrontSeq]
backCounter = 0
for BackSeqBool in BackEndsListBoolean:
    if not BackSeqBool:
        PostEndsMerge_BackEndsList = PostEndsMerge_BackEndsList + [ThresholdedSequences_BackEnds_U[backCounter]]
    backCounter = backCounter + 1

Final_BothEndsList = MergedList + ThresholdedSequences_BothEnds_U

Final_FrontEndsList = []
Final_BackEndsList = []

#def testIfMatchedEndsSubset(frontList, endList):
#    for iterFull in range(len(frontList) - 1):
#        if frontList[iterFull+1] == endList[0] and len(frontList[iterFull+1:]) <= (len(endList)-1):
#            print('FrontLen = '+str(len(frontList[iterFull+1:])) + ' BackLen = ' +str(len(endList)-1))
#            for iterSub in range(len(endList) - 1):
#                if frontList[iterFull+1 + iterSub] == frontList[-1] and frontList[iterFull+1 + iterSub] == endList[iterSub]:
#                    returnList = frontList[:iterFull + 1] + endList
#                    return returnList
#                elif frontList[iterFull+1 + iterSub] == endList[iterSub]:
#                    pass
#                else:
#                    break
#    return []




print('Ends Merged done')

#print('PostEndsMerge_FrontEndsList Length = ' + str(len(PostEndsMerge_FrontEndsList)))
#print('Final_BothEndsList Length = ' + str(len(Final_BothEndsList)))
#print('PostEndsMerge_BackEndsList Length = ' + str(len(PostEndsMerge_BackEndsList)))

#print(Final_BothEndsList)

for FrontSeq in PostEndsMerge_FrontEndsList:
    Subset = False
    for BothSeq in Final_BothEndsList:
        if testIfOrderedSubset(FrontSeq, BothSeq):
            Subset = True
            break
    if Subset == False:
        Final_FrontEndsList = Final_FrontEndsList + [FrontSeq]

for BackSeq in PostEndsMerge_BackEndsList:
    Subset = False
    for BothSeq in Final_BothEndsList:
        if testIfOrderedSubset(BackSeq, BothSeq):
            Subset = True
            break
    if Subset == False:
        Final_BackEndsList = Final_BackEndsList + [BackSeq]

#print(len(ThresholdedSequences_BothEnds_U))
#print(len(ThresholdedSequences_FrontEnds_U))
#print(len(ThresholdedSequences_BackEnds_U))
#print(len(Final_BothEndsList))
#print(len(Final_FrontEndsList))
#print(len(Final_BackEndsList))

print('Loops upon loops done')

Final_BothEndsList = sorted(Final_BothEndsList)
UniqueSequences_BothEnds = [Final_BothEndsList[i] for i in range(len(Final_BothEndsList)) if i == 0 or Final_BothEndsList[i] != Final_BothEndsList[i-1]]

Final_FrontEndsList = sorted(Final_FrontEndsList)
UniqueSequences_FrontEnds = [Final_FrontEndsList[i] for i in range(len(Final_FrontEndsList)) if i == 0 or Final_FrontEndsList[i] != Final_FrontEndsList[i-1]]

Final_BackEndsList = sorted(Final_BackEndsList)
UniqueSequences_BackEnds = [Final_BackEndsList[i] for i in range(len(Final_BackEndsList)) if i == 0 or Final_BackEndsList[i] != Final_BackEndsList[i-1]]


#UniqueSequences_FrontEnds = list(Final_FrontEndsList for Final_FrontEndsList, _ in itertools.groupby(Final_FrontEndsList))
#UniqueSequences_BackEnds = list(Final_BackEndsList for Final_BackEndsList, _ in itertools.groupby(Final_BackEndsList))
#UniqueThresholdedSequences_NoEnds = list(ThresholdedSequences_NoEnds for ThresholdedSequences_NoEnds, _ in itertools.groupby(ThresholdedSequences_NoEnds))

#print('Dups Removed done')

if (RUNMODE == 'SYNTH'):
    writeSequenceToFileSynth(outDirectory + "/" + RUNMODE + "_" + SUBMODE+"_BothEnds_" + outSuffix + ".fasta",
                             UniqueSequences_BothEnds, sequenceDictionary, FinalDictionary, kmerLength)
    writeSequenceToFileSynth(outDirectory + "/" + RUNMODE + "_" + SUBMODE + "_FrontEnds_" + outSuffix + ".fasta",
                             UniqueSequences_FrontEnds, sequenceDictionary, FinalDictionary, kmerLength)
    writeSequenceToFileSynth(outDirectory + "/" + RUNMODE + "_" + SUBMODE + "_BackEnds_" + outSuffix + ".fasta",
                             UniqueSequences_BackEnds, sequenceDictionary, FinalDictionary, kmerLength)
else:
    writeSequenceToFile(outDirectory + "/" + RUNMODE + "_" + SUBMODE + "_BothEnds_" + outSuffix + ".fasta", UniqueSequences_BothEnds,
                        sequenceDictionary, kmerLength)
    writeSequenceToFile(outDirectory + "/" + RUNMODE + "_" + SUBMODE + "_FrontEnds_" + outSuffix + ".fasta", UniqueSequences_FrontEnds,
                        sequenceDictionary, kmerLength)
    writeSequenceToFile(outDirectory + "/" + RUNMODE + "_" + SUBMODE + "_BackEnds_" + outSuffix + ".fasta", UniqueSequences_BackEnds,
                        sequenceDictionary, kmerLength)


#some of the batch combinations are fucked.

#gonna need a tool to extract the tiny unique subsequence since the lengths are bullshit.

#positive or negative analysis using names only. same fastg files.