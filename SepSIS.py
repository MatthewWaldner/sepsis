from recycle_utils import *
from utils import *
import os, sys, collections, pysam, copy, itertools, argparse


############# HANDLE INPUTS START #############
parser = argparse.ArgumentParser()
parser.add_argument("--RUNMODE", "-RM", help="Sets the manner in which a subsequence is evaluated as strain-unique. Options: 'ORGANIC_Z', 'ORGANIC_P', or 'SYNTH'")
parser.add_argument("--SUBMODE", "-SM", help="Sets the graph subset to be evaluated. CYCLIC is recommended in most cases. Options: 'CYCLIC', 'ISOLATED', 'BOTH'")
parser.add_argument("--fastgFileIn", "-F", help="The path to the SPAdes output .fastg file.")
parser.add_argument("--ScoreValue", "-S", help="The thresholding value used to evaluate if a subsequence is strain-specific.")
parser.add_argument("--outDirectory", "-O", help="The path to the output directory for the output .fasta sequence files.")
parser.add_argument("--outSuffix", "-OS", help="An optional string to add to the output files.")
parser.add_argument("--bamFileIn", "-B", help="The path to the .BAM file used in the SYNTH RUNMODE.")
args = parser.parse_args()

if args.RUNMODE:
    if args.RUNMODE in ['ORGANIC_Z', 'ORGANIC_P', 'SYNTH']:
        RUNMODE = args.RUNMODE
    else:
        raise Exception('Enter a valid RUNMODE. You have entered %s' % args.RUNMODE)
else:
    raise Exception('RUNMODE not entered.')

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
if args.ScoreValue:
    if RUNMODE == 'ORGANIC_Z':
        ScoreValue = args.ScoreValue
    elif RUNMODE == 'ORGANIC_P':
        if float(args.ScoreValue) < 0 or float(args.ScoreValue) > 100:
            raise Exception('Given RUNMODE ORGANIC_P, ScoreValue must  within the range [0,100]. You have entered %s' % args.ScoreValue)
        else:
            ScoreValue = args.ScoreValue
    elif RUNMODE == 'SYNTH':
        if float(args.ScoreValue) < 0 or float(args.ScoreValue) > 1:
            raise Exception('Given RUNMODE SYNTH, ScoreValue must  within the range [0,1]. You have entered %s' % args.ScoreValue)
        else:
            ScoreValue = float(args.ScoreValue)
    else:
        raise Exception('ERROR: RUNMODE Unrecognized Late')
else:
    raise Exception('ScoreValue not entered.')

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
    raise Exception('BAM file argument empty when RUNMODE is SYNTH')


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

#print('GRAPHS done')

if (SUBMODE == 'CYCLIC'):
    for component in listOfGraphOfCSCCs_Components:
        listOfSuccessorNodes += getAllSuccessorNodes(component, GraphOfCSCCs)
        listOfPredecessorNodes += getAllPredecessorsNodes(component, GraphOfCSCCs)
    for node in listOfSuccessorNodes:
        listOfSuccPaths += traverseFromCSCCNode(node, True, 'successor', GraphOfCSCCs, listOfSuccessorNodes, listOfPredecessorNodes)
    for node in listOfPredecessorNodes:
        listOfPredPaths += traverseFromCSCCNode(node, True, 'predecessor', GraphOfCSCCs, listOfSuccessorNodes, listOfPredecessorNodes)
    listOfCSCCPaths = listOfSuccPaths + listOfPredPaths
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
    # removes duplicates
    uniqueListOfPaths = list(listOfBothPaths for listOfBothPaths, _ in itertools.groupby(listOfBothPaths))


#print('Paths done')


if (RUNMODE == 'ORGANIC_Z'):
    if (SUBMODE == 'CYCLIC'):
        cov_vals = [get_cov_from_spades_name(n) for n in GraphOfCSCCs.nodes()]
    elif (SUBMODE == 'ISOLATED'):
        cov_vals = [get_cov_from_spades_name(n) for n in GraphOfISCCs.nodes()]
    elif (SUBMODE == 'BOTH'):
        cov_vals = [get_cov_from_spades_name(n) for n in WholeGraph.nodes()]
    else:
        raise Exception("Specify either CYCLIC, ISOLATED, or BOTH for argument 2.")
    mean_cov = np.mean(cov_vals)
    std_cov = np.std(cov_vals)
    MergedSequences = mergeZTestSequences(uniqueListOfPaths, mean_cov, std_cov, ScoreValue)
    ThresholdedSequences_BothEnds = thresholdSubsequencesAfterSplit_BothEnds_ZTest(MergedSequences, mean_cov,
                                                                                   std_cov, ScoreValue)
    ThresholdedSequences_FrontEnds = thresholdSubsequencesAfterSplit_FrontEnd_ZTest(MergedSequences, mean_cov,
                                                                                    std_cov, ScoreValue)
    ThresholdedSequences_BackEnds = thresholdSubsequencesAfterSplit_BackEnd_ZTest(MergedSequences, mean_cov,
                                                                                  std_cov, ScoreValue)
elif (RUNMODE) == 'ORGANIC_P':
    if (SUBMODE == 'CYCLIC'):
        cov_vals = [get_cov_from_spades_name(n) for n in GraphOfCSCCs.nodes()]
    elif (SUBMODE == 'ISOLATED'):
        cov_vals = [get_cov_from_spades_name(n) for n in GraphOfISCCs.nodes()]
    elif (SUBMODE == 'BOTH'):
        cov_vals = [get_cov_from_spades_name(n) for n in WholeGraph.nodes()]
    else:
        raise Exception("Specify either CYCLIC, ISOLATED, or BOTH for argument 2.")
    qualityTolerance = np.percentile(cov_vals, ScoreValue)
    #print(qualityTolerance)
    MergedSequences = mergePercentileSequences(uniqueListOfPaths, qualityTolerance)
    ThresholdedSequences_BothEnds = thresholdSubsequencesAfterSplit_BothEnds_Percentile(MergedSequences,
                                                                                        qualityTolerance)
    ThresholdedSequences_FrontEnds = thresholdSubsequencesAfterSplit_FrontEnd_Percentile(MergedSequences,
                                                                                         qualityTolerance)
    ThresholdedSequences_BackEnds = thresholdSubsequencesAfterSplit_BackEnd_Percentile(MergedSequences,
                                                                                       qualityTolerance)
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
    MergedSequences = mergeReadPercentageSequences(uniqueListOfPaths, ScoreValue, Dictionary)
    ThresholdedSequences_BothEnds = thresholdSubsequencesAfterSplit_BothEnds_ReadPercentage(MergedSequences,
                                                                                            ScoreValue,
                                                                                            Dictionary)
    ThresholdedSequences_FrontEnds = thresholdSubsequencesAfterSplit_FrontEnd_ReadPercentage(MergedSequences,
                                                                                             ScoreValue,
                                                                                             Dictionary)
    ThresholdedSequences_BackEnds = thresholdSubsequencesAfterSplit_BackEnd_ReadPercentage(MergedSequences,
                                                                                           ScoreValue,
                                                                                           Dictionary)
    FinalDictionary = Dictionary
else:
    raise Exception("Specify either ORGANIC_Z, ORGANIC_P, or SYNTH for argument 1.")

#print('Merges and Thresholds done')

#raise Exception("OOF")

BackEndsListBoolean = [False for i in range(len(ThresholdedSequences_BackEnds))] #need to remove backlist from list if merged int new list
MergedList = []
PostEndsMerge_FrontEndsList = []
PostEndsMerge_BackEndsList = []
for FrontSeq in ThresholdedSequences_FrontEnds:
    FrontSeqMergedBoolean = False
    backCounter = 0
    for BackSeq in ThresholdedSequences_BackEnds:
        if (FrontSeq[-1] == BackSeq[0]):
            FrontSeqMergedBoolean = True
            mergedSequence = FrontSeq + BackSeq[1:]
            MergedList = MergedList + [mergedSequence]
            BackEndsListBoolean[backCounter] = True
        backCounter = backCounter + 1
    if not FrontSeqMergedBoolean:
        PostEndsMerge_FrontEndsList = PostEndsMerge_FrontEndsList + [FrontSeq]
backCounter = 0
for BackSeqBool in BackEndsListBoolean:
    if not BackSeqBool:
        PostEndsMerge_BackEndsList = PostEndsMerge_BackEndsList + [ThresholdedSequences_BackEnds[backCounter]]
    backCounter = backCounter + 1

Final_BothEndsList = MergedList + ThresholdedSequences_BothEnds

Final_FrontEndsList = []
Final_BackEndsList = []

#print('Ends Merged done')

for FrontSeq in PostEndsMerge_FrontEndsList:
    for BothSeq in Final_BothEndsList:
        if not testIfOrderedSubset(FrontSeq, BothSeq):
            Final_FrontEndsList = Final_FrontEndsList + [FrontSeq]

for BackSeq in PostEndsMerge_BackEndsList:
    for BothSeq in Final_BothEndsList:
        if not testIfOrderedSubset(BackSeq, BothSeq):
            Final_BackEndsList = Final_BackEndsList + [BackSeq]


UniqueSequences_BothEnds = list(Final_BothEndsList for Final_BothEndsList, _ in itertools.groupby(Final_BothEndsList))
UniqueSequences_FrontEnds = list(Final_FrontEndsList for Final_FrontEndsList, _ in itertools.groupby(Final_FrontEndsList))
UniqueSequences_BackEnds = list(Final_BackEndsList for Final_BackEndsList, _ in itertools.groupby(Final_BackEndsList))
#UniqueThresholdedSequences_NoEnds = list(ThresholdedSequences_NoEnds for ThresholdedSequences_NoEnds, _ in itertools.groupby(ThresholdedSequences_NoEnds))

#print('Dups Removed done')

if (RUNMODE == 'SYNTH'):
    writeSequenceToFileSynth(outDirectory + "/" + RUNMODE + "_" + SUBMODE+"_BothEnds_" + outSuffix,
                             UniqueSequences_BothEnds, sequenceDictionary, FinalDictionary)
    writeSequenceToFileSynth(outDirectory + "/" + RUNMODE + "_" + SUBMODE + "_FrontEnds_" + outSuffix,
                             UniqueSequences_FrontEnds, sequenceDictionary, FinalDictionary)
    writeSequenceToFileSynth(outDirectory + "/" + RUNMODE + "_" + SUBMODE + "_BackEnds_" + outSuffix,
                             UniqueSequences_BackEnds, sequenceDictionary, FinalDictionary)
    #writeSequenceToFileSynth(out_Directory + "/" + RUNMODE + "_" + SUBMODE + "_NoEnds_" + out_Suffix,
    #                         UniqueThresholdedSequences_NoEnds, sequenceDictionary, FinalDictionary)
else:
    writeSequenceToFile(outDirectory + "/" + RUNMODE + "_" + SUBMODE + "_BothEnds_" + outSuffix + ".fasta", UniqueSequences_BothEnds,
                        sequenceDictionary)
    writeSequenceToFile(outDirectory + "/" + RUNMODE + "_" + SUBMODE + "_FrontEnds_" + outSuffix + ".fasta", UniqueSequences_FrontEnds,
                        sequenceDictionary)
    writeSequenceToFile(outDirectory + "/" + RUNMODE + "_" + SUBMODE + "_BackEnds_" + outSuffix + ".fasta", UniqueSequences_BackEnds,
                        sequenceDictionary)
    #writeSequenceToFile(out_Directory + "/" + RUNMODE + "_" + SUBMODE + "_NoEnds_" + out_Suffix, UniqueThresholdedSequences_NoEnds,
    #                    sequenceDictionary)
