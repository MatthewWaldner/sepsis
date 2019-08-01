from recycle_utils import *
import collections

maxNodeStringLen = 7
#maxSequenceLen = 0, dont use this due to the variation in length between nodes


def ensureSameDirection(node1,node2):
    if (node1[-1] == "'" and node2[-1] == "'"):
        return True
    elif (node1[-1] != "'" and node2[-1] != "'"):
        return True
    else:
        return False

#retrieve a list of all nodes in a component of a graph such that the node has more than 1 successor node
def getAllSuccessorNodes(component, parentGraph):
    listOfSuccessorNodes = list()
    for node in component.nodes():
        if len(list(parentGraph.successors(node))) > 1:
            listOfSuccessorNodes.append(node)
    return listOfSuccessorNodes

#retrieve a list of all nodes in a component of a graph such that the node has more than 1 predecessor node
def getAllPredecessorsNodes(component, parentGraph):
    listOfPredecessorsNodes = list()
    for node in component.nodes():
        if len(list(parentGraph.predecessors(node))) > 1:
            listOfPredecessorsNodes.append(node)
    return listOfPredecessorsNodes

#retrieve a list of all nodes without predecessors, meaning they are the start of a sequence
def getAllNodesWithoutPredecessors(component, parentGraph):
    nodesWithoutPredecessors = list()
    for node in component.nodes():
        if len(list(parentGraph.predecessors(node))) == 0:
            nodesWithoutPredecessors.append(node)
    return nodesWithoutPredecessors

#retrieve a list of all nodes without sucessors, meaning they are the end of a sequence
def getAllNodesWithoutSuccessors(component, parentGraph):
    nodesWithoutSuccessors = list()
    for node in component.nodes():
        if len(list(parentGraph.successors(node))) == 0:
            nodesWithoutSuccessors.append(node)
    return nodesWithoutSuccessors

###adds the node to the front of a list of paths
def extendListOfPathsFront(node, listOfPaths): #listoflists
    newLop = list()
    for innerPath in listOfPaths:
        innerPath = [node] + innerPath
        newLop.append(innerPath)
    return newLop

####adds a node to the back of a list of paths
def extendListOfPathsBack(node, listOfPaths): #listoflists
    for innerPath in listOfPaths:
        innerPath.append(node)
    return listOfPaths

def getReadPercentage(node, bamFile):
    try:
        x = bamFile.fetch(node)
    except:
        node = rc_node(node)
    numOfMappedReads = 0
    IDList = list()
    ID1List = list()
    for read in bamFile.fetch(node):
        numOfMappedReads += 1
        ID1 = read.query_name.split('_')[0]
        ID2 = read.query_name.split('_')[1]
        FullId = ID1 + '_' + ID2
        ID1List.append(ID1)
        IDList.append(FullId)
    Largest_ID = "No_Reads"
    Largest_IDCount = 0
    for ID in list(set(ID1List)):
        IDCount = collections.Counter(ID1List)[ID]
        if IDCount > Largest_IDCount:
            Largest_IDCount = IDCount
            Largest_ID = ID
        # returnList.append(ID + "," + str(IDCount) + "," + str(IDCount / numOfMappedReads))
    if numOfMappedReads == 0:
        return float(0)
    else:
        return float(Largest_IDCount) / float(numOfMappedReads)

def getDominantRead(node, bamFile):
    try:
        x = bamFile.fetch(node)
    except:
        node = rc_node(node)
    numOfMappedReads = 0
    IDList = list()
    ID1List = list()
    for read in bamFile.fetch(node):
        numOfMappedReads += 1
        ID1 = read.query_name.split('_')[0]
        ID2 = read.query_name.split('_')[1]
        FullId = ID1 + '_' + ID2
        ID1List.append(ID1)
        IDList.append(FullId)
    Largest_ID = "No_Reads"
    Largest_IDCount = 0
    for ID in list(set(ID1List)):
        IDCount = collections.Counter(ID1List)[ID]
        if IDCount > Largest_IDCount:
            Largest_IDCount = IDCount
            Largest_ID = ID
        # returnList.append(ID + "," + str(IDCount) + "," + str(IDCount / numOfMappedReads))
    return(Largest_ID)

def traverseFromCSCCNode(nextNode, initiation, travType, parentGraph, listOfSuccessorNodes, listOfPredecessorNodes):
    if travType == 'successor':
        ##INSERTION TO ALLOW ADDITION OF TERMINAL NODE TO THE PATH
        # BY FINAL I MEAN AN ENCOUNTERED SUCCESSOR (or  later start or end node in ISCC)
        # This will cause duplicate paths that will be dealt with through set conversion.
        if (initiation == False and (nextNode in listOfSuccessorNodes or nextNode in listOfPredecessorNodes)):
            #print('end')
            return extendListOfPathsFront(nextNode, [[]])
        elif len(list(parentGraph.successors(nextNode))) == 1:
            #print("1 successor")
            if ensureSameDirection(nextNode, list(parentGraph.successors(nextNode))[0]):
                futureListOfPaths = traverseFromCSCCNode(list(parentGraph.successors(nextNode))[0], False, travType, parentGraph, listOfSuccessorNodes,listOfPredecessorNodes)
                #print('middle')
                return extendListOfPathsFront(nextNode, futureListOfPaths)
            else:
                #print('fail')
                return [['FAIL']]
        elif len(list(parentGraph.successors(nextNode))) > 1:
            futureListOfPaths = []
            for succNode in list(parentGraph.successors(nextNode)):
                if ensureSameDirection(nextNode,succNode):
                    newPath = extendListOfPathsFront(nextNode, traverseFromCSCCNode(succNode, False, travType, parentGraph, listOfSuccessorNodes, listOfPredecessorNodes))
                    #print(newPath)
                    if newPath[-1][-1] != 'FAIL':
                        futureListOfPaths += newPath
            return futureListOfPaths
        elif len(list(parentGraph.successors(nextNode))) == 0:
            raise Exception("Somehow there were no successors in a CSCC. Something is wrong.")
        else:
            raise Exception("traverseFromSCCNodeBlind (succ) Error, how did you get here?")
    elif travType == 'predecessor': #doublecheck spelling
        if (initiation == False and (nextNode in listOfSuccessorNodes or nextNode in listOfPredecessorNodes)):
            return extendListOfPathsBack(nextNode, [[]])
        elif len(list(parentGraph.predecessors(nextNode))) == 1:
            if ensureSameDirection(nextNode, list(parentGraph.predecessors(nextNode))[0]):
                futureListOfPaths = traverseFromCSCCNode(list(parentGraph.predecessors(nextNode))[0], False, travType, parentGraph, listOfSuccessorNodes,listOfPredecessorNodes)
                return extendListOfPathsBack(nextNode, futureListOfPaths)
            else:
                return [['FAIL']]
        elif len(list(parentGraph.predecessors(nextNode))) > 1:
            futureListOfPaths = []
            for predNode in list(parentGraph.predecessors(nextNode)):
                if ensureSameDirection(nextNode, predNode):
                    newPath = extendListOfPathsBack(nextNode, traverseFromCSCCNode(predNode, False, travType, parentGraph, listOfSuccessorNodes, listOfPredecessorNodes))
                    if newPath[0][0] != 'FAIL':
                        futureListOfPaths += newPath
            return futureListOfPaths
        elif len(list(parentGraph.predecessors(nextNode))) == 0:
            raise Exception("Somehow there were no predecessors in a CSCC. Something is wrong.")
        else:
            raise Exception("traverseFromSCCNodeBlind (pred) Error, how did you get here?")
    else: #something went wrong
        raise Exception("Something other than successor or predecessor was entered")


def traverseFromISCCNode(nextNode, initiation, travType, parentGraph, listOfSuccessorNodes, listOfPredecessorNodes, listOfStartNodes, listOfEndNodes):
    if travType == 'successor' or travType == 'start':
        if (initiation == False and (nextNode in listOfSuccessorNodes or nextNode in listOfPredecessorNodes or nextNode in listOfStartNodes or nextNode in listOfEndNodes)):
            return extendListOfPathsFront(nextNode, [[]])
        elif len(list(parentGraph.successors(nextNode))) == 1:
            if ensureSameDirection(nextNode, list(parentGraph.successors(nextNode))[0]):
                futureListOfPaths = traverseFromISCCNode(list(parentGraph.successors(nextNode))[0], False, travType, parentGraph, listOfSuccessorNodes,listOfPredecessorNodes, listOfStartNodes, listOfEndNodes)
                return extendListOfPathsFront(nextNode, futureListOfPaths)
            else:
                return [['FAIL']]
        elif len(list(parentGraph.successors(nextNode))) > 1:
            futureListOfPaths = []
            for succNode in list(parentGraph.successors(nextNode)):
                if ensureSameDirection(nextNode, succNode):
                    newPath = extendListOfPathsFront(nextNode, traverseFromISCCNode(succNode, False, travType, parentGraph, listOfSuccessorNodes, listOfPredecessorNodes, listOfStartNodes, listOfEndNodes))
                    if newPath[-1][-1] != 'FAIL':
                        futureListOfPaths += newPath
            return futureListOfPaths
        elif len(list(parentGraph.successors(nextNode))) == 0:
            return extendListOfPathsFront(nextNode, [[]])
        else:
            raise Exception("traverseFromISCCNodeBlind (succ) Error, how did you get here?")
    elif travType == 'predecessor' or travType == 'end':  # doublecheck spelling
        if (initiation == False and (nextNode in listOfSuccessorNodes or nextNode in listOfPredecessorNodes or nextNode in listOfStartNodes or nextNode in listOfEndNodes)):
            return extendListOfPathsBack(nextNode, [[]])
        elif len(list(parentGraph.predecessors(nextNode))) == 1:
            if ensureSameDirection(nextNode, list(parentGraph.predecessors(nextNode))[0]):
                futureListOfPaths = traverseFromISCCNode(list(parentGraph.predecessors(nextNode))[0], False, travType, parentGraph, listOfSuccessorNodes,listOfPredecessorNodes, listOfStartNodes, listOfEndNodes)
                return extendListOfPathsBack(nextNode, futureListOfPaths)
            else:
                return [['FAIL']]
        elif len(list(parentGraph.predecessors(nextNode))) > 1:
            futureListOfPaths = []
            for predNode in list(parentGraph.predecessors(nextNode)):
                if ensureSameDirection(nextNode, predNode):
                    newPath = extendListOfPathsBack(nextNode, traverseFromISCCNode(predNode, False, travType, parentGraph, listOfSuccessorNodes, listOfPredecessorNodes, listOfStartNodes, listOfEndNodes))
                    if newPath[0][0] != 'FAIL':
                        futureListOfPaths += newPath
            return futureListOfPaths
        elif len(list(parentGraph.predecessors(nextNode))) == 0:
            return extendListOfPathsBack(nextNode, [[]])
        else:
            raise Exception("traverseFromISCCNodeBlind (pred) Error, how did you get here?")
    else:  # something went wrong
        raise Exception("Something other than successor or predecessor or start or end was entered as trav type")


def withinPercentile(node, percentileValue):
    nodeCoverage = get_cov_from_spades_name(node)
    if (nodeCoverage <= percentileValue):
        return True
    else:
        return False

def withinPercentileBackwardFail(node, nextNode, percentileValue):
    nodeCoverage = get_cov_from_spades_name(node)
    #print(nodeCoverage)
    nextNodeCoverage = get_cov_from_spades_name(nextNode)
    #print(nextNodeCoverage)
    if (nodeCoverage > percentileValue and nextNodeCoverage <= percentileValue):
        return True
    else:
        return False

def withinPercentileForwardFail(node, nextNode, percentileValue):
    nodeCoverage = get_cov_from_spades_name(node)
    #print(nodeCoverage)
    nextNodeCoverage = get_cov_from_spades_name(nextNode)
    #print(nextNodeCoverage)
    if (nodeCoverage <= percentileValue and nextNodeCoverage > percentileValue):
        return True
    else:
        return False

def withinPercentileBothTrue(node, nextNode, percentileValue):
    nodeCoverage = get_cov_from_spades_name(node)
    #print(nodeCoverage)
    nextNodeCoverage = get_cov_from_spades_name(nextNode)
    #print(nextNodeCoverage)
    if (nodeCoverage <= percentileValue and nextNodeCoverage <= percentileValue):
        return True
    else:
        return False

def withinPercentileMerge(node, nextNode, percentileValue):
    nodeCoverage = get_cov_from_spades_name(node)
    #print(nodeCoverage)
    nextNodeCoverage = get_cov_from_spades_name(nextNode)
    #print(nextNodeCoverage)
    if not (nodeCoverage > percentileValue and nextNodeCoverage > percentileValue):
        return True
    else:
        return False

def withinZValue(coverage1, meanCoverage, stdDev, ScoreValue):
    Z1 = (coverage1 - meanCoverage) / stdDev
    if (Z1 <= ScoreValue):
        return True
    else:
        return False

def withinZValueForwardFail(coverage1, coverage2, meanCoverage, stdDev, ScoreValue):
    Z1 = (coverage1 - meanCoverage) / stdDev
    Z2 = (coverage2 - meanCoverage) / stdDev
    # print(Z)
    if (Z1 <= ScoreValue and Z2 > ScoreValue): #at least one score needs to pass to merge.
        #print(Z1)
        #print(Z2)
        return True  # both nodes reject null hyppthesis
    else:
        return False # at least 1 node accepts null hypothesis


def withinZValueBackwardFail(coverage1, coverage2, meanCoverage, stdDev, ScoreValue):
    Z1 = (coverage1 - meanCoverage) / stdDev
    Z2 = (coverage2 - meanCoverage) / stdDev
    # print(Z)
    if (Z1 > ScoreValue and Z2 <= ScoreValue): #at least one score needs to pass to merge.
        #print(Z1)
        #print(Z2)
        return True  # both nodes reject null hyppthesis
    else:
        return False # at least 1 node accepts null hypothesis


def withinZValueBothTrue(coverage1, coverage2, meanCoverage, stdDev, ScoreValue):
    Z1 = (coverage1 - meanCoverage) / stdDev
    Z2 = (coverage2 - meanCoverage) / stdDev
    # print(Z)
    if (Z1 <= ScoreValue and Z2 <= ScoreValue): #at least one score needs to pass to merge.
        #print(Z1)
        #print(Z2)
        return True  # both nodes reject null hyppthesis
    else:
        return False # at least 1 node accepts null hypothesis

def withinZValueMerge(coverage1, coverage2, meanCoverage, stdDev, ScoreValue):
    Z1 = (coverage1 - meanCoverage) / stdDev
    Z2 = (coverage2 - meanCoverage) / stdDev
    # print(Z)
    if not (Z1 > ScoreValue and Z2 > ScoreValue): #at least one score needs to pass to merge.
        #print(Z1)
        #print(Z2)
        return True  # both nodes reject null hyppthesis
    else:
        return False # at least 1 node accepts null hypothesis

def withinReadPercentage(node, readPercentageThreshold, InDictionary):
    nodeFull = InDictionary[node]
    nodeRP = nodeFull[1]
    if (nodeRP >= readPercentageThreshold):
        return True
    else:
        return False

def withinReadPercentageMerge(node, nextNode, readPercentageThreshold, InDictionary):
    nodeFull = InDictionary[node]
    nextNodeFull = InDictionary[nextNode]
    nodeRP = nodeFull[1]
    nodeDR = nodeFull[0]
    nextNodeRP = nextNodeFull[1]
    nextNodeDR = nextNodeFull[0]
    if (not (nodeRP < readPercentageThreshold and nextNodeRP < readPercentageThreshold)) and nodeDR != "No_Reads" and nextNodeDR != "No_Reads":
        return True
    else:
        return False

#there is the case where is morphs from 1 majority read section to another majority read section that takes place. imperfect but deal.

def withinReadPercentageForwardFail(node, nextNode, readPercentageThreshold, InDictionary):
    nodeFull = InDictionary[node]
    nextNodeFull = InDictionary[nextNode]
    nodeRP = nodeFull[1]
    nodeDR = nodeFull[0]
    nextNodeRP = nextNodeFull[1]
    nextNodeDR = nextNodeFull[0]
    if (nodeDR ):
        if (nodeRP >= readPercentageThreshold and nodeDR != "No_Reads" and nextNodeDR != "No_Reads" and ((nextNodeRP < readPercentageThreshold and nextNodeDR == nodeDR) or nextNodeDR != nodeDR)):
            return True
        else:
            return False
    else:
        return False

def withinReadPercentageBackwardFail(node, nextNode, readPercentageThreshold, InDictionary):
    nodeFull = InDictionary[node]
    nextNodeFull = InDictionary[nextNode]
    nodeRP = nodeFull[1]
    nodeDR = nodeFull[0]
    nextNodeRP = nextNodeFull[1]
    nextNodeDR = nextNodeFull[0]
    if (nextNodeRP >= readPercentageThreshold and nodeDR != "No_Reads" and nextNodeDR != "No_Reads" and ((nodeRP < readPercentageThreshold and nextNodeDR == nodeDR) or nextNodeDR != nodeDR)):
        return True
    else:
        return False


def withinReadPercentageBothTrue(node, nextNode, readPercentageThreshold, InDictionary):
    nodeFull = InDictionary[node]
    nextNodeFull = InDictionary[nextNode]
    nodeRP = nodeFull[1]
    nodeDR = nodeFull[0]
    nextNodeRP = nextNodeFull[1]
    nextNodeDR = nextNodeFull[0]
    #no reads condition shouldnt happen but safety first kids
    if (nodeRP >= readPercentageThreshold and nextNodeRP >= readPercentageThreshold and nextNodeDR == nodeDR and nodeDR != "No_Reads"):
        return True
    else:
        return False

def areNodesRepeated(combinedList):
    if len(combinedList) == len(set(combinedList)):
        return False
    else:
        return True

def mergeZTestSequences(listOfLists,meanCoverage,stdDev, ScoreValue):
    currentList = listOfLists
    returnList = []
    while (len(currentList) > 1):
        firstList = currentList[0]
        remainingLists = currentList[1:]
        newLists = []
        if (len(firstList) != 1 and not areNodesRepeated(firstList)):
            if (withinZValue(get_cov_from_spades_name(firstList[0]), meanCoverage, stdDev, ScoreValue) and withinZValue(get_cov_from_spades_name(firstList[-1]), meanCoverage, stdDev, ScoreValue)):
                for list in remainingLists:
                    if firstList[0] == list[-1] and firstList[-1] == list[0]:
                        pass
                    elif firstList[0] == list[-1]:
                        combinedList = list + firstList[1:]
                        if not areNodesRepeated(combinedList):
                            if len(combinedList) <= maxNodeStringLen:
                                newLists = newLists + [combinedList]
                            else:
                                returnList = returnList + [combinedList]
                    elif firstList[-1] == list[0]:
                        combinedList = firstList + list[1:]
                        if not areNodesRepeated(combinedList):
                            if len(combinedList) <= maxNodeStringLen:
                                newLists = newLists + [combinedList]
                            else:
                                returnList = returnList + [combinedList]
                if len(newLists) == 0:
                    returnList = returnList + [firstList]
            elif withinZValue(get_cov_from_spades_name(firstList[0]), meanCoverage, stdDev, ScoreValue):
                for list in remainingLists:
                    if firstList[0] == list[-1]:
                        combinedList = list + firstList[1:]
                        if not areNodesRepeated(combinedList):
                            if len(combinedList) <= maxNodeStringLen:
                                newLists = newLists + [combinedList]
                            else:
                                returnList = returnList + [combinedList]
                if len(newLists) == 0:
                    returnList = returnList + [firstList]
            elif withinZValue(get_cov_from_spades_name(firstList[-1]), meanCoverage, stdDev, ScoreValue):
                for list in remainingLists:
                    if firstList[-1] == list[0]:
                        combinedList = firstList + list[1:]
                        if not areNodesRepeated(combinedList):
                            if len(combinedList) <= maxNodeStringLen:
                                newLists = newLists + [combinedList]
                            else:
                                returnList = returnList + [combinedList]
                if len(newLists) == 0:
                    returnList = returnList + [firstList]
            currentList = newLists + remainingLists
        else:
            currentList = remainingLists
    firstList = currentList[0]
    if (len(firstList) != 1 and not areNodesRepeated(firstList)):
        returnList = returnList + currentList
    return returnList


def mergePercentileSequences(listOfLists, percentileValue):
    currentList = listOfLists
    returnList = []
    while (len(currentList) > 1):
        firstList = currentList[0]
        remainingLists = currentList[1:]
        newLists = []
        if (len(firstList) != 1 and not areNodesRepeated(firstList)):
            if (withinPercentile(firstList[0], percentileValue) and withinPercentile(firstList[-1], percentileValue)):
                for list in remainingLists:
                    if firstList[0] == list[-1] and firstList[-1] == list[0]:
                        pass
                    elif firstList[0] == list[-1]:
                        combinedList = list + firstList[1:]
                        if not areNodesRepeated(combinedList):
                            if len(combinedList) <= maxNodeStringLen:
                                newLists = newLists + [combinedList]
                            else:
                                returnList = returnList + [combinedList]
                    elif firstList[-1] == list[0]:
                        combinedList = firstList + list[1:]
                        if not areNodesRepeated(combinedList):
                            if len(combinedList) <= maxNodeStringLen:
                                newLists = newLists + [combinedList]
                            else:
                                returnList = returnList + [combinedList]
                if len(newLists) == 0:
                    returnList = returnList + [firstList]
            elif withinPercentile(firstList[0], percentileValue):
                for list in remainingLists:
                    if firstList[0] == list[-1]:
                        combinedList = list + firstList[1:]
                        if not areNodesRepeated(combinedList):
                            if len(combinedList) <= maxNodeStringLen:
                                newLists = newLists + [combinedList]
                            else:
                                returnList = returnList + [combinedList]
                if len(newLists) == 0:
                    returnList = returnList + [firstList]
            elif withinPercentile(firstList[-1], percentileValue):
                for list in remainingLists:
                    if firstList[-1] == list[0]:
                        combinedList = firstList + list[1:]
                        if not areNodesRepeated(combinedList):
                            if len(combinedList) <= maxNodeStringLen:
                                newLists = newLists + [combinedList]
                            else:
                                returnList = returnList + [combinedList]
                if len(newLists) == 0:
                    returnList = returnList + [firstList]
            currentList = newLists + remainingLists
        else:
            currentList = remainingLists
    firstList = currentList[0]
    if (len(firstList) != 1 and not areNodesRepeated(firstList)):
        returnList = returnList + currentList
    return returnList

def mergeReadPercentageSequences(listOfLists, readPercentageThreshold, InDictionary):
    currentList = listOfLists
    returnList = []
    while (len(currentList) > 1):
        firstList = currentList[0]
        remainingLists = currentList[1:]
        newLists = []
        if (len(firstList) != 1 and not areNodesRepeated(firstList)):
            if (withinReadPercentage(firstList[0], readPercentageThreshold, InDictionary) and withinReadPercentage(firstList[-1], readPercentageThreshold, InDictionary)):
                for list in remainingLists:
                    if firstList[0] == list[-1] and firstList[-1] == list[0]:
                        pass
                    elif firstList[0] == list[-1]:
                        combinedList = list + firstList[1:]
                        if not areNodesRepeated(combinedList):
                            if len(combinedList) <= maxNodeStringLen:
                                newLists = newLists + [combinedList]
                            else:
                                returnList = returnList + [combinedList]
                    elif firstList[-1] == list[0]:
                        combinedList = firstList + list[1:]
                        if not areNodesRepeated(combinedList):
                            if len(combinedList) <= maxNodeStringLen:
                                newLists = newLists + [combinedList]
                            else:
                                returnList = returnList + [combinedList]
                if len(newLists) == 0:
                    returnList = returnList + [firstList]
                #print('a')
            elif withinReadPercentage(firstList[0], readPercentageThreshold, InDictionary):
                for list in remainingLists:
                    if firstList[0] == list[-1]:
                        combinedList = list + firstList[1:]
                        if not areNodesRepeated(combinedList):
                            if len(combinedList) <= maxNodeStringLen:
                                newLists = newLists + [combinedList]
                            else:
                                returnList = returnList + [combinedList]
                if len(newLists) == 0:
                    returnList = returnList + [firstList]
                #print('b')
            elif withinReadPercentage(firstList[-1], readPercentageThreshold, InDictionary):
                for list in remainingLists:
                    if firstList[-1] == list[0]:
                        combinedList = firstList + list[1:]
                        if not areNodesRepeated(combinedList):
                            if len(combinedList) <= maxNodeStringLen:
                                newLists = newLists + [combinedList]
                            else:
                                returnList = returnList + [combinedList]
                if len(newLists) == 0:
                    returnList = returnList + [firstList]
                #print('c')
            currentList = newLists + remainingLists
        else:
            currentList = remainingLists
        #print(len(currentList))
    firstList = currentList[0]
    if (len(firstList) != 1 and not areNodesRepeated(firstList)):
        returnList = returnList + currentList
    return returnList

###Thresholding Stage Functions ORGANIC_P

def thresholdSubsequencesAfterSplit_BothEnds_Percentile(listOfLists, percentileValue):
    currentList = listOfLists
    thresholdedList = []
    while (len(currentList) > 0): #####
        listOfNodes = currentList[0]
        remainingList = currentList[1:]
        innerThresholdedList = []
        lastNode = []
        thisNode = []
        while (len(listOfNodes) > 0):
            if lastNode == []:
                if (len(listOfNodes) == 1):
                    remainingNodes = []
                    listOfNodes = remainingNodes
                elif (len(listOfNodes) == 2):
                    remainingNodes = []
                    listOfNodes = remainingNodes
                else:   #the start of a new set of sequences
                    #print("The Start of a new set of nodes.")
                    lastNode = listOfNodes[0]
                    thisNode = listOfNodes[1]
                    if withinPercentileBackwardFail(lastNode, thisNode, percentileValue):
                        #print("The Start, BothEnds, combine, first node fails on purpose")
                        innerThresholdedList = [lastNode, thisNode]
                    #print("The Start, BothEnds, fail but thats okay")
                    remainingNodes = listOfNodes[2:]
                    listOfNodes = remainingNodes
                    lastNode = thisNode
            else:
                thisNode = listOfNodes[0]
                #only 1 node remaining
                if withinPercentileBothTrue(lastNode, thisNode, percentileValue) and (len(innerThresholdedList) > 1): #greater than 1 because the first has to fail
                    innerThresholdedList += [thisNode]
                    #print("Middle Node, BothEnds, Combine, with future")
                elif withinPercentileBackwardFail(lastNode, thisNode, percentileValue) and (len(innerThresholdedList) == 0): #nothing in the innerlist, need to have first fail
                    innerThresholdedList = [lastNode,thisNode]
                    #print("Middle Node, BothEnds, Start")
                elif withinPercentileForwardFail(lastNode, thisNode, percentileValue) and (len(innerThresholdedList) > 1):
                    #print("middle node, BothEnds, terminate")
                    innerThresholdedList += [thisNode]
                    thresholdedList += [innerThresholdedList]
                    innerThresholdedList = []
                else:
                    innerThresholdedList = []
                    #print('othercase, BothEnds, reset')
                if not len(listOfNodes[0]) == 1:
                    remainingNodes = listOfNodes[1:]
                else:
                    remainingNodes = []
                listOfNodes = remainingNodes
                lastNode = thisNode
        currentList = remainingList
    return thresholdedList


#need to fix, this can start in the middle
def thresholdSubsequencesAfterSplit_FrontEnd_Percentile(listOfLists, percentileValue):
    currentList = listOfLists
    thresholdedList = []
    while (len(currentList) > 0): #####
        listOfNodes = currentList[0]
        remainingList = currentList[1:]
        innerThresholdedList = []
        lastNode = []
        thisNode = []
        while (len(listOfNodes) > 0):
            if lastNode == []:
                if (len(listOfNodes) == 1):
                    remainingNodes = []
                    listOfNodes = remainingNodes
                elif (len(listOfNodes) == 2):
                    lastNode = listOfNodes[0]
                    thisNode = listOfNodes[1]
                    if (withinPercentileBackwardFail(lastNode, thisNode, percentileValue)):
                        innerThresholdedList = [lastNode, thisNode]
                        thresholdedList += [innerThresholdedList]
                        #print("FrontEnd, Two Node succeed")
                        innerThresholdedList = []
                    remainingNodes = []
                    listOfNodes = remainingNodes
                else:   #the start of a new set of sequences
                    #print("The Start of a new set of nodes.")
                    lastNode = listOfNodes[0]
                    thisNode = listOfNodes[1]
                    if withinPercentileBackwardFail(lastNode, thisNode, percentileValue):
                        #print("The Start, BothEnds, combine, first node fails on purpose")
                        innerThresholdedList = [lastNode, thisNode]
                        remainingNodes = listOfNodes[2:]
                        #print("FrontEnd, Start at start, first fails")
                    else:
                        remainingNodes = []
                        #print("FrontEnd, Start othercase")
                    listOfNodes = remainingNodes
                    lastNode = thisNode
            else:
                thisNode = listOfNodes[0]
                if (len(listOfNodes) == 1) and withinPercentileBothTrue(lastNode, thisNode, percentileValue) and (len(innerThresholdedList) > 1):
                    #last in list and fulfills criteria
                    innerThresholdedList += [thisNode]
                    thresholdedList += [innerThresholdedList]
                    #print('FrontEnd, Last Node, add to string, double True')
                    innerThresholdedList = []
                    remainingNodes = []
                elif((len(listOfNodes) == 1) and withinPercentileBackwardFail(lastNode, thisNode, percentileValue)):
                    innerThresholdedList = [lastNode, thisNode]
                    innerThresholdedList += [thisNode]
                    thresholdedList += [innerThresholdedList]
                    remainingNodes = []
                    #print('FrontEnd, Last 2 nodes, BackwardFail True, add 2')
                elif withinPercentileBackwardFail(lastNode, thisNode, percentileValue):
                    innerThresholdedList = [lastNode, thisNode]
                    remainingNodes = listOfNodes[1:]
                    #print('FrontEnd, Backwardfail true in middle, start')
                elif withinPercentileBothTrue(lastNode, thisNode, percentileValue) and (len(innerThresholdedList) > 1): #greater than 1 because the first has to fail
                    innerThresholdedList += [thisNode]
                    remainingNodes = listOfNodes[1:]
                    #print("FrontEnd, Middle Node, DoubleTrue, Backwardfail")
                else:
                    if len(listOfNodes) != 1:
                        remainingNodes = listOfNodes[1:]
                    else:
                        remainingNodes = []
                    #two fail nodes in a row or the len f-ed up
                listOfNodes = remainingNodes
                lastNode = thisNode
        currentList = remainingList
    return thresholdedList


#this cant start in the middle
def thresholdSubsequencesAfterSplit_BackEnd_Percentile(listOfLists, percentileValue):
    currentList = listOfLists
    thresholdedList = []
    while (len(currentList) > 0): #####
        listOfNodes = currentList[0]
        remainingList = currentList[1:]
        innerThresholdedList = []
        lastNode = []
        thisNode = []
        while (len(listOfNodes) > 0):
            if lastNode == []:
                if (len(listOfNodes) == 1):
                    remainingNodes = []
                    listOfNodes = remainingNodes
                elif (len(listOfNodes) == 2):
                    lastNode = listOfNodes[0]
                    thisNode = listOfNodes[1]
                    if (withinPercentileForwardFail(lastNode, thisNode, percentileValue)):
                        innerThresholdedList = [lastNode, thisNode]
                        thresholdedList += [innerThresholdedList]
                        innerThresholdedList = []
                        #print("only 2 nodes, ForwardFail, Succeed")
                    remainingNodes = []
                    listOfNodes = remainingNodes
                else:   #the start of a new set of sequences
                    #print("The Start of a new set of nodes.")
                    lastNode = listOfNodes[0]
                    thisNode = listOfNodes[1]
                    if withinPercentileBothTrue(lastNode, thisNode, percentileValue):
                        #print("ForwardFail, Start, Double True")
                        innerThresholdedList = [lastNode, thisNode]
                        remainingNodes = listOfNodes[2:]
                    else:
                        #print("ForwardFail, Start Fail")
                        remainingNodes = []
                        #print("The Start, first doesn't fail, don't need to do anything in BothEnds")
                    listOfNodes = remainingNodes
                    lastNode = thisNode
            else:
                thisNode = listOfNodes[0]
                if (len(listOfNodes) == 1) and withinPercentileForwardFail(lastNode, thisNode, percentileValue) and (len(innerThresholdedList) > 1):
                    #print("ForwardFail,,Last Node, True Fail, Success")
                    #last in list and fulfills criteria
                    innerThresholdedList += [thisNode]
                    thresholdedList += [innerThresholdedList]
                    innerThresholdedList = []
                    remainingNodes = []
                elif withinPercentileBothTrue(lastNode, thisNode, percentileValue) and (len(innerThresholdedList) > 1): #greater than 1 because the first has to fail
                    #print("ForwardFail, Middle Nodes, Add True True")
                    innerThresholdedList += [thisNode]
                    remainingNodes = listOfNodes[1:]
                    #print("Middle Node, Combine, with future")
                else:
                    #print("Failed after start. Next String of nodes")
                    remainingNodes = []
                    #two fail nodes in a row or the len f-ed up
                listOfNodes = remainingNodes
                lastNode = thisNode
        currentList = remainingList
    return thresholdedList


###Thresholding Stage Functions ORGANIC_Z

#for both ends need at least 3
def thresholdSubsequencesAfterSplit_BothEnds_ZTest(listOfLists, meanCoverage, stdDev, ScoreValue):
    currentList = listOfLists
    thresholdedList = []
    while (len(currentList) > 0): #####
        listOfNodes = currentList[0]
        remainingList = currentList[1:]
        innerThresholdedList = []
        lastNode = []
        thisNode = []
        while (len(listOfNodes) > 0):
            if lastNode == []:
                if (len(listOfNodes) == 1):
                    remainingNodes = []
                    listOfNodes = remainingNodes
                elif (len(listOfNodes) == 2):
                    remainingNodes = []
                    listOfNodes = remainingNodes
                else:   #the start of a new set of sequences
                    #print("The Start of a new set of nodes.")
                    lastNode = listOfNodes[0]
                    thisNode = listOfNodes[1]
                    if withinZValueBackwardFail(get_cov_from_spades_name(lastNode), get_cov_from_spades_name(thisNode), meanCoverage, stdDev, ScoreValue):
                        #print("The Start, BothEnds, combine, first node fails on purpose")
                        innerThresholdedList = [lastNode, thisNode]
                    #print("The Start, BothEnds, fail but thats okay")
                    remainingNodes = listOfNodes[2:]
                    listOfNodes = remainingNodes
                    lastNode = thisNode
            else:
                thisNode = listOfNodes[0]
                #only 1 node remaining
                if withinZValueBothTrue(get_cov_from_spades_name(lastNode), get_cov_from_spades_name(thisNode), meanCoverage, stdDev, ScoreValue) and (len(innerThresholdedList) > 1): #greater than 1 because the first has to fail
                    innerThresholdedList += [thisNode]
                    #print("Middle Node, BothEnds, Combine, with future")
                elif withinZValueBackwardFail(get_cov_from_spades_name(lastNode), get_cov_from_spades_name(thisNode), meanCoverage, stdDev, ScoreValue) and (len(innerThresholdedList) == 0): #nothing in the innerlist, need to have first fail
                    innerThresholdedList = [lastNode,thisNode]
                    #print("Middle Node, BothEnds, Start")
                elif withinZValueForwardFail(get_cov_from_spades_name(lastNode), get_cov_from_spades_name(thisNode), meanCoverage, stdDev, ScoreValue) and (len(innerThresholdedList) > 1):
                    #print("middle node, BothEnds, terminate")
                    innerThresholdedList += [thisNode]
                    thresholdedList += [innerThresholdedList]
                    innerThresholdedList = []
                else:
                    innerThresholdedList = []
                    #print('othercase, BothEnds, reset')
                if not len(listOfNodes[0]) == 1:
                    remainingNodes = listOfNodes[1:]
                else:
                    remainingNodes = []
                listOfNodes = remainingNodes
                lastNode = thisNode
        currentList = remainingList
    return thresholdedList


#need to fix, this can start in the middle
def thresholdSubsequencesAfterSplit_FrontEnd_ZTest(listOfLists, meanCoverage, stdDev, ScoreValue):
    currentList = listOfLists
    thresholdedList = []
    while (len(currentList) > 0): #####
        listOfNodes = currentList[0]
        remainingList = currentList[1:]
        innerThresholdedList = []
        lastNode = []
        thisNode = []
        while (len(listOfNodes) > 0):
            if lastNode == []:
                if (len(listOfNodes) == 1):
                    remainingNodes = []
                    listOfNodes = remainingNodes
                elif (len(listOfNodes) == 2):
                    lastNode = listOfNodes[0]
                    thisNode = listOfNodes[1]
                    if (withinZValueBackwardFail(get_cov_from_spades_name(lastNode), get_cov_from_spades_name(thisNode), meanCoverage, stdDev, ScoreValue)):
                        innerThresholdedList = [lastNode, thisNode]
                        thresholdedList += [innerThresholdedList]
                        #print("FrontEnd, Two Node succeed")
                        innerThresholdedList = []
                    remainingNodes = []
                    listOfNodes = remainingNodes
                else:   #the start of a new set of sequences
                    #print("The Start of a new set of nodes.")
                    lastNode = listOfNodes[0]
                    thisNode = listOfNodes[1]
                    if withinZValueBackwardFail(get_cov_from_spades_name(lastNode), get_cov_from_spades_name(thisNode), meanCoverage, stdDev, ScoreValue):
                        #print("The Start, BothEnds, combine, first node fails on purpose")
                        innerThresholdedList = [lastNode, thisNode]
                        remainingNodes = listOfNodes[2:]
                        #print("FrontEnd, Start at start, first fails")
                    else:
                        remainingNodes = []
                        #print("FrontEnd, Start othercase")
                    listOfNodes = remainingNodes
                    lastNode = thisNode
            else:
                thisNode = listOfNodes[0]
                if (len(listOfNodes) == 1) and withinZValueBothTrue(get_cov_from_spades_name(lastNode), get_cov_from_spades_name(thisNode), meanCoverage, stdDev, ScoreValue) and (len(innerThresholdedList) > 1):
                    #last in list and fulfills criteria
                    innerThresholdedList += [thisNode]
                    thresholdedList += [innerThresholdedList]
                    #print('FrontEnd, Last Node, add to string, double True')
                    innerThresholdedList = []
                    remainingNodes = []
                elif((len(listOfNodes) == 1) and withinZValueBackwardFail(get_cov_from_spades_name(lastNode), get_cov_from_spades_name(thisNode), meanCoverage, stdDev, ScoreValue)):
                    innerThresholdedList = [lastNode, thisNode]
                    innerThresholdedList += [thisNode]
                    thresholdedList += [innerThresholdedList]
                    remainingNodes = []
                    #print('FrontEnd, Last 2 nodes, BackwardFail True, add 2')
                elif withinZValueBackwardFail(get_cov_from_spades_name(lastNode), get_cov_from_spades_name(thisNode), meanCoverage, stdDev, ScoreValue):
                    innerThresholdedList = [lastNode, thisNode]
                    remainingNodes = listOfNodes[1:]
                    #print('FrontEnd, Backwardfail true in middle, start')
                elif withinZValueBothTrue(get_cov_from_spades_name(lastNode), get_cov_from_spades_name(thisNode), meanCoverage, stdDev, ScoreValue) and (len(innerThresholdedList) > 1): #greater than 1 because the first has to fail
                    innerThresholdedList += [thisNode]
                    remainingNodes = listOfNodes[1:]
                    #print("FrontEnd, Middle Node, DoubleTrue, Backwardfail")
                else:
                    if len(listOfNodes) != 1:
                        remainingNodes = listOfNodes[1:]
                    else:
                        remainingNodes = []
                    #two fail nodes in a row or the len f-ed up
                listOfNodes = remainingNodes
                lastNode = thisNode
        currentList = remainingList
    return thresholdedList


#this cant start in the middle
def thresholdSubsequencesAfterSplit_BackEnd_ZTest(listOfLists, meanCoverage, stdDev, ScoreValue):
    currentList = listOfLists
    thresholdedList = []
    while (len(currentList) > 0): #####
        listOfNodes = currentList[0]
        remainingList = currentList[1:]
        innerThresholdedList = []
        lastNode = []
        thisNode = []
        while (len(listOfNodes) > 0):
            if lastNode == []:
                if (len(listOfNodes) == 1):
                    remainingNodes = []
                    listOfNodes = remainingNodes
                elif (len(listOfNodes) == 2):
                    lastNode = listOfNodes[0]
                    thisNode = listOfNodes[1]
                    if (withinZValueForwardFail(get_cov_from_spades_name(lastNode), get_cov_from_spades_name(thisNode), meanCoverage, stdDev, ScoreValue)):
                        innerThresholdedList = [lastNode, thisNode]
                        thresholdedList += [innerThresholdedList]
                        innerThresholdedList = []
                        #print("only 2 nodes, ForwardFail, Succeed")
                    remainingNodes = []
                    listOfNodes = remainingNodes
                else:   #the start of a new set of sequences
                    #print("The Start of a new set of nodes.")
                    lastNode = listOfNodes[0]
                    thisNode = listOfNodes[1]
                    if withinZValueBothTrue(get_cov_from_spades_name(lastNode), get_cov_from_spades_name(thisNode), meanCoverage, stdDev, ScoreValue):
                        #print("ForwardFail, Start, Double True")
                        innerThresholdedList = [lastNode, thisNode]
                        remainingNodes = listOfNodes[2:]
                    else:
                        #print("ForwardFail, Start Fail")
                        remainingNodes = []
                        #print("The Start, first doesn't fail, don't need to do anything in BothEnds")
                    listOfNodes = remainingNodes
                    lastNode = thisNode
            else:
                thisNode = listOfNodes[0]
                if (len(listOfNodes) == 1) and withinZValueForwardFail(get_cov_from_spades_name(lastNode), get_cov_from_spades_name(thisNode), meanCoverage, stdDev, ScoreValue) and (len(innerThresholdedList) > 1):
                    #print("ForwardFail,,Last Node, True Fail, Success")
                    #last in list and fulfills criteria
                    innerThresholdedList += [thisNode]
                    thresholdedList += [innerThresholdedList]
                    innerThresholdedList = []
                    remainingNodes = []
                elif withinZValueBothTrue(get_cov_from_spades_name(lastNode), get_cov_from_spades_name(thisNode), meanCoverage, stdDev, ScoreValue) and (len(innerThresholdedList) > 1): #greater than 1 because the first has to fail
                    #print("ForwardFail, Middle Nodes, Add True True")
                    innerThresholdedList += [thisNode]
                    remainingNodes = listOfNodes[1:]
                    #print("Middle Node, Combine, with future")
                else:
                    #print("Failed after start. Next String of nodes")
                    remainingNodes = []
                    #two fail nodes in a row or the len f-ed up
                listOfNodes = remainingNodes
                lastNode = thisNode
        currentList = remainingList
    return thresholdedList

###Thresholding Stage Functions SYNTH

#for both ends need at least 3
def thresholdSubsequencesAfterSplit_BothEnds_ReadPercentage(listOfLists, readPercentageThreshold, InDictionary):
    currentList = listOfLists
    thresholdedList = []
    while (len(currentList) > 0): #####
        listOfNodes = currentList[0]
        remainingList = currentList[1:]
        innerThresholdedList = []
        lastNode = []
        thisNode = []
        while (len(listOfNodes) > 0):
            if lastNode == []:
                if (len(listOfNodes) == 1):
                    remainingNodes = []
                    listOfNodes = remainingNodes
                elif (len(listOfNodes) == 2):
                    remainingNodes = []
                    listOfNodes = remainingNodes
                else:   #the start of a new set of sequences
                    #print("The Start of a new set of nodes.")
                    lastNode = listOfNodes[0]
                    thisNode = listOfNodes[1]
                    if withinReadPercentageBackwardFail(lastNode, thisNode, readPercentageThreshold, InDictionary):
                        #print("The Start, BothEnds, combine, first node fails on purpose")
                        innerThresholdedList = [lastNode, thisNode]
                    #print("The Start, BothEnds, fail but thats okay")
                    remainingNodes = listOfNodes[2:]
                    listOfNodes = remainingNodes
                    lastNode = thisNode
            else:
                thisNode = listOfNodes[0]
                #only 1 node remaining
                if withinReadPercentageBothTrue(lastNode, thisNode, readPercentageThreshold, InDictionary) and (len(innerThresholdedList) > 1): #greater than 1 because the first has to fail
                    innerThresholdedList += [thisNode]
                    #print("Middle Node, BothEnds, Combine, with future")
                elif withinReadPercentageBackwardFail(lastNode, thisNode, readPercentageThreshold, InDictionary) and (len(innerThresholdedList) == 0): #nothing in the innerlist, need to have first fail
                    innerThresholdedList = [lastNode,thisNode]
                    #print("Middle Node, BothEnds, Start")
                elif withinReadPercentageForwardFail(lastNode, thisNode, readPercentageThreshold, InDictionary) and (len(innerThresholdedList) > 1):
                    #print("middle node, BothEnds, terminate")
                    innerThresholdedList += [thisNode]
                    thresholdedList += [innerThresholdedList]
                    innerThresholdedList = []
                else:
                    innerThresholdedList = []
                    #print('othercase, BothEnds, reset')
                if not len(listOfNodes[0]) == 1:
                    remainingNodes = listOfNodes[1:]
                else:
                    remainingNodes = []
                listOfNodes = remainingNodes
                lastNode = thisNode
        currentList = remainingList
    return thresholdedList


#need to fix, this can start in the middle
def thresholdSubsequencesAfterSplit_FrontEnd_ReadPercentage(listOfLists, readPercentageThreshold, InDictionary):
    currentList = listOfLists
    thresholdedList = []
    while (len(currentList) > 0): #####
        listOfNodes = currentList[0]
        remainingList = currentList[1:]
        innerThresholdedList = []
        lastNode = []
        thisNode = []
        while (len(listOfNodes) > 0):
            if lastNode == []:
                if (len(listOfNodes) == 1):
                    remainingNodes = []
                    listOfNodes = remainingNodes
                elif (len(listOfNodes) == 2):
                    lastNode = listOfNodes[0]
                    thisNode = listOfNodes[1]
                    if (withinReadPercentageBackwardFail(lastNode, thisNode, readPercentageThreshold, InDictionary)):
                        innerThresholdedList = [lastNode, thisNode]
                        thresholdedList += [innerThresholdedList]
                        #print("FrontEnd, Two Node succeed")
                        innerThresholdedList = []
                    remainingNodes = []
                    listOfNodes = remainingNodes
                else:   #the start of a new set of sequences
                    #print("The Start of a new set of nodes.")
                    lastNode = listOfNodes[0]
                    thisNode = listOfNodes[1]
                    if withinReadPercentageBackwardFail(lastNode, thisNode, readPercentageThreshold, InDictionary):
                        #print("The Start, BothEnds, combine, first node fails on purpose")
                        innerThresholdedList = [lastNode, thisNode]
                        remainingNodes = listOfNodes[2:]
                        #print("FrontEnd, Start at start, first fails")
                    else:
                        remainingNodes = []
                        #print("FrontEnd, Start othercase")
                    listOfNodes = remainingNodes
                    lastNode = thisNode
            else:
                thisNode = listOfNodes[0]
                if (len(listOfNodes) == 1) and withinReadPercentageBothTrue(lastNode, thisNode, readPercentageThreshold, InDictionary) and (len(innerThresholdedList) > 1):
                    #last in list and fulfills criteria
                    innerThresholdedList += [thisNode]
                    thresholdedList += [innerThresholdedList]
                    #print('FrontEnd, Last Node, add to string, double True')
                    innerThresholdedList = []
                    remainingNodes = []
                elif((len(listOfNodes) == 1) and withinReadPercentageBackwardFail(lastNode, thisNode, readPercentageThreshold, InDictionary)):
                    innerThresholdedList = [lastNode, thisNode]
                    innerThresholdedList += [thisNode]
                    thresholdedList += [innerThresholdedList]
                    remainingNodes = []
                    #print('FrontEnd, Last 2 nodes, BackwardFail True, add 2')
                elif withinReadPercentageBackwardFail(lastNode, thisNode, readPercentageThreshold, InDictionary):
                    innerThresholdedList = [lastNode, thisNode]
                    remainingNodes = listOfNodes[1:]
                    #print('FrontEnd, Backwardfail true in middle, start')
                elif withinReadPercentageBothTrue(lastNode, thisNode, readPercentageThreshold, InDictionary) and (len(innerThresholdedList) > 1): #greater than 1 because the first has to fail
                    innerThresholdedList += [thisNode]
                    remainingNodes = listOfNodes[1:]
                    #print("FrontEnd, Middle Node, DoubleTrue, Backwardfail")
                else:
                    if len(listOfNodes) != 1:
                        remainingNodes = listOfNodes[1:]
                    else:
                        remainingNodes = []
                    #two fail nodes in a row or the len f-ed up
                listOfNodes = remainingNodes
                lastNode = thisNode
        currentList = remainingList
    return thresholdedList


#this cant start in the middle
def thresholdSubsequencesAfterSplit_BackEnd_ReadPercentage(listOfLists, readPercentageThreshold, InDictionary):
    currentList = listOfLists
    thresholdedList = []
    while (len(currentList) > 0): #####
        listOfNodes = currentList[0]
        remainingList = currentList[1:]
        innerThresholdedList = []
        lastNode = []
        thisNode = []
        while (len(listOfNodes) > 0):
            if lastNode == []:
                if (len(listOfNodes) == 1):
                    remainingNodes = []
                    listOfNodes = remainingNodes
                elif (len(listOfNodes) == 2):
                    lastNode = listOfNodes[0]
                    thisNode = listOfNodes[1]
                    if (withinReadPercentageForwardFail(lastNode, thisNode, readPercentageThreshold, InDictionary)):
                        innerThresholdedList = [lastNode, thisNode]
                        thresholdedList += [innerThresholdedList]
                        innerThresholdedList = []
                        #print("only 2 nodes, ForwardFail, Succeed")
                    remainingNodes = []
                    listOfNodes = remainingNodes
                else:   #the start of a new set of sequences
                    #print("The Start of a new set of nodes.")
                    lastNode = listOfNodes[0]
                    thisNode = listOfNodes[1]
                    if withinReadPercentageBothTrue(lastNode, thisNode, readPercentageThreshold, InDictionary):
                        #print("ForwardFail, Start, Double True")
                        innerThresholdedList = [lastNode, thisNode]
                        remainingNodes = listOfNodes[2:]
                    else:
                        #print("ForwardFail, Start Fail")
                        remainingNodes = []
                        #print("The Start, first doesn't fail, don't need to do anything in BothEnds")
                    listOfNodes = remainingNodes
                    lastNode = thisNode
            else:
                thisNode = listOfNodes[0]
                if (len(listOfNodes) == 1) and withinReadPercentageForwardFail(lastNode, thisNode, readPercentageThreshold, InDictionary) and (len(innerThresholdedList) > 1):
                    #print("ForwardFail,,Last Node, True Fail, Success")
                    #last in list and fulfills criteria
                    innerThresholdedList += [thisNode]
                    thresholdedList += [innerThresholdedList]
                    innerThresholdedList = []
                    remainingNodes = []
                elif withinReadPercentageBothTrue(lastNode, thisNode, readPercentageThreshold, InDictionary) and (len(innerThresholdedList) > 1): #greater than 1 because the first has to fail
                    #print("ForwardFail, Middle Nodes, Add True True")
                    innerThresholdedList += [thisNode]
                    remainingNodes = listOfNodes[1:]
                    #print("Middle Node, Combine, with future")
                else:
                    #print("Failed after start. Next String of nodes")
                    remainingNodes = []
                    #two fail nodes in a row or the len f-ed up
                listOfNodes = remainingNodes
                lastNode = thisNode
        currentList = remainingList
    return thresholdedList

def testIfOrderedSubset(subsetList, fullList):
    if (len(fullList) - len(subsetList) + 1) <= 0:
        return False
    for iterFull in range(len(fullList) - len(subsetList) + 1):
        if all(subsetList[iterSub] == fullList[iterFull + iterSub] for iterSub in range(len(subsetList))):
            return True
        else:
            return False

def writeSequenceToFile(filename,sequenceList,sequenceDictionary):
    f = open(filename,"w+")
    #print(len(sequenceList))
    for list in sequenceList:
        name = ">"
        for node in list:
            name += node+"_._"
        name = name[0:-3]
        f.write(name+"\n")
        sequence = ""
        for node in list:
            subSequence = sequenceDictionary[node]
            sequence += subSequence
        f.write(sequence+"\n")
    f.close()

def writeSequenceToFileSynth(filename,sequenceList,sequenceDictionary,InDict):
    f = open(filename,"w+")
    #print(len(sequenceList))
    for list in sequenceList:
        name = ">"
        name += InDict[list[0]][0]+"___"
        for node in list:
            name += node+"_._"
        name = name[0:-3]
        f.write(name+"\n")
        sequence = ""
        for node in list:
            subSequence = sequenceDictionary[node]
            sequence += subSequence
        f.write(sequence+"\n")
    f.close()