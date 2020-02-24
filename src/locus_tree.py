from .species_tree import *
from .tree_table import *
import pprint

class LocusTree(SpeciesTree):

    def __init__(self, randomState):
        super().__init__(randomState=randomState)

    def initialize(self, nodes, skbioTree):
        self.treeTable = TreeTable()
        self.treeTable.createFromEntries(
            entries=nodes, skbioTree=skbioTree)

    # bounded coalescent for the locus tree model
    # keep trying until every gene merged in time
    def boundedCoalescent(self, distanceAboveRoot):
        coalescentProcess, genesIntoRoot = self.coalescent(distanceAboveRoot)
        if len(genesIntoRoot) == 1:
            return coalescentProcess
        else:
            return self.boundedCoalescent(distanceAboveRoot)

    # incomplete coalescent for IxDTL
    def incompleteCoalescent(self, distanceAboveRoot):
        fullCoalescentProcess, genesIntoRoot = self.coalescent(distanceAboveRoot)
        chosenGene = self.randomState.choice(genesIntoRoot)
        selectedCoalescentProcess = self.__selectCoalescentProcess(
            fullCoalescentProcess, chosenGene)
        return fullCoalescentProcess, selectedCoalescentProcess, chosenGene
    
    # seclect one haplotype tree from multiple trees
    # generated by the incomplete coalescent
    def __selectCoalescentProcess(self, coalescentProcess, chosenGene):
        selectedCoalescentProcess = defaultdict(list)
        for speciesNodeId, mergingSets in coalescentProcess.items():
            for mergingSet in mergingSets:
                distance = mergingSet['distance']
                fromSet = []
                toSet = []
                for clade in mergingSet['fromSet']:
                    if self._starInSet(target=clade, clade=chosenGene):
                        fromSet.append(clade)
                if mergingSet['toSet']:
                    for clade in mergingSet['toSet']:
                        if self._starInSet(target=clade, clade=chosenGene):
                            toSet.append(clade)
                    if toSet:
                        selectedCoalescentProcess[speciesNodeId].append({
                            'fromSet': fromSet, 
                            'toSet': toSet,
                            'distance': distance
                        })
                else:
                    if fromSet:
                        selectedCoalescentProcess[speciesNodeId].append({
                            'fromSet': fromSet, 
                            'toSet': None,
                            'distance': distance
                        })
        return selectedCoalescentProcess

    """""""""""""""""""""""""""""""""""""""""""""""""""""
    funtions for modelling (partly) linked duplictions
    """""""""""""""""""""""""""""""""""""""""""""""""""""
    
    """
    
    """
    def __filteredCoalescentWithRecombinationProcess(self, coalescentProcess, cladeSetIntoRoot):
        filteredProcess = defaultdict(list)
        filteredClades = []
        for speciesNodeId, mergingSets in coalescentProcess.items():
            distance = 0
            for mergingSet in mergingSets:
                fromSet = []
                toSet = []
                for element in mergingSet['fromSet']:
                    newClade = element.split('*')[-1]
                    if newClade:
                        newClade = self.__starReplace(newClade)
                        fromSet.append(newClade)
                if mergingSet['toSet']:
                    for element in mergingSet['toSet']:
                        newClade = element.split('*')[-1]
                        if newClade:
                            newClade = self.__starReplace(newClade)
                            toSet.append(newClade)
                else:
                    toSet = None
                if fromSet != toSet:
                    filteredProcess[speciesNodeId].append({
                    'fromSet': fromSet,
                    'toSet': toSet,
                    'distance': mergingSet['distance'] + distance
                    })
                    distance = 0
                else:
                    distance += mergingSet['distance']

        for element in cladeSetIntoRoot:
            newClade = element.split('*')[-1]
            if newClade:
                newClade = self.__starReplace(newClade)
                filteredClades.append(newClade)
        return filteredProcess, filteredClades

    """
    
    """
    def coalescentWithRecombination(self, copiedHaplotypeTree, copiedRootGene, distanceAboveRoot):
        nodes = self.getNodes()
        root = self.getRoot()
        coalescentProcess = defaultdict(list)
        cladeSetIntoRoot = None

        # leaves of the given species tree
        oldLeaves = [node.id for node in nodes if not node.children]

        # leaves set will be updated in the loop
        newLeaves = []

        # set of extant species that an ancestral gene will eventually be fixed in;
        # cladeSet[node.id] = set of the genes at node.id
        fromSets = {}
        toSets = {}
        coalSets = {}
        recomSets = {}
        
        # avoid doing repeated coalescence: 
        # a node will be lablled after finishing the coalesencent 
        labelled = {}

        # initialization: 
        # every node is labelled false
        # cladeSet[leafId] = 'leafId*' ('*' as seperater)
        for node in nodes:
            labelled[node.id] = False
            coalSets[node.id] = []
            toSets[node.id] = []
            fromSets[node.id] = \
                [str(node.id) + '*' + str(node.id) + '#'] if not node.children else []
        recomSets = fromSets.copy()
        
        while True:
            for leaf in oldLeaves:
                # coalescent finished
                if leaf == root.id:
                    cladeSetIntoRoot, _, _ = self.__speciesBranchRecurse(nodeId=root.id, 
                        branchLength=self.getNodeById(root.id).distanceToParent, 
                        coalescentProcess=coalescentProcess, copiedHaplotypeTree=copiedHaplotypeTree, 
                        fromSet=fromSets[root.id], coalSet=coalSets[root.id], recomSet=coalSets[root.id],
                        distanceAboveRoot=distanceAboveRoot)
                    break
                else:
                    if labelled[leaf]:
                        continue
                    labelled[leaf] = True

                    parent = self.getNodeById(leaf).parent
                    children = self.getNodeById(parent).children
                    # if the leaf has been labelled, skip
                    
                    # first make sure there are genes coming out of both children
                    # then do coalescent within each child branch
                    # cladeSet will be updated from the gene comming into the branch
                    # to genes comming out of the branch
                    if (len(fromSets[children[0]]) != 0 
                        and len(fromSets[children[1]]) != 0):
                        toSets[children[0]], coalSets[children[0]], recomSets[children[0]] = \
                            self.__speciesBranchRecurse(nodeId=children[0], 
                            branchLength=self.getNodeById(children[0]).distanceToParent, 
                            coalescentProcess=coalescentProcess, copiedHaplotypeTree=copiedHaplotypeTree, 
                            fromSet=fromSets[children[0]], coalSet=coalSets[children[0]], recomSet=recomSets[children[0]])
                        labelled[children[0]] = True
                        toSets[children[1]], coalSets[children[1]], recomSets[children[1]] = \
                            self.__speciesBranchRecurse(nodeId=children[1], 
                            branchLength=self.getNodeById(children[1]).distanceToParent, 
                            coalescentProcess=coalescentProcess, copiedHaplotypeTree=copiedHaplotypeTree, 
                            fromSet=fromSets[children[1]], coalSet=coalSets[children[1]], recomSet=recomSets[children[1]]) 
                        labelled[children[1]] = True
                        
                        # update cladeSet[parent] as the
                        # union of the cladeSet of its children 
                        fromSets[parent] = toSets[children[0]] + toSets[children[1]]
                        recomSets[parent] = recomSets[children[0]] + recomSets[children[1]]
                        coalSets[parent] = coalSets[children[0]]+ coalSets[children[1]]
                        

                        # if the parent is in newLeaves, do not add in any children of the parent
                        if len(newLeaves) > 0:
                            newLeaves = [e for e in newLeaves if e !=
                                         children[0] and e != children[1]]
                        newLeaves.append(parent)
                    else:
                        # updating leaves set
                        newLeaves.append(leaf)
            # coalesecent finished
            if leaf == root.id:
                break

            # delete repeated items in newLeaves 
            tempNewLeaves = [] 
            for newLeaf in newLeaves:
                if newLeaf not in tempNewLeaves:
                    tempNewLeaves.append(newLeaf)

            # re-initialization for the next recursion
            # oldLeaves <- newLeaves
            # label <- false
            oldLeaves = tempNewLeaves.copy()
            newLeaves = []
            labelled = {}
            for node in nodes:
                labelled[node.id] = False

        ancestralClades = []
        nonAncestralClades = []
        for clade in cladeSetIntoRoot:
            if '#' in clade:
                if '*' not in clade:
                    nonAncestralClades.append(clade)
                else: 
                    for e in copiedRootGene:
                        if e in clade:
                            ancestralClades.append(clade)
                            break
                # else:
                #     ancestralClades.append(clade)
        fullClades = ancestralClades + nonAncestralClades
        # print('fullclades=',fullClades)
        # print('rootclades=', cladeSetIntoRoot)
        chosenGeneName = self.randomState.choice(cladeSetIntoRoot)
        if chosenGeneName not in fullClades:
            # discad the unobservable ancestral duplication
            # print('reach here')
            return None, None, None, None, True
        ancestral = False
        geneNodeName = None
        if chosenGeneName in ancestralClades:
            ancestral = True
            geneNodeName, chosenGeneName, _ = self.__getBipartition(chosenGeneName)
        chosenGeneName = self.__starReplace(chosenGeneName)

        filteredProcess, filteredClades = self.__filteredCoalescentWithRecombinationProcess(
            coalescentProcess, cladeSetIntoRoot)

        fullProcess = filteredProcess
        
        selectedProcess = self.__selectCoalescentProcess(filteredProcess, chosenGeneName)

        return fullProcess, selectedProcess, chosenGeneName, geneNodeName, ancestral

    """
    
    """
    def __geneBranchRecurse(self, nodeId, distance, distanceToAdd, fromSet, 
            coalescentProcess, coalSet, recomSet, copiedProcess,
            initial=True):
        unitRecombinationRate = self.recombinationRate
        unitCoalescentRate = self.coalescentRate
        

        if coalSet and len(coalSet) > 1 and unitCoalescentRate > 0:
            coalescentRate  = self.binom(len(coalSet), 2)*unitCoalescentRate
            coalDistance = self.randomState.exponential(scale=1.0 / coalescentRate)
        else:
            coalDistance = float('inf')
        if recomSet and len(recomSet) > 1 and unitRecombinationRate > 0:
            recombinationRate = self.binom(len(recomSet), 2)*unitRecombinationRate
            recomDistance = self.randomState.exponential(scale=1.0/recombinationRate)
        else:
            recomDistance = float('inf')
        if coalDistance < min(recomDistance, distance):
            chosenGene = self.randomState.choice(coalSet)
            coalSet = coalSet.remove(chosenGene)
            targets = fromSet.copy()
            targets.remove(chosenGene)
            target = self.randomState.choice(targets)
            couple = ''.join([chosenGene, target])
            starString, checkString, mergedString = self.__getBipartition(couple)
            recomSet = recomSet.append(mergedString)
            toSet = fromSet.copy()
            toSet.remove(chosenGene)
            toSet.remove(target)
            toSet.append(mergedString)
            if initial:
                coalDistance += distanceToAdd
            coalescentProcess[nodeId].append({
                'fromSet': fromSet,
                'toSet': toSet,
                'distance': coalDistance
            })
            return self.__geneBranchRecurse(nodeId=nodeId, distance=distance-coalDistance, distanceToAdd=0,
                fromSet=toSet, coalescentProcess=coalescentProcess,
                coalSet=coalSet, recomSet=recomSet, 
                copiedProcess=copiedProcess,
                initial=False)

        elif recomDistance < min(coalDistance, distance):
            chosenGene = self.randomState.choice(recomSet)
            recomSet = recomSet.remove(chosenGene)
            starString, checkString, mergedString = self.__getBipartition(chosenGene)
            coalSet = coalSet.append(checkString)
            toSet = fromSet.copy()
            toSet.remove(chosenGene)
            toSet.append(starString)
            toSet.append(checkString)
            if initial:
                recomDistance += distanceToAdd
            coalescentProcess[nodeId].append({
                'fromSet': fromSet,
                'toSet': toSet,
                'distance': recomDistance
            })
            return self.__geneBranchRecurse(nodeId=nodeId, distance=distance-recomDistance, distanceToAdd=0,
                fromSet=toSet, coalescentProcess=coalescentProcess,
                coalSet=coalSet, recomSet=recomSet,
                copiedProcess=copiedProcess,
                initial=False)
        else:
            if copiedProcess:
                toSet = fromSet.copy()
                if copiedProcess['toSet']:
                    if self.__getDifference(copiedProcess['fromSet'], copiedProcess['toSet']):
                        [coupleL, coupleR] = self.__getDifference(copiedProcess['fromSet'], copiedProcess['toSet'])
                        for e in fromSet:
                            if coupleL in e:
                                coupleL = e
                            elif coupleR in e:
                                coupleR = e
                            else:
                                continue
                        couple = ''.join([coupleL, coupleR])
                        starString, checkString, mergedString = self.__getBipartition(couple)
                        toSet = fromSet.copy()
                        toSet.remove(coupleL)
                        toSet.remove(coupleR)
                        toSet.append(mergedString)

                        recomSet = []
                        coalSet = []
                        for e in toSet:
                            if '*' in e and '#' in e:
                                recomSet.append(e)
                            elif '#' in e:
                                coalSet.append(e) 

                        if initial:
                            distance += distanceToAdd
                        coalescentProcess[nodeId].append({
                            'fromSet': fromSet,
                            'toSet': toSet,
                            'distance': distance
                        })
                        distanceToAdd = 0
                    else:
                        toSet = fromSet.copy()
                        coalescentProcess[nodeId].append({
                            'fromSet': fromSet,
                            'toSet': None,
                            'distance': distance
                        })
                        distanceToAdd += distance
                        coalSet = coalSet
                        recomSet = recomSet
                else:
                    toSet = fromSet.copy()
                    coalescentProcess[nodeId].append({
                        'fromSet': fromSet,
                        'toSet': None,
                        'distance': distance
                    })
                    distanceToAdd += distance
                    if not coalSet:
                        coalSet = []
                    else: 
                        coalSet = coalSet
                    if not recomSet:
                        recomSet = []
                    else:
                        recomSet = recomSet
            return toSet, coalSet, recomSet, distanceToAdd

    """
    utility functions for coalescentWithRecombination
    """
    def __starReplace(self, string):
        newString = ''
        for e in string:
            if e == '#':
                newString += '*'
            else:
                newString += e
        return newString

    def __starSorted(self, couple):
        """
        1*4* + 2*3* -> 1*4*2*3* -> 1*2*3*4*
        """
        string = ''
        for e in couple:
            string += e
        splited = string.split('*')[:-1]
        splited = sorted([int(e) for e in splited])
        return [str(e) + '*' for e in splited]

    def __checkSorted(self, couple):
        """
        1#4# + 2#3# -> 1#4#2#3# -> 1#2#3#4#
        """
        string = ''
        for e in couple:
            string += e
        splited = string.split('#')[:-1]
        splited = sorted([int(e) for e in splited])
        return [str(e) + '#' for e in splited]

    def __getBipartition(self, name):
        starString=''
        checkString=''
        tempString=''
        for i in range(len(name)):
            if name[i] == '*':
                starString += tempString+'*'
                tempString=''
            elif name[i] == '#':
                checkString += tempString+'#'
                tempString=''
            else:
                tempString += name[i]
        starString = ''.join(self.__starSorted(starString))
        checkString = ''.join(self.__checkSorted(checkString))
        mergedString = starString + checkString
        return starString, checkString, mergedString

    def __getDifference(self, fromSet, toSet):
        temp = []
        for e in fromSet:
            if e not in toSet:
                temp.append(e)
            else:
                continue
        return temp

    """
    
    """
    def __speciesBranchRecurse(self, nodeId, branchLength, coalescentProcess, 
        copiedHaplotypeTree, fromSet, coalSet, recomSet, distanceAboveRoot=None):
        withinBranchProcess = copiedHaplotypeTree[nodeId]
        distanceToAdd = 0
        tempFromSet = fromSet.copy()
        toSet = fromSet.copy()
        if distanceAboveRoot:
            accumulatedDistance = 0
            if withinBranchProcess:
                for i in range(len(withinBranchProcess)):
                    copiedProcess = withinBranchProcess[i]
                    distance = withinBranchProcess[i]['distance']
                    tempDistance = accumulatedDistance
                    accumulatedDistance += distance
                    if accumulatedDistance < distanceAboveRoot:
                        toSet, coalSet, recomSet, _ = \
                            self.__geneBranchRecurse(nodeId=nodeId, distance=distance,
                            distanceToAdd=distanceToAdd, fromSet=tempFromSet,
                            coalescentProcess=coalescentProcess,
                            coalSet=coalSet, recomSet=recomSet, 
                            copiedProcess=copiedProcess)
                        tempFromSet = toSet.copy()
                    else:
                        coalescentProcess[nodeId].append({
                            'fromSet': tempFromSet,
                            'toSet': None,
                            'distance': distanceAboveRoot - tempDistance
                        })
                        break
        else:
            if withinBranchProcess:
                for i in range(len(withinBranchProcess)):
                    copiedProcess = withinBranchProcess[i]
                    distance = withinBranchProcess[i]['distance']
                    toSet, coalSet, recomSet, _ = \
                        self.__geneBranchRecurse(nodeId=nodeId, distance=distance,
                        distanceToAdd=distanceToAdd, fromSet=tempFromSet,
                        coalescentProcess=coalescentProcess,
                        coalSet=coalSet, recomSet=recomSet, 
                        copiedProcess=copiedProcess)
                    tempFromSet = toSet.copy()

        return toSet, coalSet, recomSet
