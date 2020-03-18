import copy
from .species_tree import *
from .tree_table import *

class LocusTree(SpeciesTree):

    def __init__(self, randomState):
        super().__init__(randomState=randomState)

    def initialize(self, nodes, skbioTree):
        self.treeTable = TreeTable()
        self.treeTable.createFromEntries(
            entries=nodes, skbioTree=skbioTree)
    



    
    """""""""""""""
    main functions
    """""""""""""""
    
    """
    bounded coalescent for the locus tree model
    keep trying until all genes coalesce in time
    """
    def boundedCoalescent(self, distanceAboveRoot):
        coalescentProcess, genesIntoRoot = self.coalescent(distanceAboveRoot)
        if len(genesIntoRoot) == 1:
            return coalescentProcess
        else:
            return self.boundedCoalescent(distanceAboveRoot)

    """
    incomplete coalescent:
    generate a haplotype forest and a selected haplotype tree;
    here we don't store the haplotype forest as multiple trees,
    instead we augment the incomplete process to form a fullCoalescentProcess;
    the process exceed the height of the locus tree will be discarded later.
    """
    def incompleteCoalescent(self, distanceAboveRoot):
        root = self.getRoot()
        fullCoalescentProcess, genesIntoRoot = self.coalescent(distanceAboveRoot)

        chosenGene = self.randomState.choice(genesIntoRoot)
        selectedCoalescentProcess = self.__selectCoalescentProcess(
            fullCoalescentProcess, chosenGene)

        fromSet = fullCoalescentProcess[root.id][-1]['fromSet']
        distance = fullCoalescentProcess[root.id][-1]['distance']
        fullCoalescentProcess[root.id].pop()
        while len(fromSet) >= 2:
            temp_set = fromSet
            couple = self.randomState.choice(
                fromSet, size=2, replace=False)
            toSet = [''.join(self.sorted(couple, seperater='*'))] \
                + [e for e in fromSet if e not in couple]
            # save process
            coalescentRate = self.binom(len(temp_set),2) * self.coalescentRate
            coalDistance = self.randomState.exponential(scale=1.0/coalescentRate)
            fullCoalescentProcess[root.id].append({
                'fromSet': temp_set,
                'toSet': toSet,
                'distance': coalDistance+distance
            })
            distance = 0
            fromSet = toSet
        fullCoalescentProcess[root.id].append({
                'fromSet': fromSet,
                'toSet': None,
                'distance': float('inf')
        })
        return fullCoalescentProcess, selectedCoalescentProcess, chosenGene
    
    """
    seclect a haplotype tree from the haplotype forest
    """
    def __selectCoalescentProcess(self, coalescentProcess, chosenGene):
        selectedCoalescentProcess = defaultdict(list)
        for speciesNodeId, mergingSets in coalescentProcess.items():
            for mergingSet in mergingSets:
                distance = mergingSet['distance']
                fromSet = []
                toSet = []
                for clade in mergingSet['fromSet']:
                    if self.starInSet(target=clade, clade=chosenGene):
                        fromSet.append(clade)
                if mergingSet['toSet']:
                    for clade in mergingSet['toSet']:
                        if self.starInSet(target=clade, clade=chosenGene):
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



    """""""""""""""""""""""""""""""""""""""""
    funtions for modelling linked duplictions
    """""""""""""""""""""""""""""""""""""""""

    """
    linked coalescent:
    the copied gene is represented by id*,
    while the new gene is represented by id#.
    """
    def linkedCoalescent(self, copiedHaplotypeTree, copiedRootGene, distanceAboveRoot):
        nodes = self.getNodes()
        root = self.getRoot()
        coalescentProcess = defaultdict(list)
        cladeSetIntoRoot = None

        # leaves of the given species tree
        currentLeaves = [node.id for node in nodes if not node.children]
        # leaves set will be updated in the loop
        newLeaves = []

        fromSets = {}
        toSets = {}
        coalSets = {}
        recomSets = {}
        
        # avoid doing repeated coalescence: 
        # a node will be lablled after finishing the coalesencent 
        labelled = {}

        # initialization: 
        # every node is labelled false
        for node in nodes:
            labelled[node.id] = False
            coalSets[node.id] = []
            toSets[node.id] = []
            fromSets[node.id] = \
                [str(node.id) + '*' + str(node.id) + '#'] if not node.children else []
        recomSets = fromSets.copy()
        
        while True:
            for leaf in currentLeaves:
                # coalescent finished
                if leaf == root.id:
                    cladeSetIntoRoot, _, _ = self.__speciesBranchRecurse(
                        nodeId=root.id, 
                        branchLength=self.getNodeById(root.id).distanceToParent, 
                        coalescentProcess=coalescentProcess, 
                        copiedHaplotypeTree=copiedHaplotypeTree, 
                        fromSet=fromSets[root.id], 
                        coalSet=coalSets[root.id], 
                        recomSet=coalSets[root.id],
                        distanceAboveRoot=distanceAboveRoot)
                    break
                else:
                    # if the leaf has been labelled, skip
                    if labelled[leaf]:
                        continue
                    labelled[leaf] = True

                    parent = self.getNodeById(leaf).parent
                    children = self.getNodeById(parent).children
                    
                    # first make sure there are genes coming out of both children
                    # then do coalescent within each child branch
                    # cladeSet will be updated from the gene comming into the branch
                    # to genes comming out of the branch
                    if (len(fromSets[children[0]]) != 0 
                        and len(fromSets[children[1]]) != 0):
                        toSets[children[0]], coalSets[children[0]], recomSets[children[0]] = \
                            self.__speciesBranchRecurse(
                                nodeId=children[0], 
                                branchLength=self.getNodeById(children[0]).distanceToParent, 
                                coalescentProcess=coalescentProcess, 
                                copiedHaplotypeTree=copiedHaplotypeTree, 
                                fromSet=fromSets[children[0]], 
                                coalSet=coalSets[children[0]], 
                                recomSet=recomSets[children[0]])
                        labelled[children[0]] = True
                        toSets[children[1]], coalSets[children[1]], recomSets[children[1]] = \
                            self.__speciesBranchRecurse(
                                nodeId=children[1], 
                                branchLength=self.getNodeById(children[1]).distanceToParent, 
                                coalescentProcess=coalescentProcess, 
                                copiedHaplotypeTree=copiedHaplotypeTree, 
                                fromSet=fromSets[children[1]], 
                                coalSet=coalSets[children[1]], 
                                recomSet=recomSets[children[1]]) 
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
            # tempNewLeaves = [] 
            # for newLeaf in newLeaves:
            #     if newLeaf not in tempNewLeaves:
            #         tempNewLeaves.append(newLeaf)

            # re-initialization for the next recursion
            # currentLeaves <- newLeaves
            # label <- false
            # currentLeaves = tempNewLeaves
            currentLeaves = newLeaves
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
                elif copiedRootGene: 
                    for e in copiedRootGene:
                        if e in clade:
                            ancestralClades.append(clade)
                            break
        fullClades = ancestralClades + nonAncestralClades
        chosenGeneName = self.randomState.choice(cladeSetIntoRoot)

        filteredProcess, filteredClades = self.__filteredLinkedCoalescentProcess(
            coalescentProcess, cladeSetIntoRoot)

        fullProcess = copy.deepcopy(filteredProcess)
        fromSet = fullProcess[root.id][-1]['fromSet']
        distance = fullProcess[root.id][-1]['distance']
        fullProcess[root.id].pop()
        while len(fromSet) >= 2:
            couple = self.randomState.choice(
                fromSet, size=2, replace=False)
            toSet = [''.join(self.sorted(couple, seperater='*'))] \
                + [e for e in fromSet if e not in couple]
            # save process
            coalescentRate = self.binom(len(fromSet),2) * self.coalescentRate
            coalDistance = self.randomState.exponential(scale=1.0/coalescentRate)
            fullProcess[root.id].append({
                'fromSet': fromSet,
                'toSet': toSet,
                'distance': coalDistance+distance
            })
            distance = 0
            fromSet = toSet
        fullProcess[root.id].append({
                'fromSet': fromSet,
                'toSet': None,
                'distance': float('inf')
        })

        if chosenGeneName in fullClades:
            ancestral = False
            geneNodeName = None
            if chosenGeneName in ancestralClades:
                ancestral = True
                geneNodeName, chosenGeneName, _ = self.__getBipartition(chosenGeneName)
            chosenGeneName = self.__starReplace(chosenGeneName)
            
            selectedProcess = self.__selectCoalescentProcess(filteredProcess, chosenGeneName)

            return fullProcess, selectedProcess, chosenGeneName, geneNodeName, ancestral
        else: 
            # discad the unobservable ancestral duplication
            return fullProcess, None, None, None, True

    """
    linked coalescent within a gene branch
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
            return self.__geneBranchRecurse(
                nodeId=nodeId, distance=distance-coalDistance, 
                distanceToAdd=0, fromSet=toSet, 
                coalescentProcess=coalescentProcess,
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
            return self.__geneBranchRecurse(
                nodeId=nodeId, distance=distance-recomDistance, 
                distanceToAdd=0, fromSet=toSet, 
                coalescentProcess=coalescentProcess,
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
                        coalescentProcess[nodeId].append({
                            'fromSet': fromSet,
                            'toSet': None,
                            'distance': distance
                        })
                        distanceToAdd += distance
                        coalSet = coalSet
                        recomSet = recomSet
                else:
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
    linked coalescent within a species branch
    """
    def __speciesBranchRecurse(self, nodeId, branchLength, coalescentProcess, 
        copiedHaplotypeTree, fromSet, coalSet, recomSet, distanceAboveRoot=None):
        withinBranchProcess = copiedHaplotypeTree[nodeId]
        distanceToAdd = 0
        tempFromSet = fromSet
        toSet = fromSet
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
                            self.__geneBranchRecurse(
                                nodeId=nodeId, distance=distance,
                                distanceToAdd=distanceToAdd, fromSet=tempFromSet,
                                coalescentProcess=coalescentProcess,
                                coalSet=coalSet, recomSet=recomSet, 
                                copiedProcess=copiedProcess)
                        tempFromSet = toSet
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
                        self.__geneBranchRecurse(
                            nodeId=nodeId, distance=distance,
                            distanceToAdd=distanceToAdd, fromSet=tempFromSet,
                            coalescentProcess=coalescentProcess,
                            coalSet=coalSet, recomSet=recomSet, 
                            copiedProcess=copiedProcess)
                    tempFromSet = toSet

        return toSet, coalSet, recomSet

    """
    filter out the copied genes in the coalescent process
    and change the seperator of the genes from '#' to '*'
    """
    def __filteredLinkedCoalescentProcess(self, coalescentProcess, cladeSetIntoRoot):
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
                if not self.__identicalSets(fromSet, toSet):
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

    """""""""""""""""""""""""""""""""""""""
    utility functions for linkedCoalescent
    """""""""""""""""""""""""""""""""""""""
    def __identicalSets(self, set1, set2):
        if not set2:
            return False
        elif len(set1) != len(set2):
            return False
        else:
            for e in set1:
                if e not in set2:
                    return False
            return True

    def __starReplace(self, string):
        newString = ''
        for e in string:
            if e == '#':
                newString += '*'
            else:
                newString += e
        return newString

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
        starString = ''.join(self.sorted(starString, seperater='*'))
        checkString = ''.join(self.sorted(checkString, seperater='#'))
        mergedString = starString + checkString
        return starString, checkString, mergedString

    def __getDifference(self, fromSet, toSet):
        diff = []
        for e in fromSet:
            if e not in toSet:
                diff.append(e)
            else:
                continue
        return diff
