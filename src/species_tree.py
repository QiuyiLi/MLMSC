import skbio
import numpy as np
from collections import defaultdict
from statistics import mean
from .tree_table import *


class SpeciesTree:
    """
    Nodes are represented in Tree Table which is introduced in tree_table.py
    """

    def __init__(self, randomState):
        self.__randomState = randomState

        self.__treeTable = None
        self.__coalescentRate = None

    def __repr__(self):
        return str(self.__treeTable)

    def __str__(self):
        return str(self.__treeTable)

    @property
    def randomState(self):
        return self.__randomState

    @property
    def treeTable(self):
        return self.__treeTable
    @treeTable.setter
    def treeTable(self, treeTable):
        self.__treeTable = treeTable

    @property
    def coalescentRate(self):
        return self.__coalescentRate
    @coalescentRate.setter
    def coalescentRate(self, coalescentRate):
        self.__coalescentRate = coalescentRate

    def setCoalescentRate(self, coalescentPrmt):
        if ('const' not in coalescentPrmt):
            self.__coalescentRate = self.randomState.gamma(
            shape=coalescentPrmt['shape'], scale=coalescentPrmt['scale'],
            size=len(self.getLeaves()))
        else:
            self.__coalescentRate = np.repeat(coalescentPrmt['const'], 
                len(self.getLeaves()))

    def getSkbioTree(self):
        return self.__treeTable.skbioTree

    def getNodeById(self, id):
        return self.__treeTable.getEntryById(id)

    def getNodes(self):
        return self.__treeTable.table

    def getNodeByName(self, name):
        return self.__treeTable.getEntryByName(name)

    def getRoot(self):
        return self.__treeTable.root

    def getLeaves(self):
        return self.__treeTable.leaves

    def getTreeHeight(self):
        return self.__treeTable.treeHeight

    def getDistanceToLeaf(self, nodeId, branchDistance=0):
        return self.__treeTable.distanceToLeaf(nodeId, branchDistance)

    def initialize(self, path):
        self.__treeTable = TreeTable()
        self.__treeTable.createFromNewickFile(path)
        for speciesNode in self.getNodes():
            for i in range(len(speciesNode.name)):  
                leafName = speciesNode.name[i]
                leafId = self.getNodeByName(leafName).id
                speciesNode.clades.append(leafId)
            if (speciesNode.children and not speciesNode.splits):
                for i in range(len(speciesNode.children)):
                    childName = self.getNodeById(speciesNode.children[i]).name
                    split = []
                    for j in range(len(childName)):  
                        leafName = childName[j]
                        leafId = self.getNodeByName(leafName).id
                        split.append(leafId)
                    speciesNode.splits.append(split)

    def coalescent(self, distanceAboveRoot):
        """
        the main multi-species coalecent function
        """
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
        cladeSet = {}

        # avoid doing repeated coalescence: 
        # a node will be lablled after finishing the coalesencent 
        labelled = {}

        # initialization: 
        # every node is labelled false
        # cladeSet[leafId] = 'leafId*' ('*' as seperater)
        for node in nodes:
            labelled[node.id] = False
            cladeSet[node.id] = \
                [str(node.id) + '*'] if not node.children else []

        while True:
            for leaf in oldLeaves:
                # coalescent finished
                if leaf == root.id:
                    cladeSetIntoRoot = self.__coalescentRecurse(
                        nodeId=root.id, branchLength=distanceAboveRoot,
                        cladeSet=cladeSet, coalescentProcess=coalescentProcess)
                    break
                else:
                    parent = self.getNodeById(leaf).parent
                    children = self.getNodeById(parent).children
                    # if the leaf has been labelled, skip
                    if labelled[leaf]:
                        continue
                    labelled[leaf] = True

                    # first make sure there are genes coming out of both children
                    # then do coalescent within each child branch
                    # cladeSet will be updated from the gene comming into the branch
                    # to genes comming out of the branch
                    if (len(cladeSet[children[0]]) != 0 
                        and len(cladeSet[children[1]]) != 0):
                        cladeSet[children[0]] = self.__coalescentRecurse(
                                                nodeId=children[0], 
                                                branchLength=self.getNodeById(
                                                    children[0]).distanceToParent,
                                                cladeSet=cladeSet, 
                                                coalescentProcess=coalescentProcess)
                        labelled[children[0]] = True
                        cladeSet[children[1]] = self.__coalescentRecurse(
                                                nodeId=children[1],
                                                branchLength=self.getNodeById(
                                                    children[1]).distanceToParent,
                                                cladeSet=cladeSet, 
                                                coalescentProcess=coalescentProcess)  
                        labelled[children[1]] = True
                        # update cladeSet[parent] as the
                        # union of the cladeSet of its children 
                        cladeSet[parent] = list(set().union(
                            cladeSet[children[0]], cladeSet[children[1]]))
                        
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
            # print('old',oldLeaves)
            # print('new',tempNewLeaves)
            oldLeaves = tempNewLeaves.copy()
            newLeaves = []
            labelled = {}
            for node in nodes:
                labelled[node.id] = False

        return coalescentProcess, cladeSetIntoRoot

    def __coalescentRecurse(self, nodeId, branchLength, cladeSet, coalescentProcess):
        """
        This is the recursive part of the multispecies coalescent process:
        Given a set of n genes gathering into a branch in the species tree 
        from the bottom, whenever we come across a point of coalescence, 
        we randomly merge 2 elements in the gene sets, and record the set 
        before the new coalescence, named "from_set", and the set after 
        the coalescence, named "to_set", and the distance from the last 
        coalescent event or the bottom of the branch.
        """
        if len(cladeSet[nodeId]) <= 1:
            coalescentProcess[nodeId].append({
                    'fromSet': cladeSet[nodeId],
                    'toSet': None,
                    'distance': branchLength
                })
            return cladeSet[nodeId]
        else:
            # rate of coalescence in an ancestral branch = mean of each child branch
            coalescentRate = len(cladeSet[nodeId]) \
                      * self.__getCoalescentRateInAncestralBranch(cladeSet[nodeId])
            fakeDistance = min(self.randomState.exponential(
                scale=1.0 / coalescentRate, size=len(cladeSet[nodeId])-1))

            # no coalescent event anymore in this branch
            if branchLength < fakeDistance:
                temp_set = sorted(cladeSet[nodeId])
                coalescentProcess[nodeId].append({
                    'fromSet': temp_set,
                    'toSet': None,
                    'distance': branchLength
                })
                # print(nodeId, coalescentProcess[nodeId],fakeDistance)
                return cladeSet[nodeId]
            else:
                # when coalescent, randomly merge 2 elements in the gene sets
                # if there are more than one genes in the cladeSet, start merging
                if len(cladeSet[nodeId]) >= 2:
                    temp_set = sorted(cladeSet[nodeId])
                    # choose a couple, merge them, then put it back
                    couple = self.randomState.choice(
                        cladeSet[nodeId], size=2, replace=False)
                    cladeSet[nodeId] = [''.join(self.__starSorted(couple))] \
                        + [e for e in cladeSet[nodeId] if e not in couple]

                    # save process
                    coalescentProcess[nodeId].append({
                        'fromSet': temp_set,
                        'toSet': cladeSet[nodeId].copy(),
                        'distance': fakeDistance
                    })
                else:
                    # stop when cladeSet only has one element
                    return cladeSet[nodeId]

                branchLength = branchLength - fakeDistance

                # use recursion to simulate the case when there is
                # more than one coalescent events in the branch
                self.__coalescentRecurse(
                    nodeId=nodeId, branchLength=branchLength, cladeSet=cladeSet, 
                    coalescentProcess=coalescentProcess)

        return cladeSet[nodeId]

    def __getCoalescentRateInAncestralBranch(self, cladeSet):
        indices = []
        for clade in cladeSet:
            splited = clade.split('*')[:-1]
            for index in splited:
                indices.append(int(index))
        return mean(self.coalescentRate[indices])

    def _starInSet(self, target, clade):
        """
        checking whether a given clade is in the target set
        modified for the "*" representation
        """
        if len(target) <= len(clade):
            splited_target = target.split('*')[:-1]
            splited_clade = clade.split('*')[:-1]
            return set(splited_target).issubset(set(splited_clade))
        else:
            return False

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

    def getTimeSequences(self, coalescentProcess):
        """
        backward-in-time coalescent process modified data structure 
        for constructing the coalescent tree in newick format
        e.g.,
        {   
            0: [('0*1*2*', 3.49629952440356)], 
            1: [('1*2*', 2.2767776526850145), ('0*1*2*', 3.49629952440356)], 
            2: [('1*2*', 2.2767776526850145), ('0*1*2*', 3.49629952440356)]
        }
        """
        timeSequences = {}
        geneNodeName = None
        for leaf in self.getLeaves():
            timeSequences[leaf.id] = self.__findAncestors(
                leafName=str(leaf.id) + '*', 
                coalescentProcess=coalescentProcess)
        emptySeq = True
        for leaf in self.getLeaves():
            if timeSequences[leaf.id]:
                emptySeq = False
        if emptySeq:
            timeSequences = {}
            sequence = []
            for speciesNodeId, mergingSets in coalescentProcess.items():
                for mergingSet in mergingSets:
                    if mergingSet['fromSet']:
                        geneNodeName = mergingSet['fromSet'][0]
                        break
                    break
        return timeSequences, geneNodeName

    def __findAncestors(self, leafName, coalescentProcess):
        """
        find the ancestors of the given leaf in reverse time order
        """
        sequence = []
        for speciesNodeId, mergingSets in coalescentProcess.items():
            branchDistance = 0.0
            for mergingSet in mergingSets:
                if mergingSet['toSet']:
                    branchDistance += mergingSet['distance']
                    if (leafName in mergingSet['fromSet'] 
                        and leafName not in mergingSet['toSet']):
                        for element in mergingSet['toSet']:
                            if (len(leafName) < len(element) 
                                and self._starInSet(leafName, element)):
                                coalescentHeight = self.getDistanceToLeaf(
                                    nodeId=speciesNodeId, 
                                    branchDistance=branchDistance)
                                pair = (element, coalescentHeight)
                                sequence.append(pair)
                                sequence += self.__findAncestors(
                                    leafName=element, 
                                    coalescentProcess=coalescentProcess)
        return sequence

    def boundedCoalescent(self, distanceAboveRoot):
        """
        abstract method
        implementation in locus_tree.py
        """
        pass
    
    def incompleteCoalescent(self, distanceAboveRoot):
        """
        abstract method
        implementation in locus_tree.py
        """
        pass







        

    