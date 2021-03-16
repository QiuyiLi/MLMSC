import skbio
import numpy as np
import copy
from collections import defaultdict
from statistics import mean
from .tree_table import *
from .locus_tree import *

class HaplotypeTree:
    """
    Nodes are represented in Tree Table which is introduced in tree_table.py
    """
    def __init__(self, randomState, speciesTree, locusTree):
        self.__randomState: RandomState = randomState
        self.__speciesTree: SpeciesTree = speciesTree
        self.__locusTree: LocusTree = locusTree

        self.__coalescentProcess: defaultdict = None
        self.__fullCoalescentProcess: defaultdict  = None
        # self.__treeTable = None
        self.__treeTable: TreeTable = None
        self.__parameters = {}

    def __repr__(self):
        return str(self.__treeTable)

    def __str__(self):
        return str(self.__treeTable)

    @property
    def randomState(self):
        return self.__randomState

    @property
    def speciesTree(self):
        return self.__speciesTree

    @property
    def locusTree(self):
        return self.__locusTree

    @property
    def coalescentProcess(self):
        return self.__coalescentProcess

    @property
    def fullCoalescentProcess(self):
        return self.__fullCoalescentProcess

    @property
    def treeTable(self):
        return self.__treeTable

    @property
    def parameters(self):
        return self.__parameters
    @parameters.setter
    def parameters(self, parameters):
        self.__parameters = parameters

    def getSkbioTree(self):
        return self.__treeTable.skbioTree
    def setSkbioTree(self, skbioTree):
        self.__treeTable.skbioTree = skbioTree

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

    def initialize(self, locusTree, coalescentProcess=None, fullCoalescentProcess=None, rename=True, event=None):
        if not coalescentProcess:
            # ordinary coalesent for constructing the original haplotype tree
            coalescentProcess, cladeSetIntoRoot = locusTree.coalescent(
                distanceAboveRoot=float('inf'))
            fullCoalescentProcess = coalescentProcess
        self.__coalescentProcess = coalescentProcess
        self.__fullCoalescentProcess = fullCoalescentProcess
        if event == 'transfer':
            timeSequences, geneNodeName = self.speciesTree.getTimeSequences(
                coalescentProcess=coalescentProcess)
        else:
            timeSequences, geneNodeName = locusTree.getTimeSequences(
                coalescentProcess=coalescentProcess)
        skbioTree = self.createSkbioTree(timeSequences, geneNodeName)
        self.readFromSkbioTree(skbioTree, rename)

        # for geneNode in self.getNodes():
        #     clades = geneNode.name.split('*')[:-1]
        #     clades = [int(j) for j in clades]
        #     geneNode.clades = clades
        #     if (geneNode.children and not geneNode.splits):
        #         for i in range(len(geneNode.children)):
        #             geneChildName = self.getNodeById(geneNode.children[i]).name 
        #             childClades = geneChildName.split('*')[:-1]
        #             childClades = [int(j) for j in childClades]
        #             geneNode.splits.append(childClades)

    def setParameters(self, duplicationPrmt, transferPrmt, lossPrmt, unlinkProb, hemiplasy, verbose):
        self.__parameters['d'] = duplicationPrmt
        self.__parameters['t'] = transferPrmt
        self.__parameters['l'] = lossPrmt
        self.__parameters['u'] = unlinkProb
        self.__parameters['h'] = hemiplasy
        self.__parameters['v'] = verbose

    def createSkbioTree(self, timeSequences, geneNodeName):
        """
        creat a tree structure consistent with the package skbio
        using timeSequence
        """
        skbioTree = skbio.tree.TreeNode()
        tempTimeSequences = copy.deepcopy(timeSequences)
        
        if geneNodeName: 
            # empty timeSequence {speicesId1: [], speicesId2: [], speicesId3: []}
            skbioTree.name = geneNodeName # speciesid*
            skbioTree.length = None
            return skbioTree
        # remove empty entry to see if there is coalescent or not
        for k, v in timeSequences.items():
            if not v: del tempTimeSequences[k]
        if len(tempTimeSequences) > 0:
            # skbioTree.name = name of the root
            # name at the end of the timeSequence
            skbioTree.name = next(iter(tempTimeSequences.values()))[-1][0]
            self.__createSkbioTreeRecurse(
                skbioTree=skbioTree, timeSequences=timeSequences)
            skbioTree.length = None
        else:
            # empty timeSequence {speicesId: []}
            skbioTree.name = str(next(iter(timeSequences))) + '*' # speciesid*
            skbioTree.length = None
        return skbioTree

    # utility function used in __createSkbioTreeRecurse
    def __distanceToParent(self, nodeName, parentName, timeSequences):
        for leaf, sequence in timeSequences.items():
            if (nodeName.count('*') == 1 and nodeName.split('*')[0] == str(leaf)):
                for pair in sequence:
                    if pair[0] == parentName:
                        return pair[1]
            else:
                prevPair = None
                for pair in sequence:
                    if (prevPair != None 
                        and prevPair[0] == nodeName 
                        and pair[0] == parentName):
                        return pair[1] - prevPair[1]
                    prevPair = pair
        return None

    # utility function used in __createSkbioTreeRecurse
    def __starReplace(self, string, substring):
        a = string.split('*')[:-1]
        b = substring.split('*')[:-1]
        diff = set(a).difference(set(b))
        return ''.join([e + '*' for e in sorted(list(diff))])

    def __createSkbioTreeRecurse(self, skbioTree, timeSequences):
        # one node (leaf) case (trivial)
        if skbioTree.name.count('*') == 1:
            skbioTree.length = self.__distanceToParent(
                nodeName=skbioTree.name, parentName=skbioTree.parent.name, 
                timeSequences=timeSequences)
            # print(skbioTree.name, skbioTree.parent.name)
            return
        # two nodes case (trivial)
        elif skbioTree.name.count('*') == 2:
            childLName = skbioTree.name.split('*')[0] + '*'
            childRName = skbioTree.name.split('*')[1] + '*'
            childL = skbio.tree.TreeNode(
                name=childLName, 
                length=self.__distanceToParent(
                    nodeName=childLName, parentName=skbioTree.name, 
                    timeSequences=timeSequences), 
                parent=skbioTree)
            childR = skbio.tree.TreeNode(
                name=childRName, 
                length=self.__distanceToParent(
                    nodeName=childRName, parentName=skbioTree.name, 
                    timeSequences=timeSequences),
                parent=skbioTree)
            skbioTree.children = [childL, childR]
            return
        # otherwise, more general cases
        else:
            isFound = False
            for _, sequence in timeSequences.items():
                prevPair = None
                for pair in sequence:
                    if (prevPair != None and skbioTree.name == pair[0]):
                        childLName = prevPair[0]
                        splited = childLName.split('*')[:-1]
                        splited = sorted([int(e) for e in splited])
                        childLName = ''.join([str(e) + '*' for e in splited])

                        childRName = self.__starReplace(skbioTree.name, prevPair[0])
                        splited = childRName.split('*')[:-1]
                        splited = sorted([int(e) for e in splited])
                        childRName = ''.join([str(e) + '*' for e in splited])
                        childL = skbio.tree.TreeNode(
                            name=childLName, 
                            length=self.__distanceToParent(
                                nodeName=childLName, 
                                parentName=skbioTree.name, 
                                timeSequences=timeSequences), 
                            parent=skbioTree)
                        childR = skbio.tree.TreeNode(
                            name=childRName, 
                            length=self.__distanceToParent(
                                nodeName=childRName, 
                                parentName=skbioTree.name, 
                                timeSequences=timeSequences),
                            parent=skbioTree)
                        self.__createSkbioTreeRecurse(
                            skbioTree=childL, timeSequences=timeSequences)
                        self.__createSkbioTreeRecurse(
                            skbioTree=childR, timeSequences=timeSequences)
                        skbioTree.children = [childL, childR]
                        isFound = True
                        break
                    prevPair = pair
                if isFound: break

    def readFromSkbioTree(self, skbioTree, rename=True):
        self.__treeTable = TreeTable()
        self.__treeTable.createFromSkbioTree(skbioTree, rename)





    
    """""""""""""""
    main functions
    """""""""""""""
    def Lprocess(self, distanceAboveRoot, event=None):
        events = []
        # trivial case
        if len(self.getNodes()) == 1:
            distanceAboveRoot = event['eventHeight']
        self.__LprocessRecurse(
            skbioTreeNode=self.getSkbioTree(), 
            distanceAboveRoot=distanceAboveRoot, events=events)
        return events

    def __LprocessRecurse(self, skbioTreeNode, distanceAboveRoot, events):
        node = self.getNodeByName(skbioTreeNode.name)
        lossRate = self.__parameters['l']
        if lossRate == 0:
            distanceL = float('inf')
        else:
            distanceL = self.randomState.exponential(
                scale=1.0 / lossRate)

        if (distanceL < distanceAboveRoot):      
            # loss happens first, the seaching process stops at the loss point
            eventHeight = self.getDistanceToLeaf(node.id) + distanceAboveRoot - distanceL
            speciesId, distanceAboveSpeciesNode = self.__mapEventToSpeciesTree(
                geneId=node.id, eventHeight=eventHeight, speciesId=None)
            events.append({
                'type': 'loss',
                'geneNodeId': node.id,      
                # closest gene node to the event from below
                'geneNodeName': node.name, 
                'distanceToGeneNode': distanceAboveRoot - distanceL,
                'eventHeight': eventHeight,
                'speciesNodeId': speciesId,     
                # closest species node to the event from below
                'distanceToSpeciesNode': distanceAboveSpeciesNode
            })
        elif (node.children):     # if children branches exist
                childL = skbioTreeNode.children[0]
                childR = skbioTreeNode.children[1]
                distanceToChildL = node.distanceToChildren[0]
                distanceToChildR = node.distanceToChildren[1]
                self.__LprocessRecurse(
                    skbioTreeNode=childL, 
                    distanceAboveRoot=distanceToChildL, events=events)
                self.__LprocessRecurse(
                    skbioTreeNode=childR, 
                    distanceAboveRoot=distanceToChildR, events=events)

    def Dprocess(self, unlinked, distanceAboveRoot, threshold=float('inf'), event=None):
        events = []
        # trivial case
        if len(self.getNodes()) == 1:
            distanceAboveRoot = event['eventHeight']
        self.__DprocessRecurse(
            unlinked = unlinked,
            skbioTreeNode=self.getSkbioTree(), 
            distanceAboveRoot=distanceAboveRoot, 
            threshold=threshold, events=events)
        return events

    def __DprocessRecurse(self, unlinked, skbioTreeNode, distanceAboveRoot, threshold, events):
        unlinkedProb = self.__parameters['u']
        node = self.getNodeByName(skbioTreeNode.name)
        if unlinked:
            duplicationRate = self.__parameters['d'] * unlinkedProb
        else:
            duplicationRate = self.__parameters['d'] * (1-unlinkedProb)
        if duplicationRate == 0:
            distanceD = float('inf')
        else:
            distanceD = self.randomState.exponential(
                scale=1.0 / duplicationRate)
        # duplication happens first
        if (distanceD < distanceAboveRoot):
            eventHeight = self.getDistanceToLeaf(node.id) + distanceAboveRoot - distanceD
            speciesId, distanceAboveSpeciesNode = self.__mapEventToSpeciesTree(
                geneId=node.id, eventHeight=eventHeight, speciesId=None)
            if eventHeight < threshold:
                events.append({
                    'type': 'duplication',
                    'eventHeight': eventHeight,
                    'speciesNodeId': speciesId,
                    'distanceToSpeciesNode': distanceAboveSpeciesNode,
                    'unlinked': unlinked
                })
            # looking for more events on the same branch
            self.__DprocessRecurse(
                unlinked=unlinked,
                skbioTreeNode=skbioTreeNode, 
                distanceAboveRoot=distanceAboveRoot - distanceD, 
                threshold=threshold, events=events)
        elif (node.children):     
            # look for events in the children branches
            childL = skbioTreeNode.children[0]
            childR = skbioTreeNode.children[1]
            distanceToChildL = node.distanceToChildren[0]
            distanceToChildR = node.distanceToChildren[1]
            self.__DprocessRecurse(
                unlinked=unlinked,
                skbioTreeNode=childL, 
                distanceAboveRoot=distanceToChildL, threshold=threshold, events=events)
            self.__DprocessRecurse(
                unlinked=unlinked,
                skbioTreeNode=childR, 
                distanceAboveRoot=distanceToChildR, threshold=threshold, events=events)

    def Tprocess(self, distanceAboveRoot, threshold=float('inf'), event=None):
        events = []
        # trivial case
        if len(self.getNodes()) == 1:
            distanceAboveRoot = event['eventHeight']
        self.__TprocessRecurse(
            skbioTreeNode=self.getSkbioTree(), 
            distanceAboveRoot=distanceAboveRoot, 
            threshold=threshold, events=events)
        return events

    def __TprocessRecurse(self, skbioTreeNode, distanceAboveRoot, threshold, events):
        node = self.getNodeByName(skbioTreeNode.name)
        transferRate = self.__parameters['t']
        if transferRate == 0:
            distanceT = float('inf')
        else:
            distanceT = self.randomState.exponential(
                scale=1.0 / transferRate)
        if distanceT < distanceAboveRoot:
            eventHeight = self.getDistanceToLeaf(node.id) + distanceAboveRoot - distanceT
            speciesTreeHeight = self.speciesTree.getTreeHeight()
            if eventHeight < speciesTreeHeight:
                targetSpeciesId, _ = self.__mapEventToSpeciesTree(
                    geneId=node.id, eventHeight=eventHeight, speciesId=None)

                originalSpeciesId, _ = self.__findTransferOrigin(
                    eventHeight=eventHeight, geneId=node.id, targetSpeciesId=targetSpeciesId)
                if originalSpeciesId and eventHeight < threshold:
                    distanceAboveSpeciesNode = eventHeight - \
                        self.speciesTree.getDistanceToLeaf(originalSpeciesId)
                    events.append({
                        'type': 'transfer',
                        'targetSpeciesId': targetSpeciesId,
                        'eventHeight': eventHeight,
                        'speciesNodeId': originalSpeciesId,      
                        # closest species node to the origin from below
                        'distanceToSpeciesNode': distanceAboveSpeciesNode,
                        'unlinked': True
                    })
            self.__TprocessRecurse(
                skbioTreeNode=skbioTreeNode, 
                distanceAboveRoot=distanceAboveRoot - distanceT, 
                threshold=threshold, events=events)
        elif (node.children):     
            # look for events in the children branches
            childL = skbioTreeNode.children[0]
            childR = skbioTreeNode.children[1]
            distanceToChildL = node.distanceToChildren[0]
            distanceToChildR = node.distanceToChildren[1]
            self.__TprocessRecurse(
                skbioTreeNode=childL, 
                distanceAboveRoot=distanceToChildL, threshold=threshold, events=events)
            self.__TprocessRecurse(
                skbioTreeNode=childR, 
                distanceAboveRoot=distanceToChildR, threshold=threshold, events=events)

    # M function: map the event point to where it occurs in the species tree
    def __mapEventToSpeciesTree(self, geneId, eventHeight, speciesId=None):
        if speciesId == None:
            speciesId = self.__mapGeneIdToSpeciesId(geneId=geneId, eventHeight=eventHeight)
        distanceAboveSpeciesNode = eventHeight - \
            self.speciesTree.getDistanceToLeaf(speciesId)
        if speciesId == self.speciesTree.getRoot().id:
            return speciesId, distanceAboveSpeciesNode
        else:
            speciesIdParent = self.speciesTree.getNodeById(speciesId).parent
            speciesDistanceParent = self.speciesTree.getDistanceToLeaf(speciesIdParent)
            if speciesDistanceParent > eventHeight:
                return speciesId, distanceAboveSpeciesNode
            else:
                return self.__mapEventToSpeciesTree(
                    geneId=geneId, eventHeight=eventHeight, speciesId=speciesIdParent)

    # map the gene node to the closest species node from below
    def __mapGeneIdToSpeciesId(self, geneId, eventHeight):
        speciesId = None
        geneName = self.getNodeById(geneId).name
        if self.coalescentProcess:
            # non-trivial case
            for speciesNodeId, mergingSets in self.coalescentProcess.items():
                for mergingSet in mergingSets:
                    if (mergingSet['toSet']
                        and geneName in mergingSet['toSet'] 
                        and geneName not in mergingSet['fromSet']):
                        speciesId = speciesNodeId
            if speciesId == None:
                speciesId = int(geneName[:-1])
                speciesNode = self.speciesTree.getNodeById(speciesId)
                while speciesNode.parent >= 0:
                    if eventHeight >= self.speciesTree.getDistanceToLeaf(speciesNode.parent):
                        speciesId = speciesNode.parent
                        speciesNode = self.speciesTree.getNodeById(speciesId)
                    else:
                        break
        return speciesId 

    def __findTransferOrigin(self, eventHeight, geneId, targetSpeciesId):    
        originNodes = self.locusTree.getNodes()
        speciesNodes = self.speciesTree.getNodes()
        originNodesList = []
        speciesNodesList = []
        for node in originNodes:
            if node.id == targetSpeciesId:
                continue
            if node.id == self.locusTree.getRoot().id:
                nodeHeight = self.locusTree.getDistanceToLeaf(node.id)
                if eventHeight > nodeHeight:
                    if eventHeight < self.locusTree.distanceAboveRoot:
                        originNodesList.append(node.id)
                    break
            else:
                parentHeight = self.locusTree.getDistanceToLeaf(
                    self.locusTree.getNodeById(node.parent).id)
            if parentHeight > eventHeight:
                nodeHeight = self.locusTree.getDistanceToLeaf(node.id)
                if nodeHeight <= eventHeight:
                    originNodesList.append(node.id)

        for node in speciesNodes:
            if node.id == targetSpeciesId:
                continue
            if node.id == self.speciesTree.getRoot().id:
                continue
            parentHeight = self.speciesTree.getDistanceToLeaf(
                self.speciesTree.getNodeById(node.parent).id)
            if parentHeight > eventHeight:
                nodeHeight = self.speciesTree.getDistanceToLeaf(node.id)
                if nodeHeight <= eventHeight:
                    speciesNodesList.append(node.id)
        
        prob = len(originNodesList)/len(speciesNodesList)
        scaleFactor = self.randomState.binomial(1, prob)
        if originNodesList:
            if scaleFactor == 1:
                return self.randomState.choice(originNodesList), targetSpeciesId
            else:
                return None, targetSpeciesId
        else:
            return None, targetSpeciesId

    def DTLprocess(self, locusTree, haplotypeTree, initial=False, event=None, rootLength=0):
        lossRate = self.__parameters['l']
        duplicationRate = self.__parameters['d']
        transferRate = self.__parameters['t']
        unlinkProb = self.__parameters['u']

        if lossRate > 0:
            haplotypeTreeEvents = haplotypeTree.Lprocess(
                    event=event, distanceAboveRoot=rootLength)
        else: 
            haplotypeTreeEvents = []

        if duplicationRate*unlinkProb > 0:
            coalescentTreeProcessUD, _ = locusTree.coalescent(
                distanceAboveRoot=float('inf'))
            coalescentTreeUD = HaplotypeTree(
                randomState=self.randomState, 
                speciesTree=self.speciesTree, 
                locusTree=locusTree)
            coalescentTreeUD.initialize(
                locusTree=locusTree, 
                coalescentProcess=coalescentTreeProcessUD, 
                fullCoalescentProcess=coalescentTreeProcessUD,
                rename=False)
            coalescentTreeUD.parameters = self.parameters
            if initial:
                rootLength = max(0, self.speciesTree.getTreeHeight()-coalescentTreeUD.getTreeHeight())
                threshold = self.speciesTree.getTreeHeight()
                newCoalescentTreeEventsUD = coalescentTreeUD.Dprocess(unlinked=True,
                    distanceAboveRoot=rootLength, threshold=threshold, event=event)
            else:
                rootLength = max(0, event['eventHeight'] - coalescentTreeUD.getTreeHeight())
                threshold = event['eventHeight']
                newCoalescentTreeEventsUD = coalescentTreeUD.Dprocess(unlinked=True,
                    distanceAboveRoot=rootLength, threshold=threshold, event=event)
        else:
            newCoalescentTreeEventsUD = []

        if duplicationRate*(1-unlinkProb) > 0: 
            coalescentTreeProcessLD = locusTree.linkedCoalescent(
                copiedHaplotypeTree=haplotypeTree.fullCoalescentProcess, 
                copiedRootGene=None,
                distanceAboveRoot=haplotypeTree.getTreeHeight()-locusTree.getTreeHeight())[0]
            locusRootId = locusTree.getRoot().id
            coalescentTreeLD = HaplotypeTree(
                randomState=self.randomState, 
                speciesTree=self.speciesTree, 
                locusTree=locusTree)
            coalescentTreeLD.initialize(
                locusTree=locusTree, 
                coalescentProcess=coalescentTreeProcessLD, 
                fullCoalescentProcess=coalescentTreeProcessLD,
                rename=False)
            coalescentTreeLD.parameters = self.parameters
            if initial:
                rootLength = max(0, self.speciesTree.getTreeHeight()-coalescentTreeLD.getTreeHeight())
                threshold = self.speciesTree.getTreeHeight()
            else:
                rootLength = max(0, event['eventHeight'] - coalescentTreeLD.getTreeHeight())
                threshold = event['eventHeight']
            newCoalescentTreeEventsLD = coalescentTreeLD.Dprocess(unlinked=False,
                distanceAboveRoot=rootLength, threshold=threshold, event=event)
        else:
            newCoalescentTreeEventsLD = []

        if transferRate > 0:     
            coalescentTreeProcessT, _ = self.speciesTree.coalescent(
                distanceAboveRoot=float('inf'))
            coalescentTreeT = HaplotypeTree(
                randomState=self.randomState, 
                speciesTree=self.speciesTree, 
                locusTree=locusTree)
            coalescentTreeT.initialize(
                locusTree=locusTree, 
                coalescentProcess=coalescentTreeProcessT, 
                fullCoalescentProcess=coalescentTreeProcessT,
                rename=False,
                event = 'transfer')
            coalescentTreeT.parameters = self.parameters
            if initial:
                rootLength = 0
                threshold = float('inf')
            else:
                rootLength = max(0, event['eventHeight'] - coalescentTreeT.getTreeHeight())
                threshold = event['eventHeight']
            newCoalescentTreeEventsT = coalescentTreeT.Tprocess(
                distanceAboveRoot=rootLength, threshold=threshold, event=event)
        else:
            newCoalescentTreeEventsT = []

        # newEvents = haplotypeTreeEvents + newCoalescentTreeEvents
        # newEvents.sort(reverse=True, key=lambda x: x['eventHeight'])
        newCoalescentTreeEvents = newCoalescentTreeEventsUD + newCoalescentTreeEventsLD + \
            newCoalescentTreeEventsT
        # newCoalescentTreeEvents = newCoalescentTreeEventsUD + \
        #     newCoalescentTreeEventsT
        newCoalescentTreeEvents.sort(reverse=True, key=lambda x: x['eventHeight'])
        newEvents = newCoalescentTreeEvents + haplotypeTreeEvents
        return newEvents
    
    def coalescentJoining(self, event, haplotypeTree):
        fullCoalescentProcess=haplotypeTree.fullCoalescentProcess
        selectedCoalescentProcess=haplotypeTree.coalescentProcess
        nodeId = event['speciesNodeId']
        distanceAboveSpeciesNode = event['distanceToSpeciesNode']
        fullCoalescentTree = HaplotypeTree(
            randomState=self.randomState, 
            speciesTree=self.speciesTree, 
            locusTree=self.speciesTree)
        fullCoalescentTree.initialize(
            locusTree=self.speciesTree, 
            coalescentProcess=fullCoalescentProcess, 
            fullCoalescentProcess=fullCoalescentProcess,
            rename=False)
        geneNodeName, distanceAboveGeneNode, branchLength = \
            self.__coalescentJoiningRecurse(
                nodeId=nodeId, event=event, 
                fullCoalescentProcess=fullCoalescentProcess, 
                fullCoalescentTree=fullCoalescentTree, 
                selectedCoalescentProcess=selectedCoalescentProcess,
                haplotypeTree=haplotypeTree, 
                branchLength=0, passed=False)
        return geneNodeName, distanceAboveGeneNode, branchLength

    def __coalescentJoiningRecurse(self, nodeId, event, 
        fullCoalescentProcess, fullCoalescentTree, 
        selectedCoalescentProcess, haplotypeTree, 
        branchLength, passed):
        unitCoalescentRate = self.speciesTree.coalescentRate
        lossRate = self.__parameters['l']
        fullProcess = fullCoalescentProcess[nodeId]
        selectedProcess = selectedCoalescentProcess[nodeId]
        distanceAboveSpeciesNode = event['distanceToSpeciesNode']
        for e in fullProcess:
            if not passed:
                distanceAboveSpeciesNode = distanceAboveSpeciesNode - e['distance']
                if distanceAboveSpeciesNode < 0:
                    passed = True
                    if len(e['fromSet']) >= 1:
                        coalescentRate = len(e['fromSet'])*unitCoalescentRate
                        coalDistance = self.randomState.exponential(scale=1.0 / coalescentRate)
                    else: 
                        coalDistance = float('inf')
                    if lossRate == 0:
                        lossDistance = float('inf')
                    else:
                        lossDistance = self.randomState.exponential(scale=1.0 / lossRate)
                    if coalDistance <  min(abs(distanceAboveSpeciesNode), lossDistance):
                        geneNodeName = self.randomState.choice(e['fromSet'])
                        geneNodeId = fullCoalescentTree.getNodeByName(geneNodeName).id
                        branchLength += coalDistance
                        if branchLength + event['eventHeight'] > self.locusTree.getTreeHeight():
                            # exceed the locus tree height
                            return None, None, None
                        distanceAboveGeneNode = branchLength + event['eventHeight'] \
                            - fullCoalescentTree.getDistanceToLeaf(geneNodeId)
                        for element in selectedProcess:
                            if geneNodeName in element['fromSet']:
                                return geneNodeName, distanceAboveGeneNode, branchLength
                        # coalesecent with unselected tree
                        return None, None, None
                    elif lossDistance < min(abs(distanceAboveSpeciesNode), coalDistance):
                        # get lost when joining back
                        return None, None, None
                    else:
                        branchLength += abs(distanceAboveSpeciesNode)
                        continue    
            else:
                if len(e['fromSet']) >= 1:
                        coalescentRate = len(e['fromSet'])*unitCoalescentRate
                        coalDistance = self.randomState.exponential(scale=1.0 / coalescentRate)
                else: 
                    coalDistance = float('inf')
                if lossRate == 0:
                    lossDistance = float('inf')
                else:
                    lossDistance = self.randomState.exponential(scale=1.0 / lossRate)
                if coalDistance <  min(e['distance'], lossDistance):
                    geneNodeName = self.randomState.choice(e['fromSet'])
                    geneNodeId = fullCoalescentTree.getNodeByName(geneNodeName).id
                    branchLength += coalDistance
                    if branchLength + event['eventHeight'] > self.locusTree.getTreeHeight():
                        # exceed the locus tree height
                        return None, None, None
                    distanceAboveGeneNode = branchLength + event['eventHeight'] \
                        - fullCoalescentTree.getDistanceToLeaf(geneNodeId)
                    for element in selectedProcess:
                        if geneNodeName in element['fromSet']:
                            return geneNodeName, distanceAboveGeneNode, branchLength
                    # coalesecent with unselected tree
                    return None, None, None
                elif lossDistance < min(e['distance'], coalDistance):
                    return None, None, None
                else:
                    branchLength += e['distance']
                    continue
        parentId = haplotypeTree.speciesTree.getNodeById(nodeId).parent
        if (parentId 
            and parentId != -1 
            and haplotypeTree.speciesTree.getNodeById(parentId) in haplotypeTree.locusTree.getNodes()):   
            return self.__coalescentJoiningRecurse(nodeId=parentId, event=event, 
                fullCoalescentProcess=fullCoalescentProcess, 
                fullCoalescentTree=fullCoalescentTree,
                selectedCoalescentProcess=selectedCoalescentProcess, 
                haplotypeTree=haplotypeTree,  
                branchLength=branchLength, passed = True)
        else:
            # no coalescent
            return None, None, None

    def addNewLociShell(self, events, haplotypeTree, level, completeCount, incompleteCount, unlinked_d_number_total, unlinked_d_number_survived):
        global dic 
        dic = defaultdict(list)
        return self.addNewLoci(events=events, haplotypeTree=haplotypeTree, level=level, completeCount=completeCount, incompleteCount=incompleteCount, unlinked_d_number_total=unlinked_d_number_total, unlinked_d_number_survived=unlinked_d_number_survived)

    def addNewLoci(self, events, haplotypeTree, level, completeCount, incompleteCount, unlinked_d_number_total, unlinked_d_number_survived):
        """
        1. simulate events on the current locus tree 
        2. construct the corresponding new locus tree
        3. generate the new haplotype tree and forest in the new locus tree
        4. simulate events on the new locus tree 
        5. recurse
        """
        if not dic[level]:
            eventIndex = -1
        else: 
            eventIndex = max(dic[level])
        
        for event in events:
            if (event['type'] == 'loss'):
                # print(1)
                geneNodeName = event['geneNodeName']
                geneNode = haplotypeTree.getSkbioTree().find(geneNodeName)
                geneNode.name = geneNode.name + '_loss'

            elif (event['type'] == 'duplication'):
                speciesId = event['speciesNodeId']
                distanceAboveSpeciesNode = event['distanceToSpeciesNode']

                newHaplotypeTree, chosenGeneName, geneNodeName, ancestral, completeCount, incompleteCount, unlinked_d_number_total, unlinked_d_number_survived = \
                    self.__addNewLociRecurse(haplotypeTree=haplotypeTree,
                    event=event, newLocusRootId=speciesId,
                    distanceAboveRoot=distanceAboveSpeciesNode, 
                    level=level, completeCount=completeCount,
                    incompleteCount=incompleteCount,
                    unlinked_d_number_total=unlinked_d_number_total, 
                    unlinked_d_number_survived=unlinked_d_number_survived)

                if newHaplotypeTree:
                    verbose = self.parameters['v']
                    if verbose:
                        print(newHaplotypeTree)
                        print('new haplotype tree:')	
                        print(newHaplotypeTree.getSkbioTree().ascii_art())	
                        print('haplotype tree before:')	
                        print(haplotypeTree.getSkbioTree().ascii_art())	

                    if ancestral:
                        # print('ancestral')
                        eventIndex = eventIndex + 1
                        for node in newHaplotypeTree.getSkbioTree().traverse():
                            node.name = node.name + '_locus' + str(level) + '_event' + str(eventIndex)
                        geneNode = None
                        # A bug here in the extreme case: cannot find geneNodeName in self
                        for node in haplotypeTree.getSkbioTree().traverse():
                            if geneNodeName == node.name:
                                geneNode = node
                                break
                            else:
                                continue
                        if geneNode == None:
                            # print(haplotypeTree)
                            # print(geneNodeName)
                            pass
                        else:
                            # geneNode = haplotypeTree.getSkbioTree().find(geneNodeName)
                            geneNodeParent = geneNode.parent
                            # 1. create new node
                            newNode = skbio.tree.TreeNode()
                            newNode.name = 'd' + '_locus' + str(level) + '_event' + str(eventIndex)
                            # 2. change length
                            newHaplotypeTree.getSkbioTree().root().length = \
                                event['eventHeight'] - newHaplotypeTree.getTreeHeight()
                            if geneNode.id:
                                distanceAboveGeneNode = \
                                    event['eventHeight'] - haplotypeTree.getDistanceToLeaf(geneNode.id)
                            else:
                                distanceAboveGeneNode = event['eventHeight']

                            if not geneNodeParent:
                                newNode.length = 0  
                                geneNode.length = distanceAboveGeneNode
                            else:
                                newNode.length = geneNode.length - distanceAboveGeneNode
                                while newNode.length < 0:
                                    distanceAboveGeneNode = distanceAboveGeneNode - geneNode.length
                                    geneNode = geneNode.parent
                                    if geneNode.parent and geneNode.parent != -1:
                                        geneNodeParent = geneNode.parent                                   
                                    else:
                                        geneNodeParent = None
                                        newNode.length = 0 
                                        geneNode.length = distanceAboveGeneNode
                                        break
                                    newNode.length = geneNode.length - distanceAboveGeneNode       
                                geneNode.length = distanceAboveGeneNode
                            
                            # 3. change children
                            newNode.children = [geneNode, newHaplotypeTree.getSkbioTree()]
                            if geneNodeParent:
                                for i in range(len(geneNode.parent.children)):
                                    if geneNode.parent.children[i] == geneNode:
                                        geneNode.parent.children[i] = newNode
                                        break
                            # 4. change parent
                            geneNode.parent = newNode
                            newHaplotypeTree.getSkbioTree().root().parent = newNode
                            # check root
                            if not geneNodeParent:
                                haplotypeTree.setSkbioTree(newNode)
                            else:
                                newNode.parent = geneNodeParent

                            verbose = self.parameters['v']
                            if verbose:
                                print('haplotype tree after:')	
                                print(haplotypeTree.getSkbioTree().ascii_art())
                    
                    else:
                        # print('non-ancestral')
                        # non-ancestral
                        unlinked_d_number_total += 1
                        geneNodeName, distanceAboveGeneNode, branchLength = \
                            self.coalescentJoining(
                            event=event, haplotypeTree=haplotypeTree)
                        if geneNodeName:
                            # print(2)
                            eventIndex = eventIndex + 1
                            dic[level].append(eventIndex)
                            # print(dic)
                            unlinked_d_number_survived += 1 
                            # print(str(level) + '_event' + str(eventIndex))
                            
                            for node in newHaplotypeTree.getSkbioTree().traverse():
                                node.name = node.name + '_locus' + str(level) + '_event' + str(eventIndex)
                            geneNode = haplotypeTree.getSkbioTree().find(geneNodeName)
                            geneNodeParent = geneNode.parent
                            # 1. create new node
                            newNode = skbio.tree.TreeNode()
                            newNode.name = 'd' + '_locus' + str(level) + '_event' + str(eventIndex)
                            # 2. change length
                            newHaplotypeTree.getSkbioTree().root().length = \
                                event['eventHeight'] - newHaplotypeTree.getTreeHeight() + branchLength
                            if not geneNodeParent or geneNodeParent == -1:
                                newNode.length = 0 
                                geneNode.length = distanceAboveGeneNode
                            else:
                                newNode.length = geneNode.length - distanceAboveGeneNode
                                while newNode.length < 0:
                                    distanceAboveGeneNode = distanceAboveGeneNode - geneNode.length
                                    geneNode = geneNode.parent
                                    if geneNode.parent and geneNode.parent != -1:
                                        geneNodeParent = geneNode.parent                                   
                                    else:
                                        geneNodeParent = None
                                        newNode.length = 0 
                                        geneNode.length = distanceAboveGeneNode
                                        break
                                    newNode.length = geneNode.length - distanceAboveGeneNode       
                                geneNode.length = distanceAboveGeneNode
                            # 3. change children
                            newNode.children = [geneNode, newHaplotypeTree.getSkbioTree().root()]
                            if geneNodeParent:
                                for i in range(len(geneNode.parent.children)):
                                    if geneNode.parent.children[i] == geneNode:
                                        geneNode.parent.children[i] = newNode
                                        break
                            # 4. change parent
                            geneNode.parent = newNode
                            newHaplotypeTree.getSkbioTree().root().parent = newNode

                            # check root
                            if not geneNodeParent:
                                haplotypeTree.setSkbioTree(newNode)
                            else:
                                newNode.parent = geneNodeParent
                            verbose = self.parameters['v']
                            if verbose:
                                print('haplotype tree after:')	
                                print(haplotypeTree.getSkbioTree().ascii_art())
            elif (event['type'] == 'transfer'):
                transferTargetId = event['targetSpeciesId']
                targetHeight = self.speciesTree.getDistanceToLeaf(transferTargetId)
                distanceAboveTarget = event['eventHeight'] - targetHeight
                newHaplotypeTree, chosenGeneName, geneNodeName, ancestral, completeCount, incompleteCount= \
                    self.__addNewLociRecurse(haplotypeTree=haplotypeTree,
                    event=event, newLocusRootId=transferTargetId,
                    distanceAboveRoot=distanceAboveTarget, 
                    level=level, completeCount=completeCount,
                    incompleteCount=incompleteCount)
                
                if newHaplotypeTree:

                    verbose = self.parameters['v']
                    if verbose:
                        print(newHaplotypeTree)
                        print('new haplotype tree:')	
                        print(newHaplotypeTree.getSkbioTree().ascii_art())	
                        print('haplotype tree before:')	
                        print(haplotypeTree.getSkbioTree().ascii_art())	

                    # find gene node using coal with recom
                    geneNodeName, distanceAboveGeneNode, branchLength = self.coalescentJoining(
                        event=event, haplotypeTree=haplotypeTree)
                    if geneNodeName:
                        eventIndex = eventIndex + 1
                        for node in newHaplotypeTree.getSkbioTree().traverse():
                            node.name = node.name + '_locus' + str(level) + '_event' + str(eventIndex)
                        geneNode = haplotypeTree.getSkbioTree().find(geneNodeName)
                        geneNodeParent = geneNode.parent
                        # 1. create new node
                        newNode = skbio.tree.TreeNode()
                        newNode.name = 't' + '_locus' + str(level) + '_event' + str(eventIndex)
                        # 2. change length
                        newHaplotypeTree.getSkbioTree().root().length = \
                            event['eventHeight'] - newHaplotypeTree.getTreeHeight() + branchLength
                        if not geneNodeParent or geneNodeParent == -1:
                            newNode.length = 0 
                            geneNode.length = distanceAboveGeneNode
                        else:
                            newNode.length = geneNode.length - distanceAboveGeneNode
                            while newNode.length < 0:
                                distanceAboveGeneNode = distanceAboveGeneNode - geneNode.length
                                geneNode = geneNode.parent
                                if geneNode.parent and geneNode.parent != -1:
                                    geneNodeParent = geneNode.parent                                   
                                else:
                                    geneNodeParent = None
                                    newNode.length = 0 
                                    geneNode.length = distanceAboveGeneNode
                                    break
                                newNode.length = geneNode.length - distanceAboveGeneNode       
                            geneNode.length = distanceAboveGeneNode
                        # 3. change children
                        newNode.children = [geneNode, newHaplotypeTree.getSkbioTree()]
                        if geneNodeParent:
                            for i in range(len(geneNode.parent.children)):
                                if geneNode.parent.children[i] == geneNode:
                                    geneNode.parent.children[i] = newNode
                                    break
                        # 4. change parent
                        geneNode.parent = newNode
                        newHaplotypeTree.getSkbioTree().root().parent = newNode

                        # check root
                        if not geneNodeParent:
                            haplotypeTree.setSkbioTree(newNode)
                        else:
                            newNode.parent = geneNodeParent
                        verbose = self.parameters['v']
                        if verbose:
                            print('haplotype tree after:')	
                            print(haplotypeTree.getSkbioTree().ascii_art())
        return haplotypeTree, completeCount, incompleteCount, unlinked_d_number_total, unlinked_d_number_survived

    def __addNewLociRecurse(self, haplotypeTree, event, newLocusRootId, distanceAboveRoot, level, completeCount, incompleteCount, unlinked_d_number_total, unlinked_d_number_survived):
        copiedFullProcess = self.fullCoalescentProcess

        if (event['type'] == 'duplication' or event['type'] == 'transfer'): 
            speciesSkbioTree = self.speciesTree.getSkbioTree()
            newLocusRootName = self.speciesTree.getNodeById(newLocusRootId).name

            newLocusSkbioTree = speciesSkbioTree.find(newLocusRootName).deepcopy()
            newLocusTreeNames = [node.name for node in newLocusSkbioTree.traverse()]
            newLocusTreeNodes = [node for node in self.speciesTree.getNodes() 
                if node.name in newLocusTreeNames]
            # if len(newLocusTreeNodes) == 1 and event['type'] == 'transfer':
            #     print('find a transfer target at leaf!')
            newLocusTree = LocusTree(randomState=self.randomState)
            newLocusTree.initialize(nodes=newLocusTreeNodes, skbioTree=newLocusSkbioTree)
            newLocusTree.coalescentRate = self.speciesTree.coalescentRate
            newLocusTree.recombinationRate = self.speciesTree.recombinationRate
            newLocusTree.distanceAboveRoot = distanceAboveRoot     
            
            unlinked = event['unlinked']
            fullCoalescentProcess = None
            selectedCoalescentProcess = None
            chosenGeneName = None
            geneNodeName = None
            if unlinked:
                ancestral = False
                hemiplasy = self.__parameters['h']
                if hemiplasy == 1:
                    fullCoalescentProcess, selectedCoalescentProcess, chosenGeneName, incomplete = \
                        newLocusTree.incompleteCoalescent(distanceAboveRoot)
                    completeCount += 1
                    if incomplete:
                        incompleteCount += 1
                elif hemiplasy == 0:
                    selectedCoalescentProcess = \
                        newLocusTree.boundedCoalescent(distanceAboveRoot)
                    fullCoalescentProcess = selectedCoalescentProcess
            else:
                copiedSelectedProcess = self.coalescentProcess
                copiedRootProcess = copiedSelectedProcess[newLocusRootId]
                copiedRootGene = []
                accumulatedDistance = 0
                for e in copiedRootProcess:
                    accumulatedDistance += e['distance']
                    if accumulatedDistance < distanceAboveRoot:
                        continue
                    else:
                        copiedRootGene = e['fromSet']
                fullCoalescentProcess, selectedCoalescentProcess, chosenGeneName, geneNodeName, ancestral, incomplete = \
                    newLocusTree.linkedCoalescent(
                    copiedHaplotypeTree=copiedFullProcess, 
                    copiedRootGene=copiedRootGene,
                    distanceAboveRoot=distanceAboveRoot)
                completeCount += 1
                if incomplete:
                    incompleteCount += 1
                if not selectedCoalescentProcess:
                    return None, None, None, ancestral, completeCount, incompleteCount
            newHaplotypeTree = HaplotypeTree(
                randomState=self.randomState, 
                speciesTree=self.speciesTree, 
                locusTree=newLocusTree)
            newHaplotypeTree.initialize(
                locusTree=newLocusTree, 
                coalescentProcess=selectedCoalescentProcess, 
                fullCoalescentProcess=fullCoalescentProcess,
                rename=False)
            newHaplotypeTree.parameters = self.parameters
            rootLength = event['eventHeight'] - newHaplotypeTree.getTreeHeight()
            newEvents = self.DTLprocess(
                locusTree=newLocusTree, 
                haplotypeTree=newHaplotypeTree, 
                event=event,
                rootLength=rootLength)
            _, completeCount, incompleteCount, unlinked_d_number_total, unlinked_d_number_survived = \
                newHaplotypeTree.addNewLoci(events=newEvents, haplotypeTree=newHaplotypeTree, 
                level=level+1, completeCount=completeCount, incompleteCount=incompleteCount, 
                unlinked_d_number_total=unlinked_d_number_total, unlinked_d_number_survived=unlinked_d_number_survived)
            return newHaplotypeTree, chosenGeneName, geneNodeName, ancestral, completeCount, incompleteCount, unlinked_d_number_total, unlinked_d_number_survived

    """
    unused functions
    """
    def find_ils(self, path):
        for i in range(len(self.nodes)):
            # searching in backward order
            j = len(self.nodes)-1-i 
            geneNode = self.getNodeById(j)
            clades = geneNode.name.split('*')[:-1]
            clades = [int(j) for j in clades]
            geneNode.clades = clades

            if (geneNode.children and not geneNode.splits):
                for i in range(len(geneNode.children)):
                    geneChildName = self.getNodeById(geneNode.children[i]).name 
                    childClades = geneChildName.split('*')[:-1]
                    childClades = [int(j) for j in childClades]
                    geneNode.splits.append(childClades)

            foundSpeciesNode = False
            # find lowest possibile species node, then trace back
            for speciesNode in self.speciesTree.getNodes: 
                if (set(node.clade).issuperset(set(gene_clade))):
                    foundSpeciesNode = True
                    species_node = node
                    break
            if not foundSpeciesNode:
                continue
            species_clade = species_node.clade
            species_splits = species_node.clade_split
            find_ils = False
            if (gene_splits):
                gene_split_0 = set(gene_splits[0])
                gene_split_1 = set(gene_splits[1])
                for species_split in species_splits:
                    if (set(species_split).intersection(gene_split_0) and set(species_split).intersection(gene_split_1)):
                        find_ils = True
                        break
            if (find_ils):
                species_id = self.map_gene_id_to_species_id(j)
                index = Utility.increment()
                self.full_events.append({
                    'type': 'ils',
                    'gene_node_id': j, 
                    'gene_node_name': node.name, 
                    'species_node_id': species_id,
                    'index': index
                })
                # Debug.event_count['I'] += 1
                file_name = 'ils_' + str(index)
                f = open(os.path.join(path, file_name), 'w')
                f.write(str(gene_node.name) + ',' + str(gene_split_0) + ' ' + str(gene_split_1))
                f.close()
            else: 
                if (gene_splits):
                    species_id = self.map_gene_id_to_species_id(j)
                    index = Utility.increment()
                    self.full_events.append({
                        'type': 's',
                        'gene_node_id': j, 
                        'gene_node_name': node.name, 
                        'species_node_id': species_id,
                        'index': index
                    })
                    file_name = 's_' + str(index)
                    f = open(os.path.join(path, file_name), 'w')
                    f.write(str(gene_node.name) + ',' + str(gene_split_0) + ' ' + str(gene_split_1))
                    f.close()

    def __findTransferTarget(self, eventHeight, geneId):
        speciesNodes = self.speciesTree.getNodes()
        originSpeciesId = self.__mapGeneIdToSpeciesId(geneId=geneId, eventHeight=eventHeight)
        nodesList = []
        for node in speciesNodes:
            if node.id == originSpeciesId:
                continue
            if node.id == self.speciesTree.getRoot().id:
                continue
            parentHeight = self.speciesTree.getDistanceToLeaf(
                self.speciesTree.getNodeById(node.parent).id)
            if parentHeight > eventHeight:
                nodeHeight = self.speciesTree.getDistanceToLeaf(node.id)
                if nodeHeight <= eventHeight:
                    nodesList.append(node.id)
        return self.randomState.choice(nodesList), originSpeciesId  