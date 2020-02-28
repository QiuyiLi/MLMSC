import skbio
import numpy as np
from collections import defaultdict
from statistics import mean
from .tree_table import *
from .locus_tree import *


class HaplotypeTree:
    """
    Nodes are represented in Tree Table which is introduced in tree_table.py
    """
    def __init__(self, randomState, speciesTree, locusTree):
        self.__randomState = randomState
        self.__speciesTree = speciesTree
        self.__locusTree = locusTree

        self.__coalescentProcess = None
        self.__fullCoalescentProcess = None
        self.__treeTable = None
        self.__parameters = {}
        # self.__unlinkProb = None
        # self.__hemiplasy = None
        # self.__verbose = None

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
    
    # @property
    # def hemiplasy(self):
    #     return self.__hemiplasy

    # @property
    # def verbose(self):
    #     return self.__verbose

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
        self.__coalescentProcess = coalescentProcess
        self.__fullCoalescentProcess = coalescentProcess
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

    def setEventRates(self, duplicationPrmt, transferPrmt, lossPrmt, unlinkProb, hemiplasy, verbose):
        self.__parameters['d'] = duplicationPrmt
        self.__parameters['t'] = transferPrmt
        self.__parameters['l'] = lossPrmt
        self.__parameters['u'] = unlinkProb
        self.__parameters['h'] = hemiplasy
        self.__parameters['v'] = verbose


    # def setUnlinkProb(self, unlinkProb):
    #     self.__unlinkProb = unlinkProb

    # def setHemiplasy(self, hemiplasy):
    #     self.__hemiplasy = hemiplasy

    # detailed output for debugging
    # def setVerbose(self, verbose):
    #     self.__verbose = verbose

    def createSkbioTree(self, timeSequences, geneNodeName):
        """
        creat a tree structure consistent with the package skbio
        using timeSequence
        """
        skbioTree = skbio.tree.TreeNode()
        tempTimeSequences = timeSequences.copy()
        
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
            if (nodeName.count('*') == 1 and nodeName[0] == str(leaf)):
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
            return
        # two nodes case (trivial)
        elif len(skbioTree.name) == 4:
            childLName = skbioTree.name[:2]
            childRName = skbioTree.name[2:]
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
                        childRName = self.__starReplace(
                            skbioTree.name, prevPair[0])
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

    def Dprocess(self, distanceAboveRoot, threshold=float('inf'), event=None):
        events = []
        # trivial case
        if len(self.getNodes()) == 1:
            distanceAboveRoot = event['eventHeight']
        self.__DprocessRecurse(
            skbioTreeNode=self.getSkbioTree(), 
            distanceAboveRoot=distanceAboveRoot, 
            threshold=threshold, events=events)
        return events

    def __DprocessRecurse(self, skbioTreeNode, distanceAboveRoot, threshold, events):
        unlinkedProb = self.__parameters['u']
        node = self.getNodeByName(skbioTreeNode.name)
        duplicationRate = self.__parameters['d']
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
            unlinked = self.randomState.binomial(1, unlinkedProb)
            if unlinked == 1:
                unlinked = True
            else:
                unlinked = False
            if eventHeight < threshold:
                events.append({
                    'type': 'duplication',
                    'eventHeight': eventHeight,
                    'speciesNodeId': speciesId,
                    'distanceToSpeciesNode': distanceAboveSpeciesNode,
                    'unlinked': unlinked
                })
            # looking for more events on the same branch
            self.__DprocessRecurse(skbioTreeNode=skbioTreeNode, 
                distanceAboveRoot=distanceAboveRoot - distanceD, 
                threshold=threshold, events=events)
        elif (node.children):     
            # look for events in the children branches
            childL = skbioTreeNode.children[0]
            childR = skbioTreeNode.children[1]
            distanceToChildL = node.distanceToChildren[0]
            distanceToChildR = node.distanceToChildren[1]
            self.__DprocessRecurse(
                skbioTreeNode=childL, 
                distanceAboveRoot=distanceToChildL, threshold=threshold, events=events)
            self.__DprocessRecurse(
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
                originalSpeciesId, targetSpeciesId = self.__findTransferOrigin(
                    eventHeight=eventHeight, geneId=node.id)
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
            speciesId = self.__mapGeneIdToSpeciesId(geneId=geneId)
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
    def __mapGeneIdToSpeciesId(self, geneId):
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
        return speciesId 

    def __findTransferOrigin(self, eventHeight, geneId):    
        targetSpeciesId = self.__mapGeneIdToSpeciesId(geneId=geneId)
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

    def binom(self, n, k):
        return math.factorial(n) // math.factorial(k) // math.factorial(n - k)

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
                    if len(e['fromSet']) >= 2:
                        coalescentRate = self.binom(len(e['fromSet']), 2)*unitCoalescentRate
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
                if len(e['fromSet']) >= 2:
                        coalescentRate = self.binom(len(e['fromSet']), 2)*unitCoalescentRate
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
            and haplotypeTree.speciesTree.getNodeById(parentId)
                in haplotypeTree.locusTree.getNodes()):       
            return self.__coalescentJoiningRecurse(nodeId=parentId, event=event, 
                        fullCoalescentProcess=fullCoalescentProcess, 
                        fullCoalescentTree=fullCoalescentTree,
                        selectedCoalescentProcess=selectedCoalescentProcess, 
                        haplotypeTree=haplotypeTree,  
                        branchLength=branchLength, passed = True)
        else:
            # no coalescent
            return None, None, None

    def addNewLoci(self, events, haplotypeTree, level):
        """
        1. simulate all the events on the haplotype tree 
        2. construct the corresponding new locus tree
        3. use incomplete coalescence to generate the new haplotype tree
        4. simulate all the events on the haplotype tree 
        5. recurse
        """
        selectedCoalescentProcess = haplotypeTree.coalescentProcess
        fullCoalescentProcess = haplotypeTree.fullCoalescentProcess
        eventIndex = -1
        
        for event in events:
            if (event['type'] == 'loss'):
                geneNodeName = event['geneNodeName']
                geneNode = haplotypeTree.getSkbioTree().find(geneNodeName)
                geneNode.name = geneNode.name + '_loss'

            elif (event['type'] == 'duplication'):
                speciesId = event['speciesNodeId']
                distanceAboveSpeciesNode = event['distanceToSpeciesNode']

                newHaplotypeTree, chosenGeneName, geneNodeName, ancestral = \
                    self.__addNewLociRecurse(
                    event=event, newLocusRootId=speciesId,
                    distanceAboveRoot=distanceAboveSpeciesNode, level=level)

                if newHaplotypeTree:
                    verbose = self.parameters['v']
                    if verbose:
                        print(newHaplotypeTree)
                        print('new haplotype tree:')	
                        print(newHaplotypeTree.getSkbioTree().ascii_art())	
                        print('haplotype tree before:')	
                        print(haplotypeTree.getSkbioTree().ascii_art())	

                    if ancestral:
                        eventIndex = eventIndex + 1
                        for node in newHaplotypeTree.getSkbioTree().traverse():
                            node.name = node.name + '_lv=' + str(level) + '_id=' + str(eventIndex)
                        geneNode = haplotypeTree.getSkbioTree().find(geneNodeName)
                        geneNodeParent = geneNode.parent
                        # 1. create new node
                        newNode = skbio.tree.TreeNode()
                        newNode.name = 'd' + '_lv=' + str(level) + '_id=' + str(eventIndex)
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
                        # non-ancestral
                        geneNodeName, distanceAboveGeneNode, branchLength = \
                            self.coalescentJoining(
                            event=event, haplotypeTree=haplotypeTree)
                        if geneNodeName:
                            eventIndex = eventIndex + 1
                            
                            for node in newHaplotypeTree.getSkbioTree().traverse():
                                node.name = node.name + '_lv=' + str(level) + '_id=' + str(eventIndex)
                            geneNode = haplotypeTree.getSkbioTree().find(geneNodeName)
                            geneNodeParent = geneNode.parent
                            # 1. create new node
                            newNode = skbio.tree.TreeNode()
                            newNode.name = 'd' + '_lv=' + str(level) + '_id=' + str(eventIndex)
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

            # no need for this right now
            elif (event['type'] == 'transfer'):
                transferTargetId = event['targetSpeciesId']
                targetHeight = self.speciesTree.getDistanceToLeaf(transferTargetId)
                distanceAboveTarget = event['eventHeight'] - targetHeight
                newHaplotypeTree, chosenGeneName, geneNodeName, ancestral = self.__addNewLociRecurse(
                    event=event, newLocusRootId=transferTargetId,
                    distanceAboveRoot=distanceAboveTarget, level=level)
                
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
                            node.name = node.name + '_lv=' + str(level) + '_id=' + str(eventIndex)
                        geneNode = haplotypeTree.getSkbioTree().find(geneNodeName)
                        geneNodeParent = geneNode.parent
                        # 1. create new node
                        newNode = skbio.tree.TreeNode()
                        newNode.name = 't' + '_lv=' + str(level) + '_id=' + str(eventIndex)
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
        return haplotypeTree

    def __addNewLociRecurse(self, event, newLocusRootId, distanceAboveRoot, level):
        copiedFullProcess = self.fullCoalescentProcess

        if (event['type'] == 'duplication' or event['type'] == 'transfer'): 
            # for transfer nodeId = target_id
            speciesSkbioTree = self.speciesTree.getSkbioTree()
            newLocusRootName = self.speciesTree.getNodeById(newLocusRootId).name

            newLocusSkbioTree = speciesSkbioTree.find(newLocusRootName).deepcopy()
            newLocusTreeNames = [node.name for node in newLocusSkbioTree.traverse()]
            newLocusTreeNodes = [node for node in self.speciesTree.getNodes() 
                if node.name in newLocusTreeNames]
            newLocusTree = LocusTree(randomState=self.randomState)
            newLocusTree.initialize(nodes=newLocusTreeNodes, skbioTree=newLocusSkbioTree)
            newLocusTree.coalescentRate = self.speciesTree.coalescentRate
            newLocusTree.recombinationRate = self.speciesTree.recombinationRate
            
            unlinked = event['unlinked']
            fullCoalescentProcess = None
            selectedCoalescentProcess = None
            chosenGeneName = None
            geneNodeName = None
            if unlinked:
                ancestral = False
                hemiplasy = self.__parameters['h']
                if hemiplasy == 1:
                    fullCoalescentProcess, selectedCoalescentProcess, chosenGeneName = \
                        newLocusTree.incompleteCoalescent(distanceAboveRoot)
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
                fullCoalescentProcess, selectedCoalescentProcess, chosenGeneName, geneNodeName, ancestral = \
                    newLocusTree.linkedCoalescent(
                    copiedHaplotypeTree=copiedFullProcess, 
                    copiedRootGene=copiedRootGene,
                    distanceAboveRoot=distanceAboveRoot)
                if not selectedCoalescentProcess:
                    return None, None, None, ancestral

            newHaplotypeTree = HaplotypeTree(
                randomState=self.randomState, speciesTree=self.speciesTree, locusTree=newLocusTree)
            newHaplotypeTree.initialize(
                locusTree=newLocusTree, 
                coalescentProcess=selectedCoalescentProcess, 
                fullCoalescentProcess=fullCoalescentProcess,
                rename=False)
            newHaplotypeTree.parameters = self.parameters
            rootLength = event['eventHeight'] - newHaplotypeTree.getTreeHeight()

            newHaplotypeTreeEvents = newHaplotypeTree.Lprocess(
                event=event, distanceAboveRoot=rootLength)

            # coalescentTreeProcess = None
            coalescentTreeProcessD, _ = newLocusTree.coalescent(
                distanceAboveRoot=float('inf'))
            coalescentTreeD = HaplotypeTree(
                randomState=self.randomState, speciesTree=self.speciesTree, locusTree=newLocusTree)
            coalescentTreeD.initialize(
                locusTree=newLocusTree, 
                coalescentProcess=coalescentTreeProcessD, 
                fullCoalescentProcess=coalescentTreeProcessD,
                rename=False)
            coalescentTreeD.parameters = self.parameters
            rootLength = max(0, event['eventHeight'] - coalescentTreeD.getTreeHeight())
            threshold = event['eventHeight']
            newCoalescentTreeEventsD = coalescentTreeD.Dprocess(
                distanceAboveRoot=rootLength, threshold=threshold, event=event)

            coalescentTreeProcessT, _ = self.speciesTree.coalescent(
                distanceAboveRoot=float('inf'))
            coalescentTreeT = HaplotypeTree(
                randomState=self.randomState, speciesTree=self.speciesTree, locusTree=newLocusTree)
            coalescentTreeT.initialize(
                locusTree=newLocusTree, 
                coalescentProcess=coalescentTreeProcessT, 
                fullCoalescentProcess=coalescentTreeProcessT,
                rename=False,
                event = 'transfer')
            coalescentTreeT.parameters = self.parameters
            rootLength = max(0, event['eventHeight'] - coalescentTreeT.getTreeHeight())
            threshold = event['eventHeight']
            newCoalescentTreeEventsT = coalescentTreeT.Tprocess(
                distanceAboveRoot=rootLength, threshold=threshold, event=event)

            # newEvents = newHaplotypeTreeEvents + newCoalescentTreeEvents
            # newEvents.sort(reverse=True, key=lambda x: x['eventHeight'])
            newCoalescentTreeEvents = newCoalescentTreeEventsD + newCoalescentTreeEventsT
            newCoalescentTreeEvents.sort(reverse=True, key=lambda x: x['eventHeight'])
            newEvents = newCoalescentTreeEvents + newHaplotypeTreeEvents

            newHaplotypeTree.addNewLoci(events=newEvents, haplotypeTree=newHaplotypeTree, 
                level=level + 1)
            return newHaplotypeTree, chosenGeneName, geneNodeName, ancestral

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
                # print('find ils at gene node ' + str(gene_node.name) + ' split: ' + str(gene_split_0) + ' ' 
                #         + str(gene_split_1) + ' ' + 'species_id: ' + str(species_id))
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
                    # print('find speciation at gene node ' + str(gene_node.name) + ' split: ' + str(gene_split_0) + ' ' 
                    #         + str(gene_split_1) + ' ' + 'species_id: ' + str(species_id))

    """
    unused function
    """
    def __findTransferTarget(self, eventHeight, geneId):
        speciesNodes = self.speciesTree.getNodes()
        originSpeciesId = self.__mapGeneIdToSpeciesId(geneId=geneId)
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