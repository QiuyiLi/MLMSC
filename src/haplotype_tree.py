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

    def __init__(self, randomState, speciesTree):
        self.__randomState = randomState
        self.__speciesTree = speciesTree

        self.__coalescentProcess = None
        self.__treeTable = None
        self.__eventRates = {}
        self.__recombination = None
        self.__hemiplasy = None
        self.__verbose = None

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
    def coalescentProcess(self):
        return self.__coalescentProcess

    @property
    def treeTable(self):
        return self.__treeTable

    @property
    def eventRates(self):
        return self.__eventRates
    @eventRates.setter
    def eventRates(self, eventRates):
        self.__eventRates = eventRates

    @property
    def recombination(self):
        return self.__recombination
    
    @property
    def hemiplasy(self):
        return self.__hemiplasy

    @property
    def verbose(self):
        return self.__verbose

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

    def getDistanceToLeaf(self, nodeId, branchDistance):
        return self.__treeTable.distanceToLeaf(nodeId, branchDistance)

    def initialize(self, locusTree, coalescentProcess=None, rename=True):
        if not coalescentProcess:
            # ordinary coalesent for constructing the original haplotype tree
            coalescentProcess, cladeSetIntoRoot = locusTree.coalescent(
                distanceAboveRoot=float('inf'))
        self.__coalescentProcess = coalescentProcess

        timeSequences = locusTree.getTimeSequences(
            coalescentProcess=coalescentProcess)
        # print('time sequences:')	
        # print(timeSequences)	
        # print()

        skbioTree = self.createSkbioTree(timeSequences)
        self.readFromSkbioTree(skbioTree, rename)

        for geneNode in self.getNodes():
            clades = geneNode.name.split('*')[:-1]
            clades = [int(j) for j in clades]
            geneNode.clades = clades
            if (geneNode.children and not geneNode.splits):
                for i in range(len(geneNode.children)):
                    geneChildName = self.getNodeById(geneNode.children[i]).name 
                    childClades = geneChildName.split('*')[:-1]
                    childClades = [int(j) for j in childClades]
                    geneNode.splits.append(childClades)

    def setEventRates(self, duplicationPrmt, transferPrmt, lossPrmt):
        if ('const' not in duplicationPrmt):
            self.__eventRates['d'] = self.randomState.gamma(
                shape=duplicationPrmt['shape'], scale=duplicationPrmt['scale'],
                size=len(self.getLeaves()))
        else:
            self.__eventRates['d'] = np.repeat(duplicationPrmt['const'], 
                len(self.getLeaves()))
        if ('const' not in transferPrmt):
            self.__eventRates['t'] = self.randomState.gamma(
                shape=transferPrmt['shape'], scale=transferPrmt['scale'],
                size=len(self.getLeaves()))
        else:
            self.__eventRates['t'] = np.repeat(transferPrmt['const'], 
                len(self.getLeaves()))
        if ('const' not in lossPrmt):
            self.__eventRates['l'] = self.randomState.gamma(
                shape=lossPrmt['shape'], scale=lossPrmt['scale'],
                size=len(self.getLeaves()))
        else:
            self.__eventRates['l'] = np.repeat(lossPrmt['const'], 
                len(self.getLeaves()))

    # If true, then new locus is indpendent of the original locus
    # (locus tree model, IxDTL model)
    # otherwise the new locus is linked with the original locus
    # (IDTL model)
    def setRecombination(self, recombination):
        self.__recombination = recombination
    
    # If true: IxDTL nodel
    # otherwise： locus tree model
    def setHemiplasy(self, hemiplasy):
        self.__hemiplasy = hemiplasy

    # detailed output for debugging
    def setVerbose(self, verbose):
        self.__verbose = verbose

    def createSkbioTree(self, timeSequences):
        """
        creat a tree structure consistent with the package skbio
        using timeSequence
        """
        skbioTree = skbio.tree.TreeNode()
        tempTimeSequences = timeSequences.copy()
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
        # otherwise 
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

    def dlProcess(self, distanceAboveRoot, event=None):
        events = []
        # trivial case
        if len(self.getNodes()) == 1:
            distanceAboveRoot = event['distanceToGeneNode']
        self.__dlProcessRecurse(
            skbioTreeNode=self.getSkbioTree(), distanceAboveRoot=distanceAboveRoot, events=events)
        return events

    def __getEventRateInAncestralBranch(self, eventType, clade):
        indices = []
        splited = clade.split('*')[:-1]
        for index in splited:
            indices.append(int(index))
        return mean(self.eventRates[eventType][indices])

    def __dlProcessRecurse(self, skbioTreeNode, distanceAboveRoot, events):
        node = self.getNodeByName(skbioTreeNode.name)

        distanceD = self.randomState.exponential(
            scale=1.0 / self.__getEventRateInAncestralBranch(
                eventType='d', clade=node.name))
        distanceL = self.randomState.exponential(
            scale=1.0 / self.__getEventRateInAncestralBranch(
                eventType='l', clade=node.name))

        # duplication happens first
        if (distanceD < distanceL and distanceD < distanceAboveRoot):
            eventHeight = self.getDistanceToLeaf(node.id, 0) + distanceAboveRoot - distanceD
            speciesId, distanceAboveSpeciesNode = self.__mapEventToSpeciesTree(
                geneId=node.id, eventHeight=eventHeight, speciesId=None)
            events.append({
                'type': 'duplication',
                'geneNodeId': node.id,      # closest gene node to the event from below
                'geneNodeName': node.name, 
                'distanceToGeneNode': distanceAboveRoot - distanceD,
                'eventHeight': eventHeight,
                'speciesNodeId': speciesId,
                'distanceToSpeciesNode': distanceAboveSpeciesNode,
                'index': -1
            })
            # looking for more events on the same branch
            self.__dlProcessRecurse(
                skbioTreeNode=skbioTreeNode, 
                distanceAboveRoot=distanceAboveRoot - distanceD, events=events)
        elif (distanceL <= distanceD and distanceL < distanceAboveRoot):      
            # loss happens first, the seaching process stops at the loss point
            eventHeight = self.getDistanceToLeaf(node.id, 0) + distanceAboveRoot - distanceL
            speciesId, distanceAboveSpeciesNode = self.__mapEventToSpeciesTree(
                geneId=node.id, eventHeight=eventHeight, speciesId=None)
            events.append({
                'type': 'loss',
                'geneNodeId': node.id,      # closest gene node to the event from below
                'geneNodeName': node.name, 
                'distanceToGeneNode': distanceAboveRoot - distanceL,
                'eventHeight': eventHeight,
                'speciesNodeId': speciesId,     # closest species node to the event from below
                'distanceToSpeciesNode': distanceAboveSpeciesNode,
                'index': -1
            })
        else:
            # reach the end the current branch, looking for events in the 2 children branches
            if (node.children):     # if children branches exist
                childL = skbioTreeNode.children[0]
                childR = skbioTreeNode.children[1]
                distanceToChildL = node.distanceToChildren[0]
                distanceToChildR = node.distanceToChildren[1]
                self.__dlProcessRecurse(
                    skbioTreeNode=childL, 
                    distanceAboveRoot=distanceToChildL, events=events)
                self.__dlProcessRecurse(
                    skbioTreeNode=childR, 
                    distanceAboveRoot=distanceToChildR, events=events)
            # else: if not exist, reach the leaves of the tree, searching process stops


    def dtProcess(self, distanceAboveRoot, event=None):
        events = []
        # trivial case
        if len(self.getNodes()) == 1:
            distanceAboveRoot = event['distanceToGeneNode']
        self.__dtProcessRecurse(
            skbioTreeNode=self.getSkbioTree(), distanceAboveRoot=distanceAboveRoot, events=events)
        return events

    def __dtProcessRecurse(self, skbioTreeNode, distanceAboveRoot, events):
        threshold = self.getTreeHeight() + distanceAboveRoot
        node = self.getNodeByName(skbioTreeNode.name)

        distanceD = self.randomState.exponential(
            scale=1.0 / self.__getEventRateInAncestralBranch(
                eventType='d', clade=node.name))
        distanceT = self.randomState.exponential(
            scale=1.0 / self.__getEventRateInAncestralBranch(
                eventType='l', clade=node.name))

        # duplication happens first
        if (distanceD < distanceT and distanceD < distanceAboveRoot):
            eventHeight = self.getDistanceToLeaf(node.id, 0) + distanceAboveRoot - distanceD
            speciesId, distanceAboveSpeciesNode = self.__mapEventToSpeciesTree(
                geneId=node.id, eventHeight=eventHeight, speciesId=None)
            if eventHeight > threshold:
                events.append({
                    'type': 'duplication',
                    'geneNodeId': -1,      # closest gene node to the event from below
                    'geneNodeName': -1, 
                    'distanceToGeneNode': -1,
                    'eventHeight': eventHeight,
                    'speciesNodeId': speciesId,
                    'distanceToSpeciesNode': distanceAboveSpeciesNode,
                    'index': -1
                })
            # looking for more events on the same branch
            self.__dtProcessRecurse(
                skbioTreeNode=skbioTreeNode, 
                distanceAboveRoot=distanceAboveRoot - distanceD, events=events)
        elif (distanceT <= distanceD and distanceT < distanceAboveRoot):
            eventHeight = self.getDistanceToLeaf(node.id, 0) + distanceAboveRoot - distanceT
            speciesTreeHeight = self.speciesTree.getTreeHeight()
            if eventHeight < speciesTreeHeight:
                target, originalSpeciesId = self.__findTransferTarget(
                    eventHeight=eventHeight, geneId=node.id)
                speciesId, distanceAboveSpeciesNode = self.__mapEventToSpeciesTree(
                    geneId=node.id, eventHeight=eventHeight, speciesId=None)
                if target and eventHeight > threshold:
                    events.append({
                        'type': 'transfer',
                        'geneNodeId': -1,      # closest gene node to the event from below
                        'geneNodeName': -1, 
                        'distanceToGeneNode': -1,
                        'targetSpeciesId': target,
                        'eventHeight': eventHeight,
                        'speciesNodeId': speciesId,      # closest species node to the event from below
                        'distanceToSpeciesNode': distanceAboveSpeciesNode,
                        'index': -1
                    })
            self.__dtProcessRecurse(
                skbioTreeNode=skbioTreeNode, 
                distanceAboveRoot=distanceAboveRoot - distanceT, events=events)
        else:
            # reach the end the current branch, looking for events in the 2 children branches
            if (node.children):     # if children branches exist
                childL = skbioTreeNode.children[0]
                childR = skbioTreeNode.children[1]
                distanceToChildL = node.distanceToChildren[0]
                distanceToChildR = node.distanceToChildren[1]
                self.__dtProcessRecurse(
                    skbioTreeNode=childL, 
                    distanceAboveRoot=distanceToChildL, events=events)
                self.__dtProcessRecurse(
                    skbioTreeNode=childR, 
                    distanceAboveRoot=distanceToChildR, events=events)
            # else: if not exist, reach the leaves of the tree, searching process stops

    # M function: map the event point to where it occurs in the species tree
    def __mapEventToSpeciesTree(self, geneId, eventHeight, speciesId=None):
        if speciesId == None:
            speciesId = self.__mapGeneIdToSpeciesId(geneId=geneId)
        distanceAboveSpeciesNode = eventHeight - self.speciesTree.getDistanceToLeaf(speciesId, 0)
        if speciesId == self.speciesTree.getRoot().id:
            return speciesId, distanceAboveSpeciesNode
        else:
            speciesIdParent = self.speciesTree.getNodeById(speciesId).parent
            speciesDistanceParent = self.speciesTree.getDistanceToLeaf(speciesIdParent, 0)
            if speciesDistanceParent > eventHeight:
                return speciesId, distanceAboveSpeciesNode
            else:
                return self.__mapEventToSpeciesTree(
                    geneId=geneId, eventHeight=eventHeight, speciesId=speciesIdParent)

    # map the gene node to the species node where the coalesent happens to give birth to itself
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
                self.speciesTree.getNodeById(node.parent).id, 0)
            if parentHeight > eventHeight:
                nodeHeight = self.speciesTree.getDistanceToLeaf(node.id, 0)
                if nodeHeight <= eventHeight:
                    nodesList.append(node.id)
        return self.randomState.choice(nodesList), originSpeciesId

    # def find_ils(self, path):
    #     for i in range(len(self.nodes)):
    #         # searching in backward order
    #         j = len(self.nodes)-1-i 
    #         geneNode = self.getNodeById(j)
    #         clades = geneNode.name.split('*')[:-1]
    #         clades = [int(j) for j in clades]
    #         geneNode.clades = clades

    #         if (geneNode.children and not geneNode.splits):
    #             for i in range(len(geneNode.children)):
    #                 geneChildName = self.getNodeById(geneNode.children[i]).name 
    #                 childClades = geneChildName.split('*')[:-1]
    #                 childClades = [int(j) for j in childClades]
    #                 geneNode.splits.append(childClades)

    #         foundSpeciesNode = False
    #         # find lowest possibile species node, then trace back
    #         for speciesNode in self.speciesTree.getNodes: 
    #             if (set(node.clade).issuperset(set(gene_clade))):
    #                 foundSpeciesNode = True
    #                 species_node = node
    #                 break
    #         if not foundSpeciesNode:
    #             continue
    #         species_clade = species_node.clade
    #         species_splits = species_node.clade_split
    #         find_ils = False
    #         if (gene_splits):
    #             gene_split_0 = set(gene_splits[0])
    #             gene_split_1 = set(gene_splits[1])
    #             for species_split in species_splits:
    #                 if (set(species_split).intersection(gene_split_0) and set(species_split).intersection(gene_split_1)):
    #                     find_ils = True
    #                     break
    #         if (find_ils):
    #             species_id = self.map_gene_id_to_species_id(j)
    #             index = Utility.increment()
    #             self.full_events.append({
    #                 'type': 'ils',
    #                 'gene_node_id': j, 
    #                 'gene_node_name': node.name, 
    #                 'species_node_id': species_id,
    #                 'index': index
    #             })
    #             # Debug.event_count['I'] += 1
    #             file_name = 'ils_' + str(index)
    #             f = open(os.path.join(path, file_name), 'w')
    #             f.write(str(gene_node.name) + ',' + str(gene_split_0) + ' ' + str(gene_split_1))
    #             f.close()
    #             # print('find ils at gene node ' + str(gene_node.name) + ' split: ' + str(gene_split_0) + ' ' 
    #             #         + str(gene_split_1) + ' ' + 'species_id: ' + str(species_id))
    #         else: 
    #             if (gene_splits):
    #                 species_id = self.map_gene_id_to_species_id(j)
    #                 index = Utility.increment()
    #                 self.full_events.append({
    #                     'type': 's',
    #                     'gene_node_id': j, 
    #                     'gene_node_name': node.name, 
    #                     'species_node_id': species_id,
    #                     'index': index
    #                 })
    #                 file_name = 's_' + str(index)
    #                 f = open(os.path.join(path, file_name), 'w')
    #                 f.write(str(gene_node.name) + ',' + str(gene_split_0) + ' ' + str(gene_split_1))
    #                 f.close()
    #                 # print('find speciation at gene node ' + str(gene_node.name) + ' split: ' + str(gene_split_0) + ' ' 
    #                 #         + str(gene_split_1) + ' ' + 'species_id: ' + str(species_id))

    def cutTree(self, coalescentProcess, unlinkpoints):
        for unlinkpoint in unlinkpoints:
            nodeId = unlinkpoint['speciesNodeId']
            unlinkNode = self.getNodeById(nodeId)
            unlinkDistance = unlinkpoint['distanceToSpeciesNode']
            coalDistance = 0
            for mergingSet in coalescentProcess[nodeId]:
                coalDistance = coalDistance + mergingSet['distance']
                if coalDistance < unlinkDistance:
                    coalescentProcess[nodeId].pop(0)
                else:
                    break
            coalescentProcess = self.__cutTreeRecurse(
                coalescentProcess=coalescentProcess, unlinkNode=unlinkNode, 
                unlinkDistance=unlinkDistance)
        return coalescentProcess

    def __cutTreeRecurse(self, coalescentProcess, unlinkNode, unlinkDistance):
        if unlinkNode.children: 
            children = unlinkNode.children
            for child in children:
                coalescentProcess[child] = []
                self.__cutTreeRecurse(coalescentProcess, child, unlinkDistance)
        return coalescentProcess

    def jointCoalescent(self, copiedHaplotypeTree, unlinkedSpecies):
        
        return

    def dtSubtree(self, coalescentProcess, events, haplotypeTree, level):
        """
        1. simulate all the events on the haplotype tree 
        2. construct the corresponding new locus tree
        3. use incomplete coalescence to generate the new haplotype tree
        4. simulate all the events on the haplotype tree 
        5. recurse
        """
        eventIndex = -1
        for event in events:
            if (event['type'] == 'duplication'):
                eventIndex = eventIndex + 1
                speciesId, distanceAboveSpeciesNode = self.__mapEventToSpeciesTree(
                    geneId=event['geneNodeId'], eventHeight=event['eventHeight'], 
                    speciesId=None)
                newHaplotypeTree = self.__dtSubtreeRecurse(
                    event=event, newLocusRootId=speciesId, 
                    distanceAboveRoot=distanceAboveSpeciesNode, level=level)

                for node in newHaplotypeTree.getSkbioTree().traverse():
                    node.name = node.name + '_lv=' + str(level) + '_id=' + str(eventIndex)
                
                if self.verbose:
                    print(newHaplotypeTree)
                    print('new haplotype tree:')	
                    print(newHaplotypeTree.getSkbioTree().ascii_art())	
                    print('haplotype tree before:')	
                    print(haplotypeTree.getSkbioTree().ascii_art())	

                geneNodeName = event['geneNodeName']
                geneNode = haplotypeTree.getSkbioTree().find(geneNodeName)
                geneNodeParent = geneNode.parent
                # 1. create new node
                newNode = skbio.tree.TreeNode()
                newNode.name = 'd' + '_lv=' + str(level) + '_id=' + str(eventIndex)
                # 2. change length
                newHaplotypeTree.getSkbioTree().length = event['eventHeight'] - newHaplotypeTree.getTreeHeight()
                newNode.length = 0 if not geneNodeParent else geneNode.length - event['distanceToGeneNode']
                geneNode.length = event['distanceToGeneNode']
                # 3. change children
                newNode.children = [geneNode, newHaplotypeTree.getSkbioTree()]
                if geneNodeParent:
                    for i in range(len(geneNode.parent.children)):
                        if geneNode.parent.children[i] == geneNode:
                            geneNode.parent.children[i] = newNode
                            break
                # 4. change parent
                geneNode.parent = newNode
                newHaplotypeTree.getSkbioTree().parent = newNode

                # check root
                if not geneNodeParent:
                    haplotypeTree.setSkbioTree(newNode)
                else:
                    newNode.parent = geneNodeParent

                if self.verbose:
                    print('haplotype tree after:')	
                    print(haplotypeTree.getSkbioTree().ascii_art())

            elif (event['type'] == 'transfer'):
                eventIndex = eventIndex + 1
                transferTargetId = event['targetSpeciesId']
                targetHeight = self.speciesTree.getDistanceToLeaf(transferTargetId, 0)
                distanceAboveTarget = event['eventHeight'] - targetHeight
                newHaplotypeTree = self.__dtSubtreeRecurse(
                    event=event, newLocusRootId=transferTargetId, 
                    distanceAboveRoot=distanceAboveTarget, level=level)
                
                for node in newHaplotypeTree.getSkbioTree().traverse():
                    node.name = node.name + '_lv=' + str(level) + '_id=' + str(eventIndex)

                if self.verbose:
                    print(newHaplotypeTree)
                    print('new haplotype tree:')	
                    print(newHaplotypeTree.getSkbioTree().ascii_art())	
                    print('haplotype tree before:')	
                    print(haplotypeTree.getSkbioTree().ascii_art())	
                
                geneNodeName = event['geneNodeName']
                geneNode = haplotypeTree.getSkbioTree().find(geneNodeName)
                geneNodeParent = geneNode.parent
                # 1. create new node
                newNode = skbio.tree.TreeNode()
                newNode.name = 't' + '_lv=' + str(level) + '_id=' + str(eventIndex)
                # 2. change length
                newHaplotypeTree.getSkbioTree().length = event['eventHeight'] - newHaplotypeTree.getTreeHeight()
                newNode.length = 0 if not geneNodeParent else geneNode.length - event['distanceToGeneNode']
                geneNode.length = event['distanceToGeneNode']
                # 3. change children
                newNode.children = [geneNode, newHaplotypeTree.getSkbioTree()]
                if geneNodeParent:
                    for i in range(len(geneNode.parent.children)):
                        if geneNode.parent.children[i] == geneNode:
                            geneNode.parent.children[i] = newNode
                            break
                # 4. change parent
                geneNode.parent = newNode
                newHaplotypeTree.getSkbioTree().parent = newNode
                
                # check root
                if not geneNodeParent:
                    haplotypeTree.setSkbioTree(newNode)
                else:
                    newNode.parent = geneNodeParent

                if self.verbose:
                    print('haplotype tree after:')	
                    print(haplotypeTree.getSkbioTree().ascii_art())

            elif (event['type'] == 'loss'):
                geneNodeName = event['geneNodeName']
                geneNode = haplotypeTree.getSkbioTree().find(geneNodeName)
                geneNode.name = geneNode.name + '_loss'

                # cut tree bug here...
                # do not cut it for now, lable the loss points and delete all of them all at once when everything done
                # haplotypeSkbioTree = haplotypeTree.getSkbioTree()
                # haplotypeSkbioTree.remove_deleted(
                #     lambda x: x.name == event['geneNodeName'])
                # haplotypeSkbioTree.prune()

        return haplotypeTree

    def __dtSubtreeRecurse(self, event, newLocusRootId, distanceAboveRoot, level):
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
            # print('newLocusTree', newLocusTree)	
            newLocusTree.coalescentRate = self.speciesTree.coalescentRate

            locusTreeCoalescentProcess = None
            if self.hemiplasy == 1:
                locusTreeCoalescentProcess, chosenGeneName = \
                    newLocusTree.incompleteCoalescent(distanceAboveRoot)
            elif self.hemiplasy == 0:
                locusTreeCoalescentProcess = \
                    newLocusTree.boundedCoalescent(distanceAboveRoot)

            newHaplotypeTree = HaplotypeTree(
                randomState=self.randomState, speciesTree=self.speciesTree)
            newHaplotypeTree.initialize(
                locusTree=newLocusTree, 
                coalescentProcess=locusTreeCoalescentProcess, rename=False)
            newHaplotypeTree.eventRates = self.eventRates

            rootLength = event['eventHeight'] - newHaplotypeTree.getTreeHeight()
            newHaplotypeTree.getSkbioTree().length = rootLength
            newHaplotypeTreeEvents = newHaplotypeTree.dlProcess(
                event=event, distanceAboveRoot=rootLength)

            coalescentTreeProcess = None
            coalescentTreeProcess, cladeSetIntoRoot = newLocusTree.coalescent(
                distanceAboveRoot=float('inf'))
            coalescentTree = HaplotypeTree(
                randomState=self.randomState, speciesTree=self.speciesTree)
            coalescentTree.initialize(
                locusTree=newLocusTree, 
                coalescentProcess=coalescentTreeProcess, rename=False)
            coalescentTree.eventRates = self.eventRates
            rootLength = event['eventHeight'] - coalescentTree.getTreeHeight()
            newCoalescentTreeEvents = coalescentTree.dtProcess(
                event=event, distanceAboveRoot=rootLength)

            newEvents = newHaplotypeTreeEvents + newCoalescentTreeEvents
            newEvents.sort(reverse=True, key=lambda x: x['eventHeight'])

            newHaplotypeTree.dtSubtree(
                coalescentProcess=locusTreeCoalescentProcess, 
                events=newEvents, haplotypeTree=newHaplotypeTree, 
                level=level + 1)
            return newHaplotypeTree