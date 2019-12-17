from .species_tree import *
from .tree_table import *

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
        coalescentProcess, genesIntoRoot = self.coalescent(distanceAboveRoot)
        chosenGene = np.random.choice(genesIntoRoot)
        selectedCoalescentProcess = self.__selectCoalescentProcess(
            coalescentProcess, chosenGene)
        return selectedCoalescentProcess, chosenGene
    
    # seclect one haplotype tree from multiple trees
    # sgenerated by the incomplete coalescent
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
                for clade in mergingSet['toSet']:
                    if self._starInSet(target=clade, clade=chosenGene):
                        toSet.append(clade)
                if toSet:
                    selectedCoalescentProcess[speciesNodeId].append({
                        'fromSet': fromSet, 
                        'toSet': toSet,
                        'distance': distance
                    })
        return selectedCoalescentProcess
