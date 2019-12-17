import numpy as np
from .species_tree import *
from .haplotype_tree import *
from .exception import *


class IxDTLModel:

    # def __init__(self, seed=0):
    #     self.__randomState = np.random.RandomState(seed)
    def __init__(self):
        self.__randomState = np.random.RandomState()

        self.__speciesTree = None
        self.__haplotypeTree = None
        self.__locusTrees = []
        self.__parameters = {}

    @property
    def speciesTree(self):
        return self.__speciesTree

    @property
    def haplotypeTree(self):
        return self.__haplotypeTree

    @property
    def locusTrees(self):
        return self.__locusTrees

    @property
    def parameters(self):
        return self.__parameters

    @property
    def randomState(self):
        return self.__randomState

    def run(self, inputFile, coalescentArgs, duplicationArgs, transferArgs, 
        lossArgs, hemiplasy, recombination, verbose):
        # set parameters
        self.setParameters(
            coalescent=coalescentArgs, 
            duplication=duplicationArgs,
            transfer=transferArgs, 
            loss=lossArgs, 
            hemiplasy=hemiplasy,
            recombination=recombination,
            verbose=verbose)

        # read a species tree from input file
        self.readSpeciesTree(inputFile)

        # construct the original haplotype tree according to the species tree
        self.constructOriginalHaplotypeTree()

        # run dtl process
        events = self.haplotypeTree.dtlProcess(distanceAboveRoot=0)
        events.sort(reverse=True, key=lambda x: x['eventHeight'])

        # run dt subtree
        geneTree = self.haplotypeTree.dtSubtree(
            coalescentProcess=self.haplotypeTree.coalescentProcess, 
            events=events, haplotypeTree=self.haplotypeTree, level=0)
        geneSkbioTree = geneTree.getSkbioTree()

        # cut the tree 
        geneTreeTruncated = geneTree
        geneSkbioTreeTruncated = geneSkbioTree.deepcopy()
        for node in geneSkbioTreeTruncated.traverse():
            if (node.children 
                and 'loss' in node.children[0].name 
                and 'loss' in node.children[1].name):
                geneSkbioTreeTruncated.remove_deleted(
                    lambda x: x.name == node.name)
        geneSkbioTreeTruncated.prune()
        for node in geneSkbioTreeTruncated.traverse():
            if 'loss' in node.name:
                geneSkbioTreeTruncated.remove_deleted(
                    lambda x: x.name == node.name)
        geneSkbioTreeTruncated.prune()

        if not geneSkbioTreeTruncated:
                print('Exception: ALL LOST')
                return
            
        if self.__parameters['verbose']:
            # visualizing the untruncated tree
            print('untruncated tree:')
            print(geneSkbioTree.ascii_art())	    
            # check time consistency 
            for node in geneSkbioTree.tips():	
                print(str(geneSkbioTree.distance(node)) + ' ' + str(node.name))
            # visualizing the truncated tree
            print('truncated tree:')
            print(geneSkbioTreeTruncated.ascii_art())
            # check time consistency
            for node in geneSkbioTreeTruncated.tips():	
                print(str(geneSkbioTreeTruncated.distance(node)) + ' ' + str(node.name))
            print(geneSkbioTreeTruncated.ascii_art())
            # final gene table
            geneTreeTruncated.readFromSkbioTree(skbioTree=geneSkbioTreeTruncated, rename=False)
            print(geneTreeTruncated)
                
        # save newick to file
        f = open('./output/gene_tree_full.newick','w')
        f.write(str(geneSkbioTree))
        f.close()

        f = open('./output/gene_tree_truncated.newick','w')
        f.write(str(geneSkbioTreeTruncated))
        f.close()

    def setParameters(self, coalescent, duplication, transfer, loss, 
        hemiplasy, recombination, verbose):
        if not coalescent:
            raise IxDTLError('missing coalescent parameter')
        self.__parameters['coalescent'] = coalescent

        if not duplication:
            raise IxDTLError('missing duplication parameter')
        self.__parameters['duplication'] = duplication

        if not transfer:
            raise IxDTLError('missing transfer parameter')
        self.__parameters['transfer'] = transfer

        if not loss:
            raise IxDTLError('missing loss parameter')
        self.__parameters['loss'] = loss

        if hemiplasy is None:
            raise IxDTLError('missing hemiplasy option')
        self.__parameters['hemiplasy'] = hemiplasy

        if recombination is None:
            raise IxDTLError('missing recombination option')
        self.__parameters['recombination'] = recombination

        if verbose is None:
            raise IxDTLError('missing verbose option')
        self.__parameters['verbose'] = verbose

    def readSpeciesTree(self, path):
        self.__speciesTree = SpeciesTree(randomState=self.randomState)

        self.speciesTree.initialize(path=path)

        self.speciesTree.setCoalescentRate(
            coalescentPrmt=self.__parameters['coalescent'])
        
        if self.__parameters['verbose']:
            print('species tree:')	
            print(self.speciesTree)	
            print(self.speciesTree.getSkbioTree().ascii_art())	
            print()
            
    def constructOriginalHaplotypeTree(self):
        self.__haplotypeTree = HaplotypeTree(
            randomState=self.randomState, speciesTree=self.speciesTree)

        self.haplotypeTree.initialize(locusTree=self.speciesTree)

        if self.__parameters['verbose']:
            print('original haplotype tree:')	
            print(self.haplotypeTree)	
            print(self.haplotypeTree.getSkbioTree().ascii_art())	
            print()

        self.haplotypeTree.setEventRates(
            duplicationPrmt=self.parameters['duplication'],
            transferPrmt=self.parameters['transfer'],
            lossPrmt=self.parameters['loss'])
        self.haplotypeTree.setRecombination(
            recombination=self.parameters['recombination'])
        self.haplotypeTree.setHemiplasy(
            hemiplasy=self.parameters['hemiplasy'])
        self.haplotypeTree.setVerbose(
            verbose=self.parameters['verbose'])
        
