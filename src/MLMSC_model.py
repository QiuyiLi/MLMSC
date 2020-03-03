import numpy as np
from .species_tree import *
from .locus_tree import *
from .haplotype_tree import *
from .exception import *

class MLMSC_Model:
    def __init__(self, seed=None):
        if seed == None:
            self.__randomState = np.random.RandomState()
        else:
            self.__randomState = np.random.RandomState(seed)
        self.__speciesTree = None
        # self.__haplotypeTree = None
        # self.__locusTree = None
        self.__parameters = {}

    @property
    def speciesTree(self):
        return self.__speciesTree

    @property
    def parameters(self):
        return self.__parameters

    @property
    def randomState(self):
        return self.__randomState

    def run(self, inputFile, seedArgs, coalescentArgs, recombinationArgs, duplicationArgs, transferArgs, 
        lossArgs, unlinkArgs, repeatNumber, hemiplasy, verbose):
        # set parameters
        self.setParameters(
            coalescent=coalescentArgs, 
            recombination=recombinationArgs,
            duplication=duplicationArgs,
            transfer=transferArgs, 
            loss=lossArgs, 
            unlink=unlinkArgs,
            repeat=repeatNumber,
            hemiplasy=hemiplasy,
            verbose=verbose)

        # read a species tree from input file
        self.readSpeciesTree(inputFile)

        f = open('./output/gene_tree_untruncated.newick','w')
        f.write('')
        f.close()

        f = open('./output/gene_tree_truncated.newick','w')
        f.write('')
        f.close()

        f = open('./output/gene_tree.newick','w')
        f.write('')
        f.close()

        for i in range(repeatNumber):
            if repeatNumber > 1:
                print('Tree ' + str(i+1) + ' of ' + str(repeatNumber) + ':')
            # construct the original haplotype tree according to the species tree
            originalLocusTree = self.constructOriginalLocusTree()
            originalHaplotypeTree = self.constructOriginalHaplotypeTree()

            events = originalHaplotypeTree.DTLprocess(
                locusTree=originalLocusTree, 
                haplotypeTree=originalHaplotypeTree, 
                initial=True)

            # add new loci
            geneTree = originalHaplotypeTree.addNewLoci(events=events, 
                haplotypeTree=originalHaplotypeTree, level=0)
            geneSkbioTree = geneTree.getSkbioTree()

            # cut the tree at losses
            geneTreeTruncated = geneTree
            geneSkbioTreeTruncated = geneSkbioTree.deepcopy()
            geneSkbioTreeTruncated = self.cutTree(geneSkbioTreeTruncated)
            for node in geneSkbioTreeTruncated.traverse():
                if 'loss' in node.name:
                    geneSkbioTreeTruncated.remove_deleted(
                        lambda x: x.name == node.name)
            geneSkbioTreeTruncated.prune()
            if not geneSkbioTreeTruncated:
                print('Exception: ALL LOST')
                # save newick to file
                f = open('./output/gene_tree_untruncated.newick','a')
                f.write(str(geneSkbioTree))
                f.close()

                f = open('./output/gene_tree_truncated.newick','a')
                f.write('')
                f.close()

                f = open('./output/gene_tree.newick','a')
                f.write('')
                f.close()
            else:
                geneSkbioTreeCleaned = geneSkbioTreeTruncated.deepcopy()
                for node in geneSkbioTreeCleaned.traverse():
                    if node in geneSkbioTreeCleaned.tips():
                        speciesId = int(node.name.split('*')[0])
                        remainder = node.name.split('*')[1]
                        speciesNode = self.speciesTree.getNodeById(speciesId)
                        node.name = speciesNode.name + remainder
                    else:
                        node.name = ''

                # visualization
                if self.__parameters['verbose']:
                    # visualizing the untruncated tree
                    print('untruncated tree:')
                    print(geneSkbioTree.ascii_art())	    
                    # check time consistency 
                    print('distances from tips to root:')
                    for node in geneSkbioTree.tips():	
                        print(str(geneSkbioTree.distance(node)) + ' ' + str(node.name))
                    # visualizing the truncated tree
                    print('truncated tree:')
                    print(geneSkbioTreeTruncated.ascii_art())
                    # check time consistency
                    print('distances from tips to root:')
                    for node in geneSkbioTreeTruncated.tips():	
                        print(str(geneSkbioTreeTruncated.distance(node)) + ' ' + str(node.name))
                    print(geneSkbioTreeTruncated.ascii_art())
                    # final gene table
                    print('gene tree table:')
                    geneTreeTruncated.readFromSkbioTree(skbioTree=geneSkbioTreeTruncated, rename=False)
                    print(geneTreeTruncated)
                    print('finished.')

                    # save newick to file
                    f = open('./output/gene_tree_untruncated.newick','a')
                    f.write(str(geneSkbioTree))
                    f.close()

                    f = open('./output/gene_tree_truncated.newick','a')
                    f.write(str(geneSkbioTreeTruncated))
                    f.close()

                    f = open('./output/gene_tree.newick','a')
                    string = str(geneSkbioTreeCleaned)
                    for char in string:
                        if char == "'":
                            continue
                        else:
                            f.write(char)
                    f.close()   
                else:
                    if repeatNumber == 1:
                        print('gene tree:')
                        print(geneSkbioTreeCleaned.ascii_art())
                    else:
                        print('finished.')
        
                    # save newick to file
                    f = open('./output/gene_tree_untruncated.newick','a')
                    f.write(str(geneSkbioTree))
                    f.close()

                    f = open('./output/gene_tree_truncated.newick','a')
                    f.write(str(geneSkbioTreeTruncated))
                    f.close()

                    f = open('./output/gene_tree.newick','a')
                    string = str(geneSkbioTreeCleaned)
                    for char in string:
                        if char == "'":
                            continue
                        else:
                            f.write(char)
                    f.close()
        
    def cutTree(self, untruncatedGeneTree):
        root = untruncatedGeneTree.root()
        self.cutTreeRecurse(root, untruncatedGeneTree)
        return untruncatedGeneTree
    
    def cutTreeRecurse(self, node, untruncatedGeneTree):
        if node.children:
            if 'loss' in node.children[0].name:
                if 'loss' in node.children[1].name:
                    findIt = 1
                else:
                    findIt = self.cutTreeRecurse(node.children[1], untruncatedGeneTree)
            else:
                if 'loss' in node.children[1].name:
                    findIt = self.cutTreeRecurse(node.children[0], untruncatedGeneTree)
                else:
                    findIt = self.cutTreeRecurse(node.children[0], untruncatedGeneTree)*\
                        self.cutTreeRecurse(node.children[1], untruncatedGeneTree)
            if findIt:
                node.name += '_loss'
            return findIt
        else:
            if 'loss' in node.name:
                return 1
            else:
                return 0

    def setParameters(self, coalescent, recombination, duplication, transfer, loss, unlink,
        repeat, hemiplasy, verbose):
        if coalescent is None:
            raise MLMSC_Error('missing coalescent parameter')
        self.__parameters['coalescent'] = coalescent

        if recombination is None:
            raise MLMSC_Error('missing recombination parameter')
        self.__parameters['recombination'] = recombination

        if duplication is None:
            raise MLMSC_Error('missing duplication parameter')
        self.__parameters['duplication'] = duplication

        if transfer is None:
            raise MLMSC_Error('missing transfer parameter')
        self.__parameters['transfer'] = transfer

        if loss is None:
            raise MLMSC_Error('missing loss parameter')
        self.__parameters['loss'] = loss

        if unlink is None:
            raise MLMSC_Error('missing unlink parameter')
        self.__parameters['unlink'] = unlink

        if repeat is None:
            raise MLMSC_Error('missing repeat times')
        self.__parameters['repeat'] = repeat

        if hemiplasy is None:
            raise MLMSC_Error('missing hemiplasy option')
        self.__parameters['hemiplasy'] = hemiplasy

        if verbose is None:
            raise MLMSC_Error('missing verbose option')
        self.__parameters['verbose'] = verbose

    def readSpeciesTree(self, path):
        self.__speciesTree = SpeciesTree(randomState=self.randomState)

        self.speciesTree.initialize(path=path)

        self.speciesTree.setCoalescentRate(
            coalescentPrmt=self.parameters['coalescent'])
        
        self.speciesTree.setRecombinationRate(
            recombinationPrmt=self.parameters['recombination'])

        if self.__parameters['verbose']:
            print('species tree:')	
            print(self.speciesTree)	
            print(self.speciesTree.getSkbioTree().ascii_art())	

    def constructOriginalLocusTree(self):
        locusTree = LocusTree(randomState=self.randomState)
        locusTree.initialize(
            nodes=self.speciesTree.getNodes(), 
            skbioTree=self.speciesTree.getSkbioTree())
        locusTree.coalescentRate = self.speciesTree.coalescentRate
        locusTree.recombinationRate = self.speciesTree.recombinationRate
        return locusTree
            
    def constructOriginalHaplotypeTree(self):
        haplotypeTree = HaplotypeTree(
            randomState=self.randomState, 
            speciesTree=self.speciesTree, 
            locusTree=self.speciesTree)
        haplotypeTree.initialize(locusTree=self.speciesTree)

        if self.__parameters['verbose']:
            print('original haplotype tree:')	
            print(haplotypeTree)	
            print(haplotypeTree.getSkbioTree().ascii_art())	

        haplotypeTree.setEventRates(
            duplicationPrmt=self.parameters['duplication'],
            transferPrmt=self.parameters['transfer'],
            lossPrmt=self.parameters['loss'],
            unlinkProb=self.parameters['unlink'],
            hemiplasy=self.parameters['hemiplasy'],
            verbose=self.parameters['verbose'])
        return haplotypeTree
        
