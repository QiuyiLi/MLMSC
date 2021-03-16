import os
import numpy as np
from .species_tree import *
from .locus_tree import *
from .haplotype_tree import *
from .exception import *
from collections import defaultdict
import itertools
import ete3

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
        # f = open(inputFile, 'r')
        # speciesTreeString = f.read()
        # f.close()
        # speciesTree_ete3 = self.TreeNode.read(StringIO(speciesTreeString))
        # print(self.speciesTree)
        d_name = str(duplicationArgs/0.001)
        while d_name[-1] == '0' or d_name[-1] == '.':
            d_name = d_name[:-1]
            if not d_name:
                d_name = '0'
                break

        l_name = str(lossArgs/0.001)
        while l_name[-1] == '0' or l_name[-1] == '.':
            l_name = l_name[:-1]
            if not l_name:
                l_name = '0'
                break
        
        name = 'N' + str(repeatNumber) + 'D' + d_name + 'L' + l_name + 'C' + str(round(1/coalescentArgs))
        # print(name)
        outputDir = './output_' + name
        if not os.path.exists(outputDir):
            os.makedirs(outputDir)

        files = os.listdir(outputDir)
        if 'gene_tree_untruncated.newick' not in files:
            f = open(outputDir + '/gene_tree_untruncated.newick','w')
            f.write('')
            f.close()

        if 'gene_tree_truncated.newick' not in files:
            f = open(outputDir + '/gene_tree_truncated.newick','w')
            f.write('')
            f.close()

        if 'gene_tree.newick' not in files:
            f = open(outputDir + '/gene_tree.newick','w')
            f.write('')
            f.close()

        # f = open(outputDir + '/average_tree_size.txt','w')
        # f.write('')
        # f.close()

        # f = open(outputDir + '/dup_numbers.txt','w')
        # f.write('')
        # f.close()

        # f = open(outputDir + '/full_name.newick','w')
        # f.write('')
        # f.close()
        
        if 'quartets_mlmsc.txt' not in files:
            f = open(outputDir + '/quartets_mlmsc.txt','w')
            f.write('n_qtree,p_0,p_1,n_d,n_genes,n_species\n')
            f.close()

        i = 1 
        nLeaves_single = 0
        nLeaves_full = 0
        nTree_single = 0
        nTree_full = 0
        sum_d_survived = 0

        while i <= repeatNumber:
            if repeatNumber > 1:
                if i%1000 == 0:
                    print('Tree ' + str(i) + ' of ' + str(repeatNumber) + ':')
            # the original locus tree is the same as the species tree
            originalLocusTree = self.constructOriginalLocusTree()
            # construct the original haplotype tree according to the species tree
            originalHaplotypeTree = self.constructOriginalHaplotypeTree()

            events = originalHaplotypeTree.DTLprocess(
                locusTree=originalLocusTree, 
                haplotypeTree=originalHaplotypeTree, 
                initial=True)

            # add new loci
            geneTree, completeCount, incompleteCount, unlinked_d_number_total, unlinked_d_number_survived = \
                originalHaplotypeTree.addNewLociShell(events=events, 
                haplotypeTree=originalHaplotypeTree, level=0, completeCount=1, incompleteCount=0, unlinked_d_number_total=0, unlinked_d_number_survived=0)
            geneSkbioTree = geneTree.getSkbioTree()
            # cut the tree at losses
            geneTreeTruncated = geneTree

            # for node in geneSkbioTree.traverse():
            #     if node not in geneSkbioTree.tips():
            #         node.name = ''
            # f = open(outputDir + '/full_name.newick','a')
            # string = str(geneSkbioTree)
            # for char in string:
            #     if char == "'":
            #         continue
            #     else:
            #         f.write(char)
            # f.close()

            survived_before_loss = []
            geneSkbioTreeTruncated = geneSkbioTree.deepcopy()
            for node in geneSkbioTreeTruncated.traverse():
                if node in geneSkbioTreeTruncated.tips():
                    speciesId = int(node.name.split('*')[0])
                    remainder = (node.name.split('*')[1])
                    if remainder:
                        remainder = remainder.split('_locus')[1:]
                        for e in remainder:
                            if e not in survived_before_loss:
                                survived_before_loss.append(e)

            geneSkbioTreeTruncated = self.cutTree(geneSkbioTreeTruncated)

            for node in geneSkbioTreeTruncated.traverse():
                if 'loss' in node.name:
                    geneSkbioTreeTruncated.remove_deleted(
                        lambda x: x.name == node.name)
            geneSkbioTreeTruncated.prune()
            if not geneSkbioTreeTruncated:
                # print('none')
                continue
            else:
                geneSkbioTreeCleaned = geneSkbioTreeTruncated.deepcopy()

                Vec = defaultdict(list)
                survived = []
                for node in geneSkbioTreeCleaned.traverse():
                    if node in geneSkbioTreeCleaned.tips():
                        speciesId = int(node.name.split('*')[0])
                        remainder = (node.name.split('*')[1])
                        if remainder:
                            remainder = remainder.split('_locus')[1:]
                            for e in remainder:
                                if e not in survived:
                                    survived.append(e)
                        speciesNode = self.speciesTree.getNodeById(speciesId)
                        node.name = speciesNode.name
                        if speciesNode.name == 'A':
                            node.name += str(len(Vec['A'])+1)
                            Vec['A'].append(node.name)
                        elif speciesNode.name == 'B':
                            node.name += str(len(Vec['B'])+1)
                            Vec['B'].append(node.name)
                        elif speciesNode.name == 'C':
                            node.name += str(len(Vec['C'])+1)
                            Vec['C'].append(node.name)
                        elif speciesNode.name == 'D':
                            node.name += str(len(Vec['D'])+1)
                            Vec['D'].append(node.name)
                        elif speciesNode.name == 'E':
                            node.name += str(len(Vec['E'])+1)
                            Vec['E'].append(node.name)
                        elif speciesNode.name == 'F':
                            node.name += str(len(Vec['F'])+1)
                            Vec['F'].append(node.name)
                        elif speciesNode.name == 'G':
                            node.name += str(len(Vec['G'])+1)
                            Vec['G'].append(node.name)
                        elif speciesNode.name == 'H':
                            node.name += str(len(Vec['H'])+1)
                            Vec['H'].append(node.name)
                        elif speciesNode.name == 'I':
                            node.name += str(len(Vec['I'])+1)
                            Vec['I'].append(node.name)
                        elif speciesNode.name == 'J':
                            node.name += str(len(Vec['J'])+1)
                            Vec['J'].append(node.name)
                        elif speciesNode.name == 'K':
                            node.name += str(len(Vec['K'])+1)
                            Vec['K'].append(node.name)
                        elif speciesNode.name == 'L':
                            node.name += str(len(Vec['L'])+1)
                            Vec['L'].append(node.name)
                        elif speciesNode.name == 'M':
                            node.name += str(len(Vec['M'])+1)
                            Vec['M'].append(node.name)
                        elif speciesNode.name == 'N':
                            node.name += str(len(Vec['N'])+1)
                            Vec['N'].append(node.name)
                        elif speciesNode.name == 'O':
                            node.name += str(len(Vec['O'])+1)
                            Vec['O'].append(node.name)
                        elif speciesNode.name == 'P':
                            node.name += str(len(Vec['P'])+1)
                            Vec['P'].append(node.name)
                    else:
                        node.name = ''
                # print(survived)
                # print(unlinked_d_number_total, unlinked_d_number_survived, len(survived))
                singleTree = None
                n_genes = 0
                names = []
                if Vec['A']:
                    nameA = self.randomState.choice(Vec['A'])
                    names.append(nameA)
                    n_genes += len(Vec['A'])
                if Vec['B']:
                    nameB = self.randomState.choice(Vec['B'])
                    names.append(nameB)
                    n_genes += len(Vec['B'])
                if Vec['C']:
                    nameC = self.randomState.choice(Vec['C'])
                    names.append(nameC)
                    n_genes += len(Vec['C'])
                if Vec['D']:
                    nameD = self.randomState.choice(Vec['D'])
                    names.append(nameD)
                    n_genes += len(Vec['D'])
                if Vec['E']:
                    nameE = self.randomState.choice(Vec['E'])
                    names.append(nameE)
                    n_genes += len(Vec['E'])
                if Vec['F']:
                    nameF = self.randomState.choice(Vec['F'])
                    names.append(nameF)
                    n_genes += len(Vec['F'])
                if Vec['G']:
                    nameG = self.randomState.choice(Vec['G'])
                    names.append(nameG)
                    n_genes += len(Vec['G'])
                if Vec['H']:
                    nameH = self.randomState.choice(Vec['H'])
                    names.append(nameH)
                    n_genes += len(Vec['H'])
                if Vec['I']:
                    nameI = self.randomState.choice(Vec['I'])
                    names.append(nameI)
                    n_genes += len(Vec['I'])
                if Vec['J']:
                    nameJ = self.randomState.choice(Vec['J'])
                    names.append(nameJ)
                    n_genes += len(Vec['J'])
                if Vec['K']:
                    nameK = self.randomState.choice(Vec['K'])
                    names.append(nameK)
                    n_genes += len(Vec['K'])
                if Vec['L']:
                    nameL = self.randomState.choice(Vec['L'])
                    names.append(nameL)
                    n_genes += len(Vec['L'])
                if Vec['M']:
                    nameM = self.randomState.choice(Vec['M'])
                    names.append(nameM)
                    n_genes += len(Vec['M'])
                if Vec['N']:
                    nameN = self.randomState.choice(Vec['N'])
                    names.append(nameN)
                    n_genes += len(Vec['N'])
                if Vec['O']:
                    nameO = self.randomState.choice(Vec['O'])
                    names.append(nameO)
                    n_genes += len(Vec['O'])
                if Vec['P']:
                    nameP = self.randomState.choice(Vec['P'])
                    names.append(nameP)
                    n_genes += len(Vec['P'])

                singleTree = geneSkbioTreeCleaned.deepcopy()
                singleTree = singleTree.shear(names)
                tipNumber = 0
                tipNumber_full = 0

                if len(names) >= 4:
                    quartets_single = itertools.combinations(names, 4)
                    count_single = 0
                    count_0 = 0
                    count_2 = 0
                    for element in quartets_single:
                        count_single += 1

                        geneSubtree = singleTree.shear(element)
                        for node in geneSubtree.traverse():
                            if node in geneSubtree.tips():
                                node.name = node.name[0]
                        for node in geneSubtree.tips():
                            name = node.name
                            geneSubtree_ete3 = ete3.Tree(str(geneSubtree))
                            geneSubtree_ete3.set_outgroup(name)
                            break
                        temp = []
                        for e in element:
                            temp.append(e[0])
                        speciesSubtree = self.speciesTree.getSkbioTree().shear(temp)
                        for node in speciesSubtree.traverse():
                            if node not in speciesSubtree.tips():
                                node.name = ''
                        speciesSubtree_ete3 = ete3.Tree(str(speciesSubtree))
                        speciesSubtree_ete3.set_outgroup(name)

                        rf_quartet = speciesSubtree_ete3.robinson_foulds(geneSubtree_ete3)[0]
                        if rf_quartet == 0:
                            count_0 += 1
                        elif rf_quartet == 2:
                            count_2 += 1
                    f = open(outputDir + '/quartets_mlmsc.txt','a')
                    f.write(str(count_single) + ',' + str(count_0/count_single) + ',' + str(count_2/count_single) + ',' + str(len(survived)) + ',' + str(n_genes) + ',' + str(len(names)) + '\n') 

                for node in singleTree.tips():
                    tipNumber += 1
                    node.name = node.name[0]
                for node in geneSkbioTreeCleaned.tips():
                    tipNumber_full += 1
                    node.name = node.name[0]

                if tipNumber < 4:
                    continue 
                else:
                    # f = open(outputDir + '/dup_numbers.txt','a')
                    # f.write(str(unlinked_d_number_total) + ',' + str(unlinked_d_number_survived) + ',' + str(len(survived_before_loss)) + ',' + str(len(survived)) + '\n')
                    # f.close()
                    nLeaves_single += tipNumber
                    nLeaves_full += tipNumber_full
                    sum_d_survived += len(survived)

                    species = []
                    for e in names:
                        species.append(e[0])
                    comb = itertools.combinations(species, 4)
                    count = 0
                    for element in comb:
                        num = 1
                        for e in element:
                            num = num * len(Vec[e])
                        count += num

                    nTree_full += count
                    nTree_single += len(list(itertools.combinations(names, 4)))

                    # visualization
                    i = i + 1
                    if self.__parameters['verbose']:
                        # visualizing the untruncated tree
                        print('untruncated tree:')
                        print(geneSkbioTree.ascii_art())	    
                        # check time consistency of the untruncated tree
                        print('distances from tips to root:')
                        for node in geneSkbioTree.tips():	
                            print(str(geneSkbioTree.distance(node)) + ' ' + str(node.name))
                        # visualizing the truncated tree
                        print('truncated tree:')
                        print(geneSkbioTreeTruncated.ascii_art())
                        # check time consistency of the truncated tree
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
                        f = open(outputDir + '/gene_tree_untruncated.newick','a')
                        f.write(str(geneSkbioTree))
                        f.close()

                        f = open(outputDir + '/gene_tree_truncated.newick','a')
                        f.write(str(geneSkbioTreeTruncated))
                        f.close()

                        f = open(outputDir + '/gene_tree.newick','a')
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
                            print('sheared tree:')
                            print(singleTree.ascii_art())
                        elif i%1000 == 0:
                            print('finished.')

                        # save newick to file
                        
                        f = open(outputDir + '/gene_tree_untruncated.newick','a')
                        string = str(singleTree)
                        for char in string:
                            if char == "'":
                                continue
                            else:
                                f.write(char)
                        f.close()

                        f = open(outputDir + '/gene_tree_truncated.newick','a')
                        string =  str(geneSkbioTreeTruncated)
                        for char in string:
                            if char == "'":
                                continue
                            else:
                                f.write(char)
                        f.close()

                        f = open(outputDir + '/gene_tree.newick','a')
                        string = str(geneSkbioTreeCleaned)
                        for char in string:
                            if char == "'":
                                continue
                            else:
                                f.write(char)
                        f.close()

        average_tree_size_single = nLeaves_single/repeatNumber
        average_tree_size_full = nLeaves_full/repeatNumber
        average_sum_d_survived = sum_d_survived/repeatNumber
        # f = open(outputDir + '/average_tree_size.txt','a')
        # f.write(str(average_tree_size_single) + ' ' + str(average_tree_size_full) + ' ' + str(nTree_single) + ' ' + str(nTree_full) + ' ' + str(average_sum_d_survived))
        # f.close()
        # print(average_tree_size_single, average_tree_size_full, nTree_single/repeatNumber, nTree_full/repeatNumber, average_sum_d_survived)

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
        # else:
        #     print('species tree:')	
        #     print(self.speciesTree.getSkbioTree().ascii_art())


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

        haplotypeTree.setParameters(
            duplicationPrmt=self.parameters['duplication'],
            transferPrmt=self.parameters['transfer'],
            lossPrmt=self.parameters['loss'],
            unlinkProb=self.parameters['unlink'],
            hemiplasy=self.parameters['hemiplasy'],
            verbose=self.parameters['verbose'])
        return haplotypeTree

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
        
