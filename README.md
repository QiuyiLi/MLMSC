# MLMSC Simulator

MLMSC Simulator is a program for the simulation of gene family evolution within a species tree based on the Multilocus Multipecies Coalescent (MLMSC) model. MLMSC model generalises the multispecies coalescent to gene families, and is designed to capture all possible scenarios that can arise through incomplete lineage sorting, gene duplication, transfer and loss, and any interaction between these processes. The MLMSC combines forward- and backward-in-time modelling in order to properly account for copy number hemiplasy and linkage between loci. 
The input for MLMSC Simulator are the simulation parameters and a pre-specified species tree in Newick format as an input file. The output is a simulated gene tree in Newick format.

## Citation
If you use the MLMSC simulator, please cite:

The Multilocus Multispecies Coalescent: A Flexible New Model of Gene Family Evolution; Qiuyi Li, Nicolas Galtier, Celine Scornavacca, Yao-Ban Chan; bioRxiv 2020.05.07.081836; doi: https://doi.org/10.1101/2020.05.07.081836.

## Obtaining the program

### Download

The source code of MLMSC Simulator is available in the GitHub repository:
```
https://github.com/QiuyiLi/MLMSC
```
You can clone the sources in your computer by executing:
```
git clone https://github.com/QiuyiLi/MLMSC.git
```

### Requirements
MLMSC Simulator is built under Python 3.7.2 and requirs the installation of the following packages: 
```
scikit-bio
numpy
statistics
```
You can install these packages by executing:
```
pip install -r requirements.txt
```
### Testing

The MLMSC Simulator should be ready to use now. You can test the the simulator by executing:
```
python3 MLMSC.py -i species_trees/species_tree_0.newick -s 20200210
```
The simulator will run with the default setting and you should be able to see the following outputs:
```
gene tree:
                              /-A
                    /--------|
                   |         |          /-B
          /--------|          \--------|
         |         |                    \-A_locus0_event2
         |         |
---------|          \-E_locus0_event1
         |
         |          /-E_locus0_event0
          \--------|
                    \-E_locus1_event0_locus0_event0
```

##  Usage

### Input file
* Species tree: -i or --inputFile; the prespecified species tree written in [Newick](https://en.wikipedia.org/wiki/Newick_format) format
```
e.g., python3 MLMSC.py -i species_trees/species_tree_0.newick
```

### Input parameters

* Coalescent rate: -c or --coalescentRate; the coalescent rate in multispecies coalescent, default: -c 1.0
```
e.g., python3 MLMSC.py -i species_trees/species_tree_0.newick -c 0.5
```
  
* Recombination rate: -r or --recombinationRate; the recombiantion rate in linked coalescent, default: -r 0.5
```
e.g., python3 MLMSC.py -i species_trees/species_tree_0.newick -r 0.3
```
  
* Duplication rate: -d or --duplicationRate; occurrence rate of gene duplication, default: -d 0.2
```
e.g., python3 MLMSC.py -i species_trees/species_tree_0.newick -d 0.1
```
  
* Transfer rate; -t or --transferRate; occurrence rate of horizontal gene transfer, default: -t 0.1
```
e.g., python3 MLMSC.py -i species_trees/species_tree_0.newick -t 0.2
```
  
* Loss rate; -l or --lossRate; occurrence rate of gene loss, default: -l 0.2
```
e.g., python3 MLMSC.py -i species_trees/species_tree_0.newick -l 0.1
```
  
* Unlink probability; -u or --unlinkProb; the probability that a duplication is unlinked, default: -u 0.5
```
e.g., python3 MLMSC.py -i species_trees/species_tree_0.newick -u 0.3
```
  
* Hemiplasy option: -h or --hemiplasy; whether or not the copy number hemiplasy is allowed, default: -h 1
```
e.g., python3 MLMSC.py -i species_trees/species_tree_0.newick -h 0
```
  
* Verbose option: -v or --verbose; detailed outputs for debugging purposes, default: -v 0
```
e.g., python3 MLMSC.py -i species_trees/species_tree_0.newick -v 1
```

* Number of repeats: -n or --numRepeats; number of gene trees simulated, defult: -n 1
```
e.g., python3 MLMSC.py -i species_trees/species_tree_0.newick -n 10
```

* Seed number: -s or --seed; set seed for reproduciblility, default -s None (random seed)
```
e.g., python3 MLMSC.py -i species_trees/species_tree_0.newick -s 0
```

### Command line outputs
* gene tree: an ascii drawing of the gene tree as a first-stage visualization, which does not cover the information of branch lengths. A much better visualization can be obtained by importing gene_tree.newick to many third party programs.
* Exception: ALL LOST: all gene lineages are lost, hence there is no gene tree generated.
  
### Output files
* gene_tree.newick: final output for the end users, the final gene with internal node names removed and leaf names converted according to the names of species.
* gene_tree_untruncated.newick: intermediate output for debugging purposes, the fully labelled gene (with internal node names) tree before cutted from the loss points.
* gene_tree_truncated.newick: intermediate output for debugging purposes, the fully labelled gene (with internal node names) tree after cutted from the loss points.

## Authors

* **Qiuyi Li** - *Initial work*
* **[Geoffrey Law](https://github.com/luojiahai)** - *Data structure*



## License

This project is licensed under the GNU GPLv3 - see [LICENSE.md](LICENSE.md) for more details.

## Acknowledgments

QL would like to thank Yupei You, Yichi Zhang, and Yiling Cao for their undivided support with programming, and Celine Scornavacca, Yao-ban Chan for assistance with debugging.
