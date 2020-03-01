# MLMSC Simulator

MLMSC Simulator is a program for the simulation of gene family evolution within a species tree based on the Multilocus Multipecies Coalescent (MLMSC) model. MLMSC model generalises the multispecies coalescent to gene families, and is designed to capture all possible scenarios that can arise through incomplete lineage sorting, gene duplication, transfer and loss, and any interaction between these processes. The MLMSC combines forward- and backward-in-time modelling in order to properly account for copy number hemiplasy and linkage between loci. 
The input for MLMSC Simulator are the simulation parameter values and a pre-specified species tree in Newick format as an input file. The output is a simulated gene tree in Newick format.

## Obtaining MLMSC Simulator

### Download

The source code of MLMSC Simulator is available in the GitHub repository: 
```
https://github.com/QiuyiLi/MLMSC
```
You can clone the sources in your computer by executing 
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
You can install these packages by executing 
```
pip install -r requirements.txt
```
### Testing

The MLMSC Simulator should be ready to use now. You can test the the simulator by executing
```
python3 MLMSC.py -i data/tree_sample_0.txt -s 2020210
```
The simulator will run with the default setting and you should be able to see the following outputs:
```
distances from tips to root:
5.992582775227923 A_lv=0_id=3
5.992582775227923 A
5.992582775227923 A_lv=0_id=1
5.992582775227923 A_lv=0_id=2
5.992582775227923 D_lv=0_id=0
5.992582775227923 D
gene tree:
                                        /-A_lv=0_id=3
                              /--------|
                             |         |          /-A
                    /--------|          \--------|
                   |         |                    \-A_lv=0_id=1
          /--------|         |
         |         |          \-A_lv=0_id=2
---------|         |
         |          \-D_lv=0_id=0
         |
          \-D

```

##  Usage

### Input file
* Species tree: -i or --inputFile; the prespecified species tree written in [Newick](https://en.wikipedia.org/wiki/Newick_format) format.
```
e.g., python3 MLMSC.py -i data/tree_sample_0.txt
```

### Input parameters

* Coalescent rate: -c or --coalescentRate; the coalesent rate in multispecies coalescent, default: -c 1.0.
```
e.g., python3 MLMSC.py -i data/tree_sample_0.txt -c 0.5
```
  
* Recombination rate: -r or --recombinationRate; the recombiantion rate in linked coalescent, default: -r 0.5.
```
e.g., python3 MLMSC.py -i data/tree_sample_0.txt -r 0.5
```
  
* Duplication rate: -d or --dulpicationRate; occurrence rate of gene duplication, default: -d 0.2.
```
e.g., python3 MLMSC.py -i data/tree_sample_0.txt -d 0.5
```
  
* Transfer rate; -t or --transferRate; occurrence rate of horizontal gene transfer, default: -t 0.1.
```
e.g., python3 MLMSC.py -i data/tree_sample_0.txt -t 0.5
```
  
* Loss rate; -l or --lossRate; occurrence rate of gene loss, default: -l 0.2.
```
e.g., python3 MLMSC.py -i data/tree_sample_0.txt -l 0.5
```
  
* Unlink probability; -u or --unlinkProb; the probability that a duplication is unlinked, default: -u 0.5.
```
e.g., python3 MLMSC.py -i data/tree_sample_0.txt -u 0.5
```
  
* Hemiplasy option: -h or --hemiplasy; whether or not the copy number hemiplasy is allowed, default: -h True.
```
e.g., python3 MLMSC.py -i data/tree_sample_0.txt -h False
```
  
* Verbose option: -v or --verbose; detailed outputs for debugging purposes, default: -v False.
```
e.g., python3 MLMSC.py -i data/tree_sample_0.txt -v True
```

* Seed number: -s or --seed; set seed for reproduciblility, default -s None (ramdon).
```
e.g., python3 MLMSC.py -i data/tree_sample_0.txt -s 0
```


### Command line outputs
* distances from tips to root:
* gene tree:
  
### Output files
* gene_tree_untruncated: intermediate output for debugging purposes, 
* gene_tree_truncated: intermediate output for debugging purposes,
* gene_tree_cleaned: final output for the end users, 

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```
## Authors

* **Qiuyi Li** - *Initial work*
* **Geoffery Law** - *Structuring*

See also the list of [contributors](https://github.com/QiuyiLi/MLMSC/graphs/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone whose code was used
* Inspiration
* etc
