# MLMSC Simulator

MLMSC Simulator is a program for the simulation of gene family evolution within a species tree based on the Multilocus Multipecies Coalescent (MLMSC) model. MLMSC model generalises the multispecies coalescent to gene families, and is designed to capture all possible scenarios that can arise through ILS, gene duplication, transfer and loss, and any interaction between these processes. The MLMSC combines forward- and backward-in-time modelling in order to properly account for copy number hemiplasy and linkage between loci. 
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

### Installing the requirements
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

The MLMSC Simulator should be ready to use. You can test the the simulator by executing
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
qiuyidembp:IxDTL2.0 qiuyi_li$ 

```

##  Usage

### Input file: species_tree.newick

### Input parameters

### Output files

### Other outputs

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







To run, type: 

python3 ixdtl.py -i data/tree_sample_0.txt

the code will then run with default settings, if you want to change parameters by yourself:

-c coalescent rate			e.g. -c 0.5 (default)
  
-r recombination        e.g. -r 0.5 (no need to use under current version)
  
-d duplication rate     e.g. -d 0.2 (default)
  
-t transfer rate      	e.g. -t 0 (default)
  
-l loss rate      			e.g. -l 0.1 (default)
  
-u unlink probability 	e.g. -u 0.8 (no need to use under current version)
  
-h hemiplasy option 		e.g. -h True (default)
  
-v verbose option 			e.g. -v False (defalt)
  
For more options, check ixdtl.py


Suggested order of reading:

speciesTree.py

haplotypeTree.py & locusTree.py

ixdtl_model.py


Tree structure:

tree_table.py


Input options:

ixdtl.py
