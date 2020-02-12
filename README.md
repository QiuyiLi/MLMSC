# IxDTL
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
