import sys
import os
import numpy as np

def readCommand(argv):
    from optparse import OptionParser
    usageStr = ''
    parser = OptionParser(usageStr, add_help_option=False)

    parser.add_option(
        '-n', '--name', dest='name')

    options, otherjunk = parser.parse_args(argv)

    if len(otherjunk) != 0:
        raise Exception('Command line input not understood: ' + str(otherjunk))

    args = options.name
    
    return args

args = readCommand(sys.argv[1:]) 

n = int(args.split('N')[1].split('D')[0])
d = float(args.split('D')[1].split('L')[0])
l = float(args.split('L')[1].split('C')[0])
c = int(args.split('C')[1].split('P')[0])

mlmsc_command = 'python3 MLMSC.py -i species_trees/fungi.newick -n ' + str(n) + ' -d ' + str(0.001*d) + ' -l ' + str(0.001*l) + ' -c ' + str(1/c) + ' -u 1'

for i in range(1):
    seed = str(np.random.random())[2:10]
    seed = str(int(seed))

    os.system(mlmsc_command + ' -s ' + seed)

print("done.")