import sys
from src.MLMSC_model import *


def default(str):
    return str + ' [Default: %default]'


# def parseDistributionArgs(str):
#     if str == None:
#         return {}
#     pieces = str.split(',')
#     opts = {}
#     for p in pieces:
#         if '=' in p:
#             key, val = p.split('=')
#         opts[key] = float(val)
#     return opts


def readCommand(argv):
    """
    Processes the command used to run MLMSC_Model from the command line.
    """
    from optparse import OptionParser
    usageStr = """
    USAGE:      python MLMSC.py <options>
    EXAMPLES:   (1) python MLMSC.py
                    - runs a model
                (2) python MLMSC.py --input data/species_tree.txt
                OR  python MLMSC.py -i data/species_tree.txt
    """
    parser = OptionParser(usageStr, add_help_option=False)

    parser.add_option(
        '--help', action='store_true', help='show this help message')

    parser.add_option(
        '-i', '--inputFile', dest='inputFile',
        help='the path to an input file of a species tree',
        metavar='INPUT_FILE')

    parser.add_option(
        '-s', '--seed', type='int', dest='seedArgs',
        help=default(
            'the seed number'
            'e.g., "-s 0"'),
        default=None)

    parser.add_option(
        '-c', '--coalescentRate', type='float', dest='coalescentArgs',
        help=default(
            'the unit rate of coalescent'
            'e.g., "-c 0.5"'),
        default=1)

    parser.add_option(
        '-r', '--recombinationRate', type='float', dest='recombinationArgs',
        help=default(
            'the unit rate of recombination'
            'e.g., "-r 0.5"'),
        default=0)

    parser.add_option(
        '-d', '--duplicationRate', type='float', dest='duplicationArgs',
        help=default(
            'the rate of occurence of duplications, '
            'e.g., "-d 0.2"'),
        default=0.1)
        
    parser.add_option(
        '-t', '--transferRate', type='float', dest='transferArgs',
        help=default(
            'the rate of occurence of transfers, '
            'e.g., "-t 0.1"'),
        default=0)

    parser.add_option(
        '-l', '--lossRate', type='float', dest='lossArgs',
        help=default(
            'the rate of occurence of losses, '
            'e.g., "-l 0.2"'),
        default=0.1)

    parser.add_option(
        '-u', '--unlinkProb', type='float', dest='unlinkArgs',
        help=default('probability for a duplication to be unlinked'),
        default=1)

    parser.add_option(
        '-n', '--numRepeats', type='int', dest='repeatNumber',
        help=default(
            'number of repeats, '
            'e.g., "-n 10"'),
        default=1)

    parser.add_option(
        '-h', '--hemiplasy', type='int', dest='hemiplasy',
        help=default('hemiplasy option, 0 or 1'), metavar='HEMIPLASY',
        default=1)

    parser.add_option(
        '-v', '--verbose', type='int', dest='verbose',
        help=default('verbose option, 0 or 1'), metavar='VERBOSE',
        default=0)

    options, otherjunk = parser.parse_args(argv)
    if len(otherjunk) != 0:
        raise Exception('Command line input not understood: ' + str(otherjunk))
    args = dict()

    if options.help:
        parser.print_help()
        sys.exit()

    # input file (a species tree in newick format)
    if not options.inputFile:
        parser.error('The input filename not given')
    args['inputFile'] = options.inputFile

    # distribution arguments
    args['seedArgs'] = options.seedArgs
    args['coalescentArgs'] = options.coalescentArgs
    args['recombinationArgs'] = options.recombinationArgs
    args['duplicationArgs'] = options.duplicationArgs
    args['transferArgs'] = options.transferArgs
    args['lossArgs'] = options.lossArgs
    args['unlinkArgs'] = options.unlinkArgs
    args['repeatNumber'] = options.repeatNumber

    # hemiplasy option
    if options.hemiplasy != 0 and options.hemiplasy != 1:
        parser.error('Invalid hemiplasy option: ' + str(options.hemiplasy))
    args['hemiplasy'] = True if options.hemiplasy == 1 else False

    # verbose option
    if options.verbose != 0 and options.verbose != 1:
        parser.error('Invalid verbose option: ' +
                     str(options.verbose))
    args['verbose'] = True if options.verbose == 1 else False
    return args


def runModel(**args):
    model = MLMSC_Model(seed=args['seedArgs'])
    model.run(**args)


if __name__ == '__main__':
    """
    The MLMSC function called when MLMSC.py is run
    from the command line:

    > python MLMSC.py

    See the usage string for more details.

    > python MLMSC.py --help
    """
    args= readCommand(sys.argv[1:])  # Get MLMSC components based on input
    runModel(**args)

    pass
