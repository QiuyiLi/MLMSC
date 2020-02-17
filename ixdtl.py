import sys
from src.ixdtl_model import *


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
    Processes the command used to run IxDTLModel from the command line.
    """
    from optparse import OptionParser
    usageStr = """
    USAGE:      python ixdtl.py <options>
    EXAMPLES:   (1) python ixdtl.py
                    - runs a model
                (2) python ixdtl.py --input data/species_tree.txt
                OR  python ixdtl.py -i data/species_tree.txt
    """
    parser = OptionParser(usageStr, add_help_option=False)

    parser.add_option(
        '--help', action='store_true', help='show this help message')

    parser.add_option(
        '-i', '--inputFile', dest='inputFile',
        help='the path to an input file of a species tree',
        metavar='INPUT_FILE')

    parser.add_option(
        '-c', '--coalescentArgs', type='float', dest='coalescentArgs',
        help=default(
            'the unit rate of coalescent'
            'e.g., "-c 0.5"'),
        default=1)

    parser.add_option(
        '-r', '--recombinationRate', type='float', dest='recombinationArgs',
        help=default(
            'the unit rate of recombination'
            'e.g., "-c 0.5"'),
        default=0.5)

    parser.add_option(
        '-d', '--duplicationRate', type='float', dest='duplicationArgs',
        help=default(
            'the rate of occurence of duplications, '
            'e.g., "-d 0.2"'),
        default=0.2)
        
    parser.add_option(
        '-t', '--transferRate', type='float', dest='transferArgs',
        help=default(
            'the rate of occurence of transfers, '
            'e.g., "-t 0.1"'),
        default=0.1)

    parser.add_option(
        '-l', '--lossRate', type='float', dest='lossArgs',
        help=default(
            'the rate of occurence of losses, '
            'e.g., "-l 0.2"'),
        default=0.2)

    parser.add_option(
        '-u', '--unlinkRate', type='float', dest='unlinkArgs',
        help=default('probability for a duplication to be unlinked'),
        default=0.5)

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
    args['coalescentArgs'] = options.coalescentArgs
    args['recombinationArgs'] = options.recombinationArgs
    args['duplicationArgs'] = options.duplicationArgs
    args['transferArgs'] = options.transferArgs
    args['lossArgs'] = options.lossArgs
    args['unlinkArgs'] = options.unlinkArgs

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
    model = IxDTLModel()
    # model = IxDTLModel(seed=14)
    model.run(**args)


if __name__ == '__main__':
    """
    The ixdtl function called when ixdtl.py is run
    from the command line:

    > python ixdtl.py

    See the usage string for more details.

    > python ixdtl.py --help
    """
    args = readCommand(sys.argv[1:])  # Get ixdtl components based on input
    runModel(**args)

    pass
