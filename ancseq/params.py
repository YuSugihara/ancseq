import sys
import argparse
from ancseq.__init__ import __version__


class Params(object):

    def __init__(self, program_name):
        self.program_name = program_name

    def set_options(self):
        if self.program_name == 'ancseq':
            parser = self.ancseq_options()

        if len(sys.argv) == 1:
            args = parser.parse_args(['-h'])
        else:
            args = parser.parse_args()
        return args

    def ancseq_options(self):
        parser = argparse.ArgumentParser(description='ancseq version {}'.format(__version__),
                                         formatter_class=argparse.RawTextHelpFormatter)
        parser.usage = 'ancseq -s <ALIGNED_FASTA> -m <MODE> -o <OUT_DIR> [-t <INT>]'

        # set options
        parser.add_argument('-s',
                            '--seq',
                            action='store',
                            required=True,
                            type=str,
                            help='Sequence alignment in FASTA format.',
                            metavar='')

        parser.add_argument('-m',
                            '--mode',
                            required=True,
                            type=str,
                            help='Sequence type. [NT/AA/CODON]',
                            choices=['DNA','AA', 'CODON'],
                            metavar='')

        parser.add_argument('-o',
                            '--out',
                            action='store',
                            required=True,
                            type=str,
                            help='Output directory. The given name must not exist.',
                            metavar='')

        parser.add_argument('-t',
                            '--threads',
                            action='store',
                            default=8,
                            type=int,
                            help='Number of threads. [8]',
                            metavar='')
        
        parser.add_argument('-b',
                            '--bootstrap',
                            action='store',
                            default=1000,
                            type=int,
                            help='Replicate for bootstrap. [1000]',
                            metavar='')
        
        parser.add_argument('--max-report',
                            action='store',
                            default=5,
                            type=int,
                            help='Maximum number of ambiguous sites to report at the same position. [5]',
                            metavar='')
        
        parser.add_argument('--min-prob',
                            action='store',
                            default=0.05,
                            type=float,
                            help='Minimum probability of being reported as an ambiguous site. [0.05]',
                            metavar='')
        
        parser.add_argument('--min-gap-prob',
                            action='store',
                            default=0.5,
                            type=float,
                            help='Minimum probability of replacing the ancestral state with a gap. [0.5]',
                            metavar='')
        
        parser.add_argument('--fast',
                            action='store_true',
                            help='Use -fast option in IQ-TREE [FLASE]')
        
        parser.add_argument('--model',
                            action='store',
                            default=None,
                            type=str,
                            help='Input substitution model for IQ-TREE. [None]',
                            metavar='')
        
        parser.add_argument('--stop-pseudo-codon',
                            action='store_true',
                            help='Stop calculation of pvalues of pseudo-codon [FLASE]')
        
        parser.add_argument('--asr-only',
                            action='store_true',
                            help='Skip building tree and reconstruct ancestral states only [FLASE]')

        # set version
        parser.add_argument('-v',
                            '--version',
                            action='version',
                            version='%(prog)s {}'.format(__version__))

        return parser

    