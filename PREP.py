import argparse
import os
import sys
from src import ui
import shutil

def main():
	
	str_prog = 'PREP'
	str_usage = 'python PREP.py <command> [options]'
	str_desc = r'''
some description
'''

	p = argparse.ArgumentParser(prog=str_prog,
                                 usage=str_usage,
                                 description=str_desc,
                                 formatter_class=argparse.RawDescriptionHelpFormatter,add_help=True)

	p.add_argument('-v', '--version', action='version', version='%(prog)s 1.0.0(dev)')

	sp = p.add_subparsers(title='Command', metavar='')

	p_RE = sp.add_parser('RE',
						description='detecting RNA editing neoantigens usingRNA sequencing data',
						usage='python PREP.py RE_PE [options]',
						help='detecting RNA editing neoantigens usingRNA sequencing data')
	ui.ParseRE(p_RE)
	if len(sys.argv) == 1:
		p.print_help()
		sys.exit(1)
	elif len(sys.argv) == 2:
		if sys.argv[1] == '-h' or sys.argv[1] == '--help':
			p.print_help()
		if sys.argv[1] == 'RE':
			p_RE.print_help()
			sys.exit(1)
	opts = p.parse_args()
	if sys.argv[1] == 'RE':
		from src.core import RE 
		RE.RNAEditing(opts)


if __name__ == '__main__':
	main()
