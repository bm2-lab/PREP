import argparse

def ParseRE(p_re):
	p_re_inp = p_re.add_argument_group('Input options')
	p_re_inp.add_argument('-i', '--config_file', dest='Config_file', required=True, help='input configure file of programme(required)')	
	
