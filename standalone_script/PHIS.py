import os,argparse
import sys
from Bio import Entrez, SeqIO
def inputinstructions():
	return """PHIS 1.1.0 arguments:

Usage: PHIS [options]
Options (x is an integer number)

--input <file name>        : Query phage file path: FASTA or Multi-Fasta or GenBank file
--out <folder name>        : Output folder in which results will be stored
--model <model name>       : Model name(default: all):
                           1.all 
                              run 5 models:CRISPR,Prophage,Protein_protein interaction,Genetic homology and Oligonucleotide profile/sequence composition to predict the hosts of query phage
                           2.crispr
                              run CRISPR model to predict the hosts of query phage by making blastn with spacer database
                           3.prophage
                              run Prophage model to predict the hosts of query phage by making blastp and blastn with prophage database
                           4.protein_protein_interaction
                              run Protein_protein interaction model to predict the hosts of query phage by finding PPIs(protein_protein interaction) and DDIs(domain_domain interaction)
                           5.blast
                              run Genetic homology model to predict the hosts of query phage by making blastp and blastn with bacteria database
--min_mis_crispr <x>       : Minimal mismatch of a Blastn hit on hit spacers(default: 2)
--min_cov_crispr <x>       : Minimal % coverage of a Blastn hit on hit spacers(default: 70)
--min_per_prophage <x>     : Minimal % percentage of hit proteins on hit prophage region(default:30)
--min_id_prophage <x>      : Minimal % identity of hit region on hit prophage region by making blastn(default:70)
--min_cov_prophage <x>     : Minimal % coverage of hit region on hit prophage region by making blastn(default:30)
--min_PPI <x>              : Minimal PPI number of a pair of phage-host pair(default:1)
--min_DDI <x>              : Minimal DDI number of a pair of phage-host pair(default:5)
--min_per_blast <x>        : Minimal % percentage of hit proteins on query phage and host by making blastp(default:10)
--min_id_blast <x>         : Minimal % identity of hit region on query phage and host by making blastn(default:70)
--min_cov_blast <x>        : Minimal % coverage of hit region on query phage and host by making blastn(default:10)
"""

def get_root_path():
	crispr_script_path = os.path.split(os.path.realpath(__file__))[0]
	root_path = os.path.dirname(crispr_script_path)
	return root_path

def predict_on_all(file,outdir,parameter_dict):
	script = os.path.join(root_path,'script',"overall_phage_to_host_domain.py")
	python_path = os.path.join(root_path,'software/python/miniconda3/bin','python')
	if not os.path.exists(python_path):
		python_path = "python"
	command = python_path+" "+script+" --input "+file+" --output "+outdir
	for par_name in parameter_dict.keys():
		if par_name!='model':
			command = command+" --"+par_name+" "+str(parameter_dict[par_name])
	os.system(command)

def predict_on_crispr(file,outdir,min_mis_crispr=2,min_cov_crispr=70):
	script = os.path.join(root_path,'script','PHIS_crispr.py')
	python_path = os.path.join(root_path,'software/python/miniconda3/bin','python')
	if not os.path.exists(python_path):
		python_path = "python"
	command = python_path+" "+script+" --input "+file+" --output "+outdir+" --min_mis_crispr "+str(min_mis_crispr)+" --min_cov_crispr "+min_cov_crispr
	os.system(command)

def predict_on_prophage(file,outdir,min_per_prophage=30,min_id_prophage=70,min_cov_prophage=30):
	script = os.path.join(root_path,'script','PHIS_prophage.py')
	python_path = os.path.join(root_path,'software/python/miniconda3/bin','python')
	if not os.path.exists(python_path):
		python_path = "python"
	command = python_path+" "+script+" --input_file "+file+" --output "+outdir+" --min_per_prophage "+str(min_per_prophage)+" --min_id_prophage "+str(min_id_prophage)+" --min_cov_prophage "+str(min_cov_prophage)
	os.system(command)

def predict_on_blast(file,outdir,min_per_blast=10,min_id_blast=70,min_cov_blast=10):
	script = os.path.join(root_path,'script','PHIS_blast.py')
	python_path = os.path.join(root_path,'software/python/miniconda3/bin','python')
	if not os.path.exists(python_path):
		python_path = "python"
	command = python_path+" "+script+" --input_file "+file+" --output "+outdir+" --min_per_blast "+str(min_per_blast)+" --min_id_blast "+str(min_id_blast)+" --min_cov_blast "+str(min_cov_blast)
	os.system(command)

def predict_on_pro_int(file,outdir,min_PPI=1,min_DDI=5):
	script = os.path.join(root_path,'script','PHIS_protein_protein_int.py')
	python_path = os.path.join(root_path,'software/python/miniconda3/bin','python')
	if not os.path.exists(python_path):
		python_path = "python"
	command = python_path+" "+script+" --input_file "+file+" --output "+outdir+" --min_PPI "+str(min_PPI)+" --min_DDI "+str(min_DDI)
	os.system(command)

def set_default_parameter():
	min_mis_crispr = 2
	min_cov_crispr = 70
	min_per_prophage = 30
	min_id_prophage = 70
	min_cov_prophage = 30
	min_PPI = 1
	min_DDI = 5
	min_per_blast = 10
	min_id_blast = 70
	min_cov_blast = 10
	model = 'all'
	default_parameter_dict = {'min_mis_crispr':min_mis_crispr,
	'min_cov_crispr':min_cov_crispr,
	'min_per_prophage':min_per_prophage,
	'min_id_prophage':min_id_prophage,
	'min_cov_prophage':min_cov_prophage,
	'min_PPI':min_PPI,
	'min_DDI':min_DDI,
	'min_per_blast':min_per_blast,
	'min_id_blast':min_id_blast,
	'min_cov_blast':min_cov_blast,
	'model':model}
	return default_parameter_dict

def check_file(file):
	with open(file) as f:
		contents = f.read()
	if (contents[0]!='>'):
		try:
			SeqIO.parse(file,"genbank")
		except:
			print('phage file format error!please input genbank or fasta file')
			sys.exit(-1)

def check_outdir(outdir):
	if not os.path.exists(outdir):
		print('the output directory does not exists!')
		sys.exit(-1)

def check_parameter(args):
	default_parameter_dict = set_default_parameter()
	variable_dict = vars(args)
	for name in variable_dict.keys():
		if name=='input':
			continue
		if name=='output':
			continue
		if name=='model':
			if args.model:
				input_model = args.model
				if input_model not in ['all','crispr','prophage','protein_protein_interaction','blast','opc']:
					print('there no exists %s'%input_model)
					sys.exit(-1)
		else:
			if variable_dict[name]:
				try:
					cutoff = int(variable_dict[name])
				except:
					print('%s needs integer!%s is not integer!'%(name,variable_dict[name]))
					sys.exit(-1)
		if variable_dict[name]:
			default_parameter_dict[name] = variable_dict[name]
	return default_parameter_dict

def run(input_file,outdir,parameter_dict):
	if parameter_dict['model'] == 'all':
		predict_on_all(input_file,outdir,parameter_dict)	
	if parameter_dict['model'] == 'crispr':
		predict_on_crispr(input_file,outdir,parameter_dict['min_mis_crispr'],parameter_dict['min_cov_crispr'])
	if parameter_dict['model'] == 'prophage':
		predict_on_prophage(input_file,outdir,parameter_dict['min_per_prophage'],parameter_dict['min_id_prophage'],parameter_dict['min_cov_prophage'])
	if parameter_dict['model'] == 'protein_protein_interaction':
		predict_on_pro_int(input_file,outdir,parameter_dict['min_PPI'],parameter_dict['min_DDI'])
	if parameter_dict['model'] == 'blast':
		predict_on_blast(input_file,outdir,parameter_dict['min_per_blast'],parameter_dict['min_id_blast'],parameter_dict['min_cov_blast'])

if __name__=='__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--input',help='''Path of the input file in genbank or fasta format for bacteria or phages
		please ensure your file name consist of word,number,_,if not,the program will probably run error!''')
	parser.add_argument('--output', help='Path of the output directory.\n')
	parser.add_argument('--model', help='model name(default: all)\n')	
	parser.add_argument('--min_mis_crispr', help='Minimal mismatch of a Blastn hit on hit spacers(default: 2)\n')
	parser.add_argument('--min_cov_crispr', help='Minimal % coverage of a Blastn hit on hit spacers(default: 90)\n')
	parser.add_argument('--min_per_prophage', help='Minimal % percentage of hit proteins on hit prophage region(default:30)\n')	
	parser.add_argument('--min_id_prophage', help='Minimal % identity of hit region on hit prophage region by making blastn(default:70)\n')	
	parser.add_argument('--min_cov_prophage', help='Minimal % coverage of hit region on hit prophage region by making blastn(default:30)')	
	parser.add_argument('--min_PPI', help='Minimal PPI number of a pair of phage-host pair(default:1)\n')	
	parser.add_argument('--min_DDI', help='Minimal DDI number of a pair of phage-host pair(default:5)\n')	
	parser.add_argument('--min_per_blast', help='Minimal % percentage of hit proteins on query phage and host by making blastp(default:10)\n')	
	parser.add_argument('--min_id_blast', help='Minimal % identity of hit region on query phage and host by making blastn(default:70)\n')		
	parser.add_argument('--min_cov_blast', help='Minimal % coverage of hit region on query phage and host by making blastn(default:10)\n')		
	parser.add_argument('--h',help='''print help information''')
	
	args = parser.parse_args()
	global type,root_path
	if args.h:
		print(inputinstructions())
		sys.exit(-1)
	if not args.input:
		print('Error:please input phage file!')
		print(inputinstructions())
		sys.exit(-1)
	if not args.output:
		print('Error:please input output directory!')
		print(inputinstructions())
		sys.exit(-1)
	if args.input:
		input_file = args.input
	if args.output:
		outdir = args.output
	root_path = get_root_path()
	print(root_path)
	check_file(input_file)
	check_outdir(outdir)
	parameter_dict = check_parameter(args)
	run(input_file,outdir,parameter_dict)