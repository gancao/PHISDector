import numpy as np
import os,re,argparse
from Bio import Entrez, SeqIO

def get_root_path():
	crispr_script_path = os.path.split(os.path.realpath(__file__))[0]
	root_path = os.path.dirname(crispr_script_path)
	return root_path

def mkdir(dirName):
	if not os.path.exists(dirName):
		cmd_mkdir = 'mkdir -p %s' % dirName
		os.system(cmd_mkdir)

def build_model(bacDir,modelDir):
	wish_script = os.path.join(root_path,'script/wish','WIsH')
	cmd_mkModel = wish_script+' -c build -g %s -m %s' % (bacDir,modelDir)
	os.system(cmd_mkModel)

def predict_wish(phage_dir,modelDir,outdir):
	wish_script = os.path.join(root_path,'script/wish','WIsH')
	cmd_predict = wish_script+' -c predict -g %s -m %s -r %s -b'%(phage_dir,modelDir,outdir)
	os.system(cmd_predict)

def parse_wish(resultfile,wish_dict,kind):
	with open(resultfile) as f:
		contents = f.readlines()
	#row:bacteria, column:phage
	phages = contents[0].strip().split('\t')
	for line in contents[1:]:
		line = line.strip().split('\t')
		bac_id = line[0].split('.')[0]
		values = line[1:]
		flag = 0
		for phage_id in phages:
			phage_id = phage_id.split('.')[0]
			if kind == 'p':
				if phage_id not in wish_dict.keys():
					wish_dict.update({phage_id:{}})
				if bac_id not in wish_dict[phage_id].keys():
					wish_dict[phage_id].update({bac_id:values[flag]})
			else:
				if bac_id not in wish_dict.keys():
					wish_dict.update({bac_id:{}})
				if phage_id not in wish_dict[bac_id].keys():
					wish_dict[bac_id].update({phage_id:values[flag]})
			flag = flag + 1

def phage_to_hosts(input_dir,hits_dict,outdir,kind):
	bac_model_dir = os.path.join(root_path,'database/profiles/wish/modelDir')
	model_dir = os.path.join(outdir,'model_dir')
	result_dict_file = os.path.join(outdir,'wish_result_dict')
	outdir = os.path.join(outdir,'results')
	wish_result_dict = {}
	for phage_id in hits_dict.keys():
		mkdir(model_dir)
		outdir_now = os.path.join(outdir,phage_id)
		mkdir(outdir_now)
		phage_dir = os.path.join(input_dir,phage_id)
		for host_id in hits_dict[phage_id]:
			host_model_file = os.path.join(bac_model_dir,host_id+'.mm')
			command = "cp "+host_model_file+" "+model_dir
			os.system(command)
		predict_wish(phage_dir,model_dir,outdir_now)
		command = "rm -rf "+model_dir
		os.system(command)
		outfile = os.path.join(outdir_now,'llikelihood.matrix')
		if os.path.exists(outfile):
			parse_wish(outfile,wish_result_dict,kind)
	np.save(result_dict_file,wish_result_dict,kind)

def phage_to_hosts_all(input_dir,outdir,kind):
	bac_model_dir = os.path.join(root_path,'database/profiles/wish/modelDir')
	model_dir = os.path.join(outdir,'model_dir')
	result_dict_file = os.path.join(outdir,'wish_result_dict')
	outdir = os.path.join(outdir,'results')
	wish_result_dict = {}
	phage_dirs = os.listdir(input_dir)
	for phage_id in phage_dirs:
		outdir_now = os.path.join(outdir,phage_id)
		mkdir(outdir_now)
		phage_dir = os.path.join(input_dir,phage_id)
		predict_wish(phage_dir,bac_model_dir,outdir_now)
		outfile = os.path.join(outdir_now,'llikelihood.matrix')
		if os.path.exists(outfile):
			parse_wish(outfile,wish_result_dict,kind)
	np.save(result_dict_file,wish_result_dict,kind)

if __name__=='__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--input_dir',help='Path of the input sequence,one dir for a fasta file.\n')
	parser.add_argument('--output', help='Path of the output file.\n')
	parser.add_argument('--out_dict', help='Path of the output result file in first step.\n')
	args = parser.parse_args()
	global type,strain_id,strain_def,root_path
	if args.input_dir:
		input_dir = args.input_dir
	if args.output:
		outdir = args.output
	if args.out_dict:
		output_dict_file = args.out_dict
	else:
		output_dict_file = "no"
	root_path = get_root_path()
	kind = 'p'
	if output_dict_file!='no':
		hosts_dict = np.load(output_dict_file).item()
		phage_to_hosts(input_dir,hosts_dict,outdir,kind)
	else:
		phage_to_hosts_all(input_dir,outdir,kind)