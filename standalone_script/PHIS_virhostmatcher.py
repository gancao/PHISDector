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

def run_virhostmatcher(sequence_dir,kmer_list_file,taxo_path,outdir,kind):
	script = os.path.join(root_path,'script/virhostmatcher','vhm_modify_input_model.py')
	python_path = os.path.join(root_path,'software/python/miniconda3/bin','python')
	if not os.path.exists(python_path):
		python_path = "python"
	if kind == 'p':
		command = python_path+" "+script+" -f "+kmer_list_file+" -s "+sequence_dir+" -o "+outdir+" -t "+taxo_path+" -d 1 -k p"
		print(command)
		os.system(command)
	else:
		command = python_path+" "+script+" -f "+kmer_list_file+" -s "+sequence_dir+" -o "+outdir+" -d 1 -k b"
		os.system(command)

def parse_virhostmatcher(resultfile,virhostmatcher_dict,kind):
	with open(resultfile) as f:
		contents = f.readlines()
	#row:phage, column:bacteria
	hosts = contents[0].strip().strip(',').split(',')[1:]
	for line in contents[1:]:
		line = line.strip().strip(',').split(',')
		phage_id = line[0].split('.')[0]
		values = line[1:]
		flag = 0
		for value in values:
			value = 1-2*float(value)
			bac_id = hosts[flag].split('.')[0]
			if kind == 'p':
				if phage_id not in virhostmatcher_dict.keys():
					virhostmatcher_dict.update({phage_id:{}})
				if bac_id not in virhostmatcher_dict[phage_id].keys():
					virhostmatcher_dict[phage_id].update({bac_id:str(value)})
			else:
				if bac_id not in virhostmatcher_dict.keys():
					virhostmatcher_dict.update({bac_id:{}})
				if phage_id not in virhostmatcher_dict[bac_id].keys():
					virhostmatcher_dict[bac_id].update({phage_id:str(value)})
			flag = flag+1

def get_taxa(hosts,taxa_file):
	with open(taxa_file,'w') as f:
		f.write('hostNCBIName\thostName\thostSuperkingdom\thostPhylum\thostClass\thostOrder\thostFamily\thostGenus\thostSpecies\n')
		for host in hosts:
			f.write(host+'.fasta'+'\t'+'\t'.join(['NA']*8)+'\n')

def phage_to_hosts(input_dir,hits_dict,outdir,kind):
	bac_kmer_dir = os.path.join(root_path,"database/profiles/virhostmatcher/KC")
	virhostmatcher_result_dict = {}
	i=0
	for phage_id in hits_dict.keys():
		i = i+1
		#outdir_now = os.path.join(outdir,phage_id)
		outdir_now = os.path.join(outdir,'contig_'+str(i))
		mkdir(outdir_now)
		phage_dir = os.path.join(input_dir,phage_id)
		kmer_list_file = os.path.join(outdir_now,'bac_kmer_list')
		f = open(kmer_list_file,'w')
		for host_id in hits_dict[phage_id]:
			host_kmer_path = os.path.join(bac_kmer_dir,host_id+'.fasta')
			f.write(host_id+'.fasta '+host_kmer_path+' 2\n')
			f.flush()
		f.close()
		taxa_file = os.path.join(outdir_now,'taxa.txt')
		get_taxa(hits_dict[phage_id],taxa_file)
		run_virhostmatcher(phage_dir,kmer_list_file,taxa_file,outdir_now,kind)
		outfile = os.path.join(outdir_now,'d2star_k6.csv')
		if os.path.exists(outfile):
			parse_virhostmatcher(outfile,virhostmatcher_result_dict,kind)
	virhostmatcher_result_dict_file = os.path.join(outdir,'virhostmatcher_result_dict')
	np.save(virhostmatcher_result_dict_file,virhostmatcher_result_dict)


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
		output_dict_file = 'no'
	root_path = get_root_path()
	kind='p'
	if os.path.exists(output_dict_file):
		hosts_dict = np.load(output_dict_file).item()
		phage_to_hosts(input_dir,hosts_dict,outdir,kind)
	else:
		print('error:%s does not exists!'%output_dict_file)