import os,threading
import numpy as np
import os,re,argparse
from Bio import Entrez, SeqIO
from collections import Counter
import json
from Bio import SeqIO

codons = [
    'AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT',
    'AGA', 'AGC', 'AGG', 'AGT', 'ATA', 'ATC', 'ATG', 'ATT',
    'CAA', 'CAC', 'CAG', 'CAT', 'CCA', 'CCC', 'CCG', 'CCT',
    'CGA', 'CGC', 'CGG', 'CGT', 'CTA', 'CTC', 'CTG', 'CTT',
    'GAA', 'GAC', 'GAG', 'GAT', 'GCA', 'GCC', 'GCG', 'GCT',
    'GGA', 'GGC', 'GGG', 'GGT', 'GTA', 'GTC', 'GTG', 'GTT',
    'TAA', 'TAC', 'TAG', 'TAT', 'TCA', 'TCC', 'TCG', 'TCT',
    'TGA', 'TGC', 'TGG', 'TGT', 'TTA', 'TTC', 'TTG', 'TTT'
]
codons.sort()
start_condons = ['ATG','GTG']
end_condons = ['TAA','TAG','TGA']

def get_root_path():
	crispr_script_path = os.path.split(os.path.realpath(__file__))[0]
	root_path = os.path.dirname(crispr_script_path)
	return root_path

def pred_orf(fasta_file,faa_prefix):
	script = os.path.join(root_path,'software','FragGeneScan','run_FragGeneScan.pl')
	cmd_fragGeneScan = script+' -genome %s -out %s -complete=1 -train=complete -thread=20' % (fasta_file, faa_prefix)
	os.system(cmd_fragGeneScan)

def get_inf(file):
	with open(file) as f:
		contents = f.read().strip()
	if contents[0]=='>':
		type ='fasta'
		try:
			record = SeqIO.read(file,"fasta")
			strain_id = record.id
			strain_def = record.description
		except:
			strain_def = 'unknown'
			strain_id = 'unknown'
	else:
		type="genbank"
		record = SeqIO.read(file,"genbank")
		strain_id = record.id
		strain_def = record.description
	return strain_id,strain_def,type

def getFaaFromGB(input_file,outdir):   #parse protein from genbank files in phaster
	special_pros = ['capsid','head','plate','tail','coat','portal','holin','integrase','transposase','terminase','protease','lysis','bacteriocin','tRNA']
	records = SeqIO.parse(input_file, "gb")
	file_name = os.path.basename(input_file).split()[0].split('.')[0]
	counter = 0
	outFileName = os.path.join(outdir,file_name+'.faa')
	outfiledef = os.path.join(outdir,file_name+'_protein_def')
	savefile = open(outFileName, 'w')
	savefile_protein = open(outfiledef,'w')
	for record in records:
		for feature in record.features:
			if feature.type == 'CDS':
				location = feature.location
				if str(location).find('+') != -1:
					direction = '+'
				elif str(location).find('-') != -1:
					direction = '-'
				if '<' in str(location):
					location = str(location).replace('<','')
				if '>' in str(location):
					location = str(location).replace('>','')
				locations = re.findall("\d+\.?\d*",str(location))
				min_start = locations[0]
				max_end = locations[1]
				for loc in locations:
					if int(loc)<int(min_start):
						min_start = loc
					if int(loc)>int(max_end):
						max_end = loc
				location = min_start+'_'+max_end+'_'+direction
				counter = counter+1
				if 'product' in feature.qualifiers:
					product = feature.qualifiers['product'][0]	  
					if 'protein_id' in feature.qualifiers:
						proteinId = feature.qualifiers['protein_id'][0]
					else:
						if 'inference' in feature.qualifiers:
							strInference = str(feature.qualifiers['inference'])
							if 'RefSeq' in strInference:
								proteinId = strInference.split('RefSeq:')[1].rstrip(']').rstrip('\'')
							elif 'SwissProt' in strInference:
								proteinId = strInference.split('SwissProt:')[1].rstrip(']').rstrip('\'')
							else:
								proteinId = 'unknown'
						else:
							proteinId = 'unknown'
					if 'translation' in feature.qualifiers:
						translation = feature.qualifiers['translation'][0]
						savefile_protein.write(proteinId+'\t'+product+'\n')
						savefile.write('>' +strain_id+'|' + str(location)+ '|' + str(proteinId)+'\n')					
					# savefile.write(">"+fileID+ '_'+str(counter)+'\n')                   
						if translation[-1] == '\n':
							savefile.write(translation)
						else:
							savefile.write(translation + '\n')
	savefile.close()
	savefile_protein.close()

def GetFnaSequence(fileName,outFileName):#get fasta sequence from GenBank file
	handle = open(fileName)
	SeqIO.convert(handle, 'genbank',outFileName, 'fasta')

def get_protein_position_genbank(pro_file):
	pro_locations = []
	with open(pro_file) as f:
		contents = f.read().strip().split('>')
	for line in contents[1:]:
		pro_start = line.split('\n')[0].split('|')[1].split('_')[-3]
		pro_end = line.split('\n')[0].split('|')[1].split('_')[-2]
		pro_id = line.split('\n')[0].split('|')[2]
		direction = line.split('\n')[0].split('|')[1].split('_')[-1]
		sequence = ''.join(line.split('\n')[1:]).strip()
		pro_location = [pro_start,pro_end,direction,len(sequence)]
		pro_locations.append(pro_location)
	return pro_locations

def get_protein_position_fasta(pro_file):
	pro_locations = []
	with open(pro_file) as f:
		contents = f.read().strip().split('>')
	for line in contents[1:]:
		pro_start = line.split('\n')[0].split('_')[-3]
		pro_end = line.split('\n')[0].split('_')[-2]
		direction = line.split('\n')[0].split('_')[-1]
		sequence = ''.join(line.split('\n')[1:]).strip()
		pro_location = [pro_start,pro_end,direction,len(sequence)]
		pro_locations.append(pro_location)
	return pro_locations

def scan_nucl(nucl_sequence,pro_length):
	i=0
	code_start = 0
	while i < len(nucl_sequence):
		condon = nucl_sequence[i:i + 3].upper()
		i += 3
		if condon not in codons:
			continue
		if condon in start_condons:
			code_start = i-3
			break
	if (len(nucl_sequence)-code_start)>((pro_length+1)*3):
		code_region_sequence = nucl_sequence[code_start:code_start+pro_length*3]
	else:
		code_region_sequence = nucl_sequence[code_start:code_start+3*int((((len(nucl_sequence)-code_start))/3))]
	return code_region_sequence
	
def translate_nucl(nucl_sequence):
	complement_nucl = ''
	dict = {'A':'T','T':'A','C':'G','G':'C'}
	for nucl in nucl_sequence:
		if nucl in dict.keys():
			complement_nucl = complement_nucl+dict[nucl.upper()]
		else:
			complement_nucl = complement_nucl+nucl
	return complement_nucl

def calculate_condon_usage(condon_usage_dict,code_region_sequence):
	i =0
	condon_num = 0
	while i<len(code_region_sequence):
		condon = code_region_sequence[i:i+3].upper()
		i = i+3
		if condon not in codons:
			continue
		condon_num = condon_num+1
		if condon not in condon_usage_dict.keys():
			condon_usage_dict.update({condon:0})
		condon_usage_dict[condon] = condon_usage_dict[condon] +1
	return condon_num

def get_protein_position_multifasta(pro_file):
	#fasta or multifasta
	pro_locations = {}
	with open(pro_file) as f:
		contents = f.read().strip().split('>')
	for line in contents[1:]:
		pro_start = line.split('\n')[0].split('_')[-3]
		pro_end = line.split('\n')[0].split('_')[-2]
		direction = line.split('\n')[0].split('_')[-1]
		strain_id = '_'.join(line.split('\n')[0].split('_')[0:-3]).split('.')[0]
		sequence = ''.join(line.split('\n')[1:]).strip()
		pro_location = [pro_start,pro_end,direction,len(sequence)]
		if strain_id not in pro_locations.keys():
			pro_locations.update({strain_id:[]})
		pro_locations[strain_id].append(pro_location)
	return pro_locations

def condon_usage(fasta_file,pro_file,type,f_out,strain_id):
	condon_us = {}
	with open(fasta_file) as f:
		contents = f.read().strip().split('\n')
	sequence = ''
	for line in contents[1:]:
		sequence = sequence+line.strip()
	if type == 'fasta':
		pro_locations = get_protein_position_fasta(pro_file)
	else:
		pro_locations = get_protein_position_genbank(pro_file)
	condon_usage_dict = {}
	condons_num = 0
	for pro_location in pro_locations:
		pro_start,pro_end,direction,pro_length = pro_location
		if int(pro_end)<len(sequence):
			nucl_sequence = sequence[int(pro_start)-1:int(pro_end)]
		else:
			nucl_sequence = sequence[int(pro_start)-1:]
		if direction=='+':
			code_region_sequence = scan_nucl(nucl_sequence,pro_length)
		else:
			nucl_sequence = nucl_sequence[::-1]
			complement_nucl = translate_nucl(nucl_sequence)
			complement_code_region_sequence = scan_nucl(complement_nucl,pro_length)
			code_region_sequence = translate_nucl(complement_code_region_sequence)[::-1]
		
		condon_num = calculate_condon_usage(condon_usage_dict,code_region_sequence)
		condons_num = condons_num+condon_num
	condon_usage_list = []
	for condon in codons:
		if condon in condon_usage_dict.keys():
			condon_us = 1.0*condon_usage_dict[condon]/condons_num
		else:
			condon_us = 0.0
		condon_usage_list.append(condon_us)
	
	f_out.write(strain_id+'\t'+'\t'.join(map(str,condon_usage_list))+'\n')
	f_out.flush()
	return condon_usage_list

def condon_usage_multifasta(fasta_file,pro_file,f_out):
	codon_us = {}
	sequence_dict = {}
	with open(fasta_file) as f:
		contents = f.read().strip()
	if '\n>' in contents:
		contents = contents.split('\n>')
		for strain in contents:
			strain_id = strain.split('\n')[0].split()[0].split('.')[0].strip('>')
			sequence = ''
			for line in strain.split('\n')[1:]:
				sequence = sequence+line.strip()
			sequence_dict.update({strain_id:sequence})
	else:
		sequence = ''
		strain_id = contents.split('\n')[0].split()[0].split('.')[0].strip('>')
		for line in contents.split('\n')[1:]:
			sequence = sequence+line.strip()
		sequence_dict.update({strain_id:sequence})
	strain_pro_locations = get_protein_position_multifasta(pro_file)
	condon_usage_dict = {}
	condons_num = 0
	for strain_id in sequence_dict.keys():
		sequence = sequence_dict[strain_id]
		codon_us.update({strain_id:[]})
		if strain_id in strain_pro_locations.keys():
			pro_locations = strain_pro_locations[strain_id]
		else:
			pro_locations = []
		for pro_location in pro_locations:
			pro_start,pro_end,direction,pro_length = pro_location
			if int(pro_end)<len(sequence):
				nucl_sequence = sequence[int(pro_start)-1:int(pro_end)]
			else:
				nucl_sequence = sequence[int(pro_start)-1:]
			if direction=='+':
				code_region_sequence = scan_nucl(nucl_sequence,pro_length)
			else:
				nucl_sequence = nucl_sequence[::-1]
				complement_nucl = translate_nucl(nucl_sequence)
				complement_code_region_sequence = scan_nucl(complement_nucl,pro_length)
				code_region_sequence = translate_nucl(complement_code_region_sequence)[::-1]
		
			condon_num = calculate_condon_usage(condon_usage_dict,code_region_sequence)
			condons_num = condons_num+condon_num
		condon_usage_list = []
		for condon in codons:
			if condon in condon_usage_dict.keys():
				condon_us = 1.0*condon_usage_dict[condon]/condons_num
			else:
				condon_us = 0.0
			condon_usage_list.append(condon_us)
	
		f_out.write(strain_id+'\t'+'\t'.join(map(str,condon_usage_list))+'\n')
		f_out.flush()
		codon_us[strain_id] = condon_usage_list
	return codon_us

def distance(x, y):
	x = np.array(x)
	y = np.array(y)
	return np.sqrt(np.sum((x-y)**2))

def get_phage_hosts_codon_distance(input_file,pro_file,hosts_dict,outdir):
	#get phage codon usage
	phage_codon_usage_file = os.path.join(outdir,'phage_codon_usage.txt')
	f_out = open(phage_codon_usage_file,'w')
	f_out.write('strain_id\t'+'\t'.join(codons)+'\n')
	f_out.flush()
	phage_codon_usage_dict = condon_usage_multifasta(input_file,pro_file,f_out)
	f_out.close()
	
	#get alternative hosts codon usage
	bac_codon_usage_dict_file = os.path.join(root_path,"database/db/profile",'bac_codon_usage_dict.npy')
	bac_codon_usage_dict = np.load(bac_codon_usage_dict_file).item()
	bac_codon_usage_file = os.path.join(outdir,'hosts_codon_usage.txt')
	#get codon usage distance dict
	phage_hosts_codon_distance_dict_file = os.path.join(outdir,'codon_result_dict')
	phage_hosts_codon_distance_dict = {}
	for phage_id in hosts_dict.keys():
		phage_hosts_codon_distance_dict.update({phage_id:{}})
		phage_codon_usage = phage_codon_usage_dict[phage_id]
		phage_codon_usage = list(map(float,phage_codon_usage))
		for bac_id in hosts_dict[phage_id]:
			bac_codon_usage = bac_codon_usage_dict[bac_id]
			bac_codon_usage = list(map(float,bac_codon_usage))
			codon_usage_distance = distance(phage_codon_usage,bac_codon_usage)
			phage_hosts_codon_distance_dict[phage_id].update({bac_id:codon_usage_distance})
	np.save(phage_hosts_codon_distance_dict_file,phage_hosts_codon_distance_dict)

def get_phage_hosts_codon_distance_all(input_file,pro_file,outdir):
	#get phage codon usage
	phage_codon_usage_file = os.path.join(outdir,'phage_codon_usage.txt')
	f_out = open(phage_codon_usage_file,'w')
	f_out.write('strain_id\t'+'\t'.join(codons)+'\n')
	f_out.flush()
	phage_codon_usage_dict = condon_usage_multifasta(input_file,pro_file,f_out)
	f_out.close()

	bac_codon_usage_dict_file = os.path.join(root_path,"database/db/profile",'bac_codon_usage_dict.npy')
	bac_codon_usage_dict = np.load(bac_codon_usage_dict_file).item()
	phage_hosts_codon_distance_dict_file = os.path.join(outdir,'codon_result_dict')
	phage_hosts_codon_distance_dict = {}
	for phage_id in phage_codon_usage_dict.keys():
		phage_hosts_codon_distance_dict.update({phage_id:{}})
		phage_codon_usage = phage_codon_usage_dict[phage_id]
		phage_codon_usage = list(map(float,phage_codon_usage))
		for bac_id in bac_codon_usage_dict.keys():
			bac_codon_usage = bac_codon_usage_dict[bac_id]
			bac_codon_usage = list(map(float,bac_codon_usage))
			codon_usage_distance = distance(phage_codon_usage,bac_codon_usage)
			phage_hosts_codon_distance_dict[phage_id].update({bac_id:codon_usage_distance})
	np.save(phage_hosts_codon_distance_dict_file,phage_hosts_codon_distance_dict)

def getFnaFromGB(fileName,outDir,fileID):
	outFileName = outDir + '/' + fileID + '.fna'
	handle = open(fileName)
	if os.path.exists(outFileName):
		os.remove(outFileName)
	SeqIO.convert(handle, 'genbank', outFileName, 'fasta')

if __name__=='__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--input_file',help='Path of the input file.\n')
	parser.add_argument('--input_protein',help='Path of the input protein file.\n')
	parser.add_argument('--output', help='Path of the output file.\n')
	parser.add_argument('--type', help='fasta or GenBank.\n')
	parser.add_argument('--out_dict', help='Path of the output result file in first step.\n')
	args = parser.parse_args()
	global type,strain_id,strain_def,root_path
	if args.input_file:
		input_file = args.input_file
	flag = 'no'
	if args.input_protein:
		pro_file = args.input_protein
		if not os.path.exists(pro_file):
			flag = 'yes'
	else:
		flag = 'yes'
	if args.output:
		outdir = args.output
	if args.out_dict:
		output_dict_file = args.out_dict
	else:
		output_dict_file = "no"
	root_path = get_root_path()
	file_name = os.path.basename(input_file).split()[0].split('.')[0]
	if args.type:
		type = args.type
		strain_id = file_name
		strain_def = file_name
	else:
		strain_id,strain_def,type = get_inf(input_file)
	if flag == 'yes':
		faa_prefix = os.path.join(outdir,file_name)
		pro_file = faa_prefix+'.faa'
		if type == 'fasta':
			pred_orf(input_file,faa_prefix)
		else:
			getFaaFromGB(input_file,outdir)
			getFnaFromGB(input_file,outdir,file_name)
			input_file = os.path.join(outdir,file_name+'.fna')
	root_path = get_root_path()
	if output_dict_file!='no':
		hosts_dict = np.load(output_dict_file).item()
		get_phage_hosts_codon_distance(input_file,pro_file,hosts_dict,outdir)
	else:
		get_phage_hosts_codon_distance_all(input_file,pro_file,outdir)
	print('results in:'+outdir)