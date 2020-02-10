import os,threading
import numpy as np
import os,re,argparse
from Bio import Entrez, SeqIO
from collections import Counter
import json
from Bio import SeqIO
import sys
def get_root_path():
	crispr_script_path = os.path.split(os.path.realpath(__file__))[0]
	root_path = os.path.dirname(crispr_script_path)
	return root_path

def diamond_blastp_pro_int(file,outfile,database,format,evalue,identity,coverage):
	num_threads = 20
	diamond_path = os.path.join(root_path,'software','diamond','diamond')
	if not os.path.exists(diamond_path):
		print('%s does not exists!'%diamond_path)
	else:
		diamond_path = "diamond"
	script = diamond_path+" blastp -d "+database+" -q "+file+" -f "+str(format)+" -e "+str(evalue)+" -o "+outfile+" -p "+str(num_threads)+" --id "+str(identity)+" --query-cover "+str(coverage)
	os.system(script)

def pred_orf(fasta_file,faa_prefix):
	script = os.path.join(root_path,'software','FragGeneScan','run_FragGeneScan.pl')
	if not os.path.exists(script):
		print('error:%s does not exists!'%script)
		sys.exit(-1)
	cmd_fragGeneScan = script+' -genome %s -out %s -complete=1 -train=complete -thread=20' % (fasta_file, faa_prefix)
	os.system(cmd_fragGeneScan)

def get_inf(file,outdir):
	with open(file) as f:
		contents = f.read().strip()
	if contents[0]=='>':
		type ='fasta'
		try:
			record = SeqIO.read(file,"fasta")
			strain_id = record.id
			strain_def = record.description
			strain_info_dict = {strain_id.split()[0].split('.')[0]:strain_def}
		except:
			strain_info_dict = get_strain_info(file,outdir)
	else:
		type="genbank"
		record = SeqIO.read(file,"genbank")
		strain_id = record.id
		strain_def = record.description
		strain_info_dict = {strain_id.split()[0].split('.')[0]:strain_def}
	return strain_info_dict,type

def get_strain_info(file,outdir):
	#fasta or multi-fasta
	strain_file = os.path.join(outdir,'strain_inf.txt')
	f_result = open(strain_file,'w')
	strain_inf_dict = {}
	with open(file) as f:
		contents = f.read().strip()
	if '\n>' in contents:		
		for strain in contents.split('\n>'):
			strain_title = strain.split('\n')[0].strip()
			strain_id = strain_title.split()[0].strip('>')
			strain_def = strain_title.strip('>')
			strain_inf_dict.update({strain_id.split('.')[0]:strain_def})
			f_result.write(strain_id+'\t'+strain_def+'\n')
			f_result.flush()
	else:
		strain_title = contents.split('\n')[0].strip()
		strain_id = strain_title.split()[0].strip('>')
		strain_def = strain_title.strip('>')
		strain_inf_dict.update({strain_id.split('.')[0]:strain_def})
		f_result.write(strain_id+'\t'+strain_def+'\n')
		f_result.flush()
	f_result.close()
	return strain_inf_dict

def filter_identity_coverage(blastp_file,filter_file):
	f_result = open(filter_file,'w')
	with open(blastp_file) as f:
		contents = f.readlines()
	for line in contents:
		line = line.strip().split('\t')
		query_pro = line[0]
		print(query_pro)
		hit_pro = line[1]
		identity = line[2]
		query_start = line[-6]
		query_end = line[-5]
		hit_length = abs(int(query_start)-int(query_end))+1
		if type=='fasta':
			pro_length = int((abs(int(query_pro.split('_')[-3])-int(query_pro.split('_')[-2]))+1)/3)
		else:
			pro_length = int((abs(int(query_pro.split('|')[1].split('_')[-3])-int(query_pro.split('|')[1].split('_')[-2]))+1)/3)
		coverage = float(hit_length)/pro_length
		if (coverage>=0.7) and ((float(identity)/100)>=0.4):
			f_result.write('\t'.join(line)+'\t'+str(coverage)+'\n')
			f_result.flush()
	f_result.close()

def phage_to_hosts(phage_pro_file,outdir):
	uniprot_int_homo_bac_dict_file = os.path.join(root_path,'database/db/profile','uniprot_int_homo_bac_dict.npy')
	uniprot_database = os.path.join(root_path,'database/db/database/uniprot_protein_small','uniprot_sprot.dmnd')
	phage_homolog_file = os.path.join(outdir,'phage_protein_homolog')
	diamond_blastp_pro_int(phage_pro_file,phage_homolog_file,uniprot_database,6,1,0.4,0.7)
	filter_file = os.path.join(outdir,'phage_protein_homolog_identity_0.4_coverage_0.7')
	filter_identity_coverage(phage_homolog_file,filter_file)
	result_dict = {}
	result_dict_file = os.path.join(outdir,'phage_int_bac_homo_result_dict')
	resufile = os.path.join(outdir,'phage_int_bac_homo_result.txt')
	f_result = open(resufile,'w')
	f_result.write('phage_id\tphage_def\tphage_pro_id\tphage_uniprot_pro_id\tphage_taxid\tbac_id\tbac_def\tbac_pro_id\tbac_uniprot_pro_id\tbac_taxid\tinteraction_method\tinteraction_type'+'\n')
	if os.path.exists(filter_file):
		with open(filter_file) as f:
			phage_homologs = f.readlines()
		if len(phage_homologs)==0:
			return 0
		else:
			uniprot_int_homo_bac_dict = np.load(uniprot_int_homo_bac_dict_file).item()
		for line in phage_homologs:
			line = line.strip().split('\t')
			phage_pro = line[0]
			if phage_id not in result_dict.keys():
				result_dict.update({phage_id:{}})
			uniprot_pro = line[1]
			uniprot_pro_id = uniprot_pro.split('|')[1]
			if uniprot_pro_id in uniprot_int_homo_bac_dict.keys():
				int_bacs = uniprot_int_homo_bac_dict[uniprot_pro_id]
				for bac_id,int_infos in int_bacs.items():
					flag = 0
					for int_info in int_infos:
						int_uniprot_id = int_info[0]
						if uniprot_pro_id==int_uniprot_id:
							continue
						flag = 1
						bac_id = int_info[2]
						bac_def = bac_inf_dict[bac_id.split('.')[0]]
						# f_result.write(phage_id+'\t'+phage_def+'\t'+phage_pro+'\t'+uniprot_pro_id+'\t'+int_info[-3]+'\t'+bac_id+'\t'+bac_def+'\t'+int_info[1]+'\t'+int_info[0]+'\t'+int_info[-2]+'\t'+int_info[-4]+'\t'+int_info[-1]+'\n')
						# f_result.flush()
					if flag==1:
						if bac_id.split('.')[0] not in result_dict[phage_id].keys():
							result_dict[phage_id].update({bac_id.split('.')[0]:[]})
						result_dict[phage_id][bac_id.split('.')[0]].append([line[0],uniprot_pro_id,line[2],int_infos])
	result_dict1 = {}
	for phage_id in result_dict.keys():
		for bac_id in result_dict[phage_id].keys():
			phage_id = phage_id.split('.')[0]
			try:
				phage_def = strain_info_dict[phage_id]
			except:
				phage_def = phage_id
			bac_id = bac_id.split('.')[0]
			bac_def = bac_inf_dict[bac_id]
			resu_line = result_dict[phage_id][bac_id]
			int_pair = []
			for resu in resu_line:
				phage_pro = resu[0]
				int_bac_pros = list(set([item[1] for item in resu[-1]]))
				int_number = int_number + len(int_bac_pros)
				for int_bac_pro in int_bac_pros:
					if [phage_pro,int_bac_pro] not in int_pair:
						int_pair.append([phage_pro,int_bac_pro])
			if len(int_pair)>=min_PPI:
				if phage_id not in result_dict1.keys():
					result_dict1.update({phage_id:{}})
				if bac_id not in result_dict1[phage_id].keys():
					result_dict1[phage_id].update({bac_id:resu_line})
				for resu in resu_line:
					phage_pro = resu[0]
					uniprot_pro_id = resu[1]
					int_infos = resu[-1]
					for int_info in int_infos:
						int_uniprot_id = int_info[0]
						f_result.write(phage_id+'\t'+phage_def+'\t'+phage_pro+'\t'+uniprot_pro_id+'\t'+int_info[-3]+'\t'+bac_id+'\t'+bac_def+'\t'+int_info[1]+'\t'+int_info[0]+'\t'+int_info[-2]+'\t'+int_info[-4]+'\t'+int_info[-1]+'\n')
						f_result.flush()
	f_result.close()
	np.save(result_dict_file,result_dict1)

def annotate_domain_best(faa_file,xml_file,outfile):
	hmmscan_path = os.path.join(root_path,'software/hmm','hmmscan')
	hmm_db = os.path.join(root_path,'database/db/database/domain_int','interacting_domain.hmm')
	cmd = hmmscan_path+' -E 1e-10 --tblout %s %s %s > /dev/null' % (xml_file,hmm_db,faa_file)	os.system(cmd)
	os.system(cmd)
	f1 = open(xml_file, 'r')
	genequery = {}
	for line in f1:
		if line[0] == '#': continue
		linetab = line.split()
		pfam_id = 'pfam'+linetab[1].split('.')[0].strip('PF')
		if linetab[2] not in genequery.keys():
			def_line = linetab[2]
			if '<' in def_line:
				def_line = def_line.replace('<', '')
			if '>' in def_line:
				def_line = def_line.replace('>', '')
			if 'join' in def_line:
				def_line = def_line.replace('join{', '')
				def_line = def_line.replace('}', '')
			genequery[linetab[2]] = []
			genequery[linetab[2]].append(linetab[0])  # 0 cas protein name
			genequery[linetab[2]].append(pfam_id)  #pfam id
			genequery[linetab[2]].append(def_line)  # 1 query name
			genequery[linetab[2]].append(linetab[4])  # 2 e-value
			genequery[linetab[2]].append(linetab[5])  # 3 score
			genequery[linetab[2]].append(linetab[6])  # 4 bias
			genequery[linetab[2]].append(' '.join(linetab[18:]))  # 5 descriprion
			continue
		src_evalue = float(genequery[linetab[2]][3])
		src_score = float(genequery[linetab[2]][4])
		src_bias = float(genequery[linetab[2]][5])
		dst_evalue = float(linetab[4])
		dst_score = float(linetab[5])
		dst_bias = float(linetab[6])
		if (dst_evalue < src_evalue):
			genequery[linetab[2]][0] = linetab[0]
			genequery[linetab[2]][1] = pfam_id
			genequery[linetab[2]][2] = linetab[2]
			genequery[linetab[2]][3] = linetab[4]
			genequery[linetab[2]][4] = linetab[5]
			genequery[linetab[2]][5] = linetab[6]
			genequery[linetab[2]][6] = ' '.join(linetab[18:])
		elif (dst_score > src_score and dst_evalue == src_evalue):
			genequery[linetab[2]][0] = linetab[0]
			genequery[linetab[2]][1] = pfam_id
			genequery[linetab[2]][2] = linetab[2]
			genequery[linetab[2]][3] = linetab[4]
			genequery[linetab[2]][4] = linetab[5]
			genequery[linetab[2]][5] = linetab[6]
			genequery[linetab[2]][6] = ' '.join(linetab[18:])
		elif (dst_bias < src_bias and dst_score == src_score and dst_evalue == src_evalue):
			genequery[linetab[2]][0] = linetab[0]
			genequery[linetab[2]][1] = pfam_id
			genequery[linetab[2]][2] = linetab[2]
			genequery[linetab[2]][3] = linetab[4]
			genequery[linetab[2]][4] = linetab[5]
			genequery[linetab[2]][5] = linetab[6]
			genequery[linetab[2]][6] = ' '.join(linetab[18:])
	f2 = open(outfile, 'w')
	# genequery = sorted(genequery.items(), key=lambda e: float(e[0].split('_')[-4]))
	for item in genequery:
		f2.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (genequery[item][2],genequery[item][1], genequery[item][0], str(genequery[item][3]), str(genequery[item][4]), str(genequery[item][5]),genequery[item][6]))
	f1.close()
	f2.close()

def annotate_domain(faa_file,xml_file,outfile):
	hmmscan_path = os.path.join(root_path,'software/hmm','hmmscan')
	hmm_db = os.path.join(root_path,'database/db/database/domain_int','interacting_domain.hmm')
	cmd = hmmscan_path+' -E 1e-10 --tblout %s %s %s > /dev/null' % (xml_file,hmm_db,faa_file)	os.system(cmd)
	os.system(cmd)
	f1 = open(xml_file, 'r')
	f2 = open(outfile, 'w')
	genequery = {}
	for line in f1:
		if line[0] == '#': continue
		linetab = line.split()
		def_line = linetab[2]
		if '<' in def_line:
			def_line = def_line.replace('<', '')
		if '>' in def_line:
			def_line = def_line.replace('>', '')
		if 'join' in def_line:
			def_line = def_line.replace('join{', '')
			def_line = def_line.replace('}', '')
		domain_name = linetab[0]  # 0 domain name
		pro_name = def_line  # 1 query name
		evalue = linetab[4]  # 2 e-value
		score = linetab[5]  # 3 score
		bias = linetab[6]  # 4 bias
		description = ' '.join(linetab[18:])  # 5 descriprion
		pfam_id = 'pfam'+linetab[1].split('.')[0].strip('PF')
		f2.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (pro_name,pfam_id,domain_name,evalue,score,bias,description))
		f2.flush()
	f1.close()
	f2.close()

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

if __name__=='__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--input_file',help='Path of the input file.\n')
	parser.add_argument('--input_protein',help='Path of the input file.\n')
	parser.add_argument('--output', help='Path of the output file.\n')
	parser.add_argument('--type', help='Path of the output file.\n')
	parser.add_argument('--min_PPI', help='Minimal PPI number of a pair of phage-host pair(default:1)\n')

	args = parser.parse_args()
	global type,strain_info_dict,bac_inf_dict,root_path,min_PPI
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
	strain_info_dict,type = get_inf(input_file,outdir)
	file_name = os.path.basename(input_file).split()[0].split('.')[0]
	if args.type:
		type = args.type
	if args.min_PPI:
		min_PPI = args.min_PPI
	else:
		min_PPI = 1
	if flag == 'yes':
		faa_prefix = os.path.join(outdir,file_name)
		pro_file = faa_prefix+'.faa'
		if type == 'fasta':
			pred_orf(input_file,faa_prefix)
		else:
			getFaaFromGB(input_file,outdir)
	root_path = get_root_path()
	bac_inf_dict_file = os.path.join(root_path,'database/db/profile','bac_inf_dict.npy')
	bac_inf_dict = np.load(bac_inf_dict_file).item()
	phage_to_hosts(pro_file,outdir)