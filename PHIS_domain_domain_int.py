import os,threading
import numpy.random.common
import numpy.random.bounded_integers
import numpy.random.entropy
import numpy as np
import os,re,argparse
from Bio import Entrez, SeqIO
from collections import Counter
import json
from Bio import SeqIO
def get_root_path():
	crispr_script_path = os.path.split(os.path.realpath(__file__))[0]
	root_path = os.path.dirname(crispr_script_path)
	return root_path

def diamond_blastp_pro_int(file,outfile,database,format,evalue,identity,coverage):
	num_threads = 20
	diamond_path = os.path.join(root_path,'software','diamond','diamond')
	if not os.path.exists(diamond_path):
		diamond_path = "diamond"
		print('%s does not exists!'%diamond_path)
	script = "diamond blastp -d "+database+" -q "+file+" -f "+str(format)+" -e "+str(evalue)+" -o "+outfile+" -p "+str(num_threads)+" --id "+str(identity)+" --query-cover "+str(coverage)
	os.system(script)

def pred_orf(fasta_file,faa_prefix):
	script = os.path.join(root_path,'software','FragGeneScan','run_FragGeneScan.pl')
	if not os.path.exists(diamond_path):
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

def annotate_domain_best(faa_file,xml_file,outfile):
	#cmd = 'hmmscan -E 1e-10 --tblout %s /zrom1/ganrui/PHIS/consensus_version/data/database/domain_domain_int/interacting_domain.hmm %s > /dev/null' % (xml_file, faa_file)
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

def phage_to_hosts_domain(phage_pro,outdir):
	pfam_int_domain_bac_dict_file = os.path.join(root_path,'database/db/profile','pfam_int_domain_bac_dict.npy')
	xml_file = os.path.join(outdir,'phage_domain.xml')
	phage_domain_file = os.path.join(outdir,'phage_domain.txt')
	resufile = os.path.join(outdir,'phage_int_bac_domain_result.txt')
	f_result = open(resufile,'w')
	f_result.write('phage_id\tphage_def\tphage_pro_id\tphage_domain_id\tphage_domain_name\tphage_domain_description\tbac_id\tbac_def\tbac_pro_id\tbac_domain_id\tbac_domain_name\tbac_domain_description'+'\n')	
	annotate_domain_best(phage_pro,xml_file,phage_domain_file)
	result_dict = {}
	result_dict_file = os.path.join(outdir,'phage_int_bac_domain_result_dict')
	if os.path.getsize(phage_domain_file)>0:
		pfam_int_domain_bac_dict = np.load(pfam_int_domain_bac_dict_file).item()
		with open(phage_domain_file) as f:
			contents = f.readlines()
		for line in contents:
			line = line.strip().split('\t')
			phage_pro = line[0]
			if type == 'fasta':
				phage_id = '_'.join(phage_pro.split('_')[0:-3]).split('.')[0]
			else:
				phage_id = phage_pro.split('|')[0].split('.')[0]
			try:
				phage_def = strain_info_dict[phage_id]
			except:
				phage_def = phage_id
			if phage_id not in result_dict.keys():
				result_dict.update({phage_id:{}})
			pfam_id = line[1]
			if pfam_id in pfam_int_domain_bac_dict.keys():
				int_bac_domain_infos = pfam_int_domain_bac_dict[pfam_id]
				for bac_id,int_bac_domains in int_bac_domain_infos.items():
					if bac_id.split('.')[0] not in result_dict[phage_id].keys():
						result_dict[phage_id].update({bac_id.split('.')[0]:[]})
					result_dict[phage_id][bac_id.split('.')[0]].append([line,int_bac_domains])
					for int_bac_domain in int_bac_domains:
						bac_id = int_bac_domain[2]
						bac_def = bac_inf_dict[bac_id.split('.')[0]]
	# 					f_result.write(phage_id+'\t'+phage_def+'\t'+'\t'.join(line[0:3])+'\t'+line[-1]+'\t'+bac_id+'\t'+bac_def+'\t'+int_bac_domain[1]+'\t'+int_bac_domain[0]+'\t'+int_bac_domain[3]+'\t'+int_bac_domain[-1]+'\n')
	# 					f_result.flush()
	# f_result.close()
	result_dict1 = {}
	for phage_id in result_dict1.keys():
		for bac_id in result_dict1[phage_id].keys():
			resu_line = result_dict1[phage_id][bac_id]
			int_bac_domains = resu_line[-1]
			line = resu_line[0]
			int_pair = []
			for resu in resu_line:
				phage_pro = resu[0][0]
				int_bac_pros = list(set([item[1] for item in resu[1]]))
				for int_bac_pro in int_bac_pros:
					if [phage_pro,int_bac_pro] not in int_pair:
						int_pair.append([phage_pro,int_bac_pro])
			if len(int_pair)>=min_DDI:
				if phage_id not in result_dict1.keys():
					result_dict1.update({phage_id:{}})
				if bac_id not in result_dict1[phage_id].keys():
					result_dict1[phage_id].update({bac_id:resu_line})
				for int_bac_domain in int_bac_domains:
					bac_id = int_bac_domain[2]
					bac_def = bac_inf_dict[bac_id]
					f_result.write(phage_id+'\t'+phage_def+'\t'+'\t'.join(line[0:3])+'\t'+line[-1]+'\t'+bac_id+'\t'+bac_def+'\t'+int_bac_domain[1]+'\t'+int_bac_domain[0]+'\t'+int_bac_domain[3]+'\t'+int_bac_domain[-1]+'\n')
					f_result.flush()
	f_result.close()	
	np.save(result_dict_file,result_dict1)

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
	parser.add_argument('--out_dict', help='Path of the output result file in first step.\n')
	parser.add_argument('--output', help='Path of the output file.\n')
	parser.add_argument('--type', help='Path of the output file.\n')
	parser.add_argument('--min_DDI', help='Minimal DDI number of a pair of phage-host pair(default:5)\n')

	args = parser.parse_args()
	global type,strain_info_dict,bac_inf_dict,root_path,min_DDI
	if args.input_file:
		input_file = args.input_file
	if args.out_dict:
		output_dict_file = args.out_dict
	else:
		output_dict_file = "no"
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
	if args.min_DDI:
		min_DDI = args.min_DDI
	else:
		min_DDI = 5
	if flag == 'yes':
		faa_prefix = os.path.join(outdir,file_name)
		pro_file = faa_prefix+'.faa'
		if type == 'fasta':
			pred_orf(input_file,faa_prefix)
		else:
			getFaaFromGB(input_file,outdir)	
	bac_inf_dict_file = os.path.join(root_path,'database/db/profile','bac_inf_dict.npy')
	bac_inf_dict = np.load(bac_inf_dict_file).item()	
	if output_dict_file=='no':
		print('error:%s does not exists!'%output_dict_file)
	else:
		hosts_dict = np.load(output_dict_file).item()
		phage_to_hosts_domain(pro_file,hosts_dict,outdir)