import numpy as np
import re,argparse,os
from collections import OrderedDict
import joblib

def get_root_path():
	crispr_script_path = os.path.split(os.path.realpath(__file__))[0]
	root_path = os.path.dirname(crispr_script_path)
	return root_path

def get_model_result(predict_data,outdir):
	root_path = get_root_path()
	bayes_bnb_model = joblib.load(os.path.join(root_path,'database/db/model/bayes_bnb.model'))
	bayes_gnb_model = joblib.load(os.path.join(root_path,'database/db/model/bayes_gnb.model'))
	randomforest_rfc_model = joblib.load(os.path.join(root_path,'database/db/model/randomforest_rfc.model'))
	randomforest_dtc_model = joblib.load(os.path.join(root_path,'database/db/model/randomforest_dtc.model'))
	svm_linear_model = joblib.load(os.path.join(root_path,'database/db/model/svm_linear.model'))
	svm_rbf_model = joblib.load(os.path.join(root_path,'database/db/model/svm_rbf.model'))
	logistic_model = joblib.load(os.path.join(root_path,'database/db/model/logistic.model'))
	bayes_bnb_result = bayes_bnb_model.predict(predict_data)
	bayes_bnb_result_probablity = bayes_bnb_model.predict_proba(predict_data)
	bayes_gnb_result = bayes_gnb_model.predict(predict_data)
	bayes_gnb_result_probablity = bayes_gnb_model.predict_proba(predict_data)
	randomforest_rfc_result = randomforest_rfc_model.predict(predict_data)
	randomforest_rfc_result_probablity = randomforest_rfc_model.predict_proba(predict_data)
	randomforest_dtc_result = randomforest_dtc_model.predict(predict_data)
	randomforest_dtc_result_probablity = randomforest_dtc_model.predict_proba(predict_data)
	svm_linear_result = svm_linear_model.predict(predict_data)
	svm_linear_result_probablity = svm_linear_model.predict_proba(predict_data)
	svm_rbf_result = svm_rbf_model.predict(predict_data)
	svm_rbf_result_probablity = svm_rbf_model.predict_proba(predict_data)
	logistic_result = logistic_model.predict(predict_data)
	logistic_result_probablity = logistic_model.predict_proba(predict_data)
	results = {'bayes_bnb':[bayes_bnb_result,bayes_bnb_result_probablity],
	'bayes_gnb':[bayes_gnb_result,bayes_gnb_result_probablity],
	'randomforest_rfc':[randomforest_rfc_result,randomforest_rfc_result_probablity],
	'randomforest_dtc':[randomforest_dtc_result,randomforest_dtc_result_probablity],
	'svm_rbf':[svm_rbf_result,svm_rbf_result_probablity],
	'svm_linear':[svm_linear_result,svm_linear_result_probablity],
	'logistic':[logistic_result,logistic_result_probablity]}
	outfile = os.path.join(outdir,'predict_result_dict')
	np.save(outfile,results)

if __name__=='__main__':	
	parser = argparse.ArgumentParser()
	parser.add_argument('--input',help='Path of the input file.\n')
	parser.add_argument('--output',help='Path of the input file.\n')
	args = parser.parse_args()
	if args.input:
		record_file = args.input
	if args.output:
		outdir = args.output
	hits_dict = np.load(record_file).item()
	crispr_file = os.path.join(outdir,'crispr_result.txt')
	prophage_file = os.path.join(outdir,'prophage_result.txt')
	pro_int_file = os.path.join(outdir,'pro_int_result.txt')
	blast_file = os.path.join(outdir,'blast_result.txt')
	genehomo_file = os.path.join(outdir,'genehomo_result.txt')

	#features_dict = OrderedDict()
	with open(crispr_file) as f:
		contents_crispr = f.readlines()
	with open(prophage_file) as f:
		contents_prophage = f.readlines()
	with open(pro_int_file) as f:
		contents_pro_int = f.readlines()
	with open(blast_file) as f:
		contents_blast = f.readlines()
	with open(genehomo_file) as f:
		contents_genehomo = f.readlines()
	feature_contents = [contents_crispr,contents_prophage,contents_pro_int,contents_blast,contents_genehomo]	
	feature_array = []
	feature_array_dict = {}
	for i in range(1,len(contents_crispr)):
		each_array = []
		signal_index = 0
		for array in feature_contents:
			line = array[i].strip().split('\t')
			query_id = line[0]
			query_def = line[1]
			hit_id = line[2]
			hit_def = line[3]
			if signal_index==2:
				#domain-domain
				if float(line[7])>5 and float(line[4])>0:
					c_array = line[4:]
				else:
					c_array = ['0']*6
			else:
				c_array = line[4:]
			each_array = each_array+[float(x) for x in c_array]
		feature_array.append(each_array)
		if query_id not in feature_array_dict.keys():
			feature_array_dict.update({query_id:{}})
		if hit_id not in feature_array_dict[query_id].keys():
			feature_array_dict[query_id].update({hit_id:each_array})
	if len(feature_array)>0:
		feature_array_dict_file = os.path.join(outdir,'feature_array_dict')
		np.save(feature_array_dict_file,feature_array_dict)
		get_model_result(feature_array,outdir)
	else:
		print('No result in crispr,prophage,pro_pro_int,blast,virhostmatcher,wish,condon!')
	# 	line = line.strip().split('\t')
	# 	query_id = line[0]
	# 	hit_id = line[2]
	# 	query_def = line[1]
	# 	if query_id not in features_dict.keys():
	# 		features_dict.update({query_id:{}})
	# 	if hit_id not in features_dict[query_id].keys():