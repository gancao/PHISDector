import webbrowser

def create_results(result_file,templeate_file,outdir):
	write_html = os.path.join(outdir,"test.html")
	f_save = open(GEN_HTML,'w')
	header='''
		'''
	message = """
	<html>
	<head></head>
	<body>
	<p>%s</p>
	<p>%s</p>
	</body>
	</html>
	"""%(str1,str2)
	f.write(message) 
	f.close()
	webbrowser.open(GEN_HTML,new = 1)