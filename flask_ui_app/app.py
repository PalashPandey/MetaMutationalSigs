import os, time, subprocess, shutil, glob
from flask import Flask, flash, request, redirect, render_template, Response, send_file, send_from_directory
from werkzeug.utils import secure_filename
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen

app=Flask(__name__)

app.secret_key = "secret key"
app.config['MAX_CONTENT_LENGTH'] = 16 * 1024 * 1024
path = os.getcwd()
UPLOAD_FOLDER = os.path.join(path, 'uploads')
if not os.path.isdir(UPLOAD_FOLDER):
	os.mkdir(UPLOAD_FOLDER)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
ALLOWED_EXTENSIONS = set(['vcf'])
mutationalPattern = "FALSE"
sigflow = "FALSE"
sigfit = "FALSE"
deconstructSigs = "FALSE"


def allowed_file(filename):
	return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

@app.route('/')
def upload_form():
	return render_template('index.html')

@app.route('/', methods=['POST'])
def upload_file():
	if request.method == 'POST':
		if 'files[]' not in request.files:
			flash('No file part')
			return redirect(request.url)
		files = request.files.getlist('files[]')
		whichtorun_list = request.form.getlist('whichToRun')
		print(whichtorun_list)
		if "runMutationalPatterns" in whichtorun_list   :
			global mutationalPattern
			mutationalPattern = "TRUE"   
		if "runSigflow" in whichtorun_list   :
			global sigflow
			sigflow = "TRUE"
		if "runSigfit" in whichtorun_list :
			global sigfit
			sigfit = "TRUE"
		if "runDeconstructSigs" in whichtorun_list :
			global deconstructSigs
			deconstructSigs = "TRUE"   
		print(mutationalPattern ,sigflow, sigfit, deconstructSigs )
		for file in files:
			if file and allowed_file(file.filename):
				filename = secure_filename(file.filename)
				file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
		# time.sleep(10) 
		return redirect('/upload')

@app.route('/upload')
def render_after_upload():
	return render_template('upload.html')


@app.route('/progress')
def progress():
	def generate():
		x = 1
		yield "data:" + str(x) + "\n\n"
		if glob.glob("./uploads/*.vcf"):
			matGen.SigProfilerMatrixGeneratorFunc("MetaMutationalSigs",'GRCh37' , "/uploads")
			x = x + 33
			print(x)
			yield "data:" + str(x) + "\n\n"
			subprocess.call(['C:\\Program Files\\R\\R-4.0.2\\bin\\Rscript', "\\metaSignatures\\meta_sig_main_flask.R", "metaSignatures\\flaskmultiplefileupload\\uploads" , "hg19" , mutationalPattern , sigflow, sigfit, deconstructSigs])
			x = x + 33
			print(x)
			yield "data:" + str(x) + "\n\n"
			subprocess.call(['C:\\Users\\pande\\Anaconda3\\python.exe ', "\\metaSignatures\\errors_pie_heatmap.py", "\\metaSignatures\\flaskmultiplefileupload\\uploads"   , mutationalPattern , sigflow, sigfit, deconstructSigs])
			x = x + 33
			print(x)
			yield "data:" + str(x) + "\n\n"
			shutil.make_archive("metaMutationalSignatures_results", 'zip', "\\metaSignatures\\flaskmultiplefileupload\\uploads")
			shutil.rmtree("C:\\Users\\pande\\OneDriveDrexelUniversity\\Documents\\Fall-2021\\Coop\\CGC\\SanjeeVCFFiles\\PLOS_review_paper\\metaSignatures\\flaskmultiplefileupload\\uploads")
			if not os.path.isdir(app.config['UPLOAD_FOLDER']):
				os.mkdir(app.config['UPLOAD_FOLDER'])
		else:
			pass
	return Response(generate(), mimetype= 'text/event-stream')

@app.route('/dummy')
def dummy():
	return render_template("dummy.html")

@app.route('/download')
def download():
	print("download")
	return send_file("C:\\Users\\pande\\OneDriveDrexelUniversity\\Documents\\Fall-2021\\Coop\\CGC\\SanjeeVCFFiles\\PLOS_review_paper\\metaSignatures\\flaskmultiplefileupload\\metaMutationalSignatures_results.zip", attachment_filename="metaMutationalSignatures_results.zip")


if __name__ == "__main__":
	app.run(host='127.0.0.1',port=5000,debug=False,threaded=True)
