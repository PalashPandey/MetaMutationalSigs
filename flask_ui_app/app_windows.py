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
		for file in files:
			if file and allowed_file(file.filename):
				filename = secure_filename(file.filename)
				file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
				print(filename)
		print("uplaod done ")
		# time.sleep(10) 
		return redirect('/upload')

@app.route('/upload')
def render_after_upload():
	return render_template('upload.html')



@app.route('/progress')
def progress():
	def generate():
		print(mutationalPattern , sigflow, sigfit, deconstructSigs)
		x = 1
		yield "data:" + str(x) + "\n\n"
		if glob.glob("./uploads/*.vcf"):
			# matGen.SigProfilerMatrixGeneratorFunc("MetaMutationalSigs",'GRCh37' , "/uploads")
			matGen.SigProfilerMatrixGeneratorFunc("MetaMutationalSigs",'GRCh37' , "C:\\Users\\pande\\OneDrive - Drexel University\\Documents\\Fall-2021\\Coop\\CGC\\SanjeeVCFFiles\\PLOS_review_paper\\metaSignatures\\flask_ui_app\\uploads")
			x = x + 33
			yield "data:" + str(x) + "\n\n"
			# subprocess.call(['Rscript', "meta_sig_main_flask.R", "./flask_ui_app/uploads" , "GRCh37" , mutationalPattern , sigflow, sigfit, deconstructSigs])
			subprocess.call(['C:\\Program Files\\R\\R-4.0.2\\bin\\Rscript', "C:\\Users\\pande\\OneDrive - Drexel University\\Documents\\Fall-2021\\Coop\\CGC\\SanjeeVCFFiles\\PLOS_review_paper\\metaSignatures\\meta_sig_main_flask_windows.R", "C:\\Users\\pande\\OneDrive - Drexel University\\Documents\\Fall-2021\\Coop\\CGC\\SanjeeVCFFiles\\PLOS_review_paper\\metaSignatures\\flask_ui_app\\uploads" , "GRCh37" , mutationalPattern , sigflow, sigfit, deconstructSigs])
			x = x + 33
			yield "data:" + str(x) + "\n\n"
			# subprocess.call(['python3.8', "./errors_pie_heatmap.py", "flask_ui_app/uploads"   , mutationalPattern , sigflow, sigfit, deconstructSigs])
			subprocess.call(['C:\\Users\\pande\\Anaconda3\\python.exe', "C:\\Users\\pande\\OneDrive - Drexel University\\Documents\\Fall-2021\\Coop\\CGC\\SanjeeVCFFiles\\PLOS_review_paper\\metaSignatures\\errors_pie_heatmap.py", "C:\\Users\\pande\\OneDrive - Drexel University\\Documents\\Fall-2021\\Coop\\CGC\\SanjeeVCFFiles\\PLOS_review_paper\\metaSignatures\\flask_ui_app\\uploads"   , mutationalPattern , sigflow, sigfit, deconstructSigs])
			x = x + 33
			yield "data:" + str(x) + "\n\n"
			# shutil.make_archive("metaMutationalSignatures_results", 'zip', "flask_ui_app/uploads")
			# shutil.rmtree("flask_ui_app/uploads")
			shutil.make_archive("metaMutationalSignatures_results", 'zip', "C:\\Users\\pande\\OneDrive - Drexel University\\Documents\\Fall-2021\\Coop\\CGC\\SanjeeVCFFiles\\PLOS_review_paper\\metaSignatures\\flask_ui_app\\uploads")
			shutil.rmtree("C:\\Users\\pande\\OneDrive - Drexel University\\Documents\\Fall-2021\\Coop\\CGC\\SanjeeVCFFiles\\PLOS_review_paper\\metaSignatures\\flask_ui_app\\uploads")
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
	return send_file("C:\\Users\\pande\\OneDrive - Drexel University\\Documents\\Fall-2021\\Coop\\CGC\\SanjeeVCFFiles\\PLOS_review_paper\\metaSignatures\\flask_ui_app\\metaMutationalSignatures_results.zip", attachment_filename="metaMutationalSignatures_results.zip")
	# return send_file("flask_ui_app/uploads" + "/metaMutationalSignatures_results.zip", attachment_filename="metaMutationalSignatures_results.zip")


if __name__ == "__main__":
	app.run(host='127.0.0.1',port=5000,debug=False,threaded=True)
