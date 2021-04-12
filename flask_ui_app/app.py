import os, time, subprocess, shutil, glob, zipfile
from flask import Flask, flash, request, redirect, render_template, Response, send_file, send_from_directory
from werkzeug.utils import secure_filename
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
from SigProfilerMatrixGenerator import install as genInstall

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
genome_ref =  "GRCh37"

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
		global genome_ref
		genome_ref = request.form["genome"]
		print(genome_ref)
		if genome_ref ==  "GRCh38":
			genInstall.install('GRCh38', rsync=False, bash=True)
		if genome_ref ==  "GRCm37":
			genInstall.install('GRCm37', rsync=False, bash=True)
		if genome_ref ==  "GRCm38.p6":
			genInstall.install('GRCm38.p6', rsync=False, bash=True)
		if genome_ref ==  "Rnor_6.0":
			genInstall.install('Rnor_6.0', rsync=False, bash=True)
		if runMutationalPatterns ==  False:
			runMutationalPatterns = "TRUE"
		if runMutationalPatterns ==  True:
			runMutationalPatterns = "FALSE"

		if "runMutationalPatterns" in whichtorun_list:
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
		# time.sleep(10) 
		print("files succcessfully uploaded ")
		return redirect('/upload')

@app.route('/upload')
def render_after_upload():
	return render_template('upload.html')


def zipdir(path, ziph):
    # ziph is zipfile handle
    for root, dirs, files in os.walk(path):
        for file in files:
            ziph.write(os.path.join(root, file))


@app.route('/progress')
def progress():
	def generate():
		print(genome_ref , mutationalPattern , sigflow, sigfit, deconstructSigs)
		x = 1
		yield "data:" + str(x) + "\n\n"
		if glob.glob("/app/flask_ui_app/uploads/*.vcf"):
			matGen.SigProfilerMatrixGeneratorFunc("MetaMutationalSigs",'GRCh37' , "/app/flask_ui_app/uploads")
			x = x + 33
			yield "data:" + str(x) + "\n\n"
			subprocess.call(['Rscript', "../meta_sig_main_flask.r", "/app/flask_ui_app/uploads" , genome_ref , mutationalPattern , sigflow, sigfit, deconstructSigs])
			x = x + 33
			yield "data:" + str(x) + "\n\n"
			subprocess.call(['python3.8', "../errors_pie_heatmap.py", "/app/flask_ui_app/uploads"   , mutationalPattern , sigflow, sigfit, deconstructSigs])

			shutil.rmtree("/app/flask_ui_app/uploads" + "/input")
			shutil.rmtree("/app/flask_ui_app/uploads" + "/logs")
			shutil.rmtree("/app/flask_ui_app/uploads" + "/output")

			files_in_directory = os.listdir("/app/flask_ui_app/uploads")

			filtered_files = [file for file in files_in_directory if file.endswith(".vcf")]

			for file in filtered_files:
				path_to_file = os.path.join("/app/flask_ui_app/uploads", file)
				os.remove(path_to_file)

			zipf = zipfile.ZipFile("/app/flask_ui_app/metaMutationalSignatures_results.zip", 'w', zipfile.ZIP_DEFLATED)
			os.chdir("/app/flask_ui_app/")
			zipdir("./uploads/", zipf)
			zipf.close()


			shutil.rmtree("/app/flask_ui_app/uploads")
			if not os.path.isdir(app.config['UPLOAD_FOLDER']):
				os.mkdir(app.config['UPLOAD_FOLDER'])
			x = x + 33
			yield "data:" + str(x) + "\n\n"

		else:
			pass
	return Response(generate(), mimetype= 'text/event-stream')

@app.route('/dummy')
def dummy():
	return render_template("dummy.html")

@app.route('/download')
def download():
	return send_file("/app/flask_ui_app/" + "/metaMutationalSignatures_results.zip", attachment_filename="metaMutationalSignatures_results.zip")


if __name__ == "__main__":
	app.run(host='0.0.0.0',port=5001,debug=False,threaded=True)
