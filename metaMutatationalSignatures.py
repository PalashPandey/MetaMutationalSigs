# import sys, os
# welcome_message_str = """
# =================================================================
# MetaMutationalSigs
# Author: Palash Pandey (pp535@drexel.edu)
# Desc:
#   extract - extract signatures by either automatic or semi-automatic way.
# 			Of note, when you use manual way, you need to run 2 times,
# 			firstly you should set --manual to get signature estimation results,
# 			and secondly you should set --manual --number N to get N signatures.
# Usage:
#   MetaMutationalSigs --inputdir=<input_directoryectory> [--genome=<genome>]
# Options:
#   -h --help     Show help message.
#   -     Show help message.
#   -i <file>, --input <file>       input VCF directory path
#   -g <genome>, --genome <genome>  can be hg19, hg38 or mm10, [default: hg19].
# =================================================================" -> doc

# =================================================================
#   __  __      _        __  __       _        _   _                   _  _____ _
#  |  \/  |    | |      |  \/  |     | |      | | (_)                 | |/ ____(_)
#  | \  / | ___| |_ __ _| \  / |_   _| |_ __ _| |_ _  ___  _ __   __ _| | (___  _  __ _ ___
#  | |\/| |/ _ \ __/ _` | |\/| | | | | __/ _` | __| |/ _ \| '_ \ / _` | |\___ \| |/ _` / __|
#  | |  | |  __/ || (_| | |  | | |_| | || (_| | |_| | (_) | | | | (_| | |____) | | (_| \__ \
#  |_|  |_|\___|\__\__,_|_|  |_|\__,_|\__\__,_|\__|_|\___/|_| |_|\__,_|_|_____/|_|\__, |___/
# 																				 __/ |
# 																				|___/
# Name   :       MetaMutationalSigs
# Link   :       https://github.com/PalashPandey/MetaMutationalSigs
# Doc    :       https://github.com/PalashPandey/MetaMutationalSigs
# 
#  START
# """





"""MetaMutationalSigs
Usage:
	metaMutatationalSignatures.py --i=<file> [--output=<outdir>] [--genome=<genome>] [--browser] [--mutationalSignatures] [--sigflow] [--sigfit] [--deconstructSigs]
	metaMutatationalSignatures.py --i=<file> [--output=<outdir>]
	metaMutatationalSignatures.py -h | --help
	metaMutatationalSignatures.py --version
Options:
	-h --help     Show this screen.
	--version     Show version.
	-i <file>      Input directory with VCF files.
	-o <outdir>, --output <outdir>  output file path [default: ./metaMutationalSignatures_results.zip/].
	-g <genome>, --genome <genome>  genome build, can be GRCh37, GRCh38 [default: GRCh37].
	--browser    Open in browser.
	--mutationalSignatures                   run mutationalSignatures [default: TRUE]. 
	--sigflow                   run sigflow [default: TRUE].  
	--sigfit                   run sigfit [default: TRUE].  
	--deconstructSigs                   run deconstructSigs [default: TRUE].  

  __  __      _        __  __       _        _   _                   _  _____ _
 |  \/  |    | |      |  \/  |     | |      | | (_)                 | |/ ____(_)
 | \  / | ___| |_ __ _| \  / |_   _| |_ __ _| |_ _  ___  _ __   __ _| | (___  _  __ _ ___
 | |\/| |/ _ \ __/ _` | |\/| | | | | __/ _` | __| |/ _ \| '_ \ / _` | |\___ \| |/ _` / __|
 | |  | |  __/ || (_| | |  | | |_| | || (_| | |_| | (_) | | | | (_| | |____) | | (_| \__ \
 |_|  |_|\___|\__\__,_|_|  |_|\__,_|\__\__,_|\__|_|\___/|_| |_|\__,_|_|_____/|_|\__, |___/
																				 __/ |
																				|___/
Name   :       MetaMutationalSigs
Link   :       https://github.com/PalashPandey/MetaMutationalSigs
Doc    :       https://github.com/PalashPandey/MetaMutationalSigs


"""
from docopt import docopt
import os, time, subprocess, shutil, glob
if __name__ == '__main__':
		arguments = docopt(__doc__, version='Naval Fate 2.0')
		print(arguments)

		if arguments["--browser"]:
			subprocess.call(['python3.8', "flask_ui_app/app.py"])
		else:
			from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
			input_dir = arguments["--i"]
			output_dir = arguments["--o"]
			genome_ref = arguments["--genome"]
			input_dir = arguments["--i"]
			runMutationalPatterns = arguments["--mutationalSignatures"]
			runsigflow = arguments["--sigflow"]
			runsigfit = arguments["--sigfit"]
			runDeconstructSigs = arguments["--deconstructSigs"]

			matGen.SigProfilerMatrixGeneratorFunc("MetaMutationalSigs", genome_ref , input_dir)
			subprocess.call(['Rscript', "meta_sig_main_flask.R", input_dir , genome_ref , runMutationalPatterns , runsigflow, runsigfit, runDeconstructSigs])
			subprocess.call(['python3.8', "errors_pie_heatmap.py", input_dir   , runMutationalPatterns , runsigflow, runsigfit, runDeconstructSigs])

			shutil.make_archive(output_dir + "/metaMutationalSignatures_results", 'zip', input_dir)
			shutil.rmtree(input_dir)
