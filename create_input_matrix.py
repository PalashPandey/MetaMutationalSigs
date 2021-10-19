# from SigProfilerMatrixGenerator import install as genInstall
# genInstall.install('GRCh37')


from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
matGen.SigProfilerMatrixGeneratorFunc("MetaMutationalSigs",'GRCh37' , "C:\\Users\\pande\\OneDriveDrexelUniversity\\Documents\\Fall-2021\\Coop\\CGC\\SanjeeVCFFiles\\PLOS_review_paper\\metaSignatures\\flaskmultiplefileupload\\uploads" + user_file)

# matrices = matGen.SigProfilerMatrixGeneratorFunc("Sigprofiler",'GRCh37' , "C:\\Users\\pande\\OneDrive - Drexel University\\Documents\\Fall-2021\\Coop\\CGC\\SanjeeVCFFiles\kidney_vcf\\indels\\plink")

