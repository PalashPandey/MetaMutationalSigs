# from SigProfilerMatrixGenerator import install as genInstall
# genInstall.install('GRCh37')

from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
matrices = matGen.SigProfilerMatrixGeneratorFunc("Sigprofiler",'GRCh37' , "C:\\Users\\pande\\OneDrive - Drexel University\\Documents\\Fall-2021\\Coop\\CGC\\SanjeeVCFFiles\kidney_vcf\\indels\\plink")



# from sigproSS import spss 
# data = spss.importdata()
# spss.single_sample("C:/Users/pande/OneDrive - Drexel University/Documents/Fall-2021/Coop/CGC/SanjeeVCFFiles/rectal_vcfs/plink_vcf", "extraction_output",  ref="GRCh37", exome=False)