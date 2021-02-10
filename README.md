## MetaMutationalSigs

Motivation:
Mutational signature analysis is very active and important area of interest. There are several packages available now for mutational signature analysis and they all use different approaches and give nontrivially different results. Because of the differences in their results, it is important for researchers to survey the available tools and make choose the one that best suits their application. There is a need for software that can aggregate the results from different packages and present them in a user friendly way so as to facilitate effective comparison. 

Because of the variation in the results from these different packages it is valuable to survey them 

Results: 
We created this package *MetaMutationalSigs* to facilitate comprehensive mutational signature analysis by creating a wrapper for different packages and providing a standard format for their outputs so that they can be effectively compared. We have also standardized the input formats accepted by various packages so ease interoperability. We also create standard visualizations for the results of all packages to ensure easy analysis. Our software is easy to install and use through Docker ,a package manager that automates the dependencies. 



Availability and implementation: https://github.com/PalashPandey/MetaMutationalSigs  
Contact: pp535@drexel.edu



Introduction: 


Input: 
VCF files/ MAF files. 

Commands:
Use MAF file as input and save the resulting data and visualizations in the output_dir directory    `` Rscript meta_sig_main.R input_data.maf hg19 ./output_dir`` 

Use VCF files as input and save the resulting data and visualizations in the output_dir directory    `` Rscript input_vcf_dir/ input_data.maf hg19 ./output_dir`` 

Use hg38 reference genome. Default is hg19. `` Rscript input_vcf_dir/ input_data.maf hg38 ./output_dir`` 
