## MetaMutationalSigs

Motivation:
Mutational signature analysis is very active and important area of interest. There are several packages available now for mutational signature analysis and they all use different approaches and give nontrivially different results. Because of the differences in their results, it is important for researchers to survey the available tools and make choose the one that best suits their application. There is a need for software that can aggregate the results from different packages and present them in a user friendly way so as to facilitate effective comparison. 

Results: 
We created this package *MetaMutationalSigs* to facilitate comprehensive mutational signature analysis by creating a wrapper for different packages and providing a standard format for their outputs so that they can be effectively compared. We have also standardized the input formats accepted by various packages so ease interoperability. We also create standard visualizations for the results of all packages to ensure easy analysis. Our software is easy to install and use through Docker ,a package manager that automates the dependencies. 



Availability and implementation: https://github.com/PalashPandey/MetaMutationalSigs  
Contact: pp535@drexel.edu



Introduction: 

The basic idea behind mutational signatures is that mutational process create a specific patterns of mutations. Thus it follows that if one can identify these patterns in a given sample then they can essentially detect the corresponding mutational processes. 

Since there are a huge variety of mutations possible, they are grouped into 6 major types based on the base where the mutation was observed. These 6 mutation types are: C>A, C>G, C>T, T>A, T>C and T>G. Now these 6 types of mutations are further divided based on their location, i.e. other bases that are in their immediate proximity.

Mutational signature analysis involves multiple steps which require different amount of time and processing power. Typically one starts with BAM files which aligned to some reference genome and then proceeds to the variant calling step which outputs VCF files. These steps are usually very resource intensive and thus don't allow for much experimentation, the downstream steps of variant filtering and annotation are much less expensive. The final step is the actual mutational signature analysis, which is the least resource intensive and thus allows for experimentation with different methods. 

Approach: 

We chose signature refitting as our primary task and implemented high performing packages available in R language. We also 



Discussion: 

The massive increase in the number of software packages has made managing dependencies quite burdensome, this coupled with the absence of any standard data formats for signature matrices can make mutational signature analysis difficult and hard to reproduce. Our package provides and easy way of performing these tasks that takes care of the more mundane aspects of the analysis while ensuring reproducibility. 

## Install Using Docker

Input: 
VCF files/ MAF files. 

Commands:
Use MAF file as input and save the resulting data and visualizations in the output_dir directory    `` Rscript meta_sig_main.R input_data.maf hg19 ./output_dir`` 

Use VCF files as input and save the resulting data and visualizations in the output_dir directory    `` Rscript input_vcf_dir/ input_data.maf hg19 ./output_dir`` 

Use hg38 reference genome. Default is hg19. `` Rscript input_vcf_dir/ input_data.maf hg38 ./output_dir`` 
