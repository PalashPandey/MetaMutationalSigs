## MetaMutationalSigs
![Logo](flask_ui_app\static\cover_photo.png)

Mutational signature analysis is very active and important area of interest. There are several packages available now for mutational signature analysis and they all use different approaches and give nontrivially different results. Because of the differences in their results, it is important for researchers to survey the available tools and make choose the one that best suits their application. There is a need for software that can aggregate the results from different packages and present them in a user friendly way so as to facilitate effective comparison. 

We created this package *MetaMutationalSigs* to facilitate comprehensive mutational signature analysis by creating a wrapper for different packages and providing a standard format for their outputs so that they can be effectively compared. We have also standardized the input formats accepted by various packages so ease interoperability. We also create standard visualizations for the results of all packages to ensure easy analysis. Our software is easy to install and use through Docker ,a package manager that automates the dependencies. 


## Install Using Docker

docker pull pp535/metamutationalsigs

docker run --rm -p 5001:5001 pp535/metamutationalsigs --browser

And then go to your browser at http://localhost:5001/ for the browser user interface.

docker run --rm -p 5001:5001 pp535/metamutationalsigs --browser

### Input: 
VCF files.



### Output: 

The output is returned as a compressed directory called `MetaMutationalResults`. Once uncompressed, this looks below. Directory `MetaMutationalResults` has the relevant results. 

![Image of Yaktocat](/markdown_images/fs_level_1.jpg)

Inside `MetaMutationalResults`, we can find a folder for each tool that was selected. We also have files 
![Image of Yaktocat](/markdown_images/fs_level_2.jpg)

Here is a summary of the files generated: 

|File Name                                      |Format|Description                                                                                          |
|-----------------------------------------------|------|-----------------------------------------------------------------------------------------------------|
|Heatmap_contributions_all_sigs_legacy.pdf      |pdf   |Contributions for all COSMIC Legacy SBS signatures to the overall signature.                         |
|Heatmap_contributions_all_sigs_SBS.pdf         |pdf   |Contributions for all COSMIC V3 SBS signatures.                                                      |
|Heatmap_COSMIC_legacy.pdf                      |pdf   |Heatmap for difference between the predicted contributions by different tools.  One for each sample. |
|Heatmap_COSMIC_V3.pdf                          |pdf   |Heatmap for difference between the predicted contributions by different tools.  One for each sample. |
|legacy_pie_charts.html                         |html  |Interactive pie charts of COSMIC legacy SBS contribution, per sample and for each tool.              |
|sbs_pie_charts.html                            |html  |Interactive pie charts of COSMIC V3 SBS signature contributions, per sample and for each tool.       |
|legacy_rmse_bar_plot.pdf                       |pdf   |Reconstruction error using COSMIC Legacy SBS signatures for each tool.                               |
|sbs_rmse_bar_plot.png                          |pdf   |Reconstruction error using COSMIC V3 SBS signatures for each tool.                                   |
|toolname_results\legacy_sample_error.csv       |csv   |Data used to create the bar plot.                                                                    |
|toolname_results\legacy_sample_contribution.csv|csv   |Data used to create heatmap and pie chart.                                                           |
|toolname_results\sbs_sample_error.csv          |csv   |Data used to create the bar plot.                                                                    |
|toolname_results\sbs_sample_contribution.csv   |csv   |Data used to create heatmap and pie chart.                                                           |



## Changelog
- V1 - COSMIC reference signatures V3.1 June 2020

## In progress
- V2 - COSMIC reference signatures updated to V3.2 March 2020