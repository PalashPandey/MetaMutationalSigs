import pandas as pd 
import numpy as np
import math
from plotly.subplots import make_subplots
import plotly.graph_objects as go
from plotly.offline import plot
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.spatial.distance import squareform, pdist
import plotly.express as px
import sys

np.random.seed(1234)


def make_piecharts(data_df, n_rows, n_cols, fig_title):
	specs_list = []
	subplot_coordinate = []
	for x in range(n_rows):
		temp_spec_list = []
		for  y in range(n_cols):
			temp_spec_list.append({"type": "pie"})
			subplot_coordinate.append([x+1, y+1])
		specs_list.append(temp_spec_list)
	fig = make_subplots(rows=n_rows, cols=n_cols, specs=specs_list, subplot_titles= list(data_df["sample"]))
	for i in range(data_df.shape[0]):
		df = pd.DataFrame(data_df.iloc[i,1:])
		sample_name = data_df.iloc[i,:][0]
		df.rename( columns={i : sample_name}, inplace = True)
		fig.add_trace(go.Pie(labels=list(df[sample_name].index), values= df[sample_name].values, name=sample_name, textinfo= "none"), subplot_coordinate[i][0] , subplot_coordinate[i][1] )
	fig.update_layout(width = 2000, height = 2000, title_text=fig_title)
	return fig


def add_scatterplot_layer(data_df, _color , fig ): 
	try:
		data_df["sample"] = [ i.split("_")[0].split(".")[0] for i in data_df["sample"]]
	except:
		pass
	for i in range(len(data_df)):
		fig.add_trace(go.Scatter( name = data_df.iloc[i]["sample"] , legendgroup=data_df.iloc[i]["sample"] ,  marker=dict(color= _color), x = data_df.columns[1:], y = data_df.iloc[i, 1:], mode='markers') )

r_output_file_dir = sys.argv[1]


siglfow_legacy_errors_df = pd.read_csv(r_output_file_dir + "/sigflow/legacy_fitting_reconstruction_errors.csv")
sigflow_sbs_errors_df = pd.read_csv(r_output_file_dir+ "/sigflow/SBS_fitting_reconstruction_errors.csv")
siglfow_id_errors_df = pd.read_csv(r_output_file_dir + "/sigflow/legacy_fitting_reconstruction_errors.cs/v")
sigflow_dbs_errors_df = pd.read_csv(r_output_file_dir+ "/sigflow/SBS_fitting_reconstruction_errors.csv")


sigfit_legacy_errors_df = pd.read_csv(r_output_file_dir + "/sigfit_results/sample_errors_legacy.csv")
sigfit_sbs_errors_df = pd.read_csv(r_output_file_dir + "/sigfit_results/sample_errors_sbs.csv")
sigfit_id_errors_df = pd.read_csv(r_output_file_dir + "/sigfit_results/indel_sample_errors.csv")
sigfit_dbs_errors_df = pd.read_csv(r_output_file_dir + "/sigfit_results/dbs_sample_errors.csv")


mutational_pattern_strict_legacy_errors_df = pd.read_csv(r_output_file_dir +"/mutational_patterns_results/legacy_strict_sample_errors.csv")
mutational_pattern_strict_sbs_errors_df = pd.read_csv(r_output_file_dir + "/mutational_patterns_results/strict_sample_errors.csv")
mutational_pattern_legacy_errors_df = pd.read_csv(r_output_file_dir + "/mutational_patterns_results/legacy_sample_errors.csv")
mutational_pattern_sbs_errors_df = pd.read_csv(r_output_file_dir + "/mutational_patterns_results/sample_errors.csv")
mutational_pattern_id_errors_df = pd.read_csv(r_output_file_dir + "/mutational_patterns_results/id_sample_errors.csv")
mutational_pattern_dbs_errors_df = pd.read_csv(r_output_file_dir + "/mutational_patterns_results/dbs_sample_errors.csv")

deconstruct_legacy_errors_df = pd.read_csv(r_output_file_dir + "/deconstructsigs_results/legacy_sample_errors.csv")
deconstruct_sbs_errors_df = pd.read_csv(r_output_file_dir +"/deconstructsigs_results/sbs_sample_errors.csv")
deconstruct_id_errors_df = pd.read_csv(r_output_file_dir + "/deconstructsigs_results/indel_sample_errors.csv")
deconstruct_dbs_errors_df = pd.read_csv(r_output_file_dir +"/deconstructsigs_results/dbs_sample_errors.csv")

legacy_error_list = [siglfow_legacy_errors_df.iloc[:,1] ,  sigfit_legacy_errors_df.iloc[:,1] , mutational_pattern_strict_legacy_errors_df.iloc[:,1] ,  mutational_pattern_legacy_errors_df.iloc[:,1] , deconstruct_legacy_errors_df.iloc[:,1] ]
legacy_rmse_name_list = ["sigflow","sigfit" , "mutationalPatterns_strict", "mutationalPatterns", "deconstructSigs" ]

sbs_error_list = [sigflow_sbs_errors_df.iloc[:,1] ,  sigfit_sbs_errors_df.iloc[:,1] , mutational_pattern_strict_sbs_errors_df.iloc[:,1] , mutational_pattern_sbs_errors_df.iloc[:,1] ,  deconstruct_sbs_errors_df.iloc[:,1] ]
sbs_rmse_name_list = ["sigflow", "sigfit", "mutationalPatterns_strict", "mutationalPatterns" , "deconstructSigs"]

id_error_list = [siglfow_id_errors_df.iloc[:,1] ,  sigfit_id_errors_df.iloc[:,1] , mutational_pattern_id_errors_df.iloc[:,1] ,  deconstruct_id_errors_df.iloc[:,1] ]
id_rmse_name_list = ["sigflow", "sigfit",  "mutationalPatterns" , "deconstructSigs"]

dbs_error_list = [sigflow_dbs_errors_df.iloc[:,1] ,  sigfit_dbs_errors_df.iloc[:,1] , mutational_pattern_dbs_errors_df.iloc[:,1] ,  deconstruct_dbs_errors_df.iloc[:,1] ]
dbs_rmse_name_list = ["sigflow", "sigfit",  "mutationalPatterns" , "deconstructSigs"]


legacy_rmse_list = []
for i in legacy_error_list:
	MSE = np.square(i).mean() 
	RMSE = math.sqrt(MSE)
	legacy_rmse_list.append(RMSE)
	
legacy_rmse_df = pd.DataFrame(legacy_rmse_list, index=legacy_rmse_name_list)
legacy_rmse_df.columns = ["RMSE"]
legacy_rmse_df["toolname"] = legacy_rmse_name_list
legacy_rmse_df.to_csv(r_output_file_dir + "/legacy_rmse_data.csv")
sbs_rmse_list = []
for i in sbs_error_list:
	MSE = np.square(i).mean() 
	RMSE = math.sqrt(MSE)
	sbs_rmse_list.append(RMSE)
	
sbs_rmse_df = pd.DataFrame(sbs_rmse_list, index=sbs_rmse_name_list)        
sbs_rmse_df.columns = ["RMSE"]
sbs_rmse_df.sort_values(by= ["RMSE"])
sbs_rmse_df["toolname"] = sbs_rmse_name_list
sbs_rmse_df.to_csv(r_output_file_dir + "/sbs_rmse_data.csv")


id_rmse_list = []
for i in sbs_error_list:
	MSE = np.square(i).mean() 
	RMSE = math.sqrt(MSE)
	id_rmse_list.append(RMSE)
	
id_rmse_df = pd.DataFrame(id_rmse_list, index=id_rmse_name_list)        
id_rmse_df.columns = ["RMSE"]
id_rmse_df.sort_values(by= ["RMSE"])
id_rmse_df["toolname"] = id_rmse_name_list
id_rmse_df.to_csv(r_output_file_dir + "/id_rmse_data.csv")


dbs_rmse_list = []
for i in sbs_error_list:
	MSE = np.square(i).mean() 
	RMSE = math.sqrt(MSE)
	dbs_rmse_list.append(RMSE)
	
dbs_rmse_df = pd.DataFrame(dbs_rmse_list, index=dbs_rmse_name_list)        
dbs_rmse_df.columns = ["RMSE"]
dbs_rmse_df.sort_values(by= ["RMSE"])
dbs_rmse_df["toolname"] = dbs_rmse_name_list
dbs_rmse_df.to_csv(r_output_file_dir + "/dbs_rmse_data.csv")

fig = px.bar(legacy_rmse_df, x= "toolname" , y='RMSE', title="RMSE of reconstruction error using COSMIC legacy signatures")
fig.write_image(r_output_file_dir +"/legacy_rmse_bar_plot.pdf")

fig = px.bar(sbs_rmse_df, x="toolname", y='RMSE', title="RMSE of reconstruction error using COSMIC SBS signatures")
fig.write_image(r_output_file_dir + "/sbs_rmse_bar_plot.pdf")

artefacts_cosmic_sigs = ["COSMIC_1", "COSMIC_5", "COSMIC_12"]

sigflow_legacy_exposure_df = pd.read_csv(r_output_file_dir + "/sigflow/legacy_fitting_relative_exposure.csv")
# sigflow_legacy_exposure_df.drop(artefacts_cosmic_sigs, inplace= True , axis = 1)
# sigflow_legacy_exposure_df.drop(artefacts_cosmic_sigs, inplace= True , axis = 1)


artefacts_sbs_sigs = ["SBS1","SBS5","SBS27","SBS43","SBS45","SBS46","SBS47","SBS48","SBS49","SBS50","SBS51","SBS52","SBS53","SBS54","SBS55","SBS56","SBS57","SBS58","SBS59","SBS60"]

sigflow_sbs_exposure_df = pd.read_csv(r_output_file_dir + "/sigflow/SBS_fitting_relative_exposure.csv")
# sigflow_sbs_exposure_df.drop(artefacts_sbs_sigs, inplace= True , axis = 1)

sigflow_id_exposure_df = pd.read_csv(r_output_file_dir + "/sigflow/SBS_fitting_relative_exposure.csv")
sigflow_dbs_exposure_df = pd.read_csv(r_output_file_dir + "/sigflow/SBS_fitting_relative_exposure.csv")


mutationalPatterns_legacy_exposure_df = pd.read_csv(r_output_file_dir + "/mutational_patterns_results/legacy_sample_exposures.csv")
mutationalPatterns_legacy_exposure_df = mutationalPatterns_legacy_exposure_df.transpose()
mutationalPatterns_legacy_exposure_df.columns = mutationalPatterns_legacy_exposure_df.iloc[0]
mutationalPatterns_legacy_exposure_df = mutationalPatterns_legacy_exposure_df[1:]
mutationalPatterns_legacy_exposure_df.reset_index(inplace= True)
mutationalPatterns_legacy_exposure_df = mutationalPatterns_legacy_exposure_df.rename(columns={"index": "sample"})
mutationalPatterns_legacy_exposure_df[['Signature 1', 'Signature 2', 'Signature 3', 'Signature 4',
       'Signature 5', 'Signature 6', 'Signature 7', 'Signature 8',
       'Signature 9', 'Signature 10', 'Signature 11', 'Signature 12',
       'Signature 13', 'Signature 14', 'Signature 15', 'Signature 16',
       'Signature 17', 'Signature 18', 'Signature 19', 'Signature 20',
       'Signature 21', 'Signature 22', 'Signature 23', 'Signature 24',
       'Signature 25', 'Signature 26', 'Signature 27', 'Signature 28',
       'Signature 29', 'Signature 30']] = mutationalPatterns_legacy_exposure_df[['Signature 1', 'Signature 2', 'Signature 3', 'Signature 4',
       'Signature 5', 'Signature 6', 'Signature 7', 'Signature 8',
       'Signature 9', 'Signature 10', 'Signature 11', 'Signature 12',
       'Signature 13', 'Signature 14', 'Signature 15', 'Signature 16',
       'Signature 17', 'Signature 18', 'Signature 19', 'Signature 20',
       'Signature 21', 'Signature 22', 'Signature 23', 'Signature 24',
       'Signature 25', 'Signature 26', 'Signature 27', 'Signature 28',
       'Signature 29', 'Signature 30']].apply(lambda x: x/x.sum(), axis=1)


mutationalPatterns_strict_legacy_exposure_df = pd.read_csv(r_output_file_dir + "/mutational_patterns_results/legacy_strict_sample_exposures.csv")
mutationalPatterns_strict_legacy_exposure_df = mutationalPatterns_strict_legacy_exposure_df.transpose()
mutationalPatterns_strict_legacy_exposure_df.columns = mutationalPatterns_strict_legacy_exposure_df.iloc[0]
mutationalPatterns_strict_legacy_exposure_df = mutationalPatterns_strict_legacy_exposure_df[1:]
mutationalPatterns_strict_legacy_exposure_df.reset_index(inplace= True)
mutationalPatterns_strict_legacy_exposure_df = mutationalPatterns_strict_legacy_exposure_df.rename(columns={"index": "sample"})
mutationalPatterns_strict_legacy_exposure_df[['Signature 1', 'Signature 2', 'Signature 3', 'Signature 4',
       'Signature 5', 'Signature 6', 'Signature 7', 'Signature 8',
       'Signature 9', 'Signature 10', 'Signature 11', 'Signature 12',
       'Signature 13', 'Signature 14', 'Signature 15', 'Signature 16',
       'Signature 17', 'Signature 18', 'Signature 19', 'Signature 20',
       'Signature 21', 'Signature 22', 'Signature 23', 'Signature 24',
       'Signature 25', 'Signature 26', 'Signature 27', 'Signature 28',
       'Signature 29', 'Signature 30']] = mutationalPatterns_strict_legacy_exposure_df[['Signature 1', 'Signature 2', 'Signature 3', 'Signature 4',
       'Signature 5', 'Signature 6', 'Signature 7', 'Signature 8',
       'Signature 9', 'Signature 10', 'Signature 11', 'Signature 12',
       'Signature 13', 'Signature 14', 'Signature 15', 'Signature 16',
       'Signature 17', 'Signature 18', 'Signature 19', 'Signature 20',
       'Signature 21', 'Signature 22', 'Signature 23', 'Signature 24',
       'Signature 25', 'Signature 26', 'Signature 27', 'Signature 28',
       'Signature 29', 'Signature 30']].apply(lambda x: x/x.sum(), axis=1)



mutationalPatterns_sbs_exposure_df = pd.read_csv(r_output_file_dir + "/mutational_patterns_results/sample_exposures.csv")
mutationalPatterns_sbs_exposure_df = mutationalPatterns_sbs_exposure_df.transpose()
mutationalPatterns_sbs_exposure_df.columns = mutationalPatterns_sbs_exposure_df.iloc[0]
mutationalPatterns_sbs_exposure_df = mutationalPatterns_sbs_exposure_df[1:]
mutationalPatterns_sbs_exposure_df.reset_index(inplace= True)
mutationalPatterns_sbs_exposure_df = mutationalPatterns_sbs_exposure_df.rename(columns={"index": "sample"})
# mutationalPatterns_sbs_exposure_df = mutationalPatterns_sbs_exposure_df.div(mutationalPatterns_sbs_exposure_df.sum(axis=0), axis=0)	
# mutationalPatterns_sbs_exposure_df.drop(["SBS1", "SBS5"], inplace= True , axis = 1)
mutationalPatterns_sbs_exposure_df[['SBS2', 'SBS3', 'SBS4', 'SBS6', 'SBS7a', 'SBS7b', 'SBS7c',
	   'SBS7d', 'SBS8', 'SBS9', 'SBS10a', 'SBS10b', 'SBS11', 'SBS12', 'SBS13',
	   'SBS14', 'SBS15', 'SBS16', 'SBS17a', 'SBS17b', 'SBS18', 'SBS19',
	   'SBS20', 'SBS21', 'SBS22', 'SBS23', 'SBS24', 'SBS25', 'SBS26', 'SBS28',
	   'SBS29', 'SBS30', 'SBS31', 'SBS32', 'SBS33', 'SBS34', 'SBS35', 'SBS36',
	   'SBS37', 'SBS38', 'SBS39', 'SBS40', 'SBS41', 'SBS42', 'SBS44', 'SBS84',
	   'SBS85', 'SBS86', 'SBS87', 'SBS88', 'SBS89', 'SBS90']] = mutationalPatterns_sbs_exposure_df[['SBS2', 'SBS3', 'SBS4', 'SBS6', 'SBS7a', 'SBS7b', 'SBS7c',
	   'SBS7d', 'SBS8', 'SBS9', 'SBS10a', 'SBS10b', 'SBS11', 'SBS12', 'SBS13',
	   'SBS14', 'SBS15', 'SBS16', 'SBS17a', 'SBS17b', 'SBS18', 'SBS19',
	   'SBS20', 'SBS21', 'SBS22', 'SBS23', 'SBS24', 'SBS25', 'SBS26', 'SBS28',
	   'SBS29', 'SBS30', 'SBS31', 'SBS32', 'SBS33', 'SBS34', 'SBS35', 'SBS36',
	   'SBS37', 'SBS38', 'SBS39', 'SBS40', 'SBS41', 'SBS42', 'SBS44', 'SBS84',
	   'SBS85', 'SBS86', 'SBS87', 'SBS88', 'SBS89', 'SBS90']].apply(lambda x: x/x.sum(), axis=1)
# mutationalPatterns_sbs_exposure_df.drop(artefacts_sbs_sigs, inplace= True , axis = 1, errors='ignore')


mutationalPatterns_strict_sbs_exposure_df = pd.read_csv(r_output_file_dir + "/mutational_patterns_results/strict_sample_exposures.csv")
mutationalPatterns_strict_sbs_exposure_df = mutationalPatterns_strict_sbs_exposure_df.transpose()
mutationalPatterns_strict_sbs_exposure_df.columns = mutationalPatterns_strict_sbs_exposure_df.iloc[0]
mutationalPatterns_strict_sbs_exposure_df = mutationalPatterns_strict_sbs_exposure_df[1:]
mutationalPatterns_strict_sbs_exposure_df.reset_index(inplace= True)
mutationalPatterns_strict_sbs_exposure_df = mutationalPatterns_strict_sbs_exposure_df.rename(columns={"index": "sample"})
# mutationalPatterns_strict_sbs_exposure_df = mutationalPatterns_strict_sbs_exposure_df.div(mutationalPatterns_strict_sbs_exposure_df.sum(axis=0), axis=0)
# mutationalPatterns_strict_sbs_exposure_df.drop(["SBS1", "SBS5"], inplace= True , axis = 1)
mutationalPatterns_strict_sbs_exposure_df[['SBS2', 'SBS3', 'SBS4', 'SBS6', 'SBS7a', 'SBS7b', 'SBS7c',
	   'SBS7d', 'SBS8', 'SBS9', 'SBS10a', 'SBS10b', 'SBS11', 'SBS12', 'SBS13',
	   'SBS14', 'SBS15', 'SBS16', 'SBS17a', 'SBS17b', 'SBS18', 'SBS19',
	   'SBS20', 'SBS21', 'SBS22', 'SBS23', 'SBS24', 'SBS25', 'SBS26', 'SBS28',
	   'SBS29', 'SBS30', 'SBS31', 'SBS32', 'SBS33', 'SBS34', 'SBS35', 'SBS36',
	   'SBS37', 'SBS38', 'SBS39', 'SBS40', 'SBS41', 'SBS42', 'SBS44', 'SBS84',
	   'SBS85', 'SBS86', 'SBS87', 'SBS88', 'SBS89', 'SBS90']] = mutationalPatterns_strict_sbs_exposure_df[['SBS2', 'SBS3', 'SBS4', 'SBS6', 'SBS7a', 'SBS7b', 'SBS7c',
	   'SBS7d', 'SBS8', 'SBS9', 'SBS10a', 'SBS10b', 'SBS11', 'SBS12', 'SBS13',
	   'SBS14', 'SBS15', 'SBS16', 'SBS17a', 'SBS17b', 'SBS18', 'SBS19',
	   'SBS20', 'SBS21', 'SBS22', 'SBS23', 'SBS24', 'SBS25', 'SBS26', 'SBS28',
	   'SBS29', 'SBS30', 'SBS31', 'SBS32', 'SBS33', 'SBS34', 'SBS35', 'SBS36',
	   'SBS37', 'SBS38', 'SBS39', 'SBS40', 'SBS41', 'SBS42', 'SBS44', 'SBS84',
	   'SBS85', 'SBS86', 'SBS87', 'SBS88', 'SBS89', 'SBS90']].apply(lambda x: x/x.sum(), axis=1)
# mutationalPatterns_strict_sbs_exposure_df.drop(artefacts_sbs_sigs, inplace= True , axis = 1, errors='ignore')

mutationalPatterns_id_exposure_df = pd.read_csv(r_output_file_dir + "/mutational_patterns_results/legacy_sample_exposures.csv")
mutationalPatterns_id_exposure_df = mutationalPatterns_id_exposure_df.transpose()
mutationalPatterns_id_exposure_df.columns = mutationalPatterns_id_exposure_df.iloc[0]
mutationalPatterns_id_exposure_df = mutationalPatterns_id_exposure_df[1:]
mutationalPatterns_id_exposure_df.reset_index(inplace= True)
mutationalPatterns_id_exposure_df = mutationalPatterns_id_exposure_df.rename(columns={"index": "sample"})
# mutationalPatterns_id_exposure_df[['Signature 1', 'Signature 2', 'Signature 3', 'Signature 4',
#        'Signature 5', 'Signature 6', 'Signature 7', 'Signature 8',
#        'Signature 9', 'Signature 10', 'Signature 11', 'Signature 12',
#        'Signature 13', 'Signature 14', 'Signature 15', 'Signature 16',
#        'Signature 17', 'Signature 18', 'Signature 19', 'Signature 20',
#        'Signature 21', 'Signature 22', 'Signature 23', 'Signature 24',
#        'Signature 25', 'Signature 26', 'Signature 27', 'Signature 28',
#        'Signature 29', 'Signature 30']] = mutationalPatterns_id_exposure_df[['Signature 1', 'Signature 2', 'Signature 3', 'Signature 4',
#        'Signature 5', 'Signature 6', 'Signature 7', 'Signature 8',
#        'Signature 9', 'Signature 10', 'Signature 11', 'Signature 12',
#        'Signature 13', 'Signature 14', 'Signature 15', 'Signature 16',
#        'Signature 17', 'Signature 18', 'Signature 19', 'Signature 20',
#        'Signature 21', 'Signature 22', 'Signature 23', 'Signature 24',
#        'Signature 25', 'Signature 26', 'Signature 27', 'Signature 28',
#        'Signature 29', 'Signature 30']].apply(lambda x: x/x.sum(), axis=1)


mutationalPatterns_dbs_exposure_df = pd.read_csv(r_output_file_dir + "/mutational_patterns_results/legacy_sample_exposures.csv")
mutationalPatterns_dbs_exposure_df = mutationalPatterns_dbs_exposure_df.transpose()
mutationalPatterns_dbs_exposure_df.columns = mutationalPatterns_dbs_exposure_df.iloc[0]
mutationalPatterns_dbs_exposure_df = mutationalPatterns_dbs_exposure_df[1:]
mutationalPatterns_dbs_exposure_df.reset_index(inplace= True)
mutationalPatterns_dbs_exposure_df = mutationalPatterns_dbs_exposure_df.rename(columns={"index": "sample"})
# mutationalPatterns_dbs_exposure_df[['Signature 1', 'Signature 2', 'Signature 3', 'Signature 4',
#        'Signature 5', 'Signature 6', 'Signature 7', 'Signature 8',
#        'Signature 9', 'Signature 10', 'Signature 11', 'Signature 12',
#        'Signature 13', 'Signature 14', 'Signature 15', 'Signature 16',
#        'Signature 17', 'Signature 18', 'Signature 19', 'Signature 20',
#        'Signature 21', 'Signature 22', 'Signature 23', 'Signature 24',
#        'Signature 25', 'Signature 26', 'Signature 27', 'Signature 28',
#        'Signature 29', 'Signature 30']] = mutationalPatterns_dbs_exposure_df[['Signature 1', 'Signature 2', 'Signature 3', 'Signature 4',
#        'Signature 5', 'Signature 6', 'Signature 7', 'Signature 8',
#        'Signature 9', 'Signature 10', 'Signature 11', 'Signature 12',
#        'Signature 13', 'Signature 14', 'Signature 15', 'Signature 16',
#        'Signature 17', 'Signature 18', 'Signature 19', 'Signature 20',
#        'Signature 21', 'Signature 22', 'Signature 23', 'Signature 24',
#        'Signature 25', 'Signature 26', 'Signature 27', 'Signature 28',
#        'Signature 29', 'Signature 30']].apply(lambda x: x/x.sum(), axis=1)


deconstructisgs_legacy_drop_list = ["Signature.1", "Signature.5", "Signature.12"]

deconstructSigs_legacy_exposure_df = pd.read_csv(r_output_file_dir + "/deconstructsigs_results/legacy_sample_exposures.csv")
deconstructSigs_legacy_exposure_df.rename( columns = {"Unnamed: 0": "sample"}, inplace= True)
# deconstructSigs_legacy_exposure_df.drop(deconstructisgs_legacy_drop_list, inplace= True , axis = 1)

deconstructSigs_sbs_exposure_df = pd.read_csv(r_output_file_dir + "/deconstructsigs_results/sbs_sample_exposures.csv")
deconstructSigs_sbs_exposure_df.rename( columns = {"Unnamed: 0": "sample"}, inplace= True)
# deconstructSigs_sbs_exposure_df.drop(artefacts_sbs_sigs, inplace= True , axis = 1)

deconstructSigs_id_exposure_df = pd.read_csv(r_output_file_dir + "/deconstructsigs_results/indel_sample_exposures.csv")
deconstructSigs_id_exposure_df.rename( columns = {"Unnamed: 0": "sample"}, inplace= True)
# deconstructSigs_sbs_exposure_df.drop(artefacts_sbs_sigs, inplace= True , axis = 1)

deconstructSigs_dbs_exposure_df = pd.read_csv(r_output_file_dir + "/deconstructsigs_results/sbs_sample_exposures.csv")
deconstructSigs_dbs_exposure_df.rename( columns = {"Unnamed: 0": "sample"}, inplace= True)
# deconstructSigs_sbs_exposure_df.drop(artefacts_sbs_sigs, inplace= True , axis = 1)

artefacts_cosmic_sigs_sigfit = ["mean.Signature.1", "mean.Signature.5", "mean.Signature.12"]

sigfit_legacy_exposure_df = pd.read_csv(r_output_file_dir + "/sigfit_results/sample_exposures_legacy.csv")
# sigfit_legacy_exposure_df.drop(artefacts_cosmic_sigs_sigfit, inplace= True , axis = 1)
sigfit_legacy_exposure_df.rename(columns={"Unnamed: 0": "sample" }, inplace= True)
sigfit_legacy_exposure_df = sigfit_legacy_exposure_df[["sample"] +[col for col in sigfit_legacy_exposure_df if col.startswith('mean')]]
sigfit_legacy_exposure_df.columns = ["sample"] +  [ "COSMIC_" + i.split(".")[-1] for i in sigfit_legacy_exposure_df.columns[1:]]


artefacts_sbs_sigs_sigfit = ["mean.SBS1","mean.SBS5","mean.SBS27","mean.SBS43","mean.SBS45","mean.SBS46","mean.SBS47","mean.SBS48","mean.SBS49","mean.SBS50","mean.SBS51","mean.SBS52","mean.SBS53","mean.SBS54","mean.SBS55","mean.SBS56","mean.SBS57","mean.SBS58","mean.SBS59","mean.SBS60"]
sigfit_sbs_exposure_df = pd.read_csv(r_output_file_dir + "/sigfit_results/sample_exposures_sbs.csv")
# sigfit_sbs_exposure_df.drop(artefacts_sbs_sigs_sigfit, inplace= True , axis = 1)
sigfit_sbs_exposure_df.rename(columns={"Unnamed: 0": "sample" }, inplace= True)
sigfit_sbs_exposure_df = sigfit_sbs_exposure_df[["sample"] +[col for col in sigfit_sbs_exposure_df if col.startswith('mean')]]
sigfit_sbs_exposure_df.columns = ["sample"] +  [i.split(".")[1] for i in sigfit_sbs_exposure_df.columns[1:]]

sigfit_id_exposure_df = pd.read_csv(r_output_file_dir + "/sigfit_results/sample_exposures_.csv")
sigfit_id_exposure_df.rename(columns={"Unnamed: 0": "sample" }, inplace= True)
sigfit_id_exposure_df = sigfit_id_exposure_df[["sample"] +[col for col in sigfit_id_exposure_df if col.startswith('mean')]]
sigfit_id_exposure_df.columns = ["sample"] +  [i.split(".")[1] for i in sigfit_id_exposure_df.columns[1:]]

sigfit_dbs_exposure_df = pd.read_csv(r_output_file_dir + "/sigfit_results/sample_exposures_.csv")
sigfit_dbs_exposure_df.rename(columns={"Unnamed: 0": "sample" }, inplace= True)
sigfit_dbs_exposure_df = sigfit_dbs_exposure_df[["sample"] +[col for col in sigfit_dbs_exposure_df if col.startswith('mean')]]
sigfit_dbs_exposure_df.columns = ["sample"] +  [i.split(".")[1] for i in sigfit_dbs_exposure_df.columns[1:]]


import math 
pie_chart_rows = math.ceil(sigflow_legacy_exposure_df.shape[0] / 3) 
pie_chart_cols = 3
sigflow_legacy_exposure_df_fig = make_piecharts(sigflow_legacy_exposure_df, pie_chart_rows, pie_chart_cols, "COSMIC V2 'Legacy' exposures data Sigflow")    
sigflow_sbs_exposure_df_fig = make_piecharts(sigflow_sbs_exposure_df, pie_chart_rows, pie_chart_cols, "COSMIC V3 'SBS' exposures data Sigflow")   

sbs_mutational_patters_fig = make_piecharts(mutationalPatterns_sbs_exposure_df, pie_chart_rows, pie_chart_cols, "COSMIC V3 'SBS' exposures  data MutationalPatterns")    
sbs_mutational_patters_strict_fig = make_piecharts(mutationalPatterns_strict_sbs_exposure_df, pie_chart_rows, pie_chart_cols, "COSMIC V3 'SBS' exposures data MutationalPatterns Strict")    

mutationalPatters_legacy_fig = make_piecharts(mutationalPatterns_legacy_exposure_df, pie_chart_rows, pie_chart_cols, "COSMIC V2 'Legacy' exposures  data MutationalPatterns")    
mutationalPatters_strict_legacy_fig = make_piecharts(mutationalPatterns_strict_legacy_exposure_df, pie_chart_rows, pie_chart_cols, "COSMIC V2 'Legacy' exposures data MutationalPatterns Strict")    


sigfit_legacy_exposure_df_fig = make_piecharts(sigfit_legacy_exposure_df, pie_chart_rows, pie_chart_cols, "COSMIC V2 'Legacy' exposures  data Sigfit")    
sigfit_sbs_exposure_df_fig = make_piecharts(sigfit_sbs_exposure_df, pie_chart_rows, pie_chart_cols, "COSMIC V3 'SBS' exposures data Sigfit")   

deconstructSigs_sbs_exposure_df_fig = make_piecharts(deconstructSigs_sbs_exposure_df, pie_chart_rows, pie_chart_cols, "COSMIC V3 'Legacy' exposures data DeconstructSigs")    
deconstructSigs_legacy_exposure_df_fig = make_piecharts(deconstructSigs_legacy_exposure_df, pie_chart_rows, pie_chart_cols, "COSMIC V2 'Legacy' exposures data DeconstructSigs")    


with open(r_output_file_dir + '/legacy_pie_charts.html', 'a') as f:
	f.write(sigflow_legacy_exposure_df_fig.to_html(full_html=False, include_plotlyjs='cdn'))
	f.write(mutationalPatters_legacy_fig.to_html(full_html=False, include_plotlyjs='cdn'))
	f.write(mutationalPatters_strict_legacy_fig.to_html(full_html=False, include_plotlyjs='cdn'))

	f.write(sigfit_legacy_exposure_df_fig.to_html(full_html=False, include_plotlyjs='cdn'))
	f.write(deconstructSigs_legacy_exposure_df_fig.to_html(full_html=False, include_plotlyjs='cdn'))

with open(r_output_file_dir + '/sbs_pie_charts.html', 'a') as f:
	f.write(sigflow_sbs_exposure_df_fig.to_html(full_html=False, include_plotlyjs='cdn'))
	f.write(sbs_mutational_patters_fig.to_html(full_html=False, include_plotlyjs='cdn'))        
	f.write(sbs_mutational_patters_strict_fig.to_html(full_html=False, include_plotlyjs='cdn'))           
	f.write(sigfit_sbs_exposure_df_fig.to_html(full_html=False, include_plotlyjs='cdn'))
	f.write(deconstructSigs_sbs_exposure_df_fig.to_html(full_html=False, include_plotlyjs='cdn'))



legacy_df_list = [ mutationalPatterns_legacy_exposure_df.sort_values(by=["sample"]) , mutationalPatterns_strict_legacy_exposure_df.sort_values(by=["sample"]) ,  sigfit_legacy_exposure_df.sort_values(by=["sample"]), sigflow_legacy_exposure_df.sort_values(by=["sample"]) ,  deconstructSigs_legacy_exposure_df.sort_values(by=["sample"]) ]
legacy_df_name_list = ["mutationalPatterns", "mutationalPatterns_strict",  "sigfit" , "sigflow" ,  "deconstructSigs"]
distance_df_list = []
for i in range(legacy_df_list[0].shape[0]):
	temp_df = []
	for index , j in enumerate(legacy_df_list):
		j.iloc[i,0] = legacy_df_name_list[index]
		temp_df.append(list(j.iloc[i,:]))
	_df = pd.DataFrame(temp_df, columns = legacy_df_list[0].columns)
	_df.fillna(0, inplace = True)
	distance_df_list.append(_df)

fig = plt.figure(figsize=(30, 30))
fig.subplots_adjust(hspace=2, wspace=2)
rows = 6
columns = 6
for i in range(len(distance_df_list)):
	df = distance_df_list[i]
	df__ = pd.DataFrame(squareform(pdist(df.iloc[:, 1:])), columns=df["sample"].unique(), index=df["sample"].unique())
	sns.set()
	ax = fig.add_subplot(rows, columns, i+1)    
	sns.heatmap(df__)
	plt.title("Heatmap legacy " + sigfit_legacy_exposure_df.iloc[i,0] )
plt.savefig(r_output_file_dir + "/Heatmap_legacy.pdf", dpi =500)



mutationalPatterns_sbs_exposure_df = mutationalPatterns_sbs_exposure_df.reindex(columns=['sample', 'SBS2', 'SBS3', 'SBS4', 'SBS6', 'SBS7a', 'SBS7b', 'SBS7c', 'SBS7d',
	   'SBS8', 'SBS9', 'SBS10a', 'SBS10b', 'SBS11', 'SBS12', 'SBS13', 'SBS14',
	   'SBS15', 'SBS16', 'SBS17a', 'SBS17b', 'SBS18', 'SBS19', 'SBS20',
	   'SBS21', 'SBS22', 'SBS23', 'SBS24', 'SBS25', 'SBS26', 'SBS28', 'SBS29',
	   'SBS30', 'SBS31', 'SBS32', 'SBS33', 'SBS34', 'SBS35', 'SBS36', 'SBS37',
	   'SBS38', 'SBS39', 'SBS40', 'SBS41', 'SBS42', 'SBS44', 'SBS84', 'SBS85',
	   'SBS86', 'SBS87', 'SBS88', 'SBS89', 'SBS90'])


sbs_df_list = [mutationalPatterns_sbs_exposure_df.sort_values(by=["sample"]), mutationalPatterns_strict_sbs_exposure_df.sort_values(by=["sample"]),  sigfit_sbs_exposure_df.sort_values(by=["sample"]), sigflow_sbs_exposure_df.sort_values(by=["sample"]), deconstructSigs_sbs_exposure_df.sort_values(by=["sample"]) ]
sbs_df_name_list = [ "mutationalPatterns", "mutationalPatterns_strict" ,  "sigfit" , "sigflow" , "deconstructSigs"]
distance_df_list = []
for i in range(sbs_df_list[0].shape[0]):
	temp_df = []
	for index , j in enumerate(sbs_df_list):
		j["sample"] = sbs_df_name_list[index]
		temp_df.append(list(j.iloc[i,:]))
	_df = pd.DataFrame(temp_df, columns = sbs_df_list[3].columns)
	_df.fillna(0, inplace= True)
	distance_df_list.append(_df)

fig = plt.figure(figsize=(30, 30))
fig.subplots_adjust(hspace=2, wspace=2)
rows = 6
columns = 6
for i in range(len(distance_df_list)):  
	df = distance_df_list[i]
	df_ = pd.DataFrame(squareform(pdist(df.iloc[:, 1:])), columns=df["sample"].unique(), index=df["sample"].unique())
	sns.set()
	ax = fig.add_subplot(rows, columns, i+1)    
	sns.heatmap(df_)
	plt.title("Heatmap SBS " + str(sigfit_sbs_exposure_df.iloc[i,0]))
plt.savefig(r_output_file_dir + "/Heatmap_SBS.pdf", dpi =500)



legacy_df_list = [ mutationalPatterns_legacy_exposure_df.sort_values(by=["sample"]) , mutationalPatterns_strict_legacy_exposure_df.sort_values(by=["sample"]) ,sigfit_legacy_exposure_df.sort_values(by=["sample"]), sigflow_legacy_exposure_df.sort_values(by=["sample"]),  deconstructSigs_legacy_exposure_df.sort_values(by=["sample"]) ]
legacy_df_name_list = [ "mutationalPatterns", "mutationalPatterns", "sigfit" , "sigflow" ,  "deconstruct"]
fig = plt.figure(figsize=( 120, 30 ))
fig.subplots_adjust(right = 0.9, wspace=.01)
rows = 6
columns = 6
for i in range(len(legacy_df_list)):
	df = legacy_df_list[i]
	sns.set()
	ax = fig.add_subplot(rows, columns, i+1)    
	sns.heatmap(df.set_index("sample") , cmap="YlGnBu" , xticklabels = True, yticklabels= True)
	sns.set(rc={'figure.figsize':(30, 30)})
	plt.title("Heatmap legacy " + legacy_df_name_list[i] )
plt.savefig(r_output_file_dir + "/Heatmap_exposures_all_sigs_legacy.pdf")




sbs_df_list = [mutationalPatterns_sbs_exposure_df.sort_values(by=["sample"]), sigfit_sbs_exposure_df.sort_values(by=["sample"]), sigflow_sbs_exposure_df.sort_values(by=["sample"]), deconstructSigs_sbs_exposure_df.sort_values(by=["sample"]) ]
sbs_df_name_list = [ "mutationalPatters",  "sigfit" , "sigflow" ,  "deconstructSigs"]
fig = plt.figure(figsize=(120, 30))
fig.subplots_adjust(right = 0.9, wspace=.01)
rows = 6
columns = 6
for i in range(len(sbs_df_list)):  
	df = sbs_df_list[i].set_index("sample").astype(float)
	sns.set()
	ax = fig.add_subplot(rows, columns, i+1)
	sns.heatmap(df , cmap="YlGnBu", xticklabels = True, yticklabels= True )
	sns.set(rc={'figure.figsize':(30, 30)})
	plt.title("Heatmap SBS " + str(sbs_df_name_list[i]))
plt.savefig(r_output_file_dir + "/Heatmap_exposures_all_sigs_SBS.pdf")
