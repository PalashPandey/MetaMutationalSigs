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
# r_output_file_dir = "MetaMutationalResults"

def run_legacy(sigfit = True, sigflow = True, deconstructSigs= True, mutationalPattern = True):

    try:
        deconstructisgs_legacy_drop_list = ["Signature.1", "Signature.5", "Signature.12"]
        artefacts_cosmic_sigs_sigfit = ["mean.Signature.1", "mean.Signature.5", "mean.Signature.12"]
        artefacts_sbs_sigs_sigfit = ["mean.SBS1","mean.SBS5","mean.SBS27","mean.SBS43","mean.SBS45","mean.SBS46","mean.SBS47","mean.SBS48","mean.SBS49","mean.SBS50","mean.SBS51","mean.SBS52","mean.SBS53","mean.SBS54","mean.SBS55","mean.SBS56","mean.SBS57","mean.SBS58","mean.SBS59","mean.SBS60"]        
        legacy_df_list = []
        legacy_df_name_list = []
        if sigflow:
            sigflow_legacy_exposure_df = pd.read_csv(r_output_file_dir + "/sigflow/legacy_fitting_relative_exposure.csv")
            pie_chart_rows = math.ceil(sigflow_legacy_exposure_df.shape[0] / 3) 
            pie_chart_cols = 3
            sigflow_legacy_exposure_df_fig = make_piecharts(sigflow_legacy_exposure_df, pie_chart_rows, pie_chart_cols, "COSMIC V2 'Legacy' exposures data Sigflow")    
            legacy_df_list.append(sigflow_legacy_exposure_df.sort_values(by=["sample"]))
            legacy_df_name_list.append("sigflow")
        if mutationalPattern:
            mutationalPatterns_legacy_exposure_df = pd.read_csv(r_output_file_dir + "/mutational_patterns_results/legacy_sample_exposures.csv")
            mutationalPatterns_legacy_exposure_df = mutationalPatterns_legacy_exposure_df.rename(columns={"Unnamed: 0": "sample"})
            mutationalPatterns_legacy_exposure_df.set_index("sample", inplace = True)
            mutationalPatterns_legacy_exposure_df = mutationalPatterns_legacy_exposure_df.div(mutationalPatterns_legacy_exposure_df.sum(axis=1), axis=0)
            mutationalPatterns_legacy_exposure_df.reset_index(inplace= True)

            mutationalPatterns_strict_legacy_exposure_df = pd.read_csv(r_output_file_dir + "/mutational_patterns_results/legacy_strict_sample_exposures.csv")
            mutationalPatterns_strict_legacy_exposure_df = mutationalPatterns_strict_legacy_exposure_df.rename(columns={"Unnamed: 0": "sample"})
            mutationalPatterns_strict_legacy_exposure_df.set_index("sample", inplace = True)
            mutationalPatterns_strict_legacy_exposure_df = mutationalPatterns_strict_legacy_exposure_df.div(mutationalPatterns_strict_legacy_exposure_df.sum(axis=1), axis=0)
            mutationalPatterns_strict_legacy_exposure_df.reset_index(inplace= True)

            mutationalPatters_legacy_fig = make_piecharts(mutationalPatterns_legacy_exposure_df, pie_chart_rows, pie_chart_cols, "COSMIC V2 'Legacy' exposures  data MutationalPatterns")    
            mutationalPatters_strict_legacy_fig = make_piecharts(mutationalPatterns_strict_legacy_exposure_df, pie_chart_rows, pie_chart_cols, "COSMIC V2 'Legacy' exposures data MutationalPatterns Strict")                
            legacy_df_list.append(mutationalPatterns_legacy_exposure_df.sort_values(by=["sample"]))
            legacy_df_name_list.append("mutationalPatterns")
            legacy_df_list.append(mutationalPatterns_strict_legacy_exposure_df.sort_values(by=["sample"]))
            legacy_df_name_list.append("mutationalPatterns_strict")
        if deconstructSigs:
            deconstructSigs_legacy_exposure_df = pd.read_csv(r_output_file_dir + "/deconstructsigs_results/legacy_sample_exposures.csv")
            deconstructSigs_legacy_exposure_df.rename( columns = {"Unnamed: 0": "sample"}, inplace= True)
            deconstructSigs_legacy_exposure_df_fig = make_piecharts(deconstructSigs_legacy_exposure_df, pie_chart_rows, pie_chart_cols, "COSMIC V2 'Legacy' exposures data DeconstructSigs")    
            legacy_df_list.append(deconstructSigs_legacy_exposure_df.sort_values(by=["sample"]))
            legacy_df_name_list.append("deconstructSigs")
        if sigfit:
            sigfit_legacy_exposure_df = pd.read_csv(r_output_file_dir + "/sigfit_results/sample_exposures_legacy.csv")
            sigfit_legacy_exposure_df.rename(columns={"Unnamed: 0": "sample" }, inplace= True)
            sigfit_legacy_exposure_df.columns = ["sample"] +  [ "COSMIC_" + i.split(".")[-1] for i in sigfit_legacy_exposure_df.columns[1:]]
            sigfit_legacy_exposure_df_fig = make_piecharts(sigfit_legacy_exposure_df, pie_chart_rows, pie_chart_cols, "COSMIC V2 'Legacy' exposures  data Sigfit")    
            legacy_df_list.append(sigfit_legacy_exposure_df.sort_values(by=["sample"]))
            legacy_df_name_list.append("sigfit")

        with open(r_output_file_dir + '/legacy_pie_charts.html', 'a') as f:
            if sigflow:
                f.write(sigflow_legacy_exposure_df_fig.to_html(full_html=False, include_plotlyjs='cdn'))
            if mutationalPattern:
                f.write(mutationalPatters_legacy_fig.to_html(full_html=False, include_plotlyjs='cdn'))
                f.write(mutationalPatters_strict_legacy_fig.to_html(full_html=False, include_plotlyjs='cdn'))
            if sigfit:
                f.write(sigfit_legacy_exposure_df_fig.to_html(full_html=False, include_plotlyjs='cdn'))
            if deconstructSigs:
                f.write(deconstructSigs_legacy_exposure_df_fig.to_html(full_html=False, include_plotlyjs='cdn'))

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
        rows = math.ceil(sigflow_legacy_exposure_df.shape[0] / 6)
        columns = 6
        for i in range(len(distance_df_list)):
            df = distance_df_list[i]
            df__ = pd.DataFrame(squareform(pdist(df.iloc[:, 1:])), columns=df["sample"].unique(), index=df["sample"].unique())
            sns.set()
            ax = fig.add_subplot(rows, columns, i+1)    
            sns.heatmap(df__)
            plt.title("Heatmap legacy " + sigflow_legacy_exposure_df.iloc[i,0] )
        plt.savefig(r_output_file_dir + "/Heatmap_legacy.pdf", dpi =500)

        fig = plt.figure(figsize=( 120, 30 ))
        fig.subplots_adjust(right = 0.9, wspace=.01)
        rows = 1
        columns = 5
        for i in range(len(legacy_df_list)):
            df = legacy_df_list[i]
            sns.set()
            ax = fig.add_subplot(rows, columns, i+1)    
            sns.heatmap(df.set_index("sample") , cmap="YlGnBu" , xticklabels = True, yticklabels= True)
            sns.set(rc={'figure.figsize':(30, 30)})
            plt.title("Heatmap legacy " + legacy_df_name_list[i] )
        plt.savefig(r_output_file_dir + "/Heatmap_exposures_all_sigs_legacy.pdf")

    except Exception as e:
        print("No legacy sigs found")    
        print(e.__doc__)
        print(e.message)
        pass    

def run_sbs(sigfit = True, sigflow = True, deconstructSigs= True, mutationalPattern = True):
    try:
        sbs_df_list = []    
        sbs_df_name_list = []
        if sigflow:
            sigflow_sbs_exposure_df = pd.read_csv(r_output_file_dir + "/sigflow/SBS_fitting_relative_exposure.csv")
            pie_chart_rows = math.ceil(sigflow_sbs_exposure_df.shape[0] / 3) 
            pie_chart_cols = 3
            sigflow_sbs_exposure_df_fig = make_piecharts(sigflow_sbs_exposure_df, pie_chart_rows, pie_chart_cols, "COSMIC V3 'SBS' exposures data Sigflow")   
            sbs_df_list.append(sigflow_sbs_exposure_df.sort_values(by=["sample"]))
            sbs_df_name_list.append("sigflow")
        if mutationalPattern:
            mutationalPatterns_sbs_exposure_df = pd.read_csv(r_output_file_dir + "/mutational_patterns_results/sample_exposures.csv")
            mutationalPatterns_sbs_exposure_df = mutationalPatterns_sbs_exposure_df
            mutationalPatterns_sbs_exposure_df.rename(columns={"Unnamed: 0": "sample"}, inplace = True)
            mutationalPatterns_sbs_exposure_df.set_index("sample", inplace = True)
            mutationalPatterns_sbs_exposure_df = mutationalPatterns_sbs_exposure_df.rename(columns={"index": "sample"})
            mutationalPatterns_sbs_exposure_df.sum(axis=1)
            mutationalPatterns_sbs_exposure_df = mutationalPatterns_sbs_exposure_df.div(mutationalPatterns_sbs_exposure_df.sum(axis=1), axis=0)
            mutationalPatterns_sbs_exposure_df.reset_index(inplace= True)

            mutationalPatterns_strict_sbs_exposure_df = pd.read_csv(r_output_file_dir + "/mutational_patterns_results/strict_sample_exposures.csv")
            mutationalPatterns_strict_sbs_exposure_df = mutationalPatterns_strict_sbs_exposure_df
            mutationalPatterns_strict_sbs_exposure_df.rename(columns={"Unnamed: 0": "sample"}, inplace = True)
            mutationalPatterns_strict_sbs_exposure_df.set_index("sample", inplace = True)
            mutationalPatterns_strict_sbs_exposure_df = mutationalPatterns_strict_sbs_exposure_df.rename(columns={"index": "sample"})
            mutationalPatterns_strict_sbs_exposure_df.sum(axis=1)
            mutationalPatterns_strict_sbs_exposure_df.div(mutationalPatterns_strict_sbs_exposure_df.sum(axis=1), axis=0)
            mutationalPatterns_strict_sbs_exposure_df.reset_index(inplace= True)
            sbs_mutational_patters_fig = make_piecharts(mutationalPatterns_sbs_exposure_df, pie_chart_rows, pie_chart_cols, "COSMIC V3 'SBS' exposures  data MutationalPatterns")    
            sbs_mutational_patters_strict_fig = make_piecharts(mutationalPatterns_strict_sbs_exposure_df, pie_chart_rows, pie_chart_cols, "COSMIC V3 'SBS' exposures data MutationalPatterns Strict")    
            sbs_df_list.append(mutationalPatterns_strict_sbs_exposure_df.sort_values(by=["sample"]))
            sbs_df_name_list.append("mutationalPatterns_strict")
            sbs_df_list.append(mutationalPatterns_sbs_exposure_df.sort_values(by=["sample"]))
            sbs_df_name_list.append("mutationalPatterns")
        if deconstructSigs:
            deconstructSigs_sbs_exposure_df = pd.read_csv(r_output_file_dir + "/deconstructsigs_results/sbs_sample_exposures.csv")
            deconstructSigs_sbs_exposure_df.rename( columns = {"Unnamed: 0": "sample"}, inplace= True)
            deconstructSigs_sbs_exposure_df_fig = make_piecharts(deconstructSigs_sbs_exposure_df, pie_chart_rows, pie_chart_cols, "COSMIC V3 'Legacy' exposures data DeconstructSigs")    
            sbs_df_list.append(deconstructSigs_sbs_exposure_df.sort_values(by=["sample"]))
            sbs_df_name_list.append("deconstructSigs")            
        if sigfit:
            sigfit_sbs_exposure_df = pd.read_csv(r_output_file_dir + "/sigfit_results/sample_exposures_sbs.csv")
            sigfit_sbs_exposure_df.rename(columns={"Unnamed: 0": "sample" }, inplace= True)
            sigfit_sbs_exposure_df_fig = make_piecharts(sigfit_sbs_exposure_df, pie_chart_rows, pie_chart_cols, "COSMIC V3 'SBS' exposures data Sigfit")   
            sbs_df_list.append(sigfit_sbs_exposure_df.sort_values(by=["sample"]))
            sbs_df_name_list.append("sigfit")        
        with open(r_output_file_dir + '/sbs_pie_charts.html', 'a') as f:
            if siglow:
                f.write(sigflow_sbs_exposure_df_fig.to_html(full_html=False, include_plotlyjs='cdn'))
            if mutationalPattern:
                f.write(sbs_mutational_patters_fig.to_html(full_html=False, include_plotlyjs='cdn'))        
                f.write(sbs_mutational_patters_strict_fig.to_html(full_html=False, include_plotlyjs='cdn'))           
            if sigfit:
                f.write(sigfit_sbs_exposure_df_fig.to_html(full_html=False, include_plotlyjs='cdn'))
            if deconstructSigs:
                f.write(deconstructSigs_sbs_exposure_df_fig.to_html(full_html=False, include_plotlyjs='cdn'))
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
        rows = math.ceil(sigflow_sbs_exposure_df.shape[0] / 6) 
        columns = 6
        print(rows, columns)
        for i in range(len(distance_df_list)):  
            df = distance_df_list[i]
            df_ = pd.DataFrame(squareform(pdist(df.iloc[:, 1:])), columns=df["sample"].unique(), index=df["sample"].unique())
            sns.set()
            ax = fig.add_subplot(rows, columns, i+1)    
            sns.heatmap(df_)
            plt.title("Heatmap SBS " + str(sigfit_sbs_exposure_df.iloc[i,0]))
        plt.savefig(r_output_file_dir + "/Heatmap_SBS.pdf", dpi =500)

        fig = plt.figure(figsize=(120, 30))
        fig.subplots_adjust(right = 0.9, wspace=.01)
        rows = 1
        columns = 4
        for i in range(len(sbs_df_list)):  
            df = sbs_df_list[i].set_index("sample").astype(float)
            sns.set()
            ax = fig.add_subplot(rows, columns, i+1)
            sns.heatmap(df , cmap="YlGnBu", xticklabels = True, yticklabels= True )
            sns.set(rc={'figure.figsize':(30, 30)})
            plt.title("Heatmap SBS " + str(sbs_df_name_list[i]))
        plt.savefig(r_output_file_dir + "/Heatmap_exposures_all_sigs_SBS.pdf")

    except Exception as e:
        print("No SBS sigs found")    
        print(e.__doc__)
        print(e.message)
        pass

def run_id(sigfit = True, sigflow = True, deconstructSigs= True, mutationalPattern = True):
    try:
        id_df_list = []
        id_df_name_list = []
        if sigflow:
            sigflow_id_exposure_df = pd.read_csv(r_output_file_dir + "/sigflow/ID_fitting_relative_exposure.csv")
            pie_chart_rows = math.ceil(sigflow_id_exposure_df.shape[0] / 3) 
            pie_chart_cols = 3
            sigflow_id_exposure_df_fig = make_piecharts(sigflow_id_exposure_df, pie_chart_rows, pie_chart_cols, "COSMIC V3 'SBS' exposures data Sigflow")   
            id_df_list.append(sigflow_id_exposure_df.sort_values(by=["sample"]))
            id_df_name_list.append("sigflow")
        if mutationalPattern:
            mutationalPatterns_id_exposure_df = pd.read_csv(r_output_file_dir + "/mutational_patterns_results/id_sample_exposures.csv")
            mutationalPatterns_id_exposure_df = mutationalPatterns_id_exposure_df.rename(columns={"Unnamed: 0": "sample"})
            mutationalPatterns_id_exposure_df.set_index("sample", inplace = True)
            mutationalPatterns_id_exposure_df = mutationalPatterns_id_exposure_df.div(mutationalPatterns_id_exposure_df.sum(axis=1), axis=0)
            mutationalPatterns_id_exposure_df.reset_index(inplace = True)
            id_mutational_patters_fig = make_piecharts(mutationalPatterns_id_exposure_df, pie_chart_rows, pie_chart_cols, "COSMIC V3 'SBS' exposures  data MutationalPatterns")    
            id_df_list.append(mutationalPatterns_id_exposure_df.sort_values(by=["sample"]))
            id_df_name_list.append("mutationalPatters")
        if deconstructSigs:
            deconstructSigs_id_exposure_df = pd.read_csv(r_output_file_dir + "/deconstructsigs_results/indel_sample_exposures.csv")
            deconstructSigs_id_exposure_df.rename( columns = {"Unnamed: 0": "sample"}, inplace= True)
            id_df_list.append(deconstructSigs_id_exposure_df.sort_values(by=["sample"]) )
            id_df_name_list.append("deconstructSigs")
            deconstructSigs_id_exposure_df_fig = make_piecharts(deconstructSigs_id_exposure_df, pie_chart_rows, pie_chart_cols, "COSMIC V3 'Legacy' exposures data DeconstructSigs")    
        if sigfit:
            sigfit_id_exposure_df = pd.read_csv(r_output_file_dir + "/sigfit_results/sample_exposures_indel.csv")
            sigfit_id_exposure_df.rename(columns={"Unnamed: 0": "sample" }, inplace= True)
            sigfit_id_exposure_df_fig = make_piecharts(sigfit_id_exposure_df, pie_chart_rows, pie_chart_cols, "COSMIC V3 'SBS' exposures data Sigfit")   
            id_df_list.append(sigfit_id_exposure_df.sort_values(by=["sample"]))
            id_df_name_list.append("sigfit")        
        with open(r_output_file_dir + '/id_pie_charts.html', 'a') as f:
            if sigflow:
                f.write(sigflow_id_exposure_df_fig.to_html(full_html=False, include_plotlyjs='cdn'))
            if mutationalPattern:
                f.write(id_mutational_patters_fig.to_html(full_html=False, include_plotlyjs='cdn'))        
            if sigfit:
                f.write(sigfit_id_exposure_df_fig.to_html(full_html=False, include_plotlyjs='cdn'))
            if deconstructSigs:
                f.write(deconstructSigs_id_exposure_df_fig.to_html(full_html=False, include_plotlyjs='cdn'))
        distance_df_list = []
        for i in range(id_df_list[0].shape[0]):
            temp_df = []
            for index , j in enumerate(id_df_list):
                j["sample"] = id_df_name_list[index]
                temp_df.append(list(j.iloc[i,:]))
            _df = pd.DataFrame(temp_df, columns = id_df_list[3].columns)
            _df.fillna(0, inplace= True)
            distance_df_list.append(_df)

        fig = plt.figure(figsize=(30, 30))
        fig.subplots_adjust(hspace=2, wspace=2)
        rows = pie_chart_rows
        columns = 6
        for i in range(len(distance_df_list)):  
            df = distance_df_list[i]
            df_ = pd.DataFrame(squareform(pdist(df.iloc[:, 1:])), columns=df["sample"].unique(), index=df["sample"].unique())
            sns.set()
            ax = fig.add_subplot(rows, columns, i+1)    
            sns.heatmap(df_)
            plt.title("Heatmap ID " + str(sigfit_id_exposure_df.iloc[i,0]))
        plt.savefig(r_output_file_dir + "/Heatmap_ID.pdf", dpi =500)
        fig = plt.figure(figsize=(120, 30))
        fig.subplots_adjust(right = 0.9, wspace=.01)
        rows = 1
        columns = 4
        for i in range(len(id_df_list)):  
            df = id_df_list[i].set_index("sample").astype(float)
            sns.set()
            ax = fig.add_subplot(rows, columns, i+1)
            sns.heatmap(df , cmap="YlGnBu", xticklabels = True, yticklabels= True )
            sns.set(rc={'figure.figsize':(30, 30)})
            plt.title("Heatmap ID " + str(id_df_name_list[i]))
        plt.savefig(r_output_file_dir + "/Heatmap_exposures_all_sigs_ID.pdf")
        
    except Exception as e:
        print("No ID sigs found")    
        print(e.__doc__)
        print(e.message)
        pass
    
def run_dbs(sigfit = True, sigflow = True, deconstructSigs= True, mutationalPattern = True):
    try:
        dbs_df_list = []
        dbs_df_name_list = []
        if sigflow:
            sigflow_dbs_exposure_df = pd.read_csv(r_output_file_dir + "/sigflow/DBS_fitting_relative_exposure.csv")
            pie_chart_rows = math.ceil(sigflow_dbs_exposure_df.shape[0] / 3) 
            pie_chart_cols = 3
            sigflow_dbs_exposure_df_fig = make_piecharts(sigflow_dbs_exposure_df, pie_chart_rows, pie_chart_cols, "COSMIC V3 'SBS' exposures data Sigflow")   
            dbs_df_list.append(sigflow_dbs_exposure_df.sort_values(by=["sample"]))
            dbs_df_name_list.append("sigflow")
        if mutationalPattern:
            mutationalPatterns_dbs_exposure_df = pd.read_csv(r_output_file_dir + "/mutational_patterns_results/dbs_sample_exposures.csv")
            mutationalPatterns_dbs_exposure_df = mutationalPatterns_dbs_exposure_df.rename(columns={"Unnamed: 0": "sample"})
            mutationalPatterns_dbs_exposure_df.set_index("sample", inplace = True)
            mutationalPatterns_dbs_exposure_df = mutationalPatterns_dbs_exposure_df.div(mutationalPatterns_dbs_exposure_df.sum(axis=1), axis=0)
            mutationalPatterns_dbs_exposure_df.reset_index(inplace = True)
            dbs_mutational_patters_fig = make_piecharts(mutationalPatterns_dbs_exposure_df, pie_chart_rows, pie_chart_cols, "COSMIC V3 'SBS' exposures  data MutationalPatterns")    
            dbs_df_list.append(mutationalPatterns_dbs_exposure_df.sort_values(by=["sample"]))
            dbs_df_name_list.append("mutationalPatters")
        if deconstructSigs:
            deconstructSigs_dbs_exposure_df = pd.read_csv(r_output_file_dir + "/deconstructsigs_results/dbs_sample_exposures.csv")
            deconstructSigs_dbs_exposure_df.rename( columns = {"Unnamed: 0": "sample"}, inplace= True)
            deconstructSigs_dbs_exposure_df_fig = make_piecharts(deconstructSigs_dbs_exposure_df, pie_chart_rows, pie_chart_cols, "COSMIC V3 'Legacy' exposures data DeconstructSigs")    
            dbs_df_list.append(deconstructSigs_dbs_exposure_df.sort_values(by=["sample"]))
            dbs_df_name_list.append("deconstructSigs")            
        if sigfit:
            sigfit_dbs_exposure_df = pd.read_csv(r_output_file_dir + "/sigfit_results/sample_exposures_dbs.csv")
            sigfit_dbs_exposure_df.rename(columns={"Unnamed: 0": "sample" }, inplace= True)
            sigfit_dbs_exposure_df_fig = make_piecharts(sigfit_dbs_exposure_df, pie_chart_rows, pie_chart_cols, "COSMIC V3 'SBS' exposures data Sigfit")   
            dbs_df_list.append(sigfit_dbs_exposure_df.sort_values(by=["sample"]))
            dbs_df_name_list.append("sigfit")
        with open(r_output_file_dir + '/sbs_pie_charts.html', 'a') as f:
            if sigflow:
                f.write(sigflow_dbs_exposure_df_fig.to_html(full_html=False, include_plotlyjs='cdn'))
            if mutationalPattern:
                f.write(dbs_mutational_patters_fig.to_html(full_html=False, include_plotlyjs='cdn'))        
            if sigfit:
                f.write(sigfit_dbs_exposure_df_fig.to_html(full_html=False, include_plotlyjs='cdn'))
            if deconstructSigs:
                f.write(deconstructSigs_dbs_exposure_df_fig.to_html(full_html=False, include_plotlyjs='cdn'))
        distance_df_list = []
        for i in range(dbs_df_list[0].shape[0]):
            temp_df = []
            for index , j in enumerate(dbs_df_list):
                j["sample"] = dbs_df_name_list[index]
                temp_df.append(list(j.iloc[i,:]))
            _df = pd.DataFrame(temp_df, columns = dbs_df_list[3].columns)
            _df.fillna(0, inplace= True)
            distance_df_list.append(_df)

        fig = plt.figure(figsize=(30, 30))
        fig.subplots_adjust(hspace=2, wspace=2)
        rows = pie_chart_rows
        columns = 6
        for i in range(len(distance_df_list)):  
            df = distance_df_list[i]
            df_ = pd.DataFrame(squareform(pdist(df.iloc[:, 1:])), columns=df["sample"].unique(), index=df["sample"].unique())
            sns.set()
            ax = fig.add_subplot(rows, columns, i+1)    
            sns.heatmap(df_)
            plt.title("Heatmap DBS " + str(sigfit_dbs_exposure_df.iloc[i,0]))
        plt.savefig(r_output_file_dir + "/Heatmap_DBS.pdf", dpi =500)

        fig = plt.figure(figsize=(120, 30))
        fig.subplots_adjust(right = 0.9, wspace=.01)
        rows = 1
        columns = 4
        for i in range(len(dbs_df_list)):  
            df = dbs_df_list[i].set_index("sample").astype(float)
            sns.set()
            ax = fig.add_subplot(rows, columns, i+1)
            sns.heatmap(df , cmap="YlGnBu", xticklabels = True, yticklabels= True )
            sns.set(rc={'figure.figsize':(30, 30)})
            plt.title("Heatmap DBS " + str(dbs_df_name_list[i]))
        plt.savefig(r_output_file_dir + "/Heatmap_exposures_all_sigs_DBS.pdf")

    except Exception as e:
        print("No DBS sigs found")    
        print(e.__doc__)
        print(e.message)
        pass

run_legacy()
run_sbs()
run_id()
run_dbs()