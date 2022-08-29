import numpy as np
import copy
import os, sys
import yaml
import pandas as pd
import matplotlib.pyplot as plt
import argparse

import tmin_plots
sys.path.insert(0,'../')
import xmgrace_parser
import settings
import utils
plt.style.use('../spectrum.mplstyle')

"""
Example input
--------------
channels: #list of the tmin plot batches to generate, ideally all that are associated with a specific scattering channel
  - name: isotriplet_nonstrange_nucleonnucleon #name of the channel, also name of the output file
    out_dir: isotriplet_nonstrange_nucleonnucleon_stability #subdirectory to place the batch of tmin plots
    max_level: 15 #max expected level number for any given basis
    graph_type: 'E' #(optional) default 'E'. can have values 'E' or 'dE'.  determines what graph type to be plotted: 
                        #the energy or the difference between the energy level and the non interacting level. 
    dir: /latticeQCD/raid3/sarahski/lqcd/D200_R000/isotriplet_nonstrange_nucleonnucleon/.sigmond/plots/spectrum/tmin_plots/isotriplet_nonstrange_nucleonnucleon/rebin20 #directory where the xmgrace tmin plots are stored
    omit: #(optional) list of omissions from the batch, can be the fit type long name (settings.py) or the basis_pivot name
      - geometric fit
      - isotriplet_S0_A1g_P0_single_pivot_n8_m8_d16_c50
    bases: #list of the name of bases to combine. Must be a common substring of the intended tmin plots to combine
      - isotriplet_S0_A1g_P0
      
      
"""
    

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument("config", help="config file")
    args = parser.parse_args()

    with open(args.config, "r") as yamlfile:
        project_info = yaml.load(yamlfile, Loader=yaml.FullLoader)
        
    for channel in project_info["channels"]:
        print("\n",channel["name"])
        
        #generate output directory
        out_dir = channel['out_dir']
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
            
        graph_type = 'E'
        if 'graph_type' in channel:
            graph_type = channel['graph_type']
        
        #reorganize plots by basis, level, fit, then pivot
        these_tmin_plots = utils.find_tmin_spectrum_files_python(channel)
        plots_by_bases = {}
        for basis in channel['bases']:
            plots_by_bases[basis] = {}
            for pivot in these_tmin_plots.keys():
                if basis in pivot:
                    for level in these_tmin_plots[pivot].keys():
                        if level not in plots_by_bases[basis].keys():
                            plots_by_bases[basis][level] = {}
                        for fit in these_tmin_plots[pivot][level].keys():
                            if fit not in plots_by_bases[basis][level].keys():
                                plots_by_bases[basis][level][fit] = {}
                            plots_by_bases[basis][level][fit][pivot] = these_tmin_plots[pivot][level][fit]
        
        #generate plots for each bases level and fit
        files_to_zip = []
        
        f = plt.figure()
        f.set_figwidth(8)
        f.set_figheight(8)
        for basis in plots_by_bases:
            for level in plots_by_bases[basis]:
                for fit in plots_by_bases[basis][level]:
                    file_stub = f"{basis}_ROT{level}_{settings.fit_nicknames[fit]}"
                    print(file_stub)
                    data = tmin_plots.retrieve_xmgrace_data_xydydy( plots_by_bases[basis][level][fit] )
                    
                    
                    for i,this_label in enumerate(data.keys()):
                        print(f"\t{this_label}")
                        legend_label = this_label.replace(basis,"")
                        legend_label = legend_label.replace("_"," ")
                        
                        plt.errorbar(np.array(data[this_label][0]),np.array(data[this_label][1]),np.concatenate([[np.array(data[this_label][3])],[np.array(data[this_label][2])]]),  capsize=5, color=settings.colors[i], marker=settings.markers[i], linewidth=0.0, elinewidth=1.5,label = legend_label)
                    
                    plt.xlabel(r"$t_{\textup{min}}/a$")
                    plt.ylabel(r"$aE_{\textup{fit}}$")
                    if graph_type == 'dE':
                        plt.ylabel(r"$adE_{\textup{fit}}$")
                    plt.legend()
                    plt.tight_layout()
                    plt.savefig(f'{os.path.join(out_dir,file_stub)}_tmin.png')
                    files_to_zip.append(f'{file_stub}_tmin.png')
                    plt.savefig(f'{os.path.join(out_dir,file_stub)}_tmin.pdf')
                    files_to_zip.append(f'{file_stub}_tmin.pdf')
                    plt.clf()
                    
        utils.zip_channel( channel["name"], files_to_zip, "", out_dir)