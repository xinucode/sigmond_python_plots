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
    combine_fit_forms: false #(optional) default False. If false, the various fit forms will be graphed separately. If 
                                #true, the various fit forms will be graphed together on one plot
    dir: /latticeQCD/raid3/sarahski/lqcd/D200_R000/isotriplet_nonstrange_nucleonnucleon/.sigmond/plots/spectrum/tmin_plots/isotriplet_nonstrange_nucleonnucleon/rebin20 #directory where the xmgrace tmin plots are stored
    omit: #(optional) list of omissions from the batch, can be the fit type long name (settings.py) or the basis_pivot name
      - geometric fit
      - isotriplet_S0_A1g_P0_single_pivot_n8_m8_d16_c50
    bases: #list of the name of bases to combine. Must be a common substring of the intended tmin plots to combine
      - isotriplet_S0_A1g_P0
    select: #(optional) dict of pivot tags and a dict corresponding to their fit form+legend label combo
                #when select is present, only the fit forms and the spicified tags are plotted, otherwise all available are plotted
      _kN_single_pivot_n4_m4_d16_c150: 
        fit form: legend label 
        ...
        single-exponential ratio fit: single-exponential $N\overline{K}$ ratio fit
      _piS_single_pivot_n4_m4_d16_c150:
        single-exponential ratio fit: single-exponential $S\pi$ ratio fit
    fit_choices: #(optional but requires 'select' to work) dict of basis tags and a list of chosen fit forms that 
                    #correspond to each level
      isosinglet_Sm1_G1g_P0:
          - single-exponential $S\pi$ ratio fit
          - ...
      isosinglet_Sm1_G1u_P0:
          - single-exponential $S\pi$ ratio fit
          - ...
      ...
"""
  
def stability_tmin_plots():
    config_file = "isosinglet_strange.yml"
    if __name__ == "__main__":
        parser = argparse.ArgumentParser()
        parser.add_argument("config", help="config file")
        args = parser.parse_args()
        config_file = args.config

    with open(config_file, "r") as yamlfile:
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

            combine_fit_forms = False
            if 'combine_fit_forms' in channel:
                combine_fit_forms = channel['combine_fit_forms']

            selections = {} #if empty, all are selected
            if 'select' in channel:
                selections = channel['select']
            selected_bases = {}
            if selections:
                for basis in channel['bases']:
                    for tag in selections.keys():
                        selected_bases[basis+tag] = tag

            fit_choices = {}
            if 'fit_choices' in channel:
                fit_choices = channel['fit_choices']

            #reorganize plots by basis, level, fit, then pivot
            these_tmin_plots = utils.find_tmin_spectrum_files_python(channel)
            
            plots_by_bases = {}
            for basis in channel['bases']:
                plots_by_bases[basis] = {}
                for pivot in these_tmin_plots.keys():
                    if basis in pivot:
                        if selected_bases:
                            if pivot not in list(selected_bases.keys()):
                                continue
                        for level in these_tmin_plots[pivot].keys():
                            if level not in plots_by_bases[basis].keys():
                                plots_by_bases[basis][level] = {}
                            for fit in these_tmin_plots[pivot][level].keys():
                                if selections:
                                    if fit not in list(selections[selected_bases[pivot]].keys()):
                                        continue
                                if fit not in plots_by_bases[basis][level].keys():
                                    plots_by_bases[basis][level][fit] = {}
                                plots_by_bases[basis][level][fit][pivot] = these_tmin_plots[pivot][level][fit]

            #generate plots for each bases level and fit
            files_to_zip = []

            if __name__ == "__main__":
                f = plt.figure()
                f.set_figwidth(8)
                f.set_figheight(8)
            for basis in plots_by_bases:
                for level in plots_by_bases[basis]:
                    file_stub = f"{basis}_ROT{level}"
                    i=0
                    dd = 0
                    for fit in plots_by_bases[basis][level]:
                        if not combine_fit_forms:
                            file_stub = f"{basis}_ROT{level}_{settings.fit_nicknames[fit]}"
                            i=0
                        data = tmin_plots.retrieve_xmgrace_data_xydydy( plots_by_bases[basis][level][fit] )
                        fits = tmin_plots.retrieve_xmgrace_data_xy( plots_by_bases[basis][level][fit] )

                        for this_label in data.keys():
                            chosen_fit = pd.DataFrame()
                            print(f"\t{this_label}")
                            if selections:
                                legend_label = selections[selected_bases[this_label]][fit]
                                if fit_choices:
                                    if legend_label==fit_choices[basis][level]:
                                        this_fit = fits[this_label]
                                        plt.axhline(this_fit['fit'],color="black",zorder=1)
#                                         plt.axhline(this_fit['err'][0],color="black",ls="--")
#                                         plt.axhline(this_fit['err'][1],color="black",ls="--")
                                        plt.axhspan(this_fit['err'][0],this_fit['err'][1],color="darkgrey",zorder=1)
                                        chosen_fit = data[this_label].loc[ data[this_label][1]==this_fit['fit'] ]
                                legend_label = rf"{legend_label}"
                            else:
                                legend_label = this_label.replace(basis,"")
                                legend_label = legend_label.replace("_"," ")
                                if combine_fit_forms:
                                    legend_label += f" {settings.fit_nicknames[fit]}"

                            plt.errorbar(np.array(data[this_label][0])+dd,np.array(data[this_label][1]),np.concatenate([[np.array(data[this_label][3])],[np.array(data[this_label][2])]]),  capsize=5, color=settings.colors[i], marker=settings.markers[i], linewidth=0.0, elinewidth=1.5,label = legend_label,zorder=2,markerfacecolor="white")
                            if not chosen_fit.empty:
                                plt.errorbar(np.array(chosen_fit[0])+dd,np.array(chosen_fit[1]),np.concatenate([[np.array(chosen_fit[3])],[np.array(chosen_fit[2])]]),  capsize=5, color=settings.colors[i], marker=settings.markers[i], linewidth=0.0, elinewidth=1.5,zorder=3)
                                plt.ylim(np.array(chosen_fit[1])[0]-4.0*np.array(chosen_fit[3])[0],np.array(chosen_fit[1])[0]+6.0*np.array(chosen_fit[2])[0])
                            if combine_fit_forms:
                                dd+=0.1
                            i+=1
                            

                        if not combine_fit_forms:
                            print(file_stub)
                            plt.xlabel(r"$t_{\textup{min}}/a$")
                            plt.ylabel(r"$aE_{\textup{fit}}$")
                            if graph_type == 'dE':
                                plt.ylabel(r"$adE_{\textup{fit}}$")
                            plt.legend()
                            if __name__ == "__main__":
                                plt.tight_layout()
                                plt.savefig(f'{os.path.join(out_dir,file_stub)}_tmin.png')
                                files_to_zip.append(f'{file_stub}_tmin.png')
                                plt.savefig(f'{os.path.join(out_dir,file_stub)}_tmin.pdf')
                                files_to_zip.append(f'{file_stub}_tmin.pdf')
                                plt.clf()

                    if combine_fit_forms and plots_by_bases[basis][level]:
                        print(file_stub)
                        plt.xlabel(r"$t_{\textup{min}}/a$")
                        plt.ylabel(r"$aE_{\textup{fit}}$")
                        if graph_type == 'dE':
                            plt.ylabel(r"$adE_{\textup{fit}}$")
                        plt.legend() #bbox_to_anchor=(1.0, 1.5)
                        
                        #add ni level S(0)pi(0)
#                         val = 0.44916829252886237
#                         err = 0.00208263848513058
#                         plt.axhspan(val-err,val+err,color="lightgray",zorder=1)
#                         plt.axhline(val,color="darkgrey",zorder=1)
# #                         plt.axhline(val-err,color="lightgray",ls="--")
# #                         plt.axhline(val+err,color="lightgray",ls="--")
#                         plt.yticks([0.42,0.43,0.44,0.45])
#                         left, right = plt.xlim()
#                         plt.text( abs(left-right)*0.75+left -2.0,val+err+0.0005,r"$\Sigma\pi$ threshold",
#                                  color="black",zorder=3,fontsize="x-small")
                        
                        if __name__ == "__main__":
                            plt.tight_layout()
                            plt.savefig(f'{os.path.join(out_dir,file_stub)}_tmin.png')
                            files_to_zip.append(f'{file_stub}_tmin.png')
                            plt.savefig(f'{os.path.join(out_dir,file_stub)}_tmin.pdf')
                            files_to_zip.append(f'{file_stub}_tmin.pdf')
                            plt.clf()

            if __name__ == "__main__":
                utils.zip_channel( channel["name"], files_to_zip, "", out_dir)

if __name__ == "__main__":
    stability_tmin_plots()
