import numpy as np
import h5py
import os, sys
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import image as mpimg
import yaml
import argparse
import itertools

sys.path.append('../')
import utils
import settings
plt.style.use('../spectrum.mplstyle')


"""
Sample Config Yaml file
----------------------
channel: isoquartet_nonstrange_fermionic #name of channel and corresponding 
                                            #subdirectory where relevant files are located
scattering_particles: [N, pi] #list of scattering particle names
rest_mass: pi #the particle used to normalize the energy values
title: $I=\sfrac{3}{2}$ #(optional) if unspecified, plot will not print a title in the legend.
    #also no title will print if no legend
thresholds: #(optional) dict of thersholds with the key being how you want the threshold 
                #to appear on the graph (replaced with latex) and the value is a list of  
                #the rest masses to be summed over for the threshold value
    "Npi": [ N, pi ]
    "Npipi": [ N, pi, pi]
used_levels: #(optional) dict of bases and lists of their levels to be filled markers in
                #leaving the rest of the values empty. If not present, then all values will
                #be filled markers
    G1g(0): [0]
    G1u(0): [0,1]
    Hg(0): [0]
    Hu(0): [0]
    G1(1): [0,1,2]
    G2(1): [0]
    G(2): [0,1,2,3,4]
    F1(3): [0]
    F2(3): [0,1]
    G(3): [0,1,2,3]
    G1(4): [0,1,2]
    G2(4): [0]
  
fig_width: 14 #(optional) sets width of all figures
fig_height: 6 #(optional) sets height of all figures

compare_spectrums: #(optional) generates graph that compares several spectrums with similar bases side by side
    spectrum_type: energy #value can be "energy" or "mom" determines how the files are parsed and how the graph 
                            #is labelled
    file_directory: final #subdirectory inside "channel" where the data can be found and the graph will be produced
    files: #dict of files with the legend labels as the keys and filenames as the values
        Pivot(n-1,8,16): energy_estimates_isoquartet_nonstrange_fermionic_colin_trimmed_rebin20_Bootstrap.csv
        Pivot(n,8,16): energy_estimates_isoquartet_nonstrange_fermionic_colin_rebin20_Bootstrap_8-16.csv
        Pivot(n,6,12): energy_estimates_isoquartet_nonstrange_fermionic_colin_rebin20_Bootstrap_6-12.csv
        Pivot(n,4,8): energy_estimates_isoquartet_nonstrange_fermionic_colin_rebin20_Bootstrap_4-8.csv
    best_legend_loc: upper left #(optional) put legend where desired. If omitted, no legend. If legend is
                                #desired but unsure where, put "best"
    do_scattering_particles: false #(optional) generate another graph that compares the elastic levels. False if
                                    #omitted
    fig_width: 14 #(optional) sets width of this figure, supercedes the setting for all figures
    fig_height: 6 #(optional) sets height of this figure, supercedes the setting for all figures
    zshift: ['Hg(0)','G1(1)'] #(optional) list of irreps(d^2) to zigzag (2 in a row)
                                #shift the levels so they don't overlap
    sshift: ['G2(1)','G(2)'] #(optional) list of irreps(d^2) to stripe (3 in a row)
                               #shift the levels so they don't overlap #takes precidence over shifts above
    s4shift: ['G(3)','F1(3)'] #(optional) list of irreps(d^2) to stripe (4 in a row)
                               #shift the levels so they don't overlap #takes precidence over shifts above
    s5shift: ['G1(4)','G2(4)'] #(optional) list of irreps(d^2) to stripe (5 in a row)
                               #shift the levels so they don't overlap #takes precidence over shifts above
    omit: [G1g(0),Hu(0)] #(optional) list of irreps(d^2) to omit; otherwise, it graphs all irreps
    xrange: [-0.3,0.705] #(optional) sets the xrange and prints out old versus new. Must be a list of two
    yrange: [0.0,2.5] #(optional) manually select the yrange, otherwise matplotlib automatically sets it
    plot_ni_levels: true #(optional) default false; plots the non interacting levels for each spectrum
    ni_width: 40 #(optional) default 80; sets the width of the non-interacting levels (if plot_ni_levels==true)
    remove_xlabel: #(optional) default false; boolean to decide if print xlabel
    graph_unused_levels: false #(optional) default true; if unused_levels is specified and this tag is set to false, 
        #unused levels will not be graphed.
      
final_spectrum: #(optional) generates graph that plots just one spectrum
    spectrum_type: energy #value can be "energy" or "mom" determines how the files are parsed and how the graph 
                            #is labelled
    file_directory: final #subdirectory inside "channel" where the data can be found and the graph will be produced
    file: energy_estimates_isoquartet_nonstrange_fermionic_colin_rebin20_Bootstrap_8-16.csv #file with data
    zshift: ['Hg(0)','G1(1)'] #(optional) list of irreps(d^2) to zigzag (2 in a row)
                                #shift the levels so they don't overlap
    sshift: ['G2(1)','G(2)'] #(optional) list of irreps(d^2) to stripe (3 in a row)
                               #shift the levels so they don't overlap #takes precidence over shifts above
    s4shift: ['G(3)','F1(3)'] #(optional) list of irreps(d^2) to stripe (4 in a row)
                               #shift the levels so they don't overlap #takes precidence over shifts above
    s5shift: ['G1(4)','G2(4)'] #(optional) list of irreps(d^2) to stripe (5 in a row)
                               #shift the levels so they don't overlap #takes precidence over shifts above
    ni_width: 40 #(optional) default 80; sets the width of the non-interacting levels
    best_legend_loc: lower left #(optional) put legend where desired. If omitted, no legend. If legend is
                                #desired but unsure where, put "best"
    fig_width: 14 #(optional) sets width of this figure, supercedes the setting for all figures
    fig_height: 6 #(optional) sets height of this figure, supercedes the setting for all figures
    omit: [G1g(0),Hu(0)] #(optional) list of irreps(d^2) to omit; otherwise, it graphs all irreps
    remove_xlabel: true #(optional) default false; boolean to decide if print xlabel
    graph_unused_levels: false #(optional) default true; if unused_levels is specified and this tag is set to false, 
        #unused levels will not be graphed.
    xrange: [-0.3,0.705] #(optional) sets the xrange and prints out old versus new. Must be a list of two
    yrange: [0.0,2.5] #(optional) manually select the yrange, otherwise matplotlib automatically sets it
"""

def compare_spectrums():
    remove_ref = False
    do_scattering_particles = False
    file1 = ""
    file2 = ""
    file3 = ""
    file4 = ""
    file5 = ""
    figlabel1 = ""
    figlabel2 = ""
    figlabel3 = ""
    figlabel4 = ""
    figlabel5 = ""

    #import info from config
    config_file = "isosinglet_strange.yml"
    if __name__=="__main__":
        parser = argparse.ArgumentParser()
        parser.add_argument("config", help="config file")
        args = parser.parse_args()
        config_file = args.config


    with open(config_file, "r") as yamlfile:
        configdata = yaml.load(yamlfile, Loader=yaml.FullLoader)

    particle_names = configdata['scattering_particles']
    # print()

    scat_momentum_tags = ['(0)_ref','(1)_ref','(2)_ref','(3)_ref','(4)_ref']
    scattering_particles = []
    for particle in particle_names:
        for tag in scat_momentum_tags:
            scattering_particles.append(particle+tag)

    max_mom = 5
    elastic_scat_keys = []
    if len(particle_names)==1:
        for i in range(max_mom):
            for j in range(max_mom):
                elastic_scat_keys.append( f"{particle_names[0]}({i}){particle_names[0]}({j})_ref" )
    elif len(particle_names)==2:
        for i in range(max_mom):
            for j in range(max_mom):
                elastic_scat_keys.append( f"{particle_names[0]}({i}){particle_names[1]}({j})_ref" )
        for i in range(max_mom):
            for j in range(max_mom):
                elastic_scat_keys.append( f"{particle_names[1]}({i}){particle_names[0]}({j})_ref" )
    else:
        for combo in itertools.product(particle_names,particle_names):
            for i in range(max_mom):
                for j in range(max_mom):
                    elastic_scat_keys.append( f"{combo[0]}({i}){combo[1]}({j})_ref" )

    channel = configdata['channel']
    rest_mass = configdata['rest_mass']

    graphs = ['compare_spectrums','final_spectrum']
    ni_width = 80.0

    for graph in graphs:
        omit = []
        files = {}
        if graph not in configdata.keys():
            continue

        if "omit" in configdata[graph].keys():
            omit = configdata[graph]["omit"]

        if ('compare_spectrums' in configdata.keys()) and (graph=='compare_spectrums'):
            file_directory = configdata['compare_spectrums']['file_directory']
            for key in configdata['compare_spectrums']['files'].keys():
                if configdata['compare_spectrums']['file_directory']:
                    files[key] = os.path.join(channel,file_directory,configdata['compare_spectrums']['files'][key])
            spectrum_type = configdata['compare_spectrums']['spectrum_type']
            if 'do_scattering_particles' in configdata['compare_spectrums'].keys():
                do_scattering_particles = bool(configdata['compare_spectrums']['do_scattering_particles'])

            if 'best_legend_loc' in configdata['compare_spectrums'].keys():
                best_legend_loc=configdata['compare_spectrums']['best_legend_loc']
            else:
                best_legend_loc = ""

            if spectrum_type=="mom":
                do_scattering_particles = False

            file_name = channel+f"-{spectrum_type}_spectrum_comparison_graph.pdf"
            if 'plot_ni_levels' in configdata['compare_spectrums'].keys():
                plot_ni_levels = configdata['compare_spectrums']['plot_ni_levels']
                if 'ni_width' in configdata['compare_spectrums'].keys():
                    ni_width = configdata['compare_spectrums']['ni_width']
            else:
                plot_ni_levels = False


        if ('final_spectrum' in configdata.keys()) and (graph=='final_spectrum'):
            file_directory = configdata['final_spectrum']['file_directory']
            best_legend_loc=configdata['final_spectrum']['best_legend_loc']
            files["this"] = os.path.join(channel,file_directory,configdata['final_spectrum']['file'])
            spectrum_type = configdata['final_spectrum']['spectrum_type']
            file_name = channel+"_spectrum_graph.pdf"
            if 'ni_width' in configdata['final_spectrum'].keys():
                ni_width = configdata['final_spectrum']['ni_width']
    #         non_interacting_levels = configdata['final_spectrum']['non_interacting_levels']

        if 'used_levels' in configdata.keys():
            used_levels=configdata['used_levels']
        else:
            used_levels = {}

        remove_xlabel = False
        if "remove_xlabel" in configdata[graph].keys():
            remove_xlabel = configdata[graph]["remove_xlabel"]

        graph_unused_levels = True
        if "graph_unused_levels" in  configdata[graph].keys():
            graph_unused_levels = configdata[graph]["graph_unused_levels"]
            
        #########################################################
        ############### import info from files ##################
        #########################################################
        
        Refs = {}
        datasets = {}

        for file in files.keys():
            datasets[file] = utils.unpack_file(files[file],spectrum_type)
            if datasets[file] is None:
                sys.exit()

        for scat in particle_names:
            value = 0.0
            for file in files.keys():
                value = utils.find_rest_mass( datasets[file], scat )
                if value:
                    Refs[scat] = value
                    break

        expected_keys = [f"PSQ{i}" for i in range(4)]
        possible_irreps = ['G1u', 'Hg', 'G1', 'G2', 'F1', 'F2', 'G','A1g','A1gm','G1g','Hu', 'T1g','A1u','A2u','Eg','Eu',
                           'T1u','T2g','T2u','A1','A2','B1','B2','E']
#         diff_keys = [f"d_ecm_{i}" for i in range(15)] #probably change this to elab
        if spectrum_type=="mom":
            energy_keys = [f'q2cm_{i}_ref' for i in range(15)]
        elif spectrum_type=="energy":
            if remove_ref:
                energy_keys = [f'ecm_{i}' for i in range(15)]
            else:
                energy_keys = [f'ecm_{i}_ref' for i in range(15)]
        else:
            print("Bad spectrum type")

        keys = {}
        vals = {}
        errs = {}
        levs = {}
        vals_used = {}
        keys_used = {}
        errs_used = {}
        levs_used = {}
        ekeys = {}
        evals = {}
        eerrs = {}
        max_level_num = 0
        for dataset in datasets.keys():
            vals[dataset] = []
            keys[dataset] = []
            errs[dataset] = []
            levs[dataset] = []
            vals_used[dataset] = []
            keys_used[dataset] = []
            errs_used[dataset] = []
            levs_used[dataset] = []
            evals[dataset] = []
            ekeys[dataset] = []
            eerrs[dataset] = []
            for i, (mom) in enumerate(expected_keys):
                for j, (irrep) in enumerate(possible_irreps):
                    for k, (energy) in enumerate(energy_keys):

                        if not datasets[dataset].empty:
                            val1, err1 = utils.select_val(datasets[dataset], mom, irrep, energy)
                        else:
                            val1, err1 = utils.select_val_ascii(files[dataset],mom,irrep,energy,spectrum_type)

                        if val1 is not None and err1 is not None:
                            irrep_key = f"{irrep}({mom.replace('PSQ','')})"
                            if irrep_key in omit:
                                continue
                            if k>max_level_num:
                                max_level_num=k
                            if "yrange" in configdata[graph].keys():
                                if (val1 < configdata[graph]["yrange"][0]) or (val1 > configdata[graph]["yrange"][1]):
                                    continue
                            used = False
                            if used_levels:
                                if used_levels[irrep_key]:
                                    if k in used_levels[irrep_key]:
                                        vals_used[dataset].append(val1)
                                        keys_used[dataset].append(irrep_key)
                                        errs_used[dataset].append(err1)
                                        levs_used[dataset].append(k)
                                        used = True
    #                                     print(irrep_key,k,val1)
                            if not used:
                                keys[dataset].append(irrep_key)
                                vals[dataset].append(val1)
                                errs[dataset].append(err1)
                                levs[dataset].append(k)


                    for k, (energy) in enumerate(elastic_scat_keys):

                        if not datasets[dataset].empty:
                            val1, err1 = utils.select_val(datasets[dataset], mom, irrep, energy)
                        else:
                            val1, err1 = utils.select_val_ascii(files[dataset],mom,irrep,energy,spectrum_type)

                        if val1 is not None and err1 is not None:
                            irrep_key = f"{irrep}({mom.replace('PSQ','')})"
                            if irrep_key in omit:
                                continue
                            ekeys[dataset].append(irrep_key)
                            evals[dataset].append(val1)
                            eerrs[dataset].append(err1)

    #     print(ekeys,evals,eerrs)
        scat_keys2 = []
        scat_vals = {}
        scat_errs = {}

        if do_scattering_particles:
            for dataset in datasets.keys():
                scat_vals[dataset] = []
                scat_errs[dataset] = []

            for level in scattering_particles:
                this_value = {}
                this_error = {}
                good_level = True
                for dataset in datasets.keys():
                    val1, err1 = utils.select_val(datasets[dataset], None, None, level)

                    if (val1 is not None and err1 is not None):
                        this_value[dataset]=val1
                        this_error[dataset]=err1
                    else:
                        this_value[dataset]=None
                        this_error[dataset]=None
                        good_level = False

                if good_level:
                    for name in particle_names:
                        if (name+'(') in level:
                            scat_keys2.append(name)
                    for dataset in datasets.keys():
                        scat_vals[dataset].append(this_value[dataset])
                        scat_errs[dataset].append(this_error[dataset])

        #organize indexes
        if do_scattering_particles:
            scat_indexes2 = np.zeros(len(scat_keys2))
            for j, key in enumerate(scat_keys2):
                for i, unique_key in enumerate(list(set(scat_keys2))):
                    if key == unique_key:
                        scat_indexes2[j] = i
                        
        
        indexes = {}
        this_list = []
        for dataset in datasets.keys():
            indexes[dataset] = np.zeros(len(keys[dataset]))
            for j, key in enumerate(keys[dataset]):
                this_list = list(set(this_list+keys[dataset]+keys_used[dataset]))
                this_list.sort(key=utils.sort_by_mom)
                for i, unique_key in enumerate(this_list):
                    if key == unique_key:
                        indexes[dataset][j] = i

        indexes_used = {}
        for dataset in datasets.keys():
            indexes_used[dataset] = np.zeros(len(keys_used[dataset]))
            for j, ukey in enumerate(keys_used[dataset]):
                this_list = list(set(this_list+keys[dataset]+keys_used[dataset]))
                this_list.sort(key=utils.sort_by_mom)
                for i, unique_key in enumerate(this_list):
                    if ukey == unique_key:
                        indexes_used[dataset][j] = i


        eindexes = {}
        for dataset in datasets.keys():
            eindexes[dataset] = np.zeros(len(ekeys[dataset]))
            for j, ukey in enumerate(ekeys[dataset]):
                this_list = list(set(keys[dataset]+keys_used[dataset]))
                this_list.sort(key=utils.sort_by_mom)
                for i, unique_key in enumerate(this_list):
                    if ukey == unique_key:
                        eindexes[dataset][j] = i

        
        if ('compare_spectrums' in configdata.keys()) and (graph=='compare_spectrums'):
            #plot scattering particles
            if do_scattering_particles:
                f2 = plt.figure()
                f2.set_figwidth(6)
                f2.set_figheight(12)
                dd = 0.1

                for i,dataset in enumerate(scat_vals.keys()):
                    plt.scatter(scat_indexes2+dd*i,scat_vals[dataset],color=settings.colors[i], marker=settings.markers[i], label = dataset)
                    plt.errorbar(scat_indexes2+dd*i, scat_vals[dataset], scat_errs[dataset], fmt='.', capsize=5,color=settings.colors[i])
                plt.ylabel(r"$E_{cm}/m_\pi$")
                plt.xlabel("Noninteracting Scattering Levels")
                plt.legend()
                plt.xticks(scat_indexes2, scat_keys2)
#                 plt.title(channel)
                plt.savefig(os.path.join(channel,file_directory,channel+"-ni_comparison_graph.pdf"))
                plt.clf()
                        
        #plot spectrum
        if __name__=="__main__":
            f = plt.figure()
            if 'fig_width' in configdata.keys():
                f.set_figwidth(configdata['fig_width'])
            if 'fig_height' in configdata.keys():
                f.set_figheight(configdata['fig_height'])
            if 'fig_width' in configdata[graph].keys():
                f.set_figwidth(configdata[graph]['fig_width'])
            if 'fig_height' in configdata[graph].keys():
                f.set_figheight(configdata[graph]['fig_height'])
            
        if "yrange" in configdata[graph].keys():
            plt.ylim( configdata[graph]["yrange"][0],configdata[graph]["yrange"][1])

        somekey = list(indexes.keys())[0]
        minx = min(list(indexes[somekey])+list(indexes_used[somekey]))
        maxx = max(list(indexes[somekey])+list(indexes_used[somekey]))
        dd=0.4/len(files.keys())
        plt.xlim(minx-0.5-dd*(len(vals.keys())/2),maxx+0.5+dd*(len(vals.keys())/2))
        if not remove_ref and spectrum_type=="energy":
            if configdata['thresholds']:
                plt.xlim(minx-0.5-dd*(len(vals.keys())/2),(maxx+dd*(len(vals.keys())/2))*1.025+0.5)
                


        if not remove_ref and spectrum_type=="energy":
            if configdata['thresholds']:
                for threshold in configdata['thresholds'].keys():
                    threshold_value = 0.0
                    threshold_label = threshold
                    for particle in configdata['thresholds'][threshold]:
                        threshold_value+=Refs[particle][0]
                        threshold_label=threshold_label.replace(settings.latex_format[particle], particle)
                        threshold_label=threshold_label.replace(particle, settings.latex_format[particle])
                    miny,maxy = plt.ylim()
                    plt.hlines(threshold_value,minx-dd*(len(vals.keys())/2),maxx+dd*(len(vals.keys())/2),color='black', linestyle="--", zorder=1) 
                    plt.text( (maxx+dd*(len(vals.keys())/2))*(1.01),threshold_value-0.015*abs(maxy-miny), threshold_label, zorder=6,size="x-small")
        if spectrum_type=="mom":
            if "yrange" in configdata[graph].keys():
                if (0.0>configdata[graph]["yrange"][0]) and (0.0<configdata[graph]["yrange"][1]):
                    plt.hlines(0.0,minx-dd*(len(vals.keys())/2),maxx+dd*(len(vals.keys())/2),color='black', linestyle="--", zorder=1, linewidth=1.0)
            else:
                plt.hlines(0.0,minx-dd*(len(vals.keys())/2),maxx+dd*(len(vals.keys())/2),color='black', linestyle="--", linewidth=1.0, zorder=1)

        for i,dataset in enumerate(vals.keys()):
            print(dataset)
            shifted_array = 0.0 
            used_shifted_array = 0.0 
#         if ('compare_spectrums' in configdata.keys()) and (graph=='compare_spectrums'):
#             if ('zshift' in configdata[graph].keys() or 'sshift' in configdata[graph].keys() or 's4shift' in configdata[graph].keys() or 's5shift' in configdata[graph].keys()): # and dataset!="fit"
#                 shifted_array = []
#                 for j,lev in enumerate(levs[dataset]):
#                     this_shift = 0.0
#                     if 'zshift' in configdata[graph].keys():
#                         if keys[dataset][j] in configdata[graph]['zshift']:
#                             this_shift = 0.1*dd*settings.zigzag_shifts[lev]
#                     if 'sshift' in configdata[graph].keys():
#                         if keys[dataset][j] in configdata[graph]['sshift']:
#                             this_shift = 0.2*dd*settings.stripe_shifts[lev]
#                     if 's4shift' in configdata[graph].keys():
#                         if keys[dataset][j] in configdata[graph]['s4shift']:
#                             this_shift = 0.3*dd*settings.stripe4_shifts[lev]
#                     if 's5shift' in configdata[graph].keys():
#                         if keys[dataset][j] in configdata[graph]['s5shift']:
#                             this_shift = 0.2*dd*settings.stripe5_shifts[lev]
#                     shifted_array.append(this_shift)

#                 used_shifted_array = []
#                 for j,lev in enumerate(levs_used[dataset]):
#                     this_used_shift = 0.0
#                     if 'zshift' in configdata[graph].keys():
#                         if keys_used[dataset][j] in configdata[graph]['zshift']:
#                             this_used_shift = 0.1*dd*settings.zigzag_shifts[lev]
#                     if 'sshift' in configdata[graph].keys():
#                         if keys_used[dataset][j] in configdata[graph]['sshift']:
#                             this_used_shift = 0.2*dd*settings.stripe_shifts[lev]
#                     if 's4shift' in configdata[graph].keys():
#                         if keys_used[dataset][j] in configdata[graph]['s4shift']:
#                             this_used_shift = 0.3*dd*settings.stripe4_shifts[lev]
#                     if 's5shift' in configdata[graph].keys():
#                         if keys_used[dataset][j] in configdata[graph]['s5shift']:
#                             this_used_shift = 0.2*dd*settings.stripe5_shifts[lev]
#                     used_shifted_array.append(this_used_shift)

#                 shifted_array = np.array(shifted_array)
#                 used_shifted_array = np.array(used_shifted_array)

            # utils.shift_levels(np.array(indexes[dataset]),levs[dataset],vals[dataset],errs[dataset])
            used_shifted_array = 0.2*dd*utils.shift_levels(np.array(indexes_used[dataset]),levs_used[dataset],vals_used[dataset],errs_used[dataset])
            

            if used_levels:
                marker_color = 'white'
            else:
                marker_color = settings.colors[i]

            if used_levels:
                unused_label = None
            else:
                unused_label = dataset

            splitting_factor = (i-((len(vals.keys())-1)/2))

            if (not used_levels) or (used_levels and graph_unused_levels):
                if len(np.nonzero(errs[dataset])[0]):
                    plt.errorbar(indexes[dataset]+dd*splitting_factor+shifted_array, vals[dataset], np.array(errs[dataset]),  capsize=5, color=settings.colors[i], marker=settings.markers[i], linewidth=0.0, elinewidth=1.5,mfc=marker_color,zorder=4,label=unused_label)
                else:
                    plt.scatter(indexes[dataset]+dd*splitting_factor+shifted_array, vals[dataset], color=settings.colors[i], marker=settings.markers[i],linewidth=0.0, zorder=4,label=unused_label)

            if used_levels:
                if len(np.nonzero(errs_used[dataset])[0]):
                    plt.errorbar(indexes_used[dataset]+dd*splitting_factor+used_shifted_array,vals_used[dataset], np.array(errs_used[dataset]),  capsize=5, color=settings.colors[i], marker=settings.markers[i], linewidth=0.0, elinewidth=1.5, mfc=settings.colors[i], zorder=4, label=dataset)
                else:
                    plt.scatter(indexes_used[dataset]+dd*splitting_factor+used_shifted_array,vals_used[dataset], color=settings.colors[i], marker=settings.markers[i], linewidth=0.0, zorder=4,label=dataset)

        if (graph=='final_spectrum') or plot_ni_levels:
            for i,dataset in enumerate(evals.keys()):
                plt.errorbar(eindexes[dataset]+dd*(i-((len(evals.keys())-1)/2)), evals[dataset], np.array(eerrs[dataset]), color="lightgray", marker="_",linestyle="", linewidth=0.0, elinewidth=ni_width,zorder=2)
                plt.scatter(eindexes[dataset]+dd*(i-((len(evals.keys())-1)/2)),evals[dataset],color="darkgrey", marker="_",s=ni_width*ni_width,zorder=3)


        latex_rest_mass = settings.latex_format[rest_mass].replace('$',"")
        if spectrum_type=='energy':
            plt.ylabel(rf"$E_{{\textup{{cm}}}}/m_{{{latex_rest_mass}}}$")
        else:
            plt.ylabel(rf"$q^2_{{\textup{{cm}}}}/m_{{{latex_rest_mass}}}^2$")

        if not remove_xlabel:
            plt.xlabel(r"$\Lambda(d^2)$")
        latex_keys = [settings.latex_format[key.split('(')[0]]+"("+key.split('(')[1] for key in keys[somekey]+keys_used[somekey]]
        plt.xticks(utils.unique(list(indexes[somekey])+list(indexes_used[somekey])), utils.unique(latex_keys),size="small")
#         plt.xticks([])
        if 'xrange' in configdata[graph].keys():
            print("Old xrange:", plt.xlim() )
            plt.xlim( configdata[graph]['xrange'] )
            print("New xrange:", plt.xlim() )


        if (graph=='compare_spectrums'):
            if best_legend_loc and "title" in configdata.keys():
                plt.legend(title=rf"{configdata['title']}", loc=best_legend_loc)
            elif best_legend_loc:
                plt.legend(loc=best_legend_loc)

            plt.tight_layout()
        if __name__=="__main__":
            plt.savefig(os.path.join(channel,file_directory,file_name))
            plt.clf()

if __name__=="__main__":
    compare_spectrums()