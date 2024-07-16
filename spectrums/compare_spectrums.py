import numpy as np
import h5py
import os, sys
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import image as mpimg
from matplotlib import colors as mc
import yaml
import argparse
import itertools
from timeit import default_timer as timer 

sys.path.append('../')
import utils
import cmd_plot
import settings
plt.style.use('../spectrum.mplstyle')

def load_latex_config():
    # pix_to_pt = 72/100
    # plt.rcParams['errorbar.capsize'] = 1
    # plt.rcParams['lines.linewidth'] = 2.0 * pix_to_pt
    # plt.rcParams['lines.markeredgewidth'] = 2.0 * pix_to_pt
    # plt.rcParams['savefig.bbox'] = 'tight'
    # plt.rcParams['axes.labelsize'] = 21
    # plt.rcParams['axes.labelpad'] = 15
    # plt.rcParams['axes.linewidth'] = 1.0 * pix_to_pt
    # plt.rcParams['axes.titlesize'] = 11
    # plt.rcParams['lines.markersize'] = 2
    # plt.rcParams['xtick.labelsize'] = 21
    # plt.rcParams['ytick.labelsize'] = 21
    # plt.rcParams['xtick.major.width'] = 0.5 * pix_to_pt
    # plt.rcParams['xtick.minor.width'] = 0.5 * pix_to_pt
    # plt.rcParams['ytick.major.width'] = 0.5 * pix_to_pt
    # plt.rcParams['ytick.minor.width'] = 0.5 * pix_to_pt
    # plt.rcParams['text.usetex'] = True
    # plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath,amssymb}'
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.serif'] = ['Computer Modern Roman']
    # plt.rcParams['font.size'] = 18
    plt.rcParams['legend.fancybox'] = False
    plt.rcParams['legend.shadow'] = False
    plt.rcParams['legend.framealpha'] = '0.8'
    plt.rcParams['legend.facecolor'] = 'white'
    plt.rcParams['legend.edgecolor'] = '0.5'
    # plt.rcParams['patch.linewidth'] = 0.5 * pix_to_pt
    # plt.rcParams['legend.handlelength'] = 1.0
    # plt.rcParams['legend.fontsize'] = 21
    # plt.rcParams['figure.figsize'] = (8, 6)
    
    
load_latex_config()

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

def set_figure(f, config1, config2 ):
    if 'fig_width' in config1.keys():
        f.set_figwidth(config1['fig_width'])
    if 'fig_height' in config1.keys():
        f.set_figheight(config1['fig_height'])
    if 'fig_width' in config2.keys():
        f.set_figwidth(config2['fig_width'])
    if 'fig_height' in config2.keys():
        f.set_figheight(config2['fig_height'])
        
def finalize_plot(remove_xlabel, config, graph, best_legend_loc, rest_mass, 
                    energy_key, spectrum_type, remove_ref, preliminary):
    if rest_mass:
        latex_rest_mass = settings.latex_format[rest_mass].replace('$',"")
    
    energy_tag = energy_key[1:]
    if spectrum_type=='energy':
        if energy_key=='dElab':
            if remove_ref:
                plt.ylabel(rf"$\Delta E_{{\textup{{lab}}}}$")
            else:
                plt.ylabel(rf"$\Delta E_{{\textup{{lab}}}}/m_{{{latex_rest_mass}}}$")
        else:
            if remove_ref:
                plt.ylabel(rf"$E_{{\textup{{{energy_tag}}}}}$")
            else:
                plt.ylabel(rf"$E_{{\textup{{{energy_tag}}}}}/m_{{{latex_rest_mass}}}$")
            
    else:
        plt.ylabel(rf"$q^2_{{\textup{{{energy_tag}}}}}/m_{{{latex_rest_mass}}}^2$")
            
    if not remove_xlabel:
        plt.xlabel(r"$\Lambda(d^2)$")
    
    if 'xrange' in config[graph].keys():
        print("Old xrange:", plt.xlim() )
        plt.xlim( config[graph]['xrange'] )
        print("New xrange:", plt.xlim() )
    
    if preliminary:
        plt.figtext( 0.4, 0.4, "PRELIMINARY", alpha = 0.5, color="red", fontweight="bold", rotation= 30,zorder=0)


    if (graph=='compare_spectrums'):
        if best_legend_loc and "title" in config.keys():
            plt.legend(title=rf"{config['title']}", loc=best_legend_loc).set_zorder(7)
        elif best_legend_loc:
            plt.legend(loc=best_legend_loc).set_zorder(7)

    plt.tight_layout()


def compare_spectrums(config_file = "isosinglet_strange.yml"):
    start = timer()
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
    with open(config_file, "r") as yamlfile:
        configdata = yaml.load(yamlfile, Loader=yaml.FullLoader)

    particle_names = configdata.pop('scattering_particles', [] )
    # print()

    scat_momentum_tags = ['(0)','(1)','(2)','(3)','(4)']
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
    else:
        for combo in itertools.product(particle_names,particle_names):
            for i in range(max_mom):
                for j in range(max_mom):
                    elastic_scat_keys.append( f"{combo[0]}({i}){combo[1]}({j})_ref" )

    channel = configdata['channel']
    rest_mass = configdata.pop('rest_mass',None)

    graphs = ['compare_spectrums','final_spectrum']
    ni_width = 80.0
    spacial_extent = configdata.pop('Ns',0)

    for graph in graphs:
        omit = []
        files = {}
        if graph not in configdata.keys():
            continue

        if "omit" in configdata[graph].keys():
            omit = configdata[graph]["omit"]
        file_directory = configdata[graph]['file_directory']
        spectrum_type = configdata[graph].pop('spectrum_type','energy')
        best_legend_loc = configdata[graph].pop('best_legend_loc', "")
        remove_xlabel = configdata[graph].pop( 'remove_xlabel', False )
        max_level = configdata[graph].pop( "max_level", 9 )
        graph_unused_levels = configdata[graph].pop( "graph_unused_levels", True )
        preliminary = configdata[graph].pop( 'preliminary', False )
        ni_width = configdata[graph].pop('ni_width', 80.0) 
        plot_zmags = configdata[graph].pop('plot_zmags', False)
        split_mom = configdata[graph].pop('split_mom', False)
        rotation = configdata[graph].pop('rotation', 0)
        root = configdata[graph].pop('root', None)
        
        if ('compare_spectrums' in configdata.keys()) and (graph=='compare_spectrums'):
            for key in configdata['compare_spectrums']['files'].keys():
                if configdata['compare_spectrums']['file_directory']:
                    files[key] = os.path.join(channel,file_directory,configdata['compare_spectrums']['files'][key])
            
            do_scattering_particles = configdata['compare_spectrums'].pop('do_scattering_particles', False)
            if spectrum_type=="mom":
                do_scattering_particles = False
            
            plot_ni_levels = configdata['compare_spectrums'].pop('plot_ni_levels', False)
            
            if split_mom:
                file_name = channel+f"-{spectrum_type}_spectrum_comparison_graph"
            else:
                file_name = channel+f"-{spectrum_type}_spectrum_comparison_graph.pdf"

        if ('final_spectrum' in configdata.keys()) and (graph=='final_spectrum'):
            files["this"] = os.path.join(channel,file_directory,configdata['final_spectrum']['file'])
            file_name = channel+"_spectrum_graph.pdf"

        used_levels = configdata.pop( 'used_levels', {} )
        
            
        if "remove_ref" in configdata[graph].keys():
            remove_ref = configdata[graph]["remove_ref"]
            
        if not remove_ref:
            scattering_particles = [x+"_ref" for x in scattering_particles]
            
        if spectrum_type=="energy":
            energy_key = "ecm"
            if "energy_key" in  configdata[graph].keys():
                energy_key = configdata[graph]["energy_key"]
            
        #########################################################
        ############### import info from files ##################
        #########################################################
        
        Refs = {}
        datasets = {}

        for file in files.keys():
            datasets[file] = utils.unpack_file(files[file],spectrum_type,root)
            if datasets[file] is None:
                sys.exit()

        for scat in particle_names:
            value = 0.0
            for file in files.keys():
                value = utils.find_rest_mass( datasets[file], scat, remove_ref )
                if not value:
                    continue
                if not value[0]:
                    continue
                Refs[scat] = value
                break
        if not Refs and particle_names:
            print("No values found for scattering particles:", particle_names )

        expected_keys = [f"PSQ{i}" for i in range(10)]
        possible_irreps = ['G1u', 'Hg', 'G1', 'G2', 'F1', 'F2', 'G','A1g','G1g','Hu', 'T1g','A1u','A2u','Eg','Eu',
                           'T1u','T2g','T2u','A1','A2','B1','B2','E']
        ppossible_irreps = [irrep+'p' for irrep in possible_irreps]
        npossible_irreps = [irrep+'m' for irrep in possible_irreps]
        possible_irreps+=ppossible_irreps
        possible_irreps+=npossible_irreps
        
#         diff_keys = [f"d_ecm_{i}" for i in range(15)] #probably change this to elab
        # max_level = 9
        if spectrum_type=="mom":
            energy_keys = [f'q2cm_{i}_ref' for i in range(max_level)]
        elif spectrum_type=="energy":
            if remove_ref:
                energy_keys = [f'{energy_key}_{i}' for i in range(max_level)]
            else:
                energy_keys = [f'{energy_key}_{i}_ref' for i in range(max_level)]
        else:
            print("Bad spectrum type")
            

        keys = {}
        vals = {}
        errs = {}
        levs = {}
        zmags = {}
        vals_used = {}
        keys_used = {}
        errs_used = {}
        levs_used = {}
        zmags_used = {}
        ekeys = {}
        evals = {}
        eerrs = {}
        max_level_num = 0
        for dataset in datasets.keys():
            vals[dataset] = []
            keys[dataset] = []
            errs[dataset] = []
            levs[dataset] = []
            zmags[dataset] = []
            vals_used[dataset] = []
            keys_used[dataset] = []
            errs_used[dataset] = []
            levs_used[dataset] = []
            zmags_used[dataset] = []
            evals[dataset] = []
            ekeys[dataset] = []
            eerrs[dataset] = []
            for i, (mom) in enumerate(expected_keys):
                for j, (irrep) in enumerate(possible_irreps):
                    irrep_key = f"{irrep}({mom.replace('PSQ','')})"
                    if irrep_key in omit:
                        continue
                    for k, (energy) in enumerate(energy_keys):

                        if not datasets[dataset].empty:
                            val1, err1, zmag1 = utils.select_val(datasets[dataset], mom, irrep, energy)
                        else:
                            val1, err1 = utils.select_val_ascii(files[dataset],mom,irrep,energy,spectrum_type)
                            zmag1 = None

                        if val1 is not None and err1 is not None:
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
                                        zmags_used[dataset].append(zmag1)
                                        used = True
    #                                     print(irrep_key,k,val1)
                            if not used:
                                keys[dataset].append(irrep_key)
                                vals[dataset].append(val1)
                                errs[dataset].append(err1)
                                levs[dataset].append(k)
                                zmags[dataset].append(zmag1)

                    if plot_ni_levels:
                        #get free level
                        if files[dataset].endswith("hdf5"):
                            free_levels = np.array([])
                            with h5py.File(files[dataset]) as f:
                                if f"{mom}/{irrep}" in f.keys():
                                    free_levels = f[f"{mom}/{irrep}"].attrs['free_levels'][()]
                                if free_levels.any():
                                    for nilevel in free_levels:
                                        new_nilevel = []
                                        for i in range(len(nilevel)):
                                            if type(nilevel[i])!=str:
                                                new_nilevel.append(f"{rest_mass}({nilevel[i]})")
                                            else:
                                                new_nilevel.append(nilevel[i])
                                                
                                        key = new_nilevel[0]
                                        level = f[f"single_hadrons/{new_nilevel[0]}"][()]
                                        for sh in new_nilevel[1:]:
                                            key+=sh
                                            level += f[f"single_hadrons/{sh}"][()]
                                        if mom!='PSQ0' and spacial_extent:
                                            momi = int(mom[3:])
                                            level = np.sqrt(level*level-momi*(np.pi*2.0/spacial_extent)**2)
                                        ref = f["single_hadrons/ref"][()]
                                        level /= ref
                                        # print(key, level[0], utils.bootstrap_error_by_array(level) )
                                        ekeys[dataset].append(irrep_key)
                                        evals[dataset].append(level[0])
                                        eerrs[dataset].append(utils.bootstrap_error_by_array(level))
                        else:
                            for k, (energy) in enumerate(elastic_scat_keys):

                                if not datasets[dataset].empty:
                                    val1, err1, zmag1 = utils.select_val(datasets[dataset], mom, irrep, energy)
                                else:
                                    val1, err1 = utils.select_val_ascii(files[dataset],mom,irrep,energy,spectrum_type)
                                    zmag1 = None

                                if val1 is not None and err1 is not None:
                                    if "yrange" in configdata[graph].keys():
                                        if (val1 < configdata[graph]["yrange"][0]) or (val1 > configdata[graph]["yrange"][1]):
                                            continue
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
                if remove_ref:
                    plt.ylabel(r"$E_{cm}$")
                else:
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
            set_figure(f, configdata, configdata[graph])
        elif split_mom:
            print("Cannot set 'split_mom' to True when importing function")
            return 0
            
        if "yrange" in configdata[graph].keys():
            plt.ylim( configdata[graph]["yrange"][0],configdata[graph]["yrange"][1])

        if energy_key != "dElab":
            somekey = list(indexes.keys())[0]
            try:
                minx = min(list(indexes[somekey])+list(indexes_used[somekey]))
            except ValueError as err:
                print("Missing data for plot. Check energy_key is correct, check that there is data in the file, or try omitting the yrange.")
            maxx = max(list(indexes[somekey])+list(indexes_used[somekey]))
            dd=0.4/len(files.keys())
            plt.xlim(minx-0.5-dd*(len(vals.keys())/2),maxx+0.5+dd*(len(vals.keys())/2))
            if spectrum_type=="energy":
                if 'thresholds' in configdata:
                    if configdata['thresholds']:
                        plt.xlim(minx-0.5-dd*len(vals.keys()),(maxx+dd*len(vals.keys()))*1.05+1.0)

            if spectrum_type=="energy":
                if 'thresholds' in configdata:
                    if configdata['thresholds']:
                        for threshold in configdata['thresholds'].keys():
                            threshold_value = 0.0
                            threshold_label = threshold
                            for particle in configdata['thresholds'][threshold]:
                                if particle in Refs:
                                    threshold_value+=Refs[particle][0]
                                    threshold_label=threshold_label.replace(settings.latex_format[particle], particle)
                                    threshold_label=threshold_label.replace(particle, settings.latex_format[particle])
                                else:
                                    print(f"{particle} is not listed in scattering_particles.")
                                    
                            miny,maxy = plt.ylim()
                            plt.hlines(threshold_value,minx-dd*len(vals.keys()),maxx+dd*len(vals.keys()),color='black', linestyle="--", zorder=4) 
                            plt.text( (maxx+dd*len(vals.keys()))*(1.005),threshold_value-0.0015*abs(maxy-miny), threshold_label, zorder=6,size="x-small")
            if spectrum_type=="mom":
                if "yrange" in configdata[graph].keys():
                    if (0.0>configdata[graph]["yrange"][0]) and (0.0<configdata[graph]["yrange"][1]):
                        plt.hlines(0.0,minx-dd*(len(vals.keys())/2),maxx+dd*(len(vals.keys())/2),color='black', linestyle="--", zorder=4, linewidth=1.0)
                else:
                    plt.hlines(0.0,minx-dd*(len(vals.keys())/2),maxx+dd*(len(vals.keys())/2),color='black', linestyle="--", linewidth=1.0, zorder=4)

            for i,dataset in enumerate(vals.keys()):
                shifted_array = 0.0 
                used_shifted_array = 0.0 
                if indexes_used[dataset]:
                    used_shifted_array = 0.5*dd*utils.shift_levels(np.array(indexes_used[dataset]),levs_used[dataset],vals_used[dataset],errs_used[dataset],np.zeros(len(indexes_used[dataset])))
                if not used_levels:
                    shifted_array = 0.5*dd*utils.shift_levels(np.array(indexes[dataset]),levs[dataset],vals[dataset],errs[dataset],np.zeros(len(indexes[dataset])))
                
                if used_levels:
                    marker_color = 'white'
                else:
                    marker_color = settings.colors[i]

                if used_levels:
                    unused_label = None
                else:
                    unused_label = dataset

                splitting_factor = 1.6*(i-((len(vals.keys())-1)/2))

                if (not used_levels) or (used_levels and graph_unused_levels):
                    colors = [settings.colors[i]]*len(vals[dataset])
                    mfcolor = settings.colors[i]
                    if zmags[dataset][0] is not None and plot_zmags:
                        maxz = max(zmags[dataset])
                        colors = [mc.to_rgba(c,z/maxz) for c,z in zip(colors,zmags[dataset])]
                    elif plot_zmags:
                        colors = [mc.to_rgba('white',0.0)]*len(vals[dataset])
                        mfcolor = "white"
                    
                    if len(np.nonzero(errs[dataset])[0]) and zmags[dataset][0] is not None and plot_zmags:
                        # for index, val, err, color in zip(indexes[dataset],vals[dataset],errs[dataset],colors)
                        plt.errorbar(indexes[dataset]+dd*splitting_factor+shifted_array, vals[dataset], np.array(errs[dataset]),  capsize=5, color=settings.colors[i], marker=settings.markers[i], linewidth=0.0, mfc='white',elinewidth=1.5,zorder=5,label=unused_label)
                        plt.scatter(indexes[dataset]+dd*splitting_factor+shifted_array,vals[dataset], color=colors, marker=settings.markers[i], linewidth=0.0, zorder=6)
                    elif len(np.nonzero(errs[dataset])[0]):
                        plt.errorbar(indexes[dataset]+dd*splitting_factor+shifted_array, vals[dataset], np.array(errs[dataset]),  capsize=5, color=settings.colors[i], mfc=mfcolor, marker=settings.markers[i], linewidth=0.0, elinewidth=1.5, zorder=5, label=dataset)
                    else:
                        plt.scatter(indexes[dataset]+dd*splitting_factor+shifted_array, vals[dataset], color=colors, marker=settings.markers[i],linewidth=0.0, zorder=5,label=unused_label)

                if used_levels:
                    colors = [settings.colors[i]]*len(vals_used[dataset])
                    if zmags[dataset][0] is not None and plot_zmags:
                        colors = [mc.to_rgba(c,z) for c,z in zip(colors,zmags_used[dataset])]
                    elif plot_zmags:
                        colors = ['white']*len(vals_used[dataset])
                    if len(np.nonzero(errs_used[dataset])[0]) and zmags[dataset][0] is not None:
                        plt.errorbar(indexes_used[dataset]+dd*splitting_factor+used_shifted_array,vals_used[dataset], np.array(errs_used[dataset]),  capsize=5, color=settings.colors[i], mfc='white', marker=settings.markers[i], linewidth=0.0, elinewidth=1.5, zorder=5, label=dataset)
                        plt.scatter(indexes_used[dataset]+dd*splitting_factor+used_shifted_array,vals_used[dataset], color=colors, marker=settings.markers[i], linewidth=0.0, zorder=6)
                    elif len(np.nonzero(errs_used[dataset])[0]):
                        plt.errorbar(indexes_used[dataset]+dd*splitting_factor+used_shifted_array,vals_used[dataset], np.array(errs_used[dataset]),  capsize=5, color=settings.colors[i], mfc=settings.colors[i], marker=settings.markers[i], linewidth=0.0, elinewidth=1.5, zorder=5, label=dataset)
                    else:
                        plt.scatter(indexes_used[dataset]+dd*splitting_factor+used_shifted_array,vals_used[dataset], color=colors, marker=settings.markers[i], linewidth=0.0, zorder=5,label=dataset)

            if (graph=='final_spectrum') or plot_ni_levels:
                for i,dataset in enumerate(evals.keys()):
                    plt.errorbar(eindexes[dataset]+dd*(i-((len(evals.keys())-1)/2)), evals[dataset], np.array(eerrs[dataset]), color="lightgray", marker="_",linestyle="", linewidth=0.0, elinewidth=ni_width,zorder=2)
                    plt.scatter(eindexes[dataset]+dd*(i-((len(evals.keys())-1)/2)),evals[dataset],color="darkgrey", marker="_",s=ni_width*ni_width,zorder=3)
            
            latex_keys = [settings.latex_format[key.split('(')[0]]+"("+key.split('(')[1] for key in keys[somekey]+keys_used[somekey]]
            plt.xticks(utils.unique(list(indexes[somekey])+list(indexes_used[somekey])), utils.unique(latex_keys),size="small", rotation=rotation)

        else:
            ii = 0
            labeled = False
            xlabels = []
            for i, (mom) in enumerate(expected_keys):
                if split_mom:
                    plt.close()
                    f = plt.figure()
                    set_figure(f, configdata, configdata[graph])
                    ii=0
                    xlabels = []
                for j, (irrep) in enumerate(possible_irreps):
                    if f"{irrep}({mom[-1]})" not in omit:
                        for k, (energy) in enumerate(energy_keys):
                            vals = []
                            errs = []
                            for dataset in datasets.keys():
                                if not datasets[dataset].empty:
                                    val1, err1, zmag1 = utils.select_val(datasets[dataset], mom, irrep, energy)
                                else:
                                    val1, err1 = utils.select_val_ascii(files[dataset],mom,irrep,energy,spectrum_type)
                                    zmag1 = None
                                if not val1:
                                    continue
                                vals.append(val1)
                                errs.append(err1)
                            # print( mom, irrep, energy, vals, errs)
                            if len(vals)==len(datasets):
                                for d,dataset in enumerate(datasets.keys()):
                                    # print(mom, irrep, energy, dataset, vals[d], errs[d])
                                    if not labeled:
                                        plt.errorbar( ii+d*0.1, vals[d], errs[d],  capsize=5, color=settings.colors[d], mfc=settings.colors[d], marker=settings.markers[d], linewidth=0.0, elinewidth=1.5, zorder=5, label=dataset)
                                    else:
                                        plt.errorbar( ii+d*0.1, vals[d], errs[d],  capsize=5, color=settings.colors[d], mfc=settings.colors[d], marker=settings.markers[d], linewidth=0.0, elinewidth=1.5, zorder=5)
                                
                                labeled=True
                                ii+=1
                                xlabels.append( f"{settings.latex_format[irrep]}({mom[-1]}) {k}")
                if split_mom and ii:
                    labeled = False
                    plt.xticks(range(ii), xlabels, rotation=90)
                    finalize_plot(remove_xlabel, configdata, graph, best_legend_loc, rest_mass, 
                                    energy_key, spectrum_type, remove_ref, preliminary)
                    plt.savefig(os.path.join(channel,file_directory,file_name+f"-{mom}.pdf"), transparent=True)
                    plt.show()
                    plt.clf()
                            
            if not split_mom:
                plt.xticks(range(ii), xlabels, rotation=90)
        

        if not split_mom:
            finalize_plot(remove_xlabel, configdata, graph, best_legend_loc, rest_mass, 
                            energy_key, spectrum_type, remove_ref, preliminary)
                            
        print("Total Time:", timer()-start)
        if __name__=="__main__":
            plt.savefig(os.path.join(channel,file_directory,file_name), transparent=True)
            plt.show()
            plt.clf()
            

if __name__=="__main__":
    config_file = cmd_plot.pickup_config()
    compare_spectrums(config_file)