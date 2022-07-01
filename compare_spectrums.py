import numpy as np
import h5py
import os, sys
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import image as mpimg
import utils
import yaml
import argparse

titles = {"isodoublet_nonstrange":r"$I=\nicefrac{1}{2}$","isoquartet_nonstrange_fermionic":r"$I=\nicefrac{3}{2}$"}


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
parser = argparse.ArgumentParser()
parser.add_argument("config", help="config file")
args = parser.parse_args()

with open(args.config, "r") as yamlfile:
    configdata = yaml.load(yamlfile, Loader=yaml.FullLoader)
    
particle_names = configdata['scattering_particles']

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
if len(particle_names)==2:
    for i in range(max_mom):
        for j in range(max_mom):
            elastic_scat_keys.append( f"{particle_names[0]}({i}){particle_names[1]}({j})_ref" )
    for i in range(max_mom):
        for j in range(max_mom):
            elastic_scat_keys.append( f"{particle_names[1]}({i}){particle_names[0]}({j})_ref" )

channel = configdata['channel']
rest_mass = configdata['rest_mass']

graphs = ['comparison','final_spectrum']

for graph in graphs:
    files = {}
    if graph not in configdata.keys():
        continue
    if ('comparison' in configdata.keys()) and (graph=='comparison'):
        file_directory = configdata['comparison']['file_directory']
        for key in configdata['comparison']['files'].keys():
            if configdata['comparison']['file_directory']:
                files[key] = os.path.join(channel,file_directory,configdata['comparison']['files'][key])
        spectrum_type = configdata['comparison']['spectrum_type']
        do_scattering_particles = bool(configdata['comparison']['do_scattering_particles'])

        best_legend_loc=configdata['comparison']['best_legend_loc']

        if spectrum_type=="mom":
            do_scattering_particles = False

        file_name = channel+f"-{spectrum_type}_spectrum_comparison_graph.jpg"

    if ('final_spectrum' in configdata.keys()) and (graph=='final_spectrum'):
        file_directory = configdata['final_spectrum']['file_directory']
        best_legend_loc=configdata['final_spectrum']['best_legend_loc']
        files["this"] = os.path.join(channel,file_directory,configdata['final_spectrum']['file'])
        spectrum_type="energy"
        file_name = channel+"_spectrum_graph.jpg"
#         non_interacting_levels = configdata['final_spectrum']['non_interacting_levels']

    Refs = {}
    datasets = {}

    #import info from files
    for file in files.keys():
        datasets[file] = utils.unpack_file(files[file])

    for scat in particle_names:
        value = 0.0
        for file in files.keys():
            value = utils.find_rest_mass( datasets[file], scat )
            if value:
                Refs[scat] = value
                break

    used_levels=configdata['used_levels']

    expected_keys = ['PSQ0', 'PSQ1', 'PSQ2', 'PSQ3', 'PSQ4']
    possible_irreps = ['G1u', 'G1g', 'Hg', 'G1', 'G2', 'F1', 'F2', 'G','Hu']
    if spectrum_type=="mom":
        energy_keys = ['q2cm_0_ref', 'q2cm_1_ref', 'q2cm_2_ref', 'q2cm_3_ref', 'q2cm_4_ref', 'q2cm_5_ref', 'q2cm_6_ref', 'q2cm_7_ref', 'q2cm_8_ref', 'q2cm_9_ref', 'q2cm_10_ref']
    elif spectrum_type=="energy":
        if remove_ref:
            energy_keys = ['ecm_0', 'ecm_1', 'ecm_2', 'ecm_3', 'ecm_4', 'ecm_5', 'ecm_6', 'ecm_7', 'ecm_8', 'ecm_9', 'ecm_10']
        else:
            energy_keys = ['ecm_0_ref', 'ecm_1_ref', 'ecm_2_ref', 'ecm_3_ref', 'ecm_4_ref', 'ecm_5_ref', 'ecm_6_ref', 'ecm_7_ref', 'ecm_8_ref', 'ecm_9_ref', 'ecm_10_ref']
    else:
        print("Bad spectrum type")

    keys = {}
    vals = {}
    errs = {}
    vals_used = {}
    keys_used = {}
    ekeys = {}
    evals = {}
    eerrs = {}
    for dataset in datasets.keys():
        vals[dataset] = []
        keys[dataset] = []
        errs[dataset] = []
        vals_used[dataset] = []
        keys_used[dataset] = []
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
                        keys[dataset].append(irrep_key)
                        vals[dataset].append(val1)
                        errs[dataset].append(err1)
                        if used_levels[irrep_key]:
                            if k in used_levels[irrep_key]:
                                vals_used[dataset].append(val1)
                                keys_used[dataset].append(irrep_key)
                                
                for k, (energy) in enumerate(elastic_scat_keys):

                    if not datasets[dataset].empty:
                        val1, err1 = utils.select_val(datasets[dataset], mom, irrep, energy)
                    else:
                        val1, err1 = utils.select_val_ascii(files[dataset],mom,irrep,energy,spectrum_type)

                    if val1 is not None and err1 is not None:
                        irrep_key = f"{irrep}({mom.replace('PSQ','')})"
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
    for dataset in datasets.keys():
        indexes[dataset] = np.zeros(len(keys[dataset]))
        for j, key in enumerate(keys[dataset]):
            this_list = list(set(keys[dataset]))
            this_list.sort(key=utils.sort_by_mom)
            for i, unique_key in enumerate(this_list):
                if key == unique_key:
                    indexes[dataset][j] = i

    indexes_used = {}
    for dataset in datasets.keys():
        indexes_used[dataset] = np.zeros(len(keys_used[dataset]))
        for j, ukey in enumerate(keys_used[dataset]):
            for k, key in enumerate(keys[dataset]):
                if ukey==key:
                    indexes_used[dataset][j] = indexes[dataset][k]
                    break
    eindexes = {}
    for dataset in datasets.keys():
        eindexes[dataset] = np.zeros(len(ekeys[dataset]))
        for j, ukey in enumerate(ekeys[dataset]):
            for k, key in enumerate(keys[dataset]):
                if ukey==key:
                    eindexes[dataset][j] = indexes[dataset][k]
                    break
                    
    if ('comparison' in configdata.keys()) and (graph=='comparison'):
        #plot scattering particles
        if do_scattering_particles:
            f = plt.figure()
            f.set_figwidth(6)
            f.set_figheight(12)
            dd = 0.1

            for i,dataset in enumerate(scat_vals.keys()):
                plt.scatter(scat_indexes2+dd*i,scat_vals[dataset],color=utils.colors[i], marker=utils.markers[i], label = dataset)
                plt.errorbar(scat_indexes2+dd*i, scat_vals[dataset], scat_errs[dataset], fmt='.', capsize=5,color=utils.colors[i])
            plt.ylabel(r"$E_{cm}/m_\pi$")
            plt.xlabel("Noninteracting Scattering Levels")
            plt.legend()
            plt.xticks(scat_indexes2, scat_keys2)
            plt.title(channel)
            plt.savefig(os.path.join(channel,configdata['comparison']['file_directory'],channel+"-ni_comparison_graph-"+'_'.join(datasets.keys()).replace(" ","-")+".jpg"))
            
            
    #plot spectrum
    f = plt.figure()
    f.set_figwidth(15)
    f.set_figheight(10)
    dd=0.1
    plt.style.use('spectrum.mplstyle')

    somekey = list(indexes.keys())[0]
    if not remove_ref:
        if configdata['thresholds']:
            for threshold in configdata['thresholds'].keys():
                threshold_value = 0.0
                threshold_label = threshold
                for particle in configdata['thresholds'][threshold]:
                    threshold_value+=Refs[particle][0]
                    threshold_label=threshold_label.replace(utils.latex_format[particle], particle)
                    threshold_label=threshold_label.replace(particle, utils.latex_format[particle])

                plt.hlines(threshold_value,min(indexes[somekey]-dd*(len(vals.keys())/2)),max(indexes[somekey]+dd*(len(vals.keys())/2)),color='black', linestyle="solid") 
                plt.text( max(indexes[somekey])*(1-0.125),threshold_value*0.985, threshold_label)

    for i,dataset in enumerate(vals.keys()):
        plt.errorbar(indexes[dataset]+dd*(i-((len(vals.keys())-1)/2)), vals[dataset], np.array(errs[dataset]),  capsize=5, color=utils.colors[i], marker=utils.markers[i],linestyle="", linewidth=0.0, elinewidth=2.0,mfc='white',zorder=3)
        plt.scatter(indexes_used[dataset]+dd*(i-((len(vals.keys())-1)/2)),vals_used[dataset],color=utils.colors[i], marker=utils.markers[i], label = dataset, zorder=4)
        
    if ('final_spectrum' in configdata.keys()) and (graph=='final_spectrum'):
        for i,dataset in enumerate(evals.keys()):
            plt.errorbar(eindexes[dataset]+dd*(i-((len(evals.keys())-1)/2)), evals[dataset], np.array(eerrs[dataset]), color="lightgray", marker="_",linestyle="", linewidth=0.0, elinewidth=20.0,zorder=1)
            plt.scatter(eindexes[dataset]+dd*(i-((len(evals.keys())-1)/2)),evals[dataset],color="darkgrey", marker="_",s=400,zorder=2)

    if spectrum_type=='energy':
        plt.ylabel(r"$E_{cm}/$"+utils.latex_format[rest_mass])
    else:
        plt.ylabel(r"$q^2_{cm}$")
    plt.xlabel(r"$\Lambda(d^2)$")
    latex_keys = [utils.latex_format[key.split('(')[0]]+"("+key.split('(')[1] for key in keys[somekey]]
    plt.xticks(indexes[somekey], latex_keys)

    if (graph=='comparison'):
        plt.legend(title=titles[channel], loc=best_legend_loc)
    plt.savefig(os.path.join(channel,file_directory,file_name))
    plt.show()
    plt.clf()