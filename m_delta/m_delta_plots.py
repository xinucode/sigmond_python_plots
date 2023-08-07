import numpy as np
import os, sys
import matplotlib.pyplot as plt
import yaml
import argparse

sys.path.append('../')
import utils
import settings
plt.style.use('../spectrum.mplstyle')

def check_handle_if_used( handles, label):
    for line in handles:
        if line.get_label()==label:
            return True
    return False
    

def plot_mpi_vary(ax, yaml_file, root, legend = False, handles = [], keep_xlabels = True):
    with open(yaml_file, 'r') as f:
        plot_info = yaml.safe_load(f)
        
    data = plot_info[root]['data']
    physical_point = plot_info[root]['physical_point']
    
    physical_point = [physical_point['mpi'][0], physical_point['mpi'][1], physical_point['datapoint'][0], physical_point['datapoint'][1], physical_point['style_index']]
    data = {item: [data[item]['mpi'][0], data[item]['mpi'][1], data[item]['datapoint'][0], data[item]['datapoint'][1], data[item]['style_index']] for item in data}
    
    linep, = ax.plot(physical_point[0], physical_point[2], color=settings.colors[physical_point[4]], lw=0.0, marker="*", label="Physical point", markersize=15)
    linep.set_label("Physical point")
    if not check_handle_if_used( handles, "Physical point"):
        handles.append(linep)

    for label in data.keys():
        linep = ax.errorbar(x=data[label][0], xerr=data[label][1], y=data[label][2], yerr=data[label][3], capsize=5, color=settings.colors[data[label][4]], marker=settings.markers[data[label][4]], linewidth=0.0, elinewidth=1.5, label=label, markersize=7)
        linep.set_label(label)
        if not check_handle_if_used( handles, label):
            handles.append(linep)

    mpis = [np.floor(data[label][0]) for label in data.keys()]
    mpis.append(np.floor(physical_point[0]))

    if keep_xlabels:
        ax.set_xlabel(r"$m_\pi$ (MeV)")
        ax.set_xticks(ticks=list(set(mpis)))
    else:
        ax.set_xticks(ticks=[])
    ax.set_ylabel(plot_info[root]['ylabel'])
    if legend:
        ax.legend()

def main():
    data_yaml_file = 'm_delta_and_g_data.yml'
    plot_yaml_file = 'm_delta_and_g_plot.yml'
    nx_plots = 1
    ny_plots = 1
    with open(plot_yaml_file, 'r') as f:
        plot_info = yaml.safe_load(f)
    
    figheight = plot_info.pop('figheight', 6)
    figwidth = plot_info.pop('figwidth', 16)
    align = plot_info.pop('align', 'horizontal')
    data_roots = plot_info.pop('data_roots',[])
    ratios = plot_info.pop('ratios', [1]*len(data_roots) )
    
    if not data_roots:
        print("No data given for plotting. Add 'data_roots: [item1, item2...]' to your config file.")
        return
        
    xratios = [1]
    yratios = [1]
    keep_xlabels = True
        
    if align!='horizontal' and align!='vertical':
        print("Invalid entry for 'align'. Must be 'horizontal' or 'vertical'")
        return
    elif align=='horizontal':
        nx_plots = len(data_roots)
        xratios = ratios
    elif align=='vertical':
        ny_plots = len(data_roots)
        yratios = ratios
        keep_xlabels = False
    
    f, axes = plt.subplots(ny_plots, nx_plots, figsize=(figwidth, figheight), 
                gridspec_kw={'width_ratios': xratios, 'height_ratios': yratios})

    handles = []
    for i,root in enumerate(data_roots):
        if i==len(data_roots)-1:
            keep_xlabels=True
        plot_mpi_vary(axes[i], data_yaml_file, root, False, handles, keep_xlabels)
    axes[0].legend(handles=handles)

    plt.tight_layout()
    plt.savefig("m_delta.pdf")
    plt.show()

if __name__ == "__main__":
    main()
