import numpy as np
import os, sys
import matplotlib.pyplot as plt
import yaml
import argparse

sys.path.append('../')
import utils
import settings
import cmd_plot
plt.style.use('../spectrum.mplstyle')

#Code used for creating plots:

# # fig 9 of https://arxiv.org/pdf/2208.03867.pdf
# # slide 12 of https://indico.fnal.gov/event/57249/contributions/268333/attachments/169682/227804/fv_spectroscopy_lattice2023_Skinner%20%286%29.pdf

# ideal for plotting random datapoints across various lattice settings

###############################################
#
# # mpi_vary_plots config example

# # list of plots to compute
# plots:
# # Output filename for the generated plot
#   - filename: m_deltag.pdf

# # Height and width of the figure (in inches)
#     figheight: 10
#     figwidth: 8

# # Alignment of subplots, can be "horizontal" or "vertical"
#     align: vertical

# # Ratios of subplot heights, used for vertical alignment
#     ratios: [2, 1]

# # Names of the data roots (sections) in the data YAML file to be plotted
#     data_roots: [mdelta, g]

# # more plots?
#   - ...

###############################################

# # mpi_vary_plots data example

# # root tag
# mdelta
# # Label for the y-axis
  # ylabel: $m_\Delta$ (MeV)

# # Data points for different sources
  # data:
    # - "Andersen et al. 2018":
      #   mpi: [280, 0.0]       # Mass of pion (MeV) and its error
      #   datapoint: [1344, 20] # Mass of Delta (MeV) and its error
      #   style_index: 1        # Style index for coloring and marking
    # - "Silvi et al. 2021":
      #   mpi: [255.4, 1.6]
      #   datapoint: [1380, 11.4]
      #   style_index: 3
    # ...

# # Physical point data
  # physical_point:
    # mpi: [139.5704, 0.002]  # Mass of pion (MeV) and its error
    # datapoint: [1232, 2]    # Mass of Delta (MeV) and its error
    # style_index: 2          # Style index for coloring and marking

###############################################

# to run

# python mpi_vary_plots.py [-h] config.yml data.yml

############ from old stuff ##################
# mu=0.211
# Fpi=93.43 #check f = F/sqrt(2)
mNsq_fpisq = 2.0*(939.5654133**2)/(130.41**2)

# def mpia_minus_0(mpi): #mpi^2/fpi^2
#     return 2.0*mpi/(8.0*np.pi*(1+mu))*(1.0+gA*gA*mu*mu/4.0/(1-(mu*mu)/4.0))
# def mpia_plus_0(mpi):
#     return -2.0*mpi*gA*gA/(16.0*np.pi*(1+mu))*mu/(1-(mu*mu)/4.0)

def mpia_minus_0(mu,gA): #mpi^2/fpi^2
    return mNsq_fpisq*mu*mu/(8.0*np.pi*(1.0+mu))*(1.0+gA*gA*mu*mu/(4-mu*mu))
def mpia_plus_0(mu,gA):
    return mNsq_fpisq*mu*mu*gA*gA*mu/(4.0*np.pi*(1.0+mu)*(4.0-mu*mu))
def mpia_32_12(mpi,gA):
    return mpia_plus_0(mpi,gA) - mpia_minus_0(mpi,gA)
def mpia_12_12(mpi,gA):
    return mpia_plus_0(mpi,gA) + 2.0*mpia_minus_0(mpi,gA)

###############################################

# Function to check if a handle (line) with a given label is already in the legend
def check_handle_if_used(handles, label):
    for line in handles:
        if line.get_label() == label:
            return True
    return False

# Function to plot data for varying m_pi
def plot_mpi_vary(ax, yaml_file, root, legend=False, handles=[], keep_xlabels=True, ga_lines = []):
    with open(yaml_file, 'r') as f:
        plot_info = yaml.safe_load(f)

    # Extract data and physical point information from the YAML file
    data = plot_info[root]['data']
    physical_point = plot_info[root]['physical_point']
    
    if ga_lines:
        xaxis = 'mu'
    else:
        xaxis = 'mpi'

    # Unpack data for physical point
    physical_point = [
        physical_point[xaxis][0],
        physical_point[xaxis][1],
        physical_point['datapoint'][0],
        physical_point['datapoint'][1],
        physical_point['style_index']
    ]

    skip = []
    for item in data:
        if xaxis not in item[list(item.keys())[0]]:
            skip.append(item)
            
    for item in skip:
        data.remove(item)
            
    # Unpack data for all other labels
    data_list = [{
        list(item.keys())[0]: [
            item[list(item.keys())[0]][xaxis][0],
            item[list(item.keys())[0]][xaxis][1],
            item[list(item.keys())[0]]['datapoint'][0],
            item[list(item.keys())[0]]['datapoint'][1],
            item[list(item.keys())[0]]['style_index']
        ]
    } for item in data ]
        

    # Plot physical point and data
    linep, = ax.plot(physical_point[0], physical_point[2], color=settings.colors[physical_point[4]], lw=0.0, marker="*", label="Physical point", markersize=15)
    linep.set_label("Physical point")
    if not check_handle_if_used(handles, "Physical point"):
        handles.append(linep)

    for data in data_list:
        for label in data:
            linep = ax.errorbar(x=data[label][0], xerr=data[label][1], y=data[label][2], yerr=data[label][3], capsize=5, color=settings.colors[data[label][4]], marker=settings.markers[data[label][4]], linewidth=0.0, elinewidth=1.5, label=label, markersize=7)
            linep.set_label(label)
            if not check_handle_if_used(handles, label):
                handles.append(linep)

    if xaxis=='mpi':
        mpis = [np.floor(item[list(item.keys())[0]][0]) for item in data_list]
        mpis.append(np.floor(physical_point[0]))
    else:
        mpis = [np.round(item[list(item.keys())[0]][0],2) for item in data_list]
        mpis.append(np.round(physical_point[0],2))

    #plot function lines
    x = np.linspace(min(mpis),max(mpis),1000)
    lss = ['solid', 'dashed', 'dashdot', 'dotted']
    for ga in ga_lines:
        ax.plot(x,globals()[plot_info[root]['function']](x,ga),color="black", ls=lss.pop(0))
        
    if keep_xlabels:
        if xaxis == 'mpi':
            ax.set_xlabel(r"$m_\pi$ (MeV)")
        else:
            ax.set_xlabel(r"$\mu$")
        ax.set_xticks(ticks=list(set(mpis)))
    else:
        ax.set_xticks(ticks=[])
    ax.set_ylabel(plot_info[root]['ylabel'])
    if legend:
        ax.legend()
        

# Main function
def main():
    plot_yaml_file, data_yaml_file = cmd_plot.pickup_config_and_info()

    with open(plot_yaml_file, 'r') as f:
        plot_info = yaml.safe_load(f)

    for plot in plot_info['plots']:
        nx_plots = 1
        ny_plots = 1
        # Extract settings from the plot YAML file
        figheight = plot.pop('figheight', 6)
        figwidth = plot.pop('figwidth', 16)
        align = plot.pop('align', 'horizontal')
        data_roots = plot.pop('data_roots', [])
        ratios = plot.pop('ratios', [1] * len(data_roots))
        filename = plot.pop('filename', "m_delta.pdf")

        # Handle missing data_roots setting
        if not data_roots:
            print("No data given for plotting. Add 'data_roots: [item1, item2...]' to your config file.")
            return

        xratios = [1]
        yratios = [1]
        keep_xlabels = True

        # Handle alignment setting
        if align != 'horizontal' and align != 'vertical':
            print("Invalid entry for 'align'. Must be 'horizontal' or 'vertical'")
            return
        elif align == 'horizontal':
            nx_plots = len(data_roots)
            xratios = ratios
        elif align == 'vertical':
            ny_plots = len(data_roots)
            yratios = ratios
            keep_xlabels = False

        # Create subplots
        f, axes = plt.subplots(ny_plots, nx_plots, figsize=(figwidth, figheight), gridspec_kw={'width_ratios': xratios, 'height_ratios': yratios})

        # Handle non-array axes
        if type(axes) != np.ndarray:
            axes = [axes]

        handles = []

        # Plot data for each root
        function_ga_input = plot.pop('function_ga_input',[])
        for i, root in enumerate(data_roots):
            if i == len(data_roots) - 1:
                keep_xlabels = True
            plot_mpi_vary(axes[i], data_yaml_file, root, False, handles, keep_xlabels, function_ga_input)
        axes[0].legend(handles=handles)

        plt.tight_layout()
        plt.savefig(filename)
        plt.show()

if __name__ == "__main__":
    main()
