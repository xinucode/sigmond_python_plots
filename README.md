# sigmond_python_plots
Python plots from sigmond data. See https://github.com/andrewhanlon/sigmond and https://github.com/andrewhanlon/sigmond_scripts .

## contents
spectrum.mplstyle - style file for matplotlib plots that is applied to all plots

settings.py - Key plot settings such as colors, latex, and data ordering

utils.py - common functions needed for any given plotting script. Input functions and
    error functions are here.

spectrums/compare_spectrums.py - python script to generate various spectrum plots, to both plot a 
    spectrum and also compare several spectrums. Requires a config file as an argument

tmins/tmin_plots.py - python script that generate rest mass or spectra tmin plots in both
    xmgrace and python formats
    
tmins/stability_tmin_plots.py - python script to combine the tmin plots of several pivots

m_delta/ - m_delta comparison plot in development

scat_length/ - scattering length plot in development

rebin/ - rebin analysis in development

## basic commands
to run example spectrum plot:
```
python spectrums/compare_spectrums.py spectrums/hexaquark.yml
```

to run the tmin plots script
```
python tmin/tmin_plots.py some_config.yml
```
In order to run, `some_config.yml` will need to be replaced by a config file using the 
specification for the tmin config file below

to run the stability tmin plots script
```
python tmin/stability_tmin_plots.py some_config.yml
```
In order to run, `some_config.yml` will need to be replaced by a config file using the 
specification for the stability tmin config file below

## Sample config files
### Sample spectrums config file
```
channel: isoquartet_nonstrange_fermionic #name of channel and corresponding 
                                            #subdirectory where relevant files are located
scattering_particles: [N, pi] #list of scattering particle names
rest_mass: pi #the particle used to normalize the energy values
title: $I=\nicefrac{3}{2}$ #(optional) if unspecified, plot will not print a title in the legend.
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
```

### Sample tmin config file
```
Sample Yaml Config file
------
rest_mass: #(optional) generates rest mass plot(s)
    out_dir: single_hadrons #output directory where the software will save the final plots
    out_type: python #options are "python" or "xmgrace". Type of plot to generate
    combine: true #(optional) if out_type==python, puts plots side by side in one file in order of input. default: false
    particles: #dict of the particles to plot Ex: pi, N, K...
        pi: #dict that specifies the 'files' and 'outfilestub' for the given particle
            files: #list of filepaths xmgrace tmin plots to be combined into one tmin plot (assumed to have been 
                generated by sigmond_scripts)
                - /latticeQCD/raid3/sarahski/lqcd/D200_R000/isoquartet_nonstrange_fermionic_bootstrap/.sigmond/plots/spectrum/tmin_plots/isoquartet_nonstrange_fermionic_colin/rebin20/single_hadrons/tmin_fit_isotriplet-S0-P000-A1um-PSS0-0_30000A1um-P[SS0]-0T2-25-20_0.agr
                - /latticeQCD/raid3/sarahski/lqcd/D200_R000/isoquartet_nonstrange_fermionic_bootstrap/.sigmond/plots/spectrum/tmin_plots/isoquartet_nonstrange_fermionic_colin/rebin20/single_hadrons/tmin_fit_isotriplet-S0-P000-A1um-PSS0-0_30000A1um-P[SS0]-0T2-25-20_4.agr
            outfilestub: pi0_tmin_all #name attached to output files

        
spectrum_tmins: #(optional) generates conbined tmin plot(s) for spectra
    plot_format: 2 #options are integers 1,2,3. Used to specify format of output plot. 1 - each level gets their own plot; 
        2 - all levels in a basis are graphed in one plot and each level is colored differently; 
        3 - all levels in a basis are graphed in one plot and each fit type is colored differently
    out_type:  #options are "python" or "xmgrace". Type of plot to generate
    channels: #list of channels to generate plots for
      - name: isodoublet_nonstrange #name of the channel (will be used in output filenames)
        out_dir: isodoublet_npi #name of output directory where all out files will be stored
        max_level: 10 #max expected level in any given basis. Can be a higher value than max. If less, then those 
            values will not be included
        dir: /latticeQCD/raid3/sarahski/lqcd/D200_R000/isodoublet_nonstrange_bootstrap/.sigmond/plots/spectrum/tmin_plots/isodoublet_nonstrange_colin/rebin20 #directory of input xmgrace files (assumed to have been generated by sigmond_scripts)
        omit: #(optional) list of fits, bases, and files to exclude from final plots (the whole file path of files
            should not be included)
          - geometric ratio fit
          - isodoublet_S0_G1_P1_single_pivot_n4_m4_d8_c50
          - tmin_fit_isodoublet-S0-PSQ2-G-DROT-3_20p2G-ROT-3T10-25-20_0_3.agr
```
### Sample stability tmin config file
```
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
```

## Input data
### spectrum data
There are several different available input filetypes, though they were developed for specific output data from
sigmond, sigmond scripts, and collaborators. 

Spectrum data, which consists of the spectrum levels and their errors, comes in the form of energy or tranverse
momentum squared (qsqr), which is why the "spectrum_type" tag is required for both plot types, with "mom" 
referring to the latter case.

Due to the nature of the data, there are also two different methods of printing out the values and error bars for
any given level. Either the value and error would be printed as columns, calculated in the analysis scripts, or
the bootstrap/jackknife samples will be printed to the file. Commonly, the latter is in the form of hdf5 files.
The functions in utils.py address the input and calculation of these errors when necessary. 

#### CSV files
At this moment in time, sigmond scripts only produces csv files with energy values in them. CSV files will have
"obs","val", and "err" labels. The "obs" column will state the label for the energy value and state the irrep, 
momentum, and possible the level of the value. "val" gives energy value. "err" gives symmetric error. 

#### HDF5 files
HDF5 files are more versatile. Sigmond/sigmond scripts will output both energy and momentum files. Currently,
hdf5 files are all samplings files, meaning that rather than printing the error, they print the samples used
to calculate this error. In order to properly calculate the error, jackknife or bootstrap needs to be stated.
Due to the strong desire across the collaboration for boostrap, the default method is boostrap, but for
verification of boostrap methods, it is sometimes necessary or preferential to use jackknife method. If a file
uses jackknife sampling, the filename is appended to the utils.jackknife_sampling_methods list (likely,
another method should be used, but I use it so infrequently, and I always label the jackknife files, so this
is what I have coded rn).

#### DAT files
Sigmond nor sigmond scripts produce dat files. The contents of the dat files are the output of my collaborator 
Fernando Romero-Lopez's phse shift analysis code which is basically a csv file with columns "irrep", "d^2",
"predicted qsqr values","input qsqr values","something else, I dunno", "(mN/mpi)^2".

#### directories
My collaborator, John Bulava,  prints his calculated files into organized directories. He gives my a directory 
with the organized files and in a given file in that directory will have one energy or qsqr value and it's error. 
The syntax of his output may change, so beware. Right now, the function utils.select_val_ascii() correctly searches 
the directory and finds the value in the file correctly.


### tmin data
#### xmgrace files
when the spectrum task of the sigmond_scripts software is ran, xmgrace plots are saved in an organized directory 
that contains these plots. These are parsed and used to generate python or combined xmgrace plots.
