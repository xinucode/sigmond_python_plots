# sigmond_python_plots
Python plots from sigmond data

## contents
utils.py - common functions needed for any given plotting script. Input functions and
    error functions are here.
compare_spectrums.py - python script to generate various spectrum plots, to both plot a 
    spectrum and also compare several spectrums. Requires a config file as an argument
rest_mass/ - rest mass tmin plots and others in development
    
## basic commands
to run example spectrum plot:
python compare_spectrums.py hexaquark.yml

## Sample config files
### Sample spectrum config file
channel: isoquartet_nonstrange_fermionic #name of channel and corresponding 
                                            #subdirectory where relevant files are located
scattering_particles: [N, pi] #list of scattering particle names
rest_mass: pi #the particle used to normalize the energy values
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
    shift: ['Hg(0)','G1(1)','G2(1)','G(2)','G(3)','F1(3)','G1(4)','G2(4)'] #(optional) list of irreps(d^2) to zigzag
                                                                            #shift the levels so they don't overlap
    omit: [G1g(0),Hu(0)] #(optional) list of irreps(d^2) to omit; otherwise, it graphs all irreps
    yrange: [0.0,2.5] #(optional) manually select the yrange, otherwise matplotlib automatically sets it
    plot_ni_levels: true #(optional) default false; plots the non interacting levels for each spectrum
    ni_width: 40 #(optional) default 80; sets the width of the non-interacting levels (if plot_ni_levels==true)
      
final_spectrum: #(optional) generates graph that plots just one spectrum
    spectrum_type: energy #value can be "energy" or "mom" determines how the files are parsed and how the graph 
                            #is labelled
    file_directory: final #subdirectory inside "channel" where the data can be found and the graph will be produced
    file: energy_estimates_isoquartet_nonstrange_fermionic_colin_rebin20_Bootstrap_8-16.csv #file with data
    shift: ['Hg(0)','G1(1)','G2(1)','G(2)','G(3)','F1(3)','G1(4)','G2(4)'] #(optional) list of irreps(d^2) to zigzag
                                                                            #shift the levels so they don't overlap
    ni_width: 40 #(optional) default 80; sets the width of the non-interacting levels
    best_legend_loc: lower left #(optional) put legend where desired. If omitted, no legend. If legend is
                                #desired but unsure where, put "best"
    fig_width: 14 #(optional) sets width of this figure, supercedes the setting for all figures
    fig_height: 6 #(optional) sets height of this figure, supercedes the setting for all figures
    omit: [G1g(0),Hu(0)] #(optional) list of irreps(d^2) to omit; otherwise, it graphs all irreps
    yrange: [0.0,2.5] #(optional) manually select the yrange, otherwise matplotlib automatically sets it
    
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
Sigmond does not produce dat files. The contents of the dat files are the output of my collaborator 
Fernando Romero-Lopez's phse shift analysis code which is basically a csv file with columns "irrep", "d^2",
"predicted qsqr values","input qsqr values","something else, I dunno", "(mN/mpi)^2".

#### directories
My collaborator, John Bulava,  prints his calculated files into organized directories. He gives my a directory 
with the organized files and in a given file in that directory will have one energy or qsqr value and it's error. 
The syntax of his output may change, so beware. Right now, the function utils.select_val_ascii() correctly searches 
the directory and finds the value in the file correctly.