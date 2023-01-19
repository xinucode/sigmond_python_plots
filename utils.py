import pandas as pd
import numpy as np
import h5py
import os, sys
from zipfile import ZipFile

import settings

def unique(an_ordered_list):
    if len(an_ordered_list)<=1:
        return an_ordered_list
    final_list = []
    element = an_ordered_list[0]
    final_list.append(element)
    for item in an_ordered_list:
        if item!=element and item not in final_list:
            element=item
            final_list.append(element)
    return final_list

# def zig_zag_shifts( irreps, levels, used_levels ):
#     shifts = np.zeros(len(levels))
#     used_shifts = np.zeros(len(used_levels))
#     zigzag_index = 0
#     last_index = -1
#     for i,index in enumerate(indices):
#         if keys[i] in irreps:
#             if last_index==index:
#                 zigzag_index+=1
#             else:
#                 zigzag_index=0
#                 last_index = index
#             shifts[i] = zigzag_shifts[zigzag_index]
                
#     zigzag_index = 0
#     last_index = -1
#     for i,index in enumerate(used_indices):
#         if used_keys[i] in irreps:
# #             print(keys.index(used_keys[i]))
#             for used_level in used_levels[used_keys[i]]:
#                 used_shifts[i] = shifts[ keys.index(used_keys[i])+used_level]
        
#     return shifts,used_shifts

#uses the irrep and momentum in form G1u(0) and orders based on momentum number + irrep's "alphabetical" value
def sort_by_mom(irrep):
    parts = irrep.split("(")
    return float(parts[1][0])+float(settings.alphabetical[parts[0]])

#goes through data file and inserts info into the dataframe (pandas)
def split_obs_col(dataset):
    obs_names = dataset['obs']
    obs_mom = []
    obs_irrep = []
    obs_level = []
    for i, (obs) in enumerate(obs_names):
        obs_list = obs.split('/')
        if obs_list[0] == settings.single_had_key:
            obs_mom.append(None)
            obs_irrep.append(None)
            obs_level.append(obs_list[1])
        else:
            obs_mom.append(obs_list[0])
            obs_irrep.append(obs_list[1])
            obs_level.append(obs_list[2])

    dataset.insert(1, "obs-mom", obs_mom)
    dataset.insert(2, "obs-irrep", obs_irrep)
    dataset.insert(3, "obs-level", obs_level)
    
#finds a specific value in the dataframe "dataset" for a given momentum, irrep, and energy level
def select_val(dataset, mom, irrep, energy):
    if mom is None and irrep is None:
        reduced_dataset = dataset[(dataset["obs-level"]==energy)]
    else:
        reduced_dataset = dataset[(dataset["obs-mom"]==mom) & (dataset["obs-irrep"]==irrep) & (dataset["obs-level"]==energy)]
    if len(reduced_dataset["val"])==0:
        return None, None
    elif len(reduced_dataset["val"])==1:
        values = np.array(reduced_dataset["val"][:])
        errors = np.array(reduced_dataset["err"][:])
        return values[0],errors[0]
    else:
        print(f"Duplicate Values: {mom}, {irrep}, {energy}")
        values = np.array(reduced_dataset["val"][:])
        errors = np.array(reduced_dataset["err"][:])
        return values[0],errors[0]
    
#takes and the csv file produced by sigmond scripts and returns a dataframe with the data
#sampling method has already been used to calculate errors in sigmond scripts
def retrieve_sigmond_script_data(file):
    dataset = pd.read_csv(file)
    split_obs_col(dataset)
    return dataset

#takes the hdf5 file and puts the data in the dataframe 
#hdf5 files from sigmond scripts are a file with all the samples, must calculate the error
    #here and know whether it's a bootstrap or jackknife sampling method to calculate errors properly
    #assumes boostrap sampling unless filename is included in the jackknife sampling list above
def retrieve_sigmond_script_data_hdf5(file):
    q2s = h5py.File(file,'r')
    dataset = pd.DataFrame(columns=['obs', 'val', 'err'])
    if file in settings.jackknife_sampling_methods:
        print("jackknife success")
    for mom in q2s.keys():
        for basis in q2s[mom].keys():
            try:
                for value in q2s[mom+'/'+basis].keys():
                    if file in settings.jackknife_sampling_methods:
                        df2 = pd.DataFrame([[mom+'/'+basis+'/'+value, q2s[mom+'/'+basis+'/'+value][0], jackknife_error_by_array(q2s[mom+'/'+basis+'/'+value])]], columns=['obs', 'val', 'err'])
                    else:
                        df2 = pd.DataFrame([[mom+'/'+basis+'/'+value, q2s[mom+'/'+basis+'/'+value][0], bootstrap_error_by_array(q2s[mom+'/'+basis+'/'+value])]], columns=['obs', 'val', 'err'])
                    dataset = dataset.append(df2)
            except AttributeError:
                pass
    split_obs_col(dataset)
    q2s.close()
    return dataset

#retrieves data from Fernando's print statements that I put into a
# dat file
def retrieve_sigmond_script_data_dat(file,spectrum_type):
    f = open(file,"r")
    line = f.readline().strip()
    f.close()
    if line!="energy" and line!="mom":
        print("ERROR: spectrum type not sepcified at top of dat file. Please add 'energy' or 'mom' at the",
              "first line of file.")
        return None
    
    dataframe_input = pd.read_csv(file,sep=" ",header=0,names=["obs-irrep","mom","qsqr","a","b","mass"])
    irrep = np.array(dataframe_input["obs-irrep"])+np.array(dataframe_input["mom"],dtype="str")
    obs_mom = ["PSQ"+str(int(x)) for x in list(dataframe_input["mom"])]
    levels = np.zeros( len(dataframe_input["obs-irrep"]) )
    mass = dataframe_input["mass"][0]
    for i,row in enumerate(irrep):
        levels[i]=np.count_nonzero(irrep[0:i] == row)
        
    if spectrum_type=="energy":
        levels2 = ["ecm_"+str(int(level))+"_ref" for level in list(levels)]
    else:
        levels2 = ["q2cm_"+str(int(level))+"_ref" for level in list(levels)]
    dataframe_input.insert(1, "obs-mom", obs_mom)
    dataframe_input.insert(2, "obs-level", levels2)
    if spectrum_type=="energy" and line=="mom":
        ecm = np.sqrt(1.0+np.array(dataframe_input["qsqr"]))+np.sqrt(mass*mass+np.array(dataframe_input["qsqr"]))
    elif (spectrum_type=="mom" and line=="mom") or (spectrum_type=="energy" and line=="energy"):
        ecm = np.array(dataframe_input["qsqr"])
    else:
        print("Code this (utils.retrieve_sigmond_script_data_dat) if possible.")
        return None
    dataframe_input.insert(3, "val", ecm)
    errs = np.zeros( len(dataframe_input["obs-irrep"]) )
    dataframe_input.insert(4, "err", errs)
    dataset = pd.DataFrame(columns=["obs-mom","obs-irrep","obs-level", 'val', 'err'],data=dataframe_input[["obs-mom","obs-irrep","obs-level","val","err"]])
    return dataset

def retrieve_barbara_data(file):
    q2s = h5py.File(file,'r')
    dataset = pd.DataFrame(columns=['obs', 'val', 'err'])
    for key in q2s.keys():
        obs_mom = key.split("_")[0]
        obs_irrep = key.split("_")[1]
        for level in q2s[key]["1exp"]["t0_12"]["Tmin"]["Correlated"]["Mean"].keys():
            index = int(level.split("_")[1])
            data = np.array(q2s[key]["1exp"]["t0_12"]["Tmin"]["Correlated"]["Mean"][level])
            chisqr = data[-1]
            tmin = data[0]
            tmax = data[1]
            dof = tmax-tmin-2
            chisqr_dof = chisqr/dof
            for i in range(len(chisqr_dof)):
                if chisqr_dof[i]<1.0:
                    chisqr_dof[i] = 1.0/chisqr_dof[i]
            minimum_index = np.where( chisqr_dof==min(chisqr_dof) )[0][0]
            value = data[3][minimum_index]
            error = data[4][minimum_index]
            operators = q2s[key]["Operators"][index]
            print(f"{obs_mom}/{obs_irrep}/ecm_{index}",value,error,operators)
    q2s.close()
    return None
    
#retrieves data from john's ascii files
def select_val_ascii(dataset, mom, irrep, energy,energy_or_mom='energy'):
    sub_file = "dsq"+mom.replace("PSQ","")+"_"+irrep
    if sub_file in os.listdir(dataset):
        if energy_or_mom=='energy':
            subsub_file = "level_"+energy.replace("ecm_","").replace("q2cm_","").replace("_ref","").replace("elab_","") +"_ecm_over_m_err.dat" #+"_deltaE_err.dat"
        else:
            subsub_file = "level_"+energy.replace("ecm_","").replace("q2cm_","").replace("_ref","").replace("elab_","") +"_pcm_over_m_sq_err.dat"
#             print(subsub_file)
        if subsub_file in os.listdir(os.path.join(dataset,sub_file)):
            f=open(os.path.join(dataset,sub_file,subsub_file),"r")
            results = [float(x) for x in f.read().strip().split("    ")[2].split("   ")]
            f.close()
            return results[0], results[1]
    return None, None

#calculates error from john's samplings files
def bootstrap_error_by_file( datafile ):
    try:
        f=open(datafile)
        #skip first line -> N_B
        line = f.readline()
        data = []

        line = f.readline()
        while line:
            data.append( float(line.strip().split(" ")[1]) )
            line = f.readline()
        
        result = bootstrap_error_by_array( data )
        f.close()
    except:
        result = -1
        
    return result
        
#calculates boostrap error by samplings array of form [average, sample, sample...]
def bootstrap_error_by_array( array ):
    average = array[0]
    samples = np.array(array[1:])
    sigmaB = np.sum( (samples-average)*(samples-average) ) 
    sigmaB = np.sqrt(sigmaB/(len(samples)-1))
    return sigmaB

#calculates jackkife error by samplings array of form [average, sample, sample...]
def jackknife_error_by_array( array ):
    average = array[0]
    samples = np.array(array[1:])
    sigmaJ = np.sum( (samples-average)*(samples-average) ) 
    sigmaJ = np.sqrt( (len(samples)-1)/len(samples) * sigmaJ )
    return sigmaJ
                
#retireves energy and error from john's data file (john uses bootstrap only)
def select_val_dat(dataset, mom, irrep, energy, spec_type):
    sub_file = "dsq"+mom.replace("PSQ","")+"_"+irrep
#     print(sub_file)
    if sub_file in os.listdir(dataset):
        if spec_type=="mom":
            subsub_file = "level_"+energy.replace("q2cm_","").replace("_ref","").replace("elab_","")+"_pcm_over_m_sq_smpls.dat"
        elif spec_type=="energy":
            subsub_file = "level_"+energy.replace("ecm_","").replace("_ref","").replace("elab_","")+"_ecm_over_m_smpls.dat"
        else:
            print("bad spectrum type")
            return None, None
        print(subsub_file)
        if subsub_file in os.listdir(os.path.join(dataset,sub_file)):
            f=open(os.path.join(dataset,sub_file,subsub_file),"r")
            line = f.readline()
#             print(line)
            results = [float(x) for x in f.readline().strip().split(" ")]
            f.close()
            return results[1], bootstrap_error_by_file(os.path.join(dataset,sub_file,subsub_file))
    return None, None

#based on filetype, chooses how to unpack the files. Also retrieves the rest masses
def unpack_file( filename, spectrum_type):
    if os.path.isfile(filename) and filename.endswith(".csv"):
        dataset1 = retrieve_sigmond_script_data(filename)
    elif os.path.isfile(filename) and filename.endswith(".hdf5"):
        dataset1 = retrieve_sigmond_script_data_hdf5(filename)
    elif os.path.isfile(filename) and filename.endswith(".dat"):
        dataset1 = retrieve_sigmond_script_data_dat(filename, spectrum_type)
    elif os.path.isfile(filename) and filename.endswith(".h5"):
        dataset1 = retrieve_barbara_data(filename)
    elif os.path.isdir(filename):
        dataset1 = pd.DataFrame()
    else:
        print("Bad filename:",filename)
        sys.exit()
    print(dataset1)
    return dataset1

#retrieve rest mass from dataset
def find_rest_mass( dataset, rest_mass_name ):
    tagname = rest_mass_name+'(0)_ref'
    return select_val(dataset, None, None, tagname)


"""find_tmin_spectrum_files(channel_name)
        channel_name - a key in the dict given in the spectrum config file
        
        for the given channel name, uses the directory path given in the spectrum
        config file to create a list of the filepaths to all of the tmin plots 
        in that given directory and label them based on their channel and basis
        This function relies of the directory tree structure produced by sigmond 
        scripts https://github.com/andrewhanlon/sigmond_scripts.
        
        if there is 'omit' tag under the channel name it will omit any basis name
        or filename in the omissions list from the plotting lists.
"""
def find_tmin_spectrum_files(channel_name): 
    this_channel = {}
        
    if 'graph_type' in channel_name.keys():
        graph_type = channel_name['graph_type']
    else:
        graph_type = 'E'
        
    if type(channel_name)==dict:
        this_directory = channel_name['dir']
        if 'omit' in channel_name.keys():
            omissions_list = channel_name['omit']
        else:
            omissions_list = []
    else:
        print("ERROR: Incorrect spectrum config channel")
        return {}
    for file in os.listdir(this_directory):
        if (file!='single_hadrons') and (file not in omissions_list):
            this_channel[file] = {}
    for i in this_channel.keys():
        if i not in omissions_list:
            this_channel[i]["singleR"] = []
            this_channel[i]["doubleR"] = []
            this_channel[i]["single"] = []
            this_channel[i]["double"] = []
            for level in range(0,20):
                for file in os.listdir(os.path.join(this_directory,i)):
                    if file not in omissions_list:
                        if graph_type=='E':
                            if file.endswith(f"R_0_{level}.agr"):
                                this_channel[i]["singleR"].append(os.path.join(this_directory,i,file))
                            elif file.endswith(f"R_4_{level}.agr"):
                                this_channel[i]["doubleR"].append(os.path.join(this_directory,i,file))
                            elif file.endswith(f"_0_{level}.agr"):
                                this_channel[i]["single"].append(os.path.join(this_directory,i,file))
                            elif file.endswith(f"_4_{level}.agr"):
                                this_channel[i]["double"].append(os.path.join(this_directory,i,file))
                        elif graph_type=='dE':
                            if file.endswith(f"R_0_D_{level}.agr"):
                                this_channel[i]["singleR"].append(os.path.join(this_directory,i,file))
                            elif file.endswith(f"R_4_D_{level}.agr"):
                                this_channel[i]["doubleR"].append(os.path.join(this_directory,i,file))
                            elif file.endswith(f"_0_D_{level}.agr"):
                                this_channel[i]["single"].append(os.path.join(this_directory,i,file))
                            elif file.endswith(f"_4_D_{level}.agr"):
                                this_channel[i]["double"].append(os.path.join(this_directory,i,file))
                        else:
                            print("invalid graph_type in find_tmin_spectrum_files(). use 'E' or 'dE'")
    return this_channel

"""find_tmin_spectrum_files_python(channel_name)
        channel_name - a dict from the yaml config file that contains the key
            'dir' with a directory path
        
        returns a dict with levels [basis][level][fit_type], with the xmgrace tmin 
        filenames from the directorys organized by the categories. if basis, fit_type, 
        or an individual name of a file are contained in omissions list, then they will
        not be included in the output dict.
"""
def find_tmin_spectrum_files_python(channel_name):
    this_channel = {}
    this_directory = channel_name['dir']
    if 'omit' in channel_name.keys():
        omissions_list = channel_name['omit']
    else:
        omissions_list = []
        
    if 'graph_type' in channel_name.keys():
        graph_type = channel_name['graph_type']
    else:
        graph_type = 'E'
        
    if graph_type=='E':
        tag_size = -5
        tmin_tags = settings.tmin_file_tags
    elif graph_type=='dE':
        tag_size = -7
        tmin_tags = settings.dtmin_file_tags
    else:
        print("invalid graph_type in find_tmin_spectrum_files(). use 'E' or 'dE'")
        
    for basis in os.listdir(this_directory):
        if (basis!='single_hadrons') and (basis not in omissions_list):
            this_channel[basis] = {}
            for level in range(0,channel_name['max_level']):
                this_channel[basis][level] = {}
                for file in os.listdir(os.path.join(this_directory,basis)):
                    if file not in omissions_list:
                        if file.endswith(f"_{level}.agr"):
                            file2 = file.replace(f"_{level}.agr",".agr")
                            if file2[tag_size:] in tmin_tags:
                                fit_type = tmin_tags[file2[tag_size:]]
                                if file2[:tag_size].endswith("R_"):
                                    fit_type = fit_type.replace("fit","ratio fit")
                                if fit_type not in omissions_list:
                                    this_channel[basis][level][fit_type] = os.path.join(this_directory,basis,file)
    return this_channel

"""zip_channel( name,this_file_list, this_type = ".svg", sub_dir = None) -> move to utils
        name - name of the ouput zip file (name.zip)
        this_file_list - dict where the keys refer to the filestubs
            of all of the svg files to be zipped (default: key.svg) 
        this_type - the file type to attach to the end of the key names
        sub_dir - if specified, looks in this directory name for files rather than the local
            directory.
            
        uses this_file_list.keys() to find all the files to be zipped and puts
        them in name.zip
"""
def zip_channel( name,this_file_list, this_type = ".svg", sub_dir = None):
    file = name+".zip"
    this_zip = ZipFile(file,'w')
    
    if type(this_file_list)==dict:
        this_list = this_file_list.keys()
    else:
        this_list = this_file_list
        
    for basis in this_list:
        if sub_dir:
            this_zip.write(os.path.join(sub_dir,basis)+this_type)
        else:
            this_zip.write(basis+this_type)
    this_zip.close()
