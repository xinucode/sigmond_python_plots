import pandas as pd
import numpy as np
import h5py
import os, sys
from zipfile import ZipFile
import gvar
from numba import njit

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

@njit(parallel=True)
def effenergy(t, C):
    return t[:-1]+0.5*(t[1]-t[0]), np.log( C[:-1]/C[1:] )/(t[1]-t[0])

@njit(parallel=True)
def effenergy2(t, C):
    return t[1:-1], np.arctanh( (C[:-2]-C[2:])/(C[:-2]+C[2:]) ) #alternate function

@njit(parallel=True)
def effenergy3(t, C):
    return t[1:-1], -2.0*( C[2:]-2.0*C[1:-1]+C[:-2] )/( C[2:]-C[:-2] ) #-C"/C'

# def effenergy4(t, C):
#     return t[1:-1],  #-C/I[C]

# def effenergy2(t, C):

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

# @njit(parallel=True)
def shift_levels( indexes, levels, vals, errors, shifted_array, index=0 ):
    if len(indexes):
        # if not shifted_array.any():
            # shifted_array = np.array([0.0]*len(indexes))
            
        if index==len(indexes):
            return shifted_array
            
        these_indexes = indexes[index:]
        this_remaining_irrep_index = these_indexes[0]
        this_remaining_irrep_indexes = np.where( these_indexes==this_remaining_irrep_index )[0]+index
        
        if len(this_remaining_irrep_indexes)==1:
            return shift_levels( indexes, levels, vals, errors, shifted_array, index=index+1 )
        else:
            shift = 1
            this_val_upper = vals[index]+errors[index]
            this_val_lower = vals[index]-errors[index]
            this_cluster = [index]
            for i in this_remaining_irrep_indexes[1:]:
                overlap = False  
                compare_upper = vals[i]+errors[i]
                compare_lower = vals[i]-errors[i]
                if shifted_array[i]!=0.0:
                    continue
                if compare_lower<=this_val_lower and compare_upper>=this_val_upper: #set new bounds on this_value
                    overlap=True
                elif compare_lower>=this_val_lower and compare_upper<=this_val_upper:
                    overlap=True
                elif compare_lower<=this_val_upper and compare_lower>=this_val_lower and compare_upper>=this_val_upper:
                    overlap=True
                elif compare_lower<=this_val_lower and compare_upper<=this_val_upper and compare_upper>=this_val_lower:
                    overlap=True
                if overlap:
                    shifted_array[i] = shift
                    this_cluster.append(i)
                    if shift>0:
                        shift = -shift
                    else:
                        shift = abs(shift)+1
            
            if len(this_cluster)%2==0 and len(this_cluster):
                for i in this_cluster:
                    shifted_array[i]-=0.5
            
            new_index = np.where(shifted_array[index+1:]==0.0)[0]
            if new_index.any():
                return shift_levels( indexes, levels, vals, errors, shifted_array, index=new_index[0]+index+1 )
            else:
                return shifted_array
                        
    return 0.0
    

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
        return None, None, None
    elif len(reduced_dataset["val"])==1:
        values = np.array(reduced_dataset["val"][:])
        errors = np.array(reduced_dataset["err"][:])
        if "zmags" in reduced_dataset:
            zmags = np.array(reduced_dataset["zmags"][:])
        else:
            zmags = [None]
        return values[0],errors[0],zmags[0]
    else:
        print(f"Duplicate Values: {mom}, {irrep}, {energy}")
        values = np.array(reduced_dataset["val"][:])
        errors = np.array(reduced_dataset["err"][:])
        if "zmags" in reduced_dataset:
            zmags = np.array(reduced_dataset["zmags"][:])
        else:
            zmags = [None]
        return values[0],errors[0],zmags[0]
    
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
def retrieve_sigmond_script_data_hdf5(file, channel = None):
    q2s = h5py.File(file,'r')
    dataset = pd.DataFrame(columns=['obs', 'val', 'err'])
    if file in settings.jackknife_sampling_methods:
        print("jackknife success")
        
    #get single hadron data
    mom = 'single_hadrons'
    for value in q2s[mom].keys():
        try:
            if file in settings.jackknife_sampling_methods:
                df2 = pd.DataFrame([[mom+'/'+value, q2s[mom+'/'+value][0], jackknife_error_by_array(np.array(q2s[mom+'/'+value]))]], columns=['obs', 'val', 'err'])
            else:
                df2 = pd.DataFrame([[mom+'/'+value, q2s[mom+'/'+value][0], bootstrap_error_by_array(np.array(q2s[mom+'/'+value]))]], columns=['obs', 'val', 'err'])
            dataset = pd.concat([dataset,df2], ignore_index=True)
        except AttributeError:
                pass
    
    if channel is not None:
        channel_group = q2s[channel]
    else:
        channel_group = q2s
    
    
    for mom in channel_group.keys():
        if mom=='Info':
            continue
        elif mom=='single_hadrons':
            continue
        else:
            for basis in channel_group[mom].keys():
                for value in channel_group[mom+'/'+basis].keys():
                    try:
                        if file in settings.jackknife_sampling_methods:
                            df2 = pd.DataFrame([[mom+'/'+basis+'/'+value, channel_group[mom+'/'+basis+'/'+value][0], jackknife_error_by_array(np.array(channel_group[mom+'/'+basis+'/'+value]))]], columns=['obs', 'val', 'err'])
                        else:
                            df2 = pd.DataFrame([[mom+'/'+basis+'/'+value, channel_group[mom+'/'+basis+'/'+value][0], bootstrap_error_by_array(np.array(channel_group[mom+'/'+basis+'/'+value]))]], columns=['obs', 'val', 'err'])
                        dataset = pd.concat([dataset,df2], ignore_index=True)
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
    
    # dataframe_input = pd.read_csv(file,sep=" ",header=0,names=["obs-irrep","mom","qsqr","a","b","mass"])
    dataframe_input = pd.read_csv(file,sep=" ",header=0,names=["mom","obs-irrep","epredlab", "qsqr", "edatalab", "edatacm"])
    irrep = np.array(dataframe_input["obs-irrep"])+np.array(dataframe_input["mom"],dtype="str")
    obs_mom = ["PSQ"+str(round(float(x))) for x in list(dataframe_input["mom"])]
    levels = np.zeros( len(dataframe_input["obs-irrep"]) )
    # mass = dataframe_input["mass"][0]
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
        data = [float(x) for x in dataframe_input["qsqr"]]
        ecm = np.array(data)
    else:
        print("Code this (utils.retrieve_sigmond_script_data_dat) if possible.")
        return None
    dataframe_input.insert(3, "val", ecm)
    # errs = dataframe_input["a"] #np.zeros( len(dataframe_input["obs-irrep"]) )
    errs = np.zeros( len(dataframe_input["obs-irrep"]) )
    dataframe_input.insert(4, "err", errs)
    dataset = pd.DataFrame(columns=["obs-mom","obs-irrep","obs-level", 'val', 'err'],data=dataframe_input[["obs-mom","obs-irrep","obs-level","val","err"]])
    return dataset

def retrieve_barbara_data(file):
    q2s = h5py.File(file,'r')
    dataset = pd.DataFrame(columns=['obs', 'val', 'err','ops'])
#     print( q2s.keys() )
    for key in q2s.keys():
#         print( q2s[key].keys() )
#         print( q2s[key]["1exp"]["t0_12"].keys() )
        obs_mom = key.split("_")[0]
        obs_irrep = key.split("_")[1]
        for level in q2s[key]["1exp"]["t0_7"]["Tmin"]["Correlated"]["Mean"].keys():
            index = int(level.split("_")[1])
            data = np.array(q2s[key]["1exp"]["t0_7"]["Tmin"]["Correlated"]["Mean"][level])
#             data = np.array(q2s[key]["1exp"]["t0_7"]["Tmin"]["Correlated"]["Resample"][level])
#             chisqr = data[-1]
#             tmin = data[0]
#             tmax = data[1]
#             dof = tmax-tmin-2
#             chisqr_dof = chisqr/dof
#             for i in range(len(chisqr_dof)):
#                 if chisqr_dof[i]<1.0:
#                     chisqr_dof[i] = 1.0/chisqr_dof[i]
#             minimum_index = np.where( chisqr_dof==min(chisqr_dof) )[0][0]
#             value = data[3][minimum_index]
#             error = data[4][minimum_index]
            value = data[3][0]
            error = data[4][0]
            operators = q2s[key]["Single_hadron_corrs"][index]
            dataset.loc[len(dataset)] = [f"{obs_mom}/{obs_irrep}/dElab_{index}",value,error,operators]
#             print(f"{obs_mom}/{obs_irrep}/d_ecm_{index}",value,error,operators)
    q2s.close()
    split_obs_col(dataset)
#     print(dataset)
#     return None
    return dataset
    
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
        
        result = bootstrap_error_by_array( np.array(data) )
        f.close()
    except:
        result = -1
        
    return result
        
#calculates boostrap error by samplings array of form [average, sample, sample...]
#@jit
@njit(parallel=True)
def bootstrap_error_by_array( array ):
    average = array[0]
    samples = array[1:]
    sigmaB = np.sum( (samples-average)*(samples-average) ) 
    sigmaB = np.sqrt(sigmaB/(len(samples)-1))
    return sigmaB

#calculates jackkife error by samplings array of form [average, sample, sample...]
@njit(parallel=True)
def jackknife_error_by_array( array ):
    average = array[0]
    samples = array[1:]
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
def unpack_file( filename, spectrum_type, root=None):
    if os.path.isfile(filename) and filename.endswith(".csv"):
        dataset1 = retrieve_sigmond_script_data(filename)
    elif os.path.isfile(filename) and filename.endswith(".hdf5"):
        dataset1 = retrieve_sigmond_script_data_hdf5(filename, root)
    elif os.path.isfile(filename) and filename.endswith(".dat"):
        dataset1 = retrieve_sigmond_script_data_dat(filename, spectrum_type)
    elif os.path.isfile(filename) and filename.endswith(".h5") and os.path.basename(filename).startswith("andre"):
        dataset1 = retrieve_andre_data(filename)
    elif os.path.isfile(filename) and filename.endswith(".h5"):
        dataset1 = retrieve_barbara_data(filename)
    # elif os.path.isfile(filename) and filename.endswith(".data"):
        # dataset1 = retrieve_andre_data(filename)
    elif os.path.isdir(filename):
        dataset1 = pd.DataFrame()
    else:
        print("Bad filename:",filename)
        sys.exit()
        
    return dataset1

def retrieve_andre_data(file):
    # dataset0 = pd.read_csv(file, header=None )
    # dataset = pd.DataFrame(columns=['obs', 'val', 'err'])
    # for i,row in dataset0.iterrows():
        # obs_mom = row[0]
        # obs_irrep = row[1]
        # index = row[2]
        # elab = row[8]
        # dataset.loc[len(dataset)] = [f"PSQ{obs_mom}/{obs_irrep}/elab_{index}",gvar.gvar(elab).mean,gvar.gvar(elab).sdev]
        # ecm = row[9]
        # dataset.loc[len(dataset)] = [f"PSQ{obs_mom}/{obs_irrep}/ecm_{index}",gvar.gvar(ecm).mean,gvar.gvar(ecm).sdev]
    # split_obs_col(dataset)
    # return dataset
    
    
    dataset = pd.DataFrame(columns=['obs', 'val', 'err'])
    f = h5py.File(file)
    for key in f.keys():
        irrep, level, mom = key.split("_")
        
        data = f[key]['dEnn'][()]
        dataset.loc[len(dataset)] = [f"{mom.upper()}/{irrep}/dElab_{level}",data[0],bootstrap_error_by_array(np.array(data))]
        data = f[key]['dEnn'][()]+f[key]['N1'][()]+f[key]['N2'][()]
        dataset.loc[len(dataset)] = [f"{mom.upper()}/{irrep}/elab_{level}",data[0],bootstrap_error_by_array(np.array(data))]
        data = (f[key]['dEnn'][()]+f[key]['N1'][()]+f[key]['N2'][()])/f[key]['mN'][()]
        dataset.loc[len(dataset)] = [f"{mom.upper()}/{irrep}/elab_{level}_ref",data[0],bootstrap_error_by_array(np.array(data))]
        
    f.close()
    
    split_obs_col(dataset)
    return dataset

#retrieve rest mass from dataset
def find_rest_mass( dataset, rest_mass_name, remove_ref ):
    if remove_ref:
        tagname = rest_mass_name+'(0)'
    else:
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
    
def collectDiagonalRealCorrelatorAtTime(data_object, corr_str, time, tag='rotated_correlators'):
    corr_key = f"<MCObservable><CorrT>GI{{{corr_str}}} GI{{{corr_str}}} time={time} HermMat<|CorrT><Arg>Re<|Arg><|MCObservable>"
    if corr_key in data_object[tag]['Values'].keys():
        return data_object[tag]['Values'][corr_key][()]
    corr_key = f"<MCObservable><CorrT>BL{{{corr_str}}} BL{{{corr_str}}} time={time} HermMat<|CorrT><Arg>Re<|Arg><|MCObservable>"
    return data_object[tag]['Values'][corr_key][()]

def collectCorrEstimates(data_object, corr_str, tag='rotated_correlators'):
    t = []
    values = []
    errs = []
    for i in range(64): #maybe set up global max
        try:
            if type(corr_str)==list: #simple ratio correlator
                this_numt = collectDiagonalRealCorrelatorAtTime(data_object,corr_str[0],i, tag)
                this_dett = collectDiagonalRealCorrelatorAtTime(data_object,corr_str[1],i, tag)
                this_corrt = this_numt/this_dett
            else:
                this_corrt = collectDiagonalRealCorrelatorAtTime(data_object,corr_str,i, tag)
            t.append(i)
            values.append(this_corrt[0])
            errs.append(bootstrap_error_by_array(this_corrt))
        except:
            pass
        
    return t, values, errs

def collectEnergyEstimates(data_object, corr_str, tag='rotated_correlators',func=1):
    t = []
    values = []
    errs = []
    for i in range(func-1,63-func):
        try:
            if type(corr_str)==list: #simple ratio correlator
                this_numt = collectDiagonalRealCorrelatorAtTime(data_object,corr_str[0],i-1, tag)
                this_dett = collectDiagonalRealCorrelatorAtTime(data_object,corr_str[1],i-1, tag)
                this_corrt0 = this_numt/this_dett
                
                this_numt = collectDiagonalRealCorrelatorAtTime(data_object,corr_str[0],i, tag)
                this_dett = collectDiagonalRealCorrelatorAtTime(data_object,corr_str[1],i, tag)
                this_corrt = this_numt/this_dett
                
                this_numt = collectDiagonalRealCorrelatorAtTime(data_object,corr_str[0],i+1, tag)
                this_dett = collectDiagonalRealCorrelatorAtTime(data_object,corr_str[1],i+1, tag)
                this_corrt2 = this_numt/this_dett
            else:
                this_corrt0 = collectDiagonalRealCorrelatorAtTime(data_object,corr_str,i-1, tag)
                this_corrt = collectDiagonalRealCorrelatorAtTime(data_object,corr_str,i, tag)
                this_corrt2 = collectDiagonalRealCorrelatorAtTime(data_object,corr_str,i+1, tag)
                
            if func==1:
                new_t, this_effE = effenergy(np.array([i,i+1]),np.array([this_corrt,this_corrt2]))
            elif func==2:
                new_t, this_effE = effenergy2(np.array([i-1,i,i+1]),np.array([this_corrt0,this_corrt,this_corrt2]))
            elif func==3:
                new_t, this_effE = effenergy3(np.array([i-1,i,i+1]),np.array([this_corrt0,this_corrt,this_corrt2]))
            if not np.isnan([this_effE[0][0]])[0] and not np.isnan([bootstrap_error_by_array(this_effE[0])])[0]:
                t.append(new_t[0])
                values.append(this_effE[0][0])
                errs.append(bootstrap_error_by_array(this_effE[0]))
        except Exception as error:
            pass #print(error)
        
    return t, values, errs

@njit(parallel=True)
def multi_exp_func(t, E0, E1, A0, A1, A2, A3, A4, A5, n2, n3, n4, n5):
    return A0*np.exp(-E0*t)*(1.0 + A1*np.exp(-E1*E1*t) + A2*np.exp(-n2*E1*E1*t)
                + A3*np.exp(-n3*E1*E1*t) + A4*np.exp(-n4*E1*E1*t) + A5*np.exp(-n5*E1*E1*t))

@njit(parallel=True)
def multi_exp_func_dt(t, E0, E1, A0, A1, A2, A3, A4, A5, n2, n3, n4, n5):
    return -E0*A0*np.exp(-E0*t)*(1.0 +A1*np.exp(-E1*E1*t) + A2*np.exp(-n2*E1*E1*t)
                + A3*np.exp(-n3*E1*E1*t) + A4*np.exp(-n4*E1*E1*t) + A5*np.exp(-n5*E1*E1*t)) - E1*E1*A0*np.exp(-E0*t)*(A1*np.exp(-E1*E1*t) +n2*A2*np.exp(-n2*E1*E1*t)
                +n3*A3*np.exp(-n3*E1*E1*t) +n4*A4*np.exp(-n4*E1*E1*t) + n5*A5*np.exp(-n5*E1*E1*t))

@njit(parallel=True)
def multi_exp_func_eff(t, E0, E1, A0, A1, A2, A3, A4, A5, n2, n3, n4, n5):
    return np.abs(multi_exp_func_dt(t, E0, E1, A0, A1, A2, A3, A4, A5, n2, n3, n4, n5)/multi_exp_func(t, E0, E1, A0, A1, A2, A3, A4, A5, n2, n3, n4, n5))

@njit(parallel=True)
def multi_exp_func(t, E0, E1, E2, E3, E4, A0, A1, A2, A3, A4):
    return A0*np.exp(-E0*t)*(1.0 + A1*np.exp(-E1*t) + A2*np.exp(-E2*t) + A3*np.exp(-E3*t) + A4*np.exp(-E4*t))

@njit(parallel=True)
def multi_exp_func_dt(t, E0, E1, E2, E3, E4, A0, A1, A2, A3, A4):
    return -E0*A0*np.exp(-E0*t)*(1.0 + A1*np.exp(-E1*t) + A2*np.exp(-E2*t) + A3*np.exp(-E3*t) + A4*np.exp(-E4*t)) + A0*np.exp(-E0*t)*(-E1*A1*np.exp(-E1*t) + -E2*A2*np.exp(-E2*t) + -E3*A3*np.exp(-E3*t) + -E4*A4*np.exp(-E4*t))

@njit(parallel=True)
def multi_exp_func_eff(t, E0, E1, E2, E3, E4, A0, A1, A2, A3, A4):
    return np.abs(multi_exp_func_dt(t, E0, E1, E2, E3, E4, A0, A1, A2, A3, A4)/multi_exp_func(t, E0, E1, E2, E3, E4, A0, A1, A2, A3, A4))