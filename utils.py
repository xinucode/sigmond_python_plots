import pandas as pd
import numpy as np
import h5py
import os, sys

single_had_key = 'single_hadrons'

#designates colors and markers for plotting
colors = ['green','orange','blue','purple','black']
markers = ['o','s','D','v','^']

#list of files that are jackknife samples rather than bootstrap
jackknife_sampling_methods = ['isoquartet_nonstrange_fermionic\\qsqr_samplings_isoquartet_nonstrange_fermionic_colin_rebin20_jackknife.hdf5','isoquartet_nonstrange_fermionic\\qsqr_samplings_isoquartet_nonstrange_fermionic_colin_rebin20_jackknife_6-12.hdf5']

#Gives a numerical order for the irreps so that I can control how they show up in the plots
alphabetical = {'F1':0.0601,'F2':0.0602,'G': 0.07,'G1': 0.0701,'G1g': 0.07107,'G1u': 0.07121,'G2': 0.0702, 'Hg': 0.0807, 'Hu': 0.0821}
#latex format of all irreps
latex_format = {'F1':r"$F_1$",'F2':r"$F_2$",'G': r'$G$','G1': r"$G_1$",'G1g': r"$G_{1g}$",'G1u': r"$G_{1u}$",'G2': r"$G_2$", 'Hg': r"$H_g$", 'Hu': r"$H_u$", 'N': r'$m_N$', 'pi':r"$m_\pi$"}

#uses the irrep and momentum in form G1u(0) and orders based on momentum number + irrep's "alphabetical" value
def sort_by_mom(irrep):
    parts = irrep.split("(")
    return float(parts[1][0])+float(alphabetical[parts[0]])

#goes through data file and inserts info into the dataframe (pandas)
def split_obs_col(dataset):
    obs_names = dataset['obs']
    obs_mom = []
    obs_irrep = []
    obs_level = []
    for i, (obs) in enumerate(obs_names):
        obs_list = obs.split('/')
        if obs_list[0] == single_had_key:
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
    if file in jackknife_sampling_methods:
        print("jackknife success")
    for mom in q2s.keys():
        for basis in q2s[mom].keys():
            try:
                for value in q2s[mom+'/'+basis].keys():
                    if file in jackknife_sampling_methods:
                        df2 = pd.DataFrame([[mom+'/'+basis+'/'+value, q2s[mom+'/'+basis+'/'+value][0], jackknife_error_by_array(q2s[mom+'/'+basis+'/'+value])]], columns=['obs', 'val', 'err'])
                    else:
                        df2 = pd.DataFrame([[mom+'/'+basis+'/'+value, q2s[mom+'/'+basis+'/'+value][0], bootstrap_error_by_array(q2s[mom+'/'+basis+'/'+value])]], columns=['obs', 'val', 'err'])
                    dataset = dataset.append(df2)
            except AttributeError:
                pass
    split_obs_col(dataset)
    q2s.close()
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
def unpack_file( filename ):
    if os.path.isfile(filename) and filename.endswith(".csv"):
        dataset1 = retrieve_sigmond_script_data(filename)
    elif os.path.isfile(filename) and filename.endswith(".hdf5"):
        dataset1 = retrieve_sigmond_script_data_hdf5(filename)

    elif os.path.isdir(filename):
        dataset1 = pd.DataFrame()
    else:
        print("Bad filename:",filename)
        sys.exit()
    return dataset1

#retrieve rest mass from dataset
def find_rest_mass( dataset, rest_mass_name ):
    tagname = rest_mass_name+'(0)_ref'
    return select_val(dataset, None, None, tagname)