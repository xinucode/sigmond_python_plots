single_had_key = 'single_hadrons'

#designates colors and markers for plotting
colors = ['green','blue','orange','purple','black','red','brown','gray']
markers = ['o','s','D','v','^','*','x','+']
zigzag_shifts = [-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1]

#list of files that are jackknife samples rather than bootstrap
jackknife_sampling_methods = ['isoquartet_nonstrange_fermionic\\qsqr_samplings_isoquartet_nonstrange_fermionic_colin_rebin20_jackknife.hdf5','isoquartet_nonstrange_fermionic\\qsqr_samplings_isoquartet_nonstrange_fermionic_colin_rebin20_jackknife_6-12.hdf5']

#Gives a numerical order for the irreps so that I can control how they show up in the plots
alphabetical = {'A1g':0.01107,'F1':0.0601,'F2':0.0602,'G': 0.07,'G1': 0.0701,'G1g': 0.07107,'G1u': 0.07121,'G2': 0.0702, 'Hg': 0.0807, 'Hu': 0.0821, 'T1g':0.20107,'A1u':0.01121,'A2u':0.01221,'Eg':0.0507,'Eu':0.0521,'T1u':0.20121,'T2g':0.20207,'T2u':0.20221}

#latex format of all irreps and particles
latex_format = {'A1g':r"$A_{1g}$", 'F1':r"$F_1$", 'F2':r"$F_2$", 'G': r'$G$', 'G1':r"$G_1$", 'G1g':r"$G_{1g}$", 'G1u':r"$G_{1u}$", 'G2':r"$G_2$", 'Hg':r"$H_g$", 'Hu':r"$H_u$", 'N':r'$N$', 'pi':r"$\pi$", 'T1g':r"$T_{1g}$", 'A1u':r"$A_{1u}$", 'A2u':r"$A_{2u}$", 'Eg':r"$E_{g}$", 'Eu':r"$E_{u}$", 'T1u':r"$T_{1u}$", 'T2g':r"$T_{2g}$", 'T2u':r"$T_{2u}$"}

#xmgrace label format for particles #insert link here to help
xmgrace_format = {'N': 'N', 'pi': '\\xp'}

#sigmond scripts tmin plot end tags
tmin_file_tags = {'0.agr':'single','4.agr':'double','8.agr':'geometric'}
