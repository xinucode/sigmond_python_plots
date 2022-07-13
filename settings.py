single_had_key = 'single_hadrons'

#designates colors and markers for plotting
colors = ['green','blue','orange','purple','black','red','brown','gray']
markers = ['o','s','D','v','^','*','x','+']
zigzag_shifts = [-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1]

#list of files that are jackknife samples rather than bootstrap
jackknife_sampling_methods = ['isoquartet_nonstrange_fermionic\\qsqr_samplings_isoquartet_nonstrange_fermionic_colin_rebin20_jackknife.hdf5','isoquartet_nonstrange_fermionic\\qsqr_samplings_isoquartet_nonstrange_fermionic_colin_rebin20_jackknife_6-12.hdf5']

#Gives a numerical order for the irreps so that I can control how they show up in the plots
alphabetical = {'A1g':1.01107,'F1':0.0601,'F2':0.0602,'G': 0.07,'G1': 0.0701,'G1g': 0.07107,'G1u': 0.07121,'G2': 0.0702, 'Hg': 0.0807, 'Hu': 0.0821, 'T1g':0.20107}

#latex format of all irreps
latex_format = {'A1g':r"I=1 $A_{1g}$",'F1':r"$F_1$",'F2':r"$F_2$",'G': r'$G$','G1': r"$G_1$",'G1g': r"$G_{1g}$",'G1u': r"$G_{1u}$",'G2': r"$G_2$", 'Hg': r"$H_g$", 'Hu': r"$H_u$", 'N': r'$N$', 'pi':r"$\pi$", 'mpi':r"$m_\pi$",'mpi2':r"$m_\pi^2$",'T1g':r"I=0 $T_{1g}$"}