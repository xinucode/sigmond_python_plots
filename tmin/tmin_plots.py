import numpy as np
import copy
import os, sys
import yaml
import pandas as pd
import matplotlib.pyplot as plt
import argparse

sys.path.insert(0,'../')
import xmgrace_parser
import settings
import utils
plt.style.use('../spectrum.mplstyle')

#xmgrace info
symbol_pattern_tag = "symbol_fill_pattern"
symbol_color_tag = "symbol_color"
line_color_tag = "errorbar_color"
legend_tag = "legend"
legend_bool_tag = "legend_on_off"
xydydy = "xydydy"
xy = "xy"
symbol_pattern_empty = ['1'] #['0']
symbol_pattern_filled = ['1']
red_index = 2
blue_index = 4
green_index = 15
purple_index = 8
yellow_index = 11
axis_tag = "axis_label"
square_index = 2
circle_index = 1
diamond_index = 3
triangle_index = 4


###################### xmgrace functions ##############################
"""rest_mass_tmin_plot(filename_1,filename_2,label,legend_labels)
       filename_1 - filepath of single exponential fit tmin plot
       filename_2 - filepath of double exponential fit tmin plot
       label - set y axis label
       legend_labels - how the data is labelled in legend
           options: 
               "regular" - labels data "{X}-exp fit"
               "ratio" - labels data "{X}-exp ratio fit"
               other - no label
       
       combines the tmin plots of single and double exponential fits
       for one level, applies relevant settings for rest masses,
       and returns the resulting xmgrace_parser.AgrFile object. 
"""
def rest_mass_tmin_plot(filename_1,filename_2,label,good_labels = [False,False,False],legend_labels="regular", filename_3 = None):
    single_exp_plot = xmgrace_parser.AgrFile(filename_1,True) 
    double_exp_plot = xmgrace_parser.AgrFile(filename_2,True)
    if filename_3:
        geometric_plot = xmgrace_parser.AgrFile(filename_3,True)
    
    #find all fit data
    se_xydydy = get_data_by_type(single_exp_plot,xydydy)
    de_xydydy = get_data_by_type(double_exp_plot,xydydy)
    if filename_3:
        geom_xydydy = get_data_by_type(geometric_plot,xydydy)
    
    if legend_labels=="ratio":
        se_legend_label = "1-exp ratio fit"
        de_legend_label = "2-exp ratio fit"
        geom_legend_label = "geom ratio fit"
    elif legend_labels=="regular":
        se_legend_label = "1-exp fit"
        de_legend_label = "2-exp fit"
        geom_legend_label = "geometric fit"
    else:
        se_legend_label = ""
        de_legend_label = ""
        geom_legend_label = ""
    
    
    good_labels[0],labelled_data = format_data(single_exp_plot,se_xydydy,blue_index,circle_index,se_legend_label,settings.legend_setting,label)
    good_labels[1],labelled_data = format_data(double_exp_plot,de_xydydy,green_index,square_index,de_legend_label)
    if filename_3:
        good_labels[2],labelled_data = format_data(geometric_plot,geom_xydydy,purple_index,diamond_index,geom_legend_label)
    
    se_xy = get_data_by_type(single_exp_plot,xy)
    if se_xy:
#         for g,s in se_xy:
        single_exp_plot.get_set(se_xy[1][0],se_xy[1][1]).update_properties(legend="chosen t\\smin")
#             break
    
    #transfer double to single
    transfer_data( double_exp_plot, de_xydydy, single_exp_plot)
    if filename_3:
        transfer_data( geometric_plot, geom_xydydy, single_exp_plot)
        
    #reformat world
    world1 = ""
    world2 = ""
    world3 = ""
    if se_xydydy:
        world1 = single_exp_plot.graphs[0].get_properties(["world"])[0]
    if de_xydydy:
        world2 = double_exp_plot.graphs[0].get_properties(["world"])[0]
    if filename_3:
        if geom_xydydy:
            world3 = geometric_plot.graphs[0].get_properties(["world"])[0]
    if world1 and world2:    
        new_world = expand_world(world1,world2)
        if world3:
            new_world_str = str(new_world[0])+", "+str(new_world[1])+", "+str(new_world[2])+", "+str(new_world[3])
            new_world = expand_world(new_world_str,world3)
        minx = round(new_world[0])
        maxx = round(new_world[2])
        miny,maxy = stretch_estimate_lines(single_exp_plot,minx,maxx)

        ypadding = (maxy-miny)#*(2.0)
        new_world[1] = miny-ypadding
        new_world[3] = maxy+2.0*ypadding

        new_world_str = str(new_world[0])+", "+str(new_world[1])+", "+str(new_world[2])+", "+str(new_world[3])
        single_exp_plot.graphs[0].update_properties(world=new_world_str)
    
        #ticks
        spacing = round((new_world[3]-new_world[1])/settings.max_number_of_ticks,3)
        single_exp_plot.graphs[0].update_properties(yaxis_tick_major=spacing)
        
    return single_exp_plot

"""spectrum_tmin_format*(tmin_files)
        tmin_files - a list of files in dict format
            { basis1: {"single":[level1_filepath,level2_filepath...],
                "singleR":[level1_filepath,level2_filepath...],
                "doubleR":[level1_filepath,level2_filepath...]},
                basis2... }
            where the basisN refer to a basis within a channel. Each basis 
            has a dict with keys "single","singleR", and "doubleR". 
                "single" - single exponential non ratio fit
                "double" - double exp non ration fit
                "singleR" - single exponential ratio fit
                "doubleR" - double exponential ratio fit
            Of course a non ratio double exponential is possible though
            non included here as the double exponental is more fickle
            than the single exponential and the non ratio single is also
            fickle so decided not to include.
            
        Format1 and format2 take all of the tmin plots (all levels and fit types)
        in a basis and merges them into a combined spectrum tmin plot. 
        Format3 only combines the different fit types for each level. Every level
        will have a different plot.
        
        Format1 and format3 color based on the fit type. Format 2 colors based on
        the energy level.
"""
def spectrum_tmin_format1(tmin_files):
    plots_to_combine = {}
    spectrum_plots = {}
    for basis in tmin_files.keys():
        plots_to_combine[basis] = []
        print(basis)
        
        if tmin_files[basis]["singleR"]:
            for level,file in enumerate(tmin_files[basis]["singleR"]):
                plots_to_combine[basis].append(xmgrace_parser.AgrFile(file,True))
                this_xydydy = get_data_by_type(plots_to_combine[basis][-1],xydydy)
                if level ==0:
                    good_label, labelled_data = format_data(plots_to_combine[basis][-1],this_xydydy,blue_index,circle_index,"1-exp ratio fit",settings.legend_setting,settings.y_label)
                    labelled_level = level
                else:
                    if good_label:
                        format_data(plots_to_combine[basis][-1],this_xydydy,blue_index,circle_index)
                    elif labelled_data:
                        plots_to_combine[basis][labelled_level].get_set(labelled_data[0],labelled_data[1]).update_properties(legend="")
                        good_label, labelled_data = format_data(plots_to_combine[basis][-1],this_xydydy,blue_index,circle_index,"1-exp ratio fit")
                        labelled_level = level
                    else:
                        good_label, labelled_data = format_data(plots_to_combine[basis][-1],this_xydydy,blue_index,circle_index,"1-exp ratio fit")
                        labelled_level = level
                        
                
        if tmin_files[basis]["doubleR"]:
            for level,file in enumerate(tmin_files[basis]["doubleR"]):
                plots_to_combine[basis].append(xmgrace_parser.AgrFile(file,True))
                this_xydydy = get_data_by_type(plots_to_combine[basis][-1],xydydy)
                if level ==0:
                    good_label, labelled_data = format_data(plots_to_combine[basis][-1],this_xydydy,green_index,square_index,"2-exp ratio fit",settings.legend_setting,settings.y_label)
                    labelled_level = 0
                else:
                    if good_label:
                        format_data(plots_to_combine[basis][-1],this_xydydy,green_index,square_index)
                    elif labelled_data:
                        plots_to_combine[basis][labelled_level].get_set(labelled_data[0],labelled_data[1]).update_properties(legend="")
                        good_label, labelled_data = format_data(plots_to_combine[basis][-1],this_xydydy,green_index,square_index,"2-exp ratio fit")
                        labelled_level = level
                    else:
                        good_label, labelled_data = format_data(plots_to_combine[basis][-1],this_xydydy,green_index,square_index,"2-exp ratio fit")
                        labelled_level = level
                
        if tmin_files[basis]["single"]:
            for level,file in enumerate(tmin_files[basis]["single"]):
                plots_to_combine[basis].append(xmgrace_parser.AgrFile(file,True))
                this_xydydy = get_data_by_type(plots_to_combine[basis][-1],xydydy)
                if level ==0:
                    good_label, labelled_data = format_data(plots_to_combine[basis][-1],this_xydydy,purple_index,diamond_index,"1-exp fit",settings.legend_setting,settings.y_label)
                    labelled_level = level
                else:
                    if good_label:
                        format_data(plots_to_combine[basis][-1],this_xydydy,purple_index,diamond_index)
                    elif labelled_data:
                        plots_to_combine[basis][labelled_level].get_set(labelled_data[0],labelled_data[1]).update_properties(legend="")
                        good_label, labelled_data = format_data(plots_to_combine[basis][-1],this_xydydy,purple_index,diamond_index,"1-exp fit")
                        labelled_level = level
                    else:
                        good_label, labelled_data = format_data(plots_to_combine[basis][-1],this_xydydy,purple_index,diamond_index,"1-exp fit")
                        labelled_level = level
                
        if tmin_files[basis]["double"]:
            for level,file in enumerate(tmin_files[basis]["double"]):
                plots_to_combine[basis].append(xmgrace_parser.AgrFile(file,True))
                this_xydydy = get_data_by_type(plots_to_combine[basis][-1],xydydy)
                if level ==0:
                    good_label, labelled_data = format_data(plots_to_combine[basis][-1],this_xydydy,yellow_index,triangle_index,"2-exp fit",settings.legend_setting,settings.y_label)
                    labelled_level = level
                else:
                    if good_label:
                        format_data(plots_to_combine[basis][-1],this_xydydy,11,4)
                    elif labelled_data:
                        plots_to_combine[basis][labelled_level].get_set(labelled_data[0],labelled_data[1]).update_properties(legend="")
                        good_label, labelled_data = format_data(plots_to_combine[basis][-1],this_xydydy,11,4,"2-exp fit")
                        labelled_level = level
                    else:
                        good_label, labelled_data = format_data(plots_to_combine[basis][-1],this_xydydy,11,4,"2-exp fit")
                        labelled_level = level
                
                
        #combine all plots
        for i,plot in enumerate(plots_to_combine[basis]):
            if i==0:
                spectrum_plots[basis] = plot
            else:
                all_xydydy = get_data_by_type(plot,xydydy)
                transfer_data( copy.deepcopy(plots_to_combine[basis][i]), all_xydydy, spectrum_plots[basis])
                all_xy = get_data_by_type(plot,"xy")
                transfer_data( copy.deepcopy(plots_to_combine[basis][i]), all_xy, spectrum_plots[basis])
                
                world1 = plot.graphs[0].get_properties(["world"])[0]
                world2 = spectrum_plots[basis].graphs[0].get_properties(["world"])[0]
                new_world = expand_world(world1,world2)

                minx = round(new_world[0])
                maxx = round(new_world[2])
                miny,maxy = stretch_estimate_lines(spectrum_plots[basis],minx,maxx)

                ypadding = (maxy-miny)/(6.0)
                new_world[1] = miny-ypadding
                new_world[3] = maxy+2.0*ypadding
                new_world_str = str(new_world[0])+", "+str(new_world[1])+", "+str(new_world[2])+", "+str(new_world[3])
                spectrum_plots[basis].graphs[0].update_properties(world=new_world_str)
                
        #ticks
        spacing = round((new_world[3]-new_world[1])/settings.max_number_of_ticks,3)
        spectrum_plots[basis].graphs[0].update_properties(yaxis_tick_major=spacing)
                
        strip_unnecessary_data(spectrum_plots[basis], (new_world[3]-new_world[1]))
                
    return spectrum_plots

def spectrum_tmin_format2(tmin_files):
    plots_to_combine = {}
    spectrum_plots = {}
    for basis in tmin_files.keys():
        plots_to_combine[basis] = []
        print(basis)
        
        if tmin_files[basis]["singleR"]:
            for level,file in enumerate(tmin_files[basis]["singleR"]):
                plots_to_combine[basis].append(xmgrace_parser.AgrFile(file,True))
                this_xydydy = get_data_by_type(plots_to_combine[basis][-1],xydydy)
                if level==0:
                    good_label, labelled_data = format_data(plots_to_combine[basis][-1],this_xydydy,settings.spectrum_color_indices[level],circle_index,"1-exp ratio fit",settings.legend_setting,settings.y_label)
                    labelled_level = level
                else:
                    if good_label:
                        format_data(plots_to_combine[basis][-1],this_xydydy,settings.spectrum_color_indices[level],circle_index)
                    elif labelled_data:
                        plots_to_combine[basis][labelled_level].get_set(labelled_data[0],labelled_data[1]).update_properties(legend="")
                        good_label, labelled_data = format_data(plots_to_combine[basis][-1],this_xydydy,settings.spectrum_color_indices[level],circle_index,"1-exp ratio fit")
                        labelled_level = level
                    else:
                        good_label, labelled_data = format_data(plots_to_combine[basis][-1],this_xydydy,settings.spectrum_color_indices[level],circle_index,"1-exp ratio fit")
                        labelled_level = level
                        
                
        if tmin_files[basis]["doubleR"]:
            for level,file in enumerate(tmin_files[basis]["doubleR"]):
                plots_to_combine[basis].append(xmgrace_parser.AgrFile(file,True))
                this_xydydy = get_data_by_type(plots_to_combine[basis][-1],xydydy)
                if level ==0:
                    good_label, labelled_data = format_data(plots_to_combine[basis][-1],this_xydydy,settings.spectrum_color_indices[level],square_index,"2-exp ratio fit",settings.legend_setting,settings.y_label)
                    labelled_level = 0
                else:
                    if good_label:
                        format_data(plots_to_combine[basis][-1],this_xydydy,settings.spectrum_color_indices[level],square_index)
                    elif labelled_data:
                        plots_to_combine[basis][labelled_level].get_set(labelled_data[0],labelled_data[1]).update_properties(legend="")
                        good_label, labelled_data = format_data(plots_to_combine[basis][-1],this_xydydy,settings.spectrum_color_indices[level],square_index,"2-exp ratio fit")
                        labelled_level = level
                    else:
                        good_label, labelled_data = format_data(plots_to_combine[basis][-1],this_xydydy,settings.spectrum_color_indices[level],square_index,"2-exp ratio fit")
                        labelled_level = level
                
        if tmin_files[basis]["single"]:
            for level,file in enumerate(tmin_files[basis]["single"]):
                plots_to_combine[basis].append(xmgrace_parser.AgrFile(file,True))
                this_xydydy = get_data_by_type(plots_to_combine[basis][-1],xydydy)
                if level ==0:
                    good_label, labelled_data = format_data(plots_to_combine[basis][-1],this_xydydy,settings.spectrum_color_indices[level],diamond_index,"1-exp fit",settings.legend_setting,settings.y_label)
                    labelled_level = level
                else:
                    if good_label:
                        format_data(plots_to_combine[basis][-1],this_xydydy,settings.spectrum_color_indices[level],diamond_index)
                    elif labelled_data:
                        plots_to_combine[basis][labelled_level].get_set(labelled_data[0],labelled_data[1]).update_properties(legend="")
                        good_label, labelled_data = format_data(plots_to_combine[basis][-1],this_xydydy,settings.spectrum_color_indices[level],diamond_index,"1-exp fit")
                        labelled_level = level
                    else:
                        good_label, labelled_data = format_data(plots_to_combine[basis][-1],this_xydydy,settings.spectrum_color_indices[level],diamond_index,"1-exp fit")
                        labelled_level = level
                
        if tmin_files[basis]["double"]:
            for level,file in enumerate(tmin_files[basis]["double"]):
                plots_to_combine[basis].append(xmgrace_parser.AgrFile(file,True))
                this_xydydy = get_data_by_type(plots_to_combine[basis][-1],xydydy)
                if level==0:
                    good_label, labelled_data = format_data(plots_to_combine[basis][-1],this_xydydy,settings.spectrum_color_indices[level],4,"2-exp fit",settings.legend_setting,settings.y_label)
                    labelled_level = level
                else:
                    if good_label:
                        format_data(plots_to_combine[basis][-1],this_xydydy,settings.spectrum_color_indices[level],4)
                    elif labelled_data:
                        plots_to_combine[basis][labelled_level].get_set(labelled_data[0],labelled_data[1]).update_properties(legend="")
                        good_label, labelled_data = format_data(plots_to_combine[basis][-1],this_xydydy,settings.spectrum_color_indices[level],4,"2-exp fit")
                        labelled_level = level
                    else:
                        good_label, labelled_data = format_data(plots_to_combine[basis][-1],this_xydydy,settings.spectrum_color_indices[level],4,"2-exp fit")
                        labelled_level = level
                
                
        #combine all plots
        for i,plot in enumerate(plots_to_combine[basis]):
            if i==0:
                spectrum_plots[basis] = plot
            else:
                all_xydydy = get_data_by_type(plot,xydydy)
                transfer_data( copy.deepcopy(plots_to_combine[basis][i]), all_xydydy, spectrum_plots[basis])
                all_xy = get_data_by_type(plot,"xy")
                transfer_data( copy.deepcopy(plots_to_combine[basis][i]), all_xy, spectrum_plots[basis])
                
                world1 = plot.graphs[0].get_properties(["world"])[0]
                world2 = spectrum_plots[basis].graphs[0].get_properties(["world"])[0]
                new_world = expand_world(world1,world2)

                minx = round(new_world[0])
                maxx = round(new_world[2])
                miny,maxy = stretch_estimate_lines(spectrum_plots[basis],minx,maxx)

                ypadding = (maxy-miny)/(6.0)
                new_world[1] = miny-ypadding
                new_world[3] = maxy+2.0*ypadding
                new_world_str = str(new_world[0])+", "+str(new_world[1])+", "+str(new_world[2])+", "+str(new_world[3])
                spectrum_plots[basis].graphs[0].update_properties(world=new_world_str)
                
        #ticks
        spacing = round((new_world[3]-new_world[1])/settings.max_number_of_ticks,3)
        spectrum_plots[basis].graphs[0].update_properties(yaxis_tick_major=spacing)
                
        strip_unnecessary_data(spectrum_plots[basis], (new_world[3]-new_world[1]))
                
    return spectrum_plots

def spectrum_tmin_format3(tmin_files):
    plots_to_combine = {}
    spectrum_plots = {}
    for basis in tmin_files.keys():
        plots_to_combine[basis] = []
        print(basis)
        
        if tmin_files[basis]["singleR"]:
            for level,file in enumerate(tmin_files[basis]["singleR"]):
                plots_to_combine[basis].append([])
                plots_to_combine[basis][level].append(xmgrace_parser.AgrFile(file,True))
                this_xydydy = get_data_by_type(plots_to_combine[basis][level][-1],xydydy)
                format_data(plots_to_combine[basis][level][-1],this_xydydy,blue_index,circle_index,"1-exp ratio fit",settings.legend_setting,settings.y_label)
                
        if tmin_files[basis]["doubleR"]:
            for level,file in enumerate(tmin_files[basis]["doubleR"]):
                plots_to_combine[basis][level].append(xmgrace_parser.AgrFile(file,True))
                this_xydydy = get_data_by_type(plots_to_combine[basis][level][-1],xydydy)
                format_data(plots_to_combine[basis][level][-1],this_xydydy,green_index,square_index,"2-exp ratio fit",settings.legend_setting,settings.y_label)
                
        if tmin_files[basis]["double"]:
            for level,file in enumerate(tmin_files[basis]["double"]):
                plots_to_combine[basis][level].append(xmgrace_parser.AgrFile(file,True))
                this_xydydy = get_data_by_type(plots_to_combine[basis][level][-1],xydydy)
                format_data(plots_to_combine[basis][level][-1],this_xydydy,11,4,"2-exp fit",settings.legend_setting,settings.y_label)
                
        if tmin_files[basis]["single"]:
            for level,file in enumerate(tmin_files[basis]["single"]):
                plots_to_combine[basis][level].append(xmgrace_parser.AgrFile(file,True))
                this_xydydy = get_data_by_type(plots_to_combine[basis][level][-1],xydydy)
                format_data(plots_to_combine[basis][level][-1],this_xydydy,purple_index,diamond_index,"1-exp fit",settings.legend_setting,settings.y_label)
                       
        #combine all plots
        for level,plots in enumerate(plots_to_combine[basis]):
            new_world = []
            for i,plot in enumerate(plots):
                plot_title = basis+"_"+str(level)
                if i==0:
                    spectrum_plots[plot_title] = plot
                else:
                    all_xydydy = get_data_by_type(plot,xydydy)
                    transfer_data( copy.deepcopy(plots_to_combine[basis][level][i]), all_xydydy, spectrum_plots[plot_title])
                    all_xy = get_data_by_type(plot,"xy")
                    transfer_data( copy.deepcopy(plots_to_combine[basis][level][i]), all_xy, spectrum_plots[plot_title])

                    world1 = plot.graphs[0].get_properties(["world"])[0]
                    world2 = spectrum_plots[plot_title].graphs[0].get_properties(["world"])[0]
                    new_world = expand_world(world1,world2)

                    minx = round(new_world[0])
                    maxx = round(new_world[2])
                    miny,maxy = stretch_estimate_lines(spectrum_plots[plot_title],minx,maxx)

#                     ypadding = (maxy-miny)#*(2.0)
#                     new_world[1] = miny-ypadding
#                     new_world[3] = maxy+3.0*ypadding
                    new_world_str = str(new_world[0])+", "+str(new_world[1])+", "+str(new_world[2])+", "+str(new_world[3])
                    spectrum_plots[plot_title].graphs[0].update_properties(world=new_world_str)

            if new_world:
                #ticks
                spacing = round((new_world[3]-new_world[1])/settings.max_number_of_ticks,3)
                spectrum_plots[plot_title].graphs[0].update_properties(yaxis_tick_major=spacing)
                strip_unnecessary_data(spectrum_plots[plot_title], (new_world[3]-new_world[1]))
                
    return spectrum_plots

"""expand_world(world_1,world_2)
        world_1 - string of the world parameters for a given xmgrace plot
            "xmin, ymin, xmax, ymax"
        world_2 - same structure as world_1
        
        returns a list of values [xmin,ymin,xmax,ymax] of a world that includes both
        the world_1 and world_2
"""
def expand_world(world_1,world_2):
    world1 = [float(i) for i in world_1.split(",")]
    world2 = [float(i) for i in world_2.split(",")]
    return [min(world1[0],world2[0]),min(world1[1],world2[1]),max(world1[2],world2[2]),max(world1[3],world2[3])]

    
"""stretch_estimate_lines( plot_handle, xmin, xmax)
        plot_handle - xmgrace parser plot object 
        xmin - new xmin for all estimate lines
        xmax - new xmax for all estimate lines
        
        sigmond scripts produces tmin plots that contain many varying fits compared
        to the currently chosen fit for the final spectrum plot. When combining the
        tmin plots, there is often different tmin ranges and so these estimate lines
        need to be extended to the new xrange. This function goes through and finds
        all of the estimate lines (type xy) and modifies their range. It returns the
        miny and maxy of the estimates.
"""
def stretch_estimate_lines( plot_handle, xmin, xmax):
    miny = 300.0
    maxy = -300.0
    xy_data = get_data_by_type(plot_handle,xy)
    for g,s in xy_data:
        this_line = plot_handle.datasets[s].get_data()
        this_line[0,0] = xmin
        this_line[0,1] = xmax
        if this_line[1,0]<miny:
            miny = this_line[1,0]
        if this_line[1,1]<miny:
            miny = this_line[1,1]
        if this_line[1,0]>maxy:
            maxy = this_line[1,0]
        if this_line[1,1]>maxy:
            maxy = this_line[1,1]
        plot_handle.datasets[s].set_data(np.transpose(this_line))
    return miny,maxy

"""strip_unnecessary_data( plot_handle, yrange)
        plot_handle - xmgrace parser plot object 
        yrange - max size for error bars 
        
        in a given plot, goes through and removes all data points with error bars 
        larger than the yrange.
"""
def strip_unnecessary_data( plot_handle, yrange):
    xydydy_data = get_data_by_type(plot_handle,xydydy)
    for g,s in reversed(xydydy_data):
        this_line = plot_handle.datasets[s].get_data()
        if this_line.shape[0]<this_line.size and this_line.any():
            bad_indices = []
            for i in range(0,this_line.shape[1]):
                if (this_line[2][i]+this_line[3][i])>=yrange:
                    bad_indices.append(i)
            #print(bad_indices)
            for i in range(len(bad_indices)):
                this_line = np.delete(this_line,bad_indices.pop(),1)
            if this_line.any():
                plot_handle.datasets[s].set_data(np.transpose(this_line))
            else:
                plot_handle.kill_set(g,s)
        elif this_line.shape[0]==this_line.size and this_line.any():
            if (this_line[2]+this_line[3])>=yrange:
                plot_handle.kill_set(g,s)

"""transfer_data( plothandle_source, datasets, plothandle_target)
        plothandle_source - xmgrace parser plot object including the datasets to be
            copied
        datasets - list of the indexes for datasets in the plothandle_source object
            to be copied ([[g0,s0],[g0,s1],...])
        plothandle_target - xmgrace parser plot object that will recieve the datasets
        
        Copies the dataset objects from the plothandle_source and inserts them into
        plothandle_target with indices that do not conflict with plothandle_target's
        datasets. 
"""
def transfer_data( plothandle_source, datasets, plothandle_target):
    for g,s in datasets:
        this_set = copy.deepcopy(plothandle_source.get_set(g,s))
        this_dataset = copy.deepcopy(plothandle_source.get_dataset(g,s))
        this_index = len(plothandle_target.graphs[g].sets)
        plothandle_target.datasets.append(this_dataset)
        plothandle_target.graphs[g].sets.append(this_set)
        plothandle_target.datasets[-1]._force_index(g, this_index)
        plothandle_target.graphs[g].sets[-1]._force_index(this_index)
        
"""get_data_by_type(plot_handle,this_type)
        plot_handle - xmgrace parser plot object 
        this_type - xmgrace data type (xy, xydydy...)
        
        returns a list of the dataset indexes in plot_handle of the data type this_type
        ([[g0,s0],[g0,s1],...])
"""
def get_data_by_type(plot_handle,this_type):
    xydydy_data = []
    for data in plot_handle.datasets:
        if data.get_type() == this_type:
            xydydy_data.append(data.get_g_s())
    return xydydy_data
    
"""format_data(plot_handle, xydydy_list, this_color_index, symbol_index, 
    legend_label="", this_legend_setting=None, this_label=None)
        plot_handle - xmgrace parser plot object 
        xydydy_list - list of dataset indices ([[g0,s0],[g0,s1],...]) to be formatted
        this_color_index - xmgrace color index that the datasets will be recolored to
        symbol_index - xmgrace symbol index that the datasets will be converted to
        legend_label - legend label for given datasets. Default is no label.
        this_legend_setting - string of "x,y" location on plot for the legend to be 
            located. Default is original sigmond script output setting
        this_label - new yaxis label. Original sigmond script output is default
        
        formats the datasets given by xydydy_list in plot_handle and returns a boolean
        on whether the labelling was successful and the indices of the dataset that
        was labelled (empty list if not labelled). (return boolean, [g,s])
"""
def format_data(plot_handle, xydydy_list, this_color_index, symbol_index, legend_label="", this_legend_setting=None, this_label=None):
    #fill all holes in single exponential plot
    for g,s in xydydy_list:
        plot_handle.get_set(g,s).update_properties(symbol_fill_pattern=symbol_pattern_filled[0])
        plot_handle.get_set(g,s).update_properties(symbol=symbol_index)
    
    bad_pvalues = []
    #find the bad pvalue data
    for g,s in xydydy_list:
        if plot_handle.get_set(g,s).get_properties([symbol_color_tag]) == [str(red_index)]:
            bad_pvalues.append([g,s])
                
    #recolor and empty holes of bad pvalue data
    for g,s in bad_pvalues:
        plot_handle.get_set(g,s).update_properties(symbol_fill_pattern=symbol_pattern_empty[0])
           
    #recolor all data to this color
    recolor_data( plot_handle, xydydy_list, this_color_index)
#     for g,s in xydydy_list:
#         plot_handle.get_set(g,s).update_properties(symbol_color=this_color_index)
#         plot_handle.get_set(g,s).update_properties(symbol_fill_color=this_color_index)
#         plot_handle.get_set(g,s).update_properties(errorbar_color=this_color_index)
                
    #label data
    labelled = False
    labelled_data = []
    for g,s in xydydy_list:
        if [g,s] not in bad_pvalues:
            plot_handle.get_set(g,s).update_properties(legend=legend_label)
            labelled=True
            labelled_data = [g,s]
            break
    if not labelled:
        if bad_pvalues:
            plot_handle.get_set(bad_pvalues[0][0],bad_pvalues[0][1]).update_properties(legend=legend_label)
            labelled_data = [bad_pvalues[0][0],bad_pvalues[0][1]]
        
    #clear extra labels
    plot_handle.drawing_objects = []
        
    #turn on and set legend
    if this_legend_setting:
        plot_handle.graphs[0].update_properties(legend_on_off="on",legend=this_legend_setting)
        
    #set axis labels
    if this_label:
        plot_handle.graphs[0].update_properties(xaxis_label="t\\s\Qmin\\N\q/a")
        plot_handle.graphs[0].update_properties(yaxis_label=this_label)
    return labelled,labelled_data

"""recolor_data( plot_handle1, this_list, color_index)
        plot_handle1 - xmgrace parser plot object 
        this_list - list of dataset indices ([[g0,s0],[g0,s1],...]) to be recolored
        color_index - xmgrace color index that the datasets will be recolored to
        
        recolors the datasets to color_index in plot_handle1
"""
def recolor_data( plot_handle1, this_list, color_index):
    #recolor all data to this color
    for g,s in this_list:
        plot_handle1.get_set(g,s).update_properties(symbol_color=color_index)
        plot_handle1.get_set(g,s).update_properties(symbol_fill_color=color_index)
        plot_handle1.get_set(g,s).update_properties(errorbar_color=color_index)
        
"""print_to_svg(plot_handle,filestub)
        plot_handle - xmgrace parser plot object
        filestub - string for naming output file (filestub.svg)
        
        prints the plot_handle to svg file (filestub.svg)
"""
def print_to_svg(plot_handle,filestub):
    this_output = plot_handle.hardcopy(filestub+".png")
    if "truncated" in this_output:
        for i in range(0,30):
            if "truncated" in this_output:
                current_size = list(plot_handle.get_size())
                current_size[0] = current_size[0]+0.2
                plot_handle.set_size(current_size[0],current_size[1])
                plot_handle.set_graph_view(0, x_max=current_size[0]-settings.standard_xmax_gap, width=settings.standard_graph_width, y_min=settings.standard_y_min, height=settings.standard_graph_height)
                this_output = plot_handle.hardcopy(filestub+".png")
            else:
                print("SUcceess")
                break
    
"""retrieve_xmgrace_data_xydydy( files )
        files - a dict of files where the values can be a single string of a 
            filename or a list of filename strings.
            
        Takes the dict of files and makes a dict of dataframes with all the same 
        keys and pulls the xydydy type data from the file(s) and puts those values
        into the dataframe with columns: time, value, error bar 1, error bar 2.
"""
def retrieve_xmgrace_data_xydydy( files ):
    datasets = {}

    for key in files.keys():
        t = np.array([])
        val = np.array([])
        err1 = np.array([])
        err2 = np.array([])
        if type(files[key])==str:
            this_plot = xmgrace_parser.AgrFile(files[key],True) 
            xydydy_data_indexes = get_data_by_type(this_plot,xydydy)
            for (g,s) in xydydy_data_indexes:
                this_t, this_val, this_err1, this_err2 = this_plot.get_dataset(g,s).get_data()
#                 print(this_t, this_val)
                if this_t.shape:
                    t = np.concatenate((t, this_t))
                    val = np.concatenate((val, this_val))
                    err1 = np.concatenate((err1, this_err1))
                    err2 = np.concatenate((err2, this_err2))
                else:
                    if len(t):
                        t = np.append(t, [this_t])
                        val = np.append(val, [this_val])
                        err1 = np.append(err1, [this_err1])
                        err2 = np.append(err2, [this_err2])
                    else:
                        t = np.array([this_t])
                        val = np.array([this_val])
                        err1 = np.array([this_err1])
                        err2 = np.array([this_err2])
#                 print(t, val)
            
        else:
            for file in files[key]:
                this_plot = xmgrace_parser.AgrFile(file,True) 
                xydydy_data_indexes = get_data_by_type(this_plot,xydydy)
                for (g,s) in xydydy_data_indexes:
                    this_t, this_val, this_err1, this_err2 = this_plot.get_dataset(g,s).get_data()
                    if this_t.shape:
                        t = np.concatenate((t, this_t))
                        val = np.concatenate((val, this_val))
                        err1 = np.concatenate((err1, this_err1))
                        err2 = np.concatenate((err2, this_err2))
                    else:
                        if len(t):
                            np.insert(t, -1, this_t)
                            np.insert(val, -1, this_val)
                            np.insert(err1, -1, this_err1)
                            np.insert(err2, -1, this_err2)
                        else:
                            t = np.array([this_t])
                            val = np.array([this_val])
                            err1 = np.array([this_err1])
                            err2 = np.array([this_err2])
                    
        datasets[key] = pd.DataFrame()
        datasets[key].insert(0,0,t)
        datasets[key].insert(1,1,val)
        datasets[key].insert(2,2,err1)
        datasets[key].insert(3,3,err2)
    return datasets

"""retrieve_xmgrace_data_xy( files )
        files - a dict of files where the values can be a single string of a 
            filename or a list of filename strings.
            
        Takes the dict of files and makes a dict of dicts with all the same 
        keys and pulls the xy type data from the file(s) and puts those values
        into the dict with keys: fits, errs.
"""
def retrieve_xmgrace_data_xy( files ):
    fits_all = {}
    for key in files.keys():
        if type(files[key])==str:
            val = np.array([])
            this_plot = xmgrace_parser.AgrFile(files[key],True) 
            xy_data_indexes = get_data_by_type(this_plot,xy)
            for (g,s) in xy_data_indexes:
                this_t, this_val = this_plot.get_dataset(g,s).get_data()
                val = np.concatenate((val, this_val))
            sorted_values = list(set(val))
            sorted_values.sort()
            fits = sorted_values[1]
            errs = [sorted_values[0],sorted_values[2]]
            
        else:
            fits = []
            errs = []
            for file in files[key]:
                val = np.array([])
                this_plot = xmgrace_parser.AgrFile(file,True) 
                xy_data_indexes = get_data_by_type(this_plot,xy)
                for (g,s) in xy_data_indexes:
                    this_t, this_val = this_plot.get_dataset(g,s).get_data()
                    val = np.concatenate((val, this_val))
                sorted_values = list(set(val))
                sorted_values.sort()
                fits.append(sorted_values[1])
                errs.append([sorted_values[0],sorted_values[2]])
        
        fits_all[key] = {"fit":fits,"err":errs}
    return fits_all



################## python functions ##########################################

"""generate_python_rest_mass_plot( fits, tmins, level=None, plot_format=1 )
        fits - a dict with keys 'fits' and 'errs' where 'errs' is list of length two
        tmins - a dict of dataframes with keys for each fit type and the columns:
            time, value, error bar 1, error bar 2.
        level - this spectrum basis level num (if it's important)
        plot_format - plot format value (1-3) used for formatting the tmin plots in the
            same manner as the xmgrace tmin plots.
        
        generates or contributes to a tmin plot by plotting the fits, tmins, and uses the
        level and plot_format to properly color and label the plot according the chosen
        plot format.
"""
def generate_python_rest_mass_plot( fits, tmins, level=None, plot_format=1 ):
    minx = 10
    maxx = 10
    for i,label in enumerate(list(tmins.keys())):
        if not tmins[label].empty:
            minx = min(min(np.array(tmins[label][0])),minx)
            maxx = max(max(np.array(tmins[label][0])),maxx)
            
                
            if plot_format==1:
                this_color=settings.colors[i]
                if level:
                    this_label = None
                else:
                    this_label = label
            else:
                this_color=settings.colors[level]
                if level:
                    this_label = None
                else:
                    this_label = None
                    plt.errorbar([],[], capsize=5, color="black", marker=settings.markers[i], linewidth=0.0, elinewidth=1.5,label=label)
            
            
            plt.errorbar(np.array(tmins[label][0]),np.array(tmins[label][1]),np.concatenate([[np.array(tmins[label][3])],[np.array(tmins[label][2])]]),  capsize=5, color=this_color, marker=settings.markers[i], linewidth=0.0, elinewidth=1.5 ,label=this_label)
    
    plt.hlines(fits['fit'],minx,maxx,color="black")
    plt.hlines(fits['err'][0],minx,maxx,color="black",ls="--")
    plt.hlines(fits['err'][1],minx,maxx,color="black",ls="--")

        
"""
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
"""        
        
#run plots
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument("config", help="config file")
    args = parser.parse_args()

    with open(args.config, "r") as yamlfile:
        project_info = yaml.load(yamlfile, Loader=yaml.FullLoader)
    
    if "rest_mass" in project_info.keys():
        out_dir = project_info["rest_mass"]['out_dir']
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
        for i,particle in enumerate(project_info["rest_mass"]['particles'].keys()):
            print(f"\nCreating {particle} rest mass tmin plot...")
            files = project_info["rest_mass"]['particles'][particle]['files']
            organized_files = {}
            for file in files:
                organized_files[settings.tmin_file_tags[file[-5:]]] = file #assumes no ratio (single hadrons should not have ratio fits)
            if project_info["rest_mass"]['out_type']=="xmgrace":
                ylabel = f'"am\\s{settings.xmgrace_format[particle]}"'
                if ('single' in organized_files.keys()) and ('double' in organized_files.keys()):
                    if 'geometric' in organized_files.keys():
                        this_rest_mass_tmin_plot = rest_mass_tmin_plot(organized_files['single'],organized_files['double'],ylabel,filename_3 = organized_files['geometric'])
                    else:
                        this_rest_mass_tmin_plot = rest_mass_tmin_plot(organized_files['single'],organized_files['double'],ylabel)
                elif (settings.fit_longnames['single'] in organized_files.keys()) and (settings.fit_longnames['double'] in organized_files.keys()):
                    if settings.fit_longnames['geometric'] in organized_files.keys():
                        this_rest_mass_tmin_plot = rest_mass_tmin_plot(organized_files[settings.fit_longnames['single']],organized_files[settings.fit_longnames['double']],ylabel,filename_3 = organized_files[settings.fit_longnames['geometric']])
                    else:
                        this_rest_mass_tmin_plot = rest_mass_tmin_plot(organized_files[settings.fit_longnames['single']],organized_files[settings.fit_longnames['double']],ylabel)
                else:
                    print(f"The xmgrace plots are only configured for 'single' and 'double' and/or 'geometric'.")

                outfilestub = os.path.join(out_dir,project_info["rest_mass"]['particles'][particle]['outfilestub'])
                this_rest_mass_tmin_plot.write(outfilestub+".agr")
                print_to_svg(this_rest_mass_tmin_plot,outfilestub)
                
            elif project_info["rest_mass"]['out_type']=="python":
                combine = False
                n_particles = len(project_info["rest_mass"]['particles'].keys())
                if ('combine' in project_info["rest_mass"].keys()) and (n_particles > 1):
                    combine = project_info["rest_mass"]['combine']
                    
                #retrieve data
                tmin_data = retrieve_xmgrace_data_xydydy(organized_files)
                fit_data = retrieve_xmgrace_data_xy(organized_files) #update the parsing for the new changes
                
                #plot data
                if (combine and not i) or not combine:
                    f = plt.figure()
                    f.set_figwidth(8*(combine+1))
                    f.set_figheight(6)
                
                if combine:
                    plt.subplot(np.ceil(n_particles/2), 2, i+1)
                    
                generate_python_rest_mass_plot( fit_data[list(fit_data.keys())[0]], tmin_data )
                
                plt.xlabel(r"$t_{\textup{min}}/a$")
                latex_rest_mass = settings.latex_format[particle].replace('$','')
                plt.ylabel(rf"$am_{{{latex_rest_mass}}}$")
                plt.legend()
#                 plt.ylim(fit_data[0]-3.0*(fit_data[1]-fit_data[0]), fit_data[2]+6.0*(fit_data[2]-fit_data[1]) )
                
                if not combine:
                    plt.tight_layout()
                    plt.savefig( os.path.join(out_dir,f'{particle}_tmin.pdf'))
                    plt.clf()
        
        if project_info["rest_mass"]['out_type']=="python" and combine:
            plt.tight_layout()
            plt.savefig( os.path.join(out_dir,'rest_mass_tmin.pdf'))
                
        
    if "spectrum_tmins" in project_info.keys():
        plot_format = project_info["spectrum_tmins"]["plot_format"]
        if project_info["spectrum_tmins"]['out_type']=="xmgrace":
            for channel in project_info["spectrum_tmins"]["channels"]:
                out_dir = channel['out_dir']
                if not os.path.exists(out_dir):
                    os.mkdir(out_dir)
                print("\n",channel["name"])

                if plot_format == 1:
                    tmin_plots = spectrum_tmin_format1(utils.find_tmin_spectrum_files(channel))
                elif plot_format == 2:
                    tmin_plots = spectrum_tmin_format2(utils.find_tmin_spectrum_files(channel))
                else:
                    tmin_plots = spectrum_tmin_format3(utils.find_tmin_spectrum_files(channel))

                for basis in tmin_plots.keys():
                    file_stub = os.path.join(out_dir,basis)
                    tmin_plots[basis].write(file_stub+".agr")
                    print_to_svg(tmin_plots[basis],file_stub)
                utils.zip_channel(channel["name"]+f'_format{plot_format}', tmin_plots, sub_dir = out_dir)
                
        if project_info["spectrum_tmins"]['out_type']=="python":
            for channel in project_info["spectrum_tmins"]["channels"]:
                files_to_zip = []
                out_dir = channel['out_dir']
                if not os.path.exists(out_dir):
                    os.mkdir(out_dir)
                print("\n",channel["name"])
                tmin_plots = utils.find_tmin_spectrum_files_python(channel)
                if (plot_format >= 1) or (plot_format <= 3):
                    f = plt.figure()
                    f.set_figwidth(8)
                    f.set_figheight(8)
                else:
                    break
                for basis in tmin_plots.keys():
                    print(basis)
                    file_stub = os.path.join(out_dir,basis)
                    for level in tmin_plots[basis].keys():
                        if tmin_plots[basis][level]:
#                             print(tmin_plots[basis][level].keys())
                            tmin_data = retrieve_xmgrace_data_xydydy(tmin_plots[basis][level])
                            fit_data = retrieve_xmgrace_data_xy(tmin_plots[basis][level])
                            if plot_format == 3:
                                generate_python_rest_mass_plot( fit_data[list(fit_data.keys())[0]], tmin_data)
                            else:
                                generate_python_rest_mass_plot( fit_data[list(fit_data.keys())[0]], tmin_data, level, plot_format)

                            if plot_format == 3:
                                plt.xlabel(r"$t_{\textup{min}}/a$")
                                plt.ylabel(r"$aE_{\textup{fit}}$")
                                plt.legend()
                                plt.tight_layout()
                                plt.savefig(f'{file_stub}_tmin_ROT{level}.pdf')
                                files_to_zip.append(f'{basis}_tmin_ROT{level}.pdf')
                                plt.savefig(f'{file_stub}_tmin_ROT{level}.png')
                                files_to_zip.append(f'{basis}_tmin_ROT{level}.png')
                                plt.clf()
                            
                    if (plot_format == 1) or (plot_format == 2):
                        plt.xlabel(r"$t_{\textup{min}}/a$")
                        plt.ylabel(r"$aE_{\textup{fit}}$")
                        plt.legend()
                        plt.tight_layout()
                        plt.savefig(f'{file_stub}_tmin.pdf')
                        files_to_zip.append(f'{basis}_tmin.pdf')
                        plt.savefig(f'{file_stub}_tmin.png')
                        files_to_zip.append(f'{basis}_tmin.png')
                        plt.clf()
            
                utils.zip_channel(channel["name"]+f'_format{plot_format}',files_to_zip,"", sub_dir = out_dir)
