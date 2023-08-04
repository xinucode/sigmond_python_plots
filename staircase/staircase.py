#!/usr/bin/env python3

# this script was edited in haste during the final weeks of my thesis. please have
# some mercy, future Danny.

################################################################################################
# This script makes a staircase plot for an input XML of the form
# <Staircase>
#     <Title>A1g</Title>
#     <YAxisLabel>$a_t E$</YAxisLabel>
#     <Category>
#         <Name>Category name</Name>
#         <Color>matplotlib compatible color string</Color> (should be unique for each category)
#         <Colors>r g b...</Colors> (if sorting by index, you can instead specify unique colors for each level)
#         <Hatched>true false false true...</Hatched> optional
#         <Index>0</Index> choose an ordering of the categories
#         <Energies>e0 e1 e2 e3...</Energies>
#         <Uncertainties>u0 u1 u2 u3...</Uncertainties>
#     </Category>
#     ...more categories
#     <Legend> (optional: if colors are manually specified, you can change the legend here)
#         <Colors>r g b...</Colors>
#         <Labels>label0 label2...</Label>
#     </Legend>
#     <LegendLoc>matplotlib compatible legend location</LegendLoc> (optional: default is lower right)
#    <NoninteractingLevels>
#        <Color>matplotlib compatible color string</Color>
#        <Energies>e0 e1 e2 e3...</Energies>
#        <Uncertainties>u0 u1 u2 u3...</Uncertainties>
#        <Labels>
#            <Label>label 0</Label>
#            <Label>label 1</Label>
#            <Label>label 2</Label>
#            ...
#        </Labels>
#        <Fontsize>9</Fontsize>
#    </NoninteractingLevels>
#    <Thresholds>
#        <ThresholdEnergies>e0 e1 e2...</ThresholdEnergies>
#        <ThresholdUncertainties>u0 u1 u2...</ThresholdUncertainties>
#    </Thresholds>
# </Staircase>
################################################################################################

from functools import reduce
from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import rcParams
import xml.etree.ElementTree as et
import numpy as np
import sys

from matplotlib import rc
rc('text', usetex=True)
rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size':14})

try:
    input_file_name = sys.argv[1]
except:
    sys.exit('Usage: staircase input.xml (optional: output filename)')

root = et.parse(input_file_name).getroot()

# count how many boxes there will be (not counting extra boxes from splitting into multiple colors)
energy_elements = root.iter('Energies')
num_boxes = 0
for el in energy_elements:
    num_boxes += len(el.text.split())

# arbitarily choose the x-axis of our plot to have the range [0,1]
box_width = 0.8 * 1.0 / num_boxes

# To sort by category, we sort by the user-defined <Index>
#
# This reduce function just makes sure that
# every <Category> must have an <Index> element.
assert(reduce(lambda x, y: x * y, [el.find('Index') is not None for el in root.findall('Category')]))

boxes = []
for category in root.findall('Category'):
    energies_str = category.find('Energies').text
    energies = np.array([float(en) for en in energies_str.split()])
    argsort_en = np.argsort(energies)
    energies = energies[argsort_en]
    uncertainties_str = category.find('Uncertainties').text
    uncertainties = np.array([float(unc) for unc in uncertainties_str.split()])
    uncertainties = uncertainties[argsort_en]
    color = category.find('Color')
    colors = category.find('Colors')
    hatched = category.find('Hatched')
    x_ind = 0
    if color is not None:
        for energy, error in zip(energies, uncertainties):
            boxes.append({'y' : energy - error,
                        'height' : 2 * error,
                        'color' : category.find('Color').text,
                        'index' : int(category.find('Index').text),
                        'x_ind' : x_ind})
            x_ind += 1
    elif colors is not None:
        colors = np.array(colors.text.split())
        colors = colors[argsort_en]
        if hatched is not None:
            hatched = np.array(hatched.text.split())
            hatched = hatched[argsort_en]
            for energy, error, color, hatch in zip(energies, uncertainties, colors, hatched):
                color_split = color.split('/')
                if len(color_split) == 1:
                    boxes.append({'y' : energy - error,
                                'height' : 2 * error,
                                'color' : color,
                                'index' : int(category.find('Index').text),
                                'x_ind' : x_ind,
                                'hatch' : False if hatch == 'false' else True})
                    x_ind += 1
                else:
                    for i, cs in enumerate(color_split):
                        boxes.append({'y' : (energy - error) + i * 2 * error / len(color_split),
                                    'height' : 2 * error / len(color_split),
                                    'color' : cs,
                                    'index' : int(category.find('Index').text),
                                    'x_ind' : x_ind,
                                    'hatch' : False if hatch == 'false' else True})
                    x_ind += 1
        else:
            for energy, error, color in zip(energies, uncertainties, colors):
                color_split = color.split('/')
                if len(color_split) == 1:
                    boxes.append({'y' : energy - error,
                                'height' : 2 * error,
                                'color' : color,
                                'index' : int(category.find('Index').text),
                                'x_ind' : x_ind})
                    x_ind += 1
                else:
                    for i, cs in enumerate(color_split):
                        boxes.append({'y' : (energy - error) + i * 2 * error / len(color_split),
                                    'height' : 2 * error / len(color_split),
                                    'color' : cs,
                                    'index' : int(category.find('Index').text),
                                    'x_ind' : x_ind})
                    x_ind += 1
    else:
        raise(RuntimeError(f'Need to specify Color or Colors for {category.text}'))

# energies = []
category_indices = []
# for box in boxes:
#     energies.append(box['y'] + 0.5 * box['height'])
#     category_indices.append(box['index'])
x_inds = []
for box in boxes:
    x_inds.append(box['x_ind'])
    category_indices.append(box['index'])

sorted_box_indices = np.lexsort((x_inds, category_indices))

boxes = [boxes[i] for i in sorted_box_indices]

# define x_ind_global to be a global index
x_ind_global = 0
boxes[0]['x_ind_global'] = 0
for i in range(1,len(boxes)):
    if boxes[i - 1]['x_ind'] != boxes[i]['x_ind']:
        x_ind_global += 1
    boxes[i]['x_ind_global'] = x_ind_global

ni_levels = root.find('NoninteractingLevels')
if ni_levels is not None:
    energies_str = ni_levels.find('Energies').text
    energies = [float(en) for en in energies_str.split()]
    uncertainties_str = ni_levels.find('Uncertainties').text
    uncertainties = [float(unc) for unc in uncertainties_str.split()]
    for level in zip(energies, uncertainties):
        x_ind_global += 1
        boxes.append({'y' : level[0] - level[1],
                      'height' : 2 * level[1],
                      'color' : ni_levels.find('Color').text,
                      'x_ind_global' : x_ind_global})
# for ind, box in enumerate(boxes):
#     box['x'] = (float(ind) + 0.1) / num_boxes
for box in boxes:
    box['x'] = (float(box['x_ind_global'] + 0.1) / (x_ind_global + 1))

fig, ax = plt.subplots(1)
rcParams['hatch.linewidth'] = 1.0 
for box in boxes:
    hatch = None
    if 'hatch' in box:
        hatch = '///' if box['hatch'] else None
    ax.add_patch(Rectangle((box['x'], box['y']), box_width, box['height'],
                           fill = True, facecolor = box['color'], edgecolor='#404040',
                           linewidth=0.5, hatch=hatch))

if ni_levels is not None:
    ni_boxes = boxes[-1*len(energies_str.split()):]
    ni_box_heights = [ni_box['height'] for ni_box in ni_boxes]
    for box, label in zip(ni_boxes, ni_levels.find('Labels')):
        # plt.text(box['x'] + box_width*1.1, box['y'] + 0.5 * box['height'], label.text,
        #     fontsize=int(ni_levels.find('Fontsize').text), horizontalalignment='left',
        #     verticalalignment='center')
        # # if label.text != r'$\overline K \left( 1 \right) K \left(1 \right)$':
        # #     plt.text(box['x'] + box_width*1.1, box['y'] + 0.5 * box['height'], label.text,
        # #         fontsize=int(ni_levels.find('Fontsize').text), horizontalalignment='left',
        # #         verticalalignment='center')
        # # else:
        # #     plt.text(box['x'] - 0.3 * box_width, box['y'] + 2.1 * box['height'], label.text,
        # #         fontsize=int(ni_levels.find('Fontsize').text), horizontalalignment='left',
        # #         verticalalignment='bottom')
        # # if label.text == r'$\pi \left( 0 \right) \eta \left(0 \right)$':
        # #     plt.text(box['x'], box['y'] + box['height'] + 0.01, label.text,
        # #         fontsize=int(ni_levels.find('Fontsize').text), horizontalalignment='left',
        # #         verticalalignment='bottom')
        if label.text == r'$\overline K \left( 1 \right) K \left(1 \right)$':
            plt.text(box['x'] + box_width, box['y'] + box['height'] + 0.01, label.text,
                fontsize=int(ni_levels.find('Fontsize').text), horizontalalignment='right',
                verticalalignment='bottom')
        elif label.text == r'$\pi \left( 0 \right) \eta^\prime \left(0 \right)$':
            plt.text(box['x'] + box_width, box['y'] + box['height'] + 0.02, label.text,
                fontsize=int(ni_levels.find('Fontsize').text), horizontalalignment='right',
                verticalalignment='bottom')
        else:
            plt.text(box['x'], box['y'] + box['height'] + 0.01, label.text,
                fontsize=int(ni_levels.find('Fontsize').text), horizontalalignment='left',
                verticalalignment='bottom')

ceilings = [box['y'] + box['height'] for box in boxes]
floors = [box['y'] for box in boxes]
ymin = min(floors) - 0.1*(max(ceilings) - min(floors))
ymax = max(ceilings) + 0.1*(max(ceilings) - min(floors))
plt.ylim([ymin, ymax])
left_walls = [box['x'] for box in boxes]
right_walls = [box['y'] + box_width for box in boxes]
xmin = 0
xmax = 1.0 if root.find('NoninteractingLevels') is None else 1.04
plt.xlim([xmin, xmax])

# Add vertical lines to separate categories
categories = root.findall('Category')
if len(categories) > 1:
    num_energies = 0
    cat_ind = [int(cat.find('Index').text) for cat in categories]
    sorted_indices = np.argsort(cat_ind)
    categories = [categories[ind] for ind in sorted_indices]
    for category in categories[:len(categories) - 1*int(root.find('NoninteractingLevels') is None)]:
        num_energies += len(category.find('Energies').text.split())
        plt.plot([(num_energies) / num_boxes, (num_energies) / num_boxes], [-100, 100], '--k')

# plt.plot([19.0 / num_boxes, 19.0 / num_boxes], [-100, 100], '--k')

# add thresholds
thresholds = root.find('Thresholds')
if thresholds is not None:
    energies = [float(e) for e in thresholds.find('ThresholdEnergies').text.split()]
    uncertainties = [float(u) for u in thresholds.find('ThresholdUncertainties').text.split()]
    assert(len(energies) == len(uncertainties))
    for e, u in zip(energies, uncertainties):
        if u/e > 1e-6:
            ax.add_patch(Rectangle((xmin, e-u), xmax-xmin, 2*u, fill=True,
                        linestyle='--', facecolor='#787878', linewidth=0.5, edgecolor='#000000'))
    else:
        plt.plot([xmin, xmax],[e,e],linestyle='--', color='k', linewidth=0.1)

plt.xticks([])
plt.title(root.find('Title').text)
plt.ylabel(root.find('YAxisLabel').text, rotation=1, fontsize=24)
plt.gca().yaxis.set_label_coords(-0.05, 1)
handles = []
legend = root.find('Legend')
if legend is not None:
    colors = legend.find('Colors').text.split()
    labels = legend.find('Labels').text.split()
    assert(len(colors) == len(labels))
    for color, label in zip(colors, labels):
        handles.append(mpatches.Patch(facecolor=color, label=label, edgecolor='#404040',
                           linewidth=0.5))
else:    
    for category in root.findall('Category'):
        handles.append(mpatches.Patch(facecolor=category.find('Color').text,
                                    label=category.find('Name').text,
                                    edgecolor='#404040',
                                    linewidth=0.5))
legendloc = root.find('LegendLoc')
if legendloc is not None:
    loc = legendloc.text
else:
    loc = 'lower right'
plt.legend(handles=handles, loc=loc, fontsize=12, frameon=False)
fig.set_size_inches(10, 5)
# fig.set_size_inches(10, 6)
try:
    file_name = sys.argv[2]
    plt.savefig(file_name, bbox_inches='tight')
except:
    plt.show()