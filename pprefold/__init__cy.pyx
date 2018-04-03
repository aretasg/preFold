#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Author: Aretas Gaspariunas

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import sys
import re
import os

# open and parse a fasta file
cpdef dict fasta_file_parsing (object fasta, int sequence_lim=80):
    cdef str sequence, key, value
    cdef dict seq_dict = {} # here all the peptide sequences are stored
    try:
        with open(fasta) as fasta_file:
            file_2_string = fasta_file.read().split('>')
            for sequence in file_2_string:
                seq_dict[sequence.split('\n')[0]] = ''.join(sequence.split
                    ('\n')[1:]).upper().replace('*', '').replace('X', '') # excludes X from the sequences
            del seq_dict['']
            # check if there a no more than 60 sequencies
            if len(seq_dict) > sequence_lim:
                sys.exit("There are too many sequencies in the '{0}'file. The "
                    "program is designed to accept up to {1} sequencies per a run."
                    .format(fasta, sequence_lim))
            else:
                pass
            # check if the sequence is a peptide
            for key, value in seq_dict.items():
                if re.search(r'S|L|I|M|F|R|E|K|D|O|H|N|Y|V|W|P', value):
                    continue
                else:
                    sys.exit("The input '{0}' sequence(s) is not a peptide sequence."
                        .format(key))
            print ('Found {0} peptide sequence(s) in the FASTA file.\n'
                .format(len(seq_dict)))
    except IOError:
        sys.exit('Could not open the file! Please make sure {0} is in fasta format'
            .format(fasta))

    return seq_dict

# pka and hydrophobicity data import
cpdef list data_import (str pka_file='PKA_DATA_VOET.DAT'):
    # retrieving data file paths
    cdef str this_dir, DATA_PATH, line, i
    cdef dict hb_data_dictionary, pka_data_dictionary
    cdef object pka_data

    this_dir = os.path.dirname(os.path.abspath(__file__))
    DATA_PATH = os.path.join(this_dir, "datasets", "data.txt")

    # open and parse normalised amino acid hydrophobicity data
    hb_data_dictionary = {}
    with open (DATA_PATH.replace('data.txt', 'HB_DATA_NORM.DAT')) as hb_data:
        for line in hb_data:
            hb_data_dictionary[line.split()[0]] = float(line.split()[1])

    # open and parse pKa data
    try:
        pka_data_dictionary = {}
        if pka_file in ['PKA_DATA_VOET.DAT', 'PKA_DATA_CRC.DAT']:
            pka_data = open(DATA_PATH.replace('data.txt', pka_file))
        else:
            pka_data = open(pka_file)
        for line in pka_data:
            pka_data_dictionary[line.split()[0]] = [float(i) for i in line.split()[1:]]
        pka_data.close()

    except IOError:
        sys.exit('Could not open the pKa table file! Please make sure {0} is '
            'in the right format.'.format(pka_file))

    return [hb_data_dictionary, pka_data_dictionary]

hb_data_dict = data_import()[0]
pka_data_dict = data_import()[1]

# an algorithm to describe peptide charge and hydrophobicity relantionship (Uversky et al, 2000)
cpdef double uversky_algorithm (double hb_mean):
    cdef double r
    r = 2.785 * hb_mean - 1.151
    return abs(r)

# a function to calculate the foldability of a peptide sequence (Prilsuky and Felder et al, 2005)
cpdef double unfoldability_algorithm (double hb_mean, double q_mean):
    cdef double index
    index = 2.785 * hb_mean - abs(q_mean) - 1.151
    return index

# a function to calculate hydrophobicity average of a peptide sequence
cpdef double hb_avg (str sequence, dict hb_dict):
    cdef double hb_sum = 0
    cdef str i
    for i in sequence:
        hb_sum += hb_dict[i]
    return hb_sum / len(sequence)

# a function to find a net charge of a sequence (Moore, 1985)
# ph value min 0.1 max 14;
def net_charge (str sequence, dict pka_table=pka_data_dict, double ph=7.4, bint ter_c=True, bint ter_n=True):

    cdef double pos_sum = 0, neg_sum = 0
    cdef str i

    def neg_charge_calc (ph_level, pka):
        cdef double z_neg
        z_neg = -1 / (10**(-1*(ph_level - pka)) + 1)
        return z_neg

    def pos_charge_calc (ph_level, pka):
        cdef double z_pos
        z_pos = 1 / (10**(ph_level - pka) + 1)
        return z_pos

    # calculating net charge of a peptide sequence
    for i in sequence:
        if i in ['D','E','C','Y']:
            neg_sum += neg_charge_calc (ph, pka_table[i][2])
        elif i in ['H','K','R']:
            pos_sum += pos_charge_calc (ph, pka_table[i][2])
        else:
            pass
    # terminal charge calculations
    if ter_n is True:
        # N terminus
        pos_sum += pos_charge_calc (ph, pka_table[sequence[0]][1])
    if ter_c is True:
        # C terminus
        neg_sum += neg_charge_calc (ph, pka_table[sequence[-1:]][0])
    else:
        pass

    z = pos_sum + neg_sum
    return z

# a function to find hydrophobicity, charge and unfoldability of a peptide sequence
# using a sliding window mechanism
# an arbitrary choice of including terminus charges in the first and last elements or
# first and last k/2 elements (currently the former is applied)
cpdef object unfoldty_sliding_win (str seq, double ph_level=7.4, int win_size=50,
    int step=1, bint n_ter=True, bint c_ter=True, dict pka_table=pka_data_dict,
    dict hb_dict=hb_data_dict):

    cdef int win_start = 0
    cdef int win_end = win_start + win_size

    cdef object unfold_array = np.zeros (len(seq))
    cdef object hb_array = np.zeros (len(seq))
    cdef object charge_array = np.zeros (len(seq))
    cdef double hb_mean, charge_mean, unfold_i
    cdef int i

    while win_end <= len(seq):
        # calculates the mean results of a sequence in the window
        hb_mean = hb_avg (seq[win_start:win_end], hb_dict)
        charge_mean = net_charge (seq[win_start:win_end], pka_table, ph_level,
            False, False) / len(seq[win_start:win_end])
        unfold_i = unfoldability_algorithm (hb_mean, charge_mean)

        if win_start == 0:
            # populates the first k/2 array elements with the first window result
            for i in range (int(win_size/2)):
                hb_array[i] = hb_mean
                if i == 0: # adds N terminus charge into the calculation
                    charge_array[i] = net_charge (seq[win_start:win_end],
                        pka_table, ph_level, False, n_ter) / len(seq[win_start:win_end])
                    unfold_array[i] = unfoldability_algorithm (hb_mean, charge_array[i])
                else:
                    charge_array[i] = charge_mean
                    unfold_array[i] = unfold_i

        # populates the last k/2 array elements with the last window result
        elif win_start + int(win_size/2 - 1) == len(seq) - win_size/2 - 1:
            for i in range(int(len(seq) - win_size/2 - 1), len(seq)):
                hb_array[i] = hb_mean
                if i == len(seq) - 1: # adds C terminus charge into the calculation
                    charge_array[i] = net_charge (seq[win_start:win_end],
                        pka_table, ph_level, c_ter, False) / len(seq[win_start:win_end])
                    unfold_array[i] = unfoldability_algorithm (hb_mean, charge_array[i])
                else:
                    charge_array[i] = charge_mean
                    unfold_array[i] = unfold_i
        else:
            unfold_array[win_start + int(win_size/2 - 1)] = unfold_i
            hb_array[win_start + int(win_size/2 - 1)] = hb_mean
            charge_array[win_start + int(win_size/2 - 1)] = charge_mean

        win_start += step
        win_end += step
    # print (unfold_array)

    cdef dict dict_all = {}
    cdef object dt
    dict_all['hydrophobicity'] = hb_array
    dict_all['charge'] = charge_array
    dict_all['unfoldability'] = unfold_array

    dt = np.round(pd.DataFrame(dict_all), decimals=6)

    return dt

# prints some info about the unfolded regions
cpdef dict print_data_info (object value, str tag, dict seq_dict, dict hb_dict, dict pka_table, double ph_level, double boundry):
    start_pos = None # the start position of a disordered region
    # cdef int end_pos = 0 # the end position of a disordered region

    longest_region = None # the longest disordered region in the sequence
    cdef int n, number_disordered_res = 0, end_pos = 0 # the number of residues in all disordered regions
    cdef double i
    cdef dict disorder_dict = {}
    list_for_region = []
    for n, i in enumerate (value):
        if i < -1*boundry:
            list_for_region.append(i)
            if start_pos is None:
                start_pos = n
        if (i > boundry or n == (len(value) - 1)) and start_pos is not None:
            end_pos = start_pos + len(list_for_region)
            if len(list_for_region) > 4:
                disorder_dict[(start_pos + 1, end_pos)] = np.asarray(list_for_region)
                # the total number of residues in all disordered regions
                number_disordered_res += len(list_for_region)
                # finding the longest disordered region
                if longest_region is None or len(list_for_region) > longest_region:
                    longest_region = len(list_for_region)
            start_pos = None
            list_for_region = []
        else:
            pass

    # check if there are any disordered regions
    if disorder_dict:
        print ("Summary of '{0}':".format(tag))

        # unfoldability, charge and phobic mean information about the input sequence
        hb_mean = hb_avg (seq_dict[tag], hb_dict)
        charge = net_charge (seq_dict[tag], pka_table, ph_level, True, True)
        unfoldability = unfoldability_algorithm (hb_mean, charge/(len(seq_dict[tag])))
        print ('{3} residues, unfoldability {2:.3f} (Charge: {1:.3f}, Phobic: {0:.3f})'
            .format(hb_mean, charge/(len(seq_dict[tag])), unfoldability, len(seq_dict[tag])))

        # information about the disordered regions in the input sequence
        print ('Number of Disordered Regions: {0}\nLongest Disordered Region: '
            '{1}\nNumber of Disordered Residues: {2}'
            .format(len(disorder_dict), longest_region, number_disordered_res))
        for key1, value1 in disorder_dict.items():
            print ('Predicted disorder segment: {0}-{1} length: {2} score: {3:.3f} ± {4:.2f}'
                .format(key1[0], key1[1], len(value1), np.mean(value1), np.std(value1)))
        print ('')
    # if there are no disordered regions
    else:
        print ('No disordered regions found.')

    return disorder_dict

# a function to highlight disordered regions using colour (green/red)
cpdef void coloured_seq (dict seq_dict, str tag, dict disorder_dict):

    string2print = list(seq_dict[tag])
    # writting colours
    cdef int counter = 0
    cdef tuple k
    cdef object v
    for k, v in disorder_dict.items():
        string2print.insert(int(k[0]) - 1 + counter, 'red')
        string2print.insert(int(k[1]) + 1 + counter, 'reset')
        string2print.insert(int(k[1]) + 2 + counter, 'green')
        string2print.insert(int(k[0]) - 1 + counter, 'reset')
        counter += 4
    # indentations
    cdef int index = 1
    for n, i in enumerate(string2print):
        if len(i) > 1 or i == ' ' or i == '\n':
            pass
        elif index % 10 == 0:
            if index % 50 == 0:
                string2print.insert(n + 1, '\n')
            else:
                string2print.insert(n + 1, ' ')
            index += 1
        else:
            index += 1

    # numbers
    string2print = 'green' + ''.join(string2print) + 'reset'
    cdef str new_string = ''
    index = 1
    cdef str stash_colour = 'reset'
    for i in string2print.split('\n'):
        # formatting new string
        new_string = ''.join([new_string, 'reset' + str("% 4d" % index) +
            stash_colour + ' ' + i + '\n'])
        # stashing colour
        i = i[::-1]
        try:
            m_obj = i.index('green'[::-1])
            m_obj2 = i.index('red'[::-1])
            if int(m_obj) < int(m_obj2):
                stash_colour = 'green'
            else:
                stash_colour = 'red'
        except ValueError:
            try:
                i.index('red'[::-1])
                stash_colour = 'red'
            except ValueError:
                try:
                    i.index('green'[::-1])
                    stash_colour = 'green'
                except ValueError:
                    pass

        index += 50

    print (new_string.replace('red', '\033[31m').replace('green', '\033[32m')
        .replace('reset', '\033[0m'))
    # print (''.join(string2print) + '\n')
    print ('\033[31m' + '(Predicted disordered segment)' + '\033[0m' + '\n')

# writes data to .csv
cpdef void write_data_2_csv (str tag, object dataframe, str seq):
    cdef object csv_df
    csv_df = dataframe.copy()
    csv_df.insert(0, 'Residue', list(seq))
    csv_df.insert(0, 'Residue number', [i for i in range(1, len(seq) + 1)])
    csv_df.set_index('Residue number')
    csv_df.to_csv('{0}.csv'.format(tag.split('|')[0]), sep=',', index=False)

# generating unfoldability figures for each sequence
cpdef void generate_figure (object y1, object y2, object y3, int win_size, str tag, int fig_counter, phobicity, charges, int dpi):
    # removing first and last k/2 values
    cdef int i
    for i in range(0, int(win_size/2 - 1)):
        y1[i] = np.nan
    for i in range(int(len(y1) - win_size/2), len(y1)):
        y1[i] = np.nan
    # the initial setup
    cdef object x_axis
    x_axis = np.arange(1, len(y1)+1)
    fig, ax = plt.subplots()
    plt.figure(fig_counter)
    plot1 = plt.plot(x_axis, y1, color='k', linewidth=1, zorder=2)
    # plotting phobicity and charge
    if phobicity is True:
        for i in range(0, int(win_size/2 - 1)):
            y2[i] = np.nan
        for i in range(int(len(y2) - win_size/2), len(y2)):
            y2[i] = np.nan
        line = plt.plot(x_axis, y2, color='blue', linewidth=1, zorder=3)
    if charges is True:
        for i in range(0, int(win_size/2 - 1)):
            y3[i] = np.nan
        for i in range(int(len(y3) - win_size/2), len(y3)):
            y3[i] = np.nan
        plt.plot(x_axis, y3, color='pink', linewidth=1, zorder=5)
    # labels and grid
    # plt.title ('{0}'.format(tag.split('|')[0]))
    plt.ylabel ("Unfoldability")
    plt.xlabel ("Residue Number")
    plt.grid(zorder=10)
    # filling neg/pos values with green/red colour
    cdef object x_axis_new, y_axis_new
    x_axis_new = x_axis[int(win_size/2 - 1):int(len(y1) - win_size/2)]
    y_axis_new = y1[int(win_size/2 - 1):int(len(y1) - win_size/2)]
    plt.fill_between(x_axis_new, 0, y_axis_new , where=y_axis_new >0,
        interpolate=True, color='g', zorder=4)
    plt.fill_between(x_axis_new, 0, y_axis_new , where=y_axis_new <0,
        interpolate=True, color='r', zorder=4)
    # legend
    # shrinking the plot to fit the legend at the bottom
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.15, box.width, box.height * 0.85])
    # creating legend markups
    red_patch = mpatches.Patch(color='r', label='unfolded')
    green_patch = mpatches.Patch(color='g', label='folded')
    blue_line = mlines.Line2D([], [], color='blue',
                  markersize=15, label='hydrophobicity')
    pink_line = mlines.Line2D([], [], color='pink',
                  markersize=15, label='charge')
    # setting legend
    if phobicity is True and charges is True:
        ax.legend(handles=[red_patch, green_patch, blue_line, pink_line],
            loc='upper center', bbox_to_anchor=(0.5, -0.17), fancybox=True, ncol=5)
    elif phobicity is True:
        ax.legend(handles=[red_patch, green_patch, blue_line], loc='upper center',
            bbox_to_anchor=(0.5, -0.17), fancybox=True, ncol=5)
    elif charges is True:
        ax.legend(handles=[red_patch, green_patch, pink_line], loc='upper center',
            bbox_to_anchor=(0.5, -0.17), fancybox=True, ncol=5)
    else:
        ax.legend(handles=[red_patch, green_patch], loc='upper center',
            bbox_to_anchor=(0.5, -0.17), fancybox=True, ncol=5)
    # x axis tick
    locs = ax.xaxis.get_ticklocs()
    ax.set_xticks(np.append(locs[1:-1], len(y1)))
    plt.xticks(rotation='vertical')
    #plt.show()
    plt.savefig('fold_predict_{0}.png'.format(tag.split('|')[0].replace('>', '')
        + '_' + str(fig_counter)), format='png', dpi=dpi, figsize=(8,4))

# CLI argument parser and I/O handler
cpdef void main():
    import argparse
    import platform

    # CLI argument parser
    parser = argparse.ArgumentParser(
        description='''A CLI tool to predict foldability of a peptide sequence.\n
            The tool is inteded to be used with Python 3.6 or 2.7.\n
            For more information and support please visit: github.com/aretas2/preFold''',
        epilog='Example usage in CLI: "prefold.py -i foo.fasta"')
    parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('-i', '--fasta',
        help='Specify the FASTA file with a peptide sequence(s) for foldability '
        'prediciton.', required=True)
    optional.add_argument('-ph', '--ph_lvl',
        help='Specify the pH level to be used for the foldability prediction '
        'calculation (0.1 < pH < 14); Default=7.4.', default=7.4, type=int)
    optional.add_argument('-s', '--step',
        help='Specify the step to be used in the sliding window approach to '
        'calculate foldability prediction; Default=1.',
            default=1, type=int)
    optional.add_argument('-k', '--window_size',
        help='Specify the window size to be used for the calculation; '
        'Default=50.', default=50, type=int)
    optional.add_argument('-csv', '--output_csv',
        help='Specify this flag if you wish the numerical output to be '
        'provided as a .csv file.',
            action='store_true')
    optional.add_argument('-z', '--plot_charge',
        help='Specify this flag if you wish the charge of the sequence to be '
        'plotted.', action='store_true')
    optional.add_argument('-hb', '--plot_hb',
        help='Specify this flag if you wish the phobicity of the sequence to be'
        ' plotted.', action='store_true')
    optional.add_argument('-b', '--boundry',
        help='Specify the boundry for calling disordered regions of peptide '
        'sequence; Default=0.005.',
        type=float, default=0.005)
    optional.add_argument('-ter', '--ter_include',
        help='Specify the flag for N and C terminal charges to be NOT '
        'included in the calculation.',
        action='store_false')
    optional.add_argument('-dpi', '--figure_dpi',
        help='Specify the dpi (resolution) of a figure; Default=200.',
        type=int, default=200)
    optional.add_argument('-pka', '--pka_table',
        help='Specify a file with amino acid residue pKa values to be used for '
        'the calculation (must be located in datasets folder).',
        default='PKA_DATA_VOET.DAT')
    args = parser.parse_args()

    if args.ph_lvl < 0.1 or args.ph_lvl > 14:
        parser.error('Please select pH value in range: 0.1 < pH < 14')

    # egg
    if args.fasta.lower() == 'satera'[::-1]:
        print ("eil a si ekac ehT"[::-1])
        sys.exit("*POOF*")
    else:
        pass

    # opening fasta file
    cdef dict sequence_dict
    sequence_dict = fasta_file_parsing(args.fasta)

    # importing data files
    cdef dict pka_data_dict
    pka_data_dict = data_import(args.pka_table)[1]

    # performing hydrophobicity, charge and unfoldability data calculation using a sliding window
    cdef dict data_dict = {}
    cdef str key, value
    for key, value in sequence_dict.copy().items():
        if len(value) > args.window_size: # checks if the sequence len is larger than win size
            data_dict[key] = unfoldty_sliding_win (value, args.ph_lvl,
                args.window_size, args.step, args.ter_include, args.ter_include,
                pka_data_dict, hb_data_dict)
                # returns pd dataframe with hb, charge and unfold columns;
                #stores in a dict under the fasta tag as a key
        else:
            print ('The peptide sequence {0} is shorter than selected window '
                'size and is excluded from the calculation.\n'
                .format(key))
            del sequence_dict[key]

    # generating statistics, general information, sequence highlights and figures
    cdef int figure_number = 1
    cdef str k
    cdef dict v
    cdef dict disorder_dict
    for k, v in data_dict.items():
        if args.output_csv is True:
            write_data_2_csv (k, v, sequence_dict[k])

        disorder_dict = print_data_info (v['unfoldability'], k,
            sequence_dict, hb_data_dict, pka_data_dict, args.ph_lvl, abs(args.boundry))
        if disorder_dict:
            # checking OS of the computer (only Windows needs colorama)
            if re.search(platform.platform(), 'Windows'):
                from colorama import Fore
                from colorama import Style
                coloured_seq (sequence_dict, k, disorder_dict)
            else:
                coloured_seq (sequence_dict, k, disorder_dict)

            generate_figure (v['unfoldability'], v['hydrophobicity'],
                v['charge'], args.window_size, key, figure_number,
                args.plot_hb, args.plot_charge, args.figure_dpi)
            figure_number += 1
        else:
            pass

if __name__ == '__main__':
    main()
