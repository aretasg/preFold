#!/usr/bin/env python

import sys
import numpy as np # 1.13.3
import pandas as pd # 0.22.0
import math
import matplotlib.pyplot as plt # 2.1.2
import matplotlib.patches as mpatches
import matplotlib.lines as mlines

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

# a function to calculate standard error
cpdef double std_err (double sd, int n):
    return sd / math.sqrt(n)

# a function to find a net charge of a sequence (Moore, 1985)
# ph value min 0.1 max 14;
def net_charge (str sequence, dict pka_table, double ph, ter_c, ter_n):

    def neg_charge_calc (double ph_level, double pka):
        cdef double z_neg
        z_neg = -1 / (10**(-1*(ph_level - pka)) + 1)
        return z_neg

    def pos_charge_calc (double ph_level, double pka):
        cdef double z_pos
        z_pos = 1 / (10**(ph_level - pka) + 1)
        return z_pos

    cdef double pos_sum = 0, neg_sum = 0
    cdef str i
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

    return pos_sum + neg_sum

# a function to find hydrophobicity, charge and unfoldability of a peptide sequence
# using a sliding window mechanism
# an arbitrary choice of including terminus charges in the first and last elements or
# first and last k/2 elements (currently the former is applied)
cpdef object unfoldty_sliding_win (double ph_level, int win_size, int step, str seq, n_ter, c_ter, dict pka_table, dict hb_dict):

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
            for i in range (int(win_size/2)): # populates the first k/2 array elements with the first window result
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
    dict_all['hydrophobicity'] = hb_array
    dict_all['charge'] = charge_array
    dict_all['unfoldability'] = unfold_array
    return pd.DataFrame(dict_all)

# prints some info about the unfolded regions
cpdef void print_data_info (object value, str tag, dict seq_dict, dict hb_dict, dict pka_table, double ph_level, double boundry):
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

    print ("Summary of '{0}':".format(tag))

    # unfoldability, charge and phobic mean information about the input sequence
    hb_mean = hb_avg (seq_dict[tag], hb_dict)
    charge = net_charge (seq_dict[tag], pka_table, ph_level, True, True)
    unfoldability = unfoldability_algorithm (hb_mean, charge/(len(seq_dict[tag])))
    print ('{3} residues, unfoldability {2:.3f} (Charge: {1:.3f}, Phobic: {0:.3f})'
        .format(hb_mean, charge/(len(seq_dict[tag])), unfoldability, len(seq_dict[tag])))

    # information about the disordered regions in the input sequence
    print ('Number of Disordered Regions: {0}\nLongest Disordered Region: {1}\nNumber of Disordered Residues: {2}'
        .format(len(disorder_dict), longest_region, number_disordered_res))
    for key1, value1 in disorder_dict.items():
        mean_ = np.mean(value1)
        std_dev = np.std(value1)
        print ('Predicted disorder segment: {0}-{1} length: {2} score: {3:.3f} ± {4:.2f}'
            .format(key1[0], key1[1], len(value1), mean_, std_dev))
    print ('')

    # printing sequence coloured by ordered/disordered regions
    try:
        from colorama import Fore # 0.3.9
        from colorama import Style
    except:
        sys.exit('Please make sure colorama module is installed for the coloured representation of the disordered regions')

    string2print = list(seq_dict[tag])
    cdef object v
    # cdef str i2
    # writting colours
    cdef int counter = 0
    for k, v in disorder_dict.items():
        string2print.insert(int(k[0]) - 1 + counter, f'{Fore.RED}')
        string2print.insert(int(k[1]) + 1 + counter, f'{Style.RESET_ALL}')
        string2print.insert(int(k[1]) + 2 + counter, f'{Fore.GREEN}')
        string2print.insert(int(k[0]) - 1 + counter, f'{Style.RESET_ALL}')
        counter += 4
    # indentations
    cdef int index = 1
    cdef str i2
    for n, i2 in enumerate(string2print):
        if len(i2) > 1 or i2 == ' ' or i2 == '\n':
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
    string2print = f'{Fore.GREEN}' + ''.join(string2print) + f'{Style.RESET_ALL}'
    cdef str new_string = ''
    cdef int index2 = 1
    cdef str i3
    for i3 in string2print.split('\n'):
        new_string = ''.join([new_string, str("% 4d" % index2) + ' ' + i3 + '\n'])
        index2 += 50
    print (new_string)
    # print (''.join(string2print) + '\n')
    print (f'{Fore.RED}' + '(Predicted disordered segment)' + f'{Style.RESET_ALL}\n')

# writes data to .csv
def write_data_2_csv (tag, dataframe):
    np.round(dataframe, decimals=5).to_csv('{0}.csv'.format(tag.split('|')[0]), sep=',', index=False)

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
    plt.title ('{0}'.format(tag.split('|')[0]))
    plt.ylabel ("Unfoldability")
    plt.xlabel ("Residue Number")
    plt.grid(zorder=10)
    # filling neg/pos values with green/red colour
    cdef object x_axis_new, y_axis_new
    x_axis_new = x_axis[int(win_size/2 - 1):int(len(y1) - win_size/2)]
    y_axis_new = y1[int(win_size/2 - 1):int(len(y1) - win_size/2)]
    plt.fill_between(x_axis_new, 0, y_axis_new , where=y_axis_new >0, interpolate=True, color='g', zorder=4)
    plt.fill_between(x_axis_new, 0, y_axis_new , where=y_axis_new <0, interpolate=True, color='r', zorder=4)
    # legend
    # shrinking the plot to fit the legend at the bottom
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.15, box.width, box.height * 0.85])
    # creating legend markupds
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
    ax.set_xticks(np.append(locs[1:], len(y1)))
    plt.xticks(rotation='vertical')
    #plt.show()
    plt.savefig('fold_predict_{0}.png'.format(tag.split('|')[0].replace('>', '') + '_' + str(fig_counter)), format='png', dpi=dpi, figsize=(8,4))

if __name__ == '__main__':
    pass
