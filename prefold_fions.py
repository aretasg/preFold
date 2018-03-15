#!/usr/bin/env python3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines

# an algorithm to describe peptide charge and hydrophobicity relantionship (Uversky et al, 2000)
def uversky_algorithm (hb_mean):
    r = 2.785 * hb_mean - 1.151
    return abs(r)

# a function to calculate the foldability of a peptide sequence (Prilsuky and Felder et al, 2005)
def unfoldability_algorithm (hb_mean, q_mean):
    index = 2.785 * hb_mean - abs(q_mean) - 1.151
    return index

# a function to calculate hydrophobicity average of a peptide sequence
def hb_avg (sequence, hb_dict):
    hb_sum = 0
    for i in sequence:
        hb_sum += hb_dict[i]
    return hb_sum / len(sequence)

# a function to calculate standard error
def std_err (sd, n):
    return sd / math.sqrt(n)

# a function to find a net charge of a sequence (Moore, 1985)
# ph value min 0.1 max 14;
def net_charge (sequence, pka_table, ph, ter_c, ter_n):

    def neg_charge_calc (ph_level, pka):
        z_neg = -1 / (10**(-1*(ph_level - pka)) + 1)
        return z_neg

    def pos_charge_calc (ph_level, pka):
        z_pos = 1 / (10**(ph_level - pka) + 1)
        return z_pos

    pos_sum = 0
    neg_sum = 0

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
def unfoldty_sliding_win (ph_level, win_size, step, seq, n_ter, c_ter, pka_table, hb_dict):

    win_start = 0
    win_end = win_start + win_size

    unfold_array = np.zeros (len(seq))
    hb_array = np.zeros (len(seq))
    charge_array = np.zeros (len(seq))

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

    dict_all = {}
    dict_all['hydrophobicity'] = hb_array
    dict_all['charge'] = charge_array
    dict_all['unfoldability'] = unfold_array

    dt = np.round(pd.DataFrame(dict_all), decimals=6)

    return dt

# prints some info about the unfolded regions
def print_data_info (value, tag, seq_dict, hb_dict, pka_table, ph_level, boundry):
    start_pos = None # the start position of a disordered region
    end_pos = 0 # the end position of a disordered region

    longest_region = None # the longest disordered region in the sequence
    number_disordered_res = 0 # the number of residues in all disordered regions

    disorder_dict = {}
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
        print ('Number of Disordered Regions: {0}\nLongest Disordered Region: {1}\nNumber of Disordered Residues: {2}'
            .format(len(disorder_dict), longest_region, number_disordered_res))
        for key1, value1 in disorder_dict.items():
            mean_ = np.mean(value1)
            std_dev = np.std(value1)
            print (u'Predicted disorder segment: {0}-{1} length: {2} score: {3:.3f} ± {4:.2f}'
                .format(key1[0], key1[1], len(value1), mean_, std_dev))
        print ('')
    # if there are no disordered regions
    else:
        print ('No disordered regions found.')

    return disorder_dict

# a function to highlight disordered regions using colour (green/red)
def coloured_seq (seq_dict, tag, disorder_dict):
    try:
        from colorama import Fore # 0.3.9
        from colorama import Style
    except:
        sys.exit('Please make sure colorama module is installed for the coloured\
         representation of the disordered regions')

    string2print = list(seq_dict[tag])
    # writting colours
    counter = 0
    for k, v in disorder_dict.items():
        string2print.insert(int(k[0]) - 1 + counter, 'red')
        string2print.insert(int(k[1]) + 1 + counter, 'reset')
        string2print.insert(int(k[1]) + 2 + counter, 'green')
        string2print.insert(int(k[0]) - 1 + counter, 'reset')
        counter += 4
    # indentations
    index = 1
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
    new_string = ''
    index = 1
    stash_colour = 'reset'
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

    print (new_string.replace('red', f'{Fore.RED}').replace('green',
        f'{Fore.GREEN}').replace('reset', f'{Style.RESET_ALL}'))
    # print (''.join(string2print) + '\n')
    print (f'{Fore.RED}' + '(Predicted disordered segment)' + f'{Style.RESET_ALL}\n')

# writes data to .csv
def write_data_2_csv (tag, dataframe, seq):
    csv_df = dataframe.copy()
    csv_df.insert(0, 'Residue', list(seq))
    csv_df.insert(0, 'Residue number', [i for i in range(1, len(seq) + 1)])
    csv_df.set_index('Residue number')
    csv_df.to_csv('{0}.csv'.format(tag.split('|')[0]), sep=',', index=False)

# generating unfoldability figures for each sequence
def generate_figure (y1, y2, y3, win_size, tag, fig_counter, phobicity, charges, dpi):
    # removing first and last k/2 values
    for i in range(0, int(win_size/2 - 1)):
        y1[i] = np.nan
    for i in range(int(len(y1) - win_size/2), len(y1)):
        y1[i] = np.nan
    # the initial setup
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

if __name__ == '__main__':
    pass
