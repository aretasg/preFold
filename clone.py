import sys
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

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

# a function to find a net charge of a sequence
# ph value min 0.1 max 14;
# algorithm from innovagen
def net_charge (sequence, pka_table, ph, ter_c, ter_n):

    def neg_charge_calc (ph_level, pka):
        z_neg = 10**ph_level / (10**ph_level + 10**pka)
        return z_neg

    def pos_charge_calc (ph_level, pka):
        z_pos = 10**pka / (10**ph_level + 10**pka)
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

    if ter_n is True:
        # N terminus
        if pka_table[sequence[0]][1] > ph:
            pos_sum += pos_charge_calc (ph, pka_table[sequence[0]][1])
        elif pka_table[sequence[0]][1] < ph:
            pass

    if ter_c is True:
        # C terminus
        if pka_table[sequence[-1:]][0] < ph:
            neg_sum += neg_charge_calc (ph, pka_table[sequence[-1:]][0])
        elif pka_table[sequence[-1:]][0] > ph:
            pass
    else:
        pass

    z = pos_sum - neg_sum
    return z

def net_charge2 (sequence, pka_table, ph):
    pass

# open and parse a fasta file
# 1. fix the creation of an empty key 2. add additional file format checks
sequence_dict = {} # here all the peptide sequence are stored
try:
    with open (sys.argv[1], 'r') as fasta_file:
        file_2_string = fasta_file.read().split('>')
        for sequence in file_2_string:
            sequence_dict[sequence.split('\n')[0]] = ''.join(sequence.split('\n')[1:]).upper().replace('*', '')
        del sequence_dict['']
        print ('Found {0} peptide sequence(s) in the FASTA file.'.format(len(sequence_dict)))
except IOError:
    print('Cannot open the file! The format must be in FASTA!')
    sys.exit()

# open and parse normalised amino acid hydrophobicity data
hb_data_dict = {}
with open ('HB_DATA_NORM.DAT', 'r') as hb_data:
    for line in hb_data:
        hb_data_dict[line.split()[0]] = float(line.split()[1])

# open and parse pka data
pka_data_dict = {}
with open ('PKA_DATA_CRC.DAT', 'r') as pka_data:
    for line in pka_data:
        pka_data_dict[line.split()[0]] = [float(i) for i in line.split()[1:]]

# calculating the foldability using a sliding window
# an arbitrary choice of including terminus charge in the first and last elements or
# first and last k/2 elements (currently the former is applied)
unfold_plot_dict = {} # fasta tag : np array
hb_plot_dict = {}
charge_plot_dict = {}

ph = 7
k = 20 # window size
step = 1

for key, value in sequence_dict.items():

    win_start = 0
    win_end = win_start + k

    unfold_array = np.zeros (len(value))
    hb_array = np.zeros (len(value))
    charge_array = np.zeros (len(value))

    while win_end <= len(value):
        # calculates the mean results of a sequence in the window
        hb_mean = hb_avg (value[win_start:win_end], hb_data_dict)
        charge_mean = net_charge (value[win_start:win_end], pka_data_dict, ph,
            False, False) / len(value[win_start:win_end])
        unfold_i = unfoldability_algorithm (hb_mean, charge_mean)

        if win_start == 0:
            for i in range (int(k/2)): # populates the first k/2 array elements with the first window result
                hb_array[i] = hb_mean
                if i == 0: # adds N terminus charge into the calculation
                    charge_array[i] = net_charge (value[win_start:win_end],
                        pka_data_dict, ph, False, True) / len(value[win_start:win_end])
                    unfold_array[i] = unfoldability_algorithm (hb_mean, charge_array[i])
                else:
                    charge_array[i] = charge_mean
                    unfold_array[i] = unfold_i

        # populates the last k/2 array elements with the last window result
        elif win_start + int(k/2 - 1) == len(value) - k/2 - 1:
            for i in range(int(len(value) - k/2 - 1), len(value)):
                hb_array[i] = hb_mean
                if i == len(value) - 1: # adds C terminus charge into the calculation
                    charge_array[i] = net_charge (value[win_start:win_end],
                        pka_data_dict, ph, True, False) / len(value[win_start:win_end])
                    unfold_array[i] = unfoldability_algorithm (hb_mean, charge_array[i])
                else:
                    charge_array[i] = charge_mean
                    unfold_array[i] = unfold_i
        else:
            unfold_array[win_start + int(k/2 - 1)] = unfold_i
            hb_array[win_start + int(k/2 - 1)] = hb_mean
            charge_array[win_start + int(k/2 - 1)] = charge_mean

        win_start += step
        win_end += step

    unfold_plot_dict[key] = unfold_array
    hb_plot_dict[key] = hb_array
    charge_plot_dict[key] = charge_array
    # print (unfold_array)

# some info about the unfolded regions
for key, value in unfold_plot_dict.items():
    start_pos = None # the start position of a disordered region
    end_pos = 0 # the end position of a disordered region

    longest_region = None # the longest disordered region in the sequence
    number_disordered_res = 0 # the number of residues in all disordered regions

    disorder_dict = {}
    list_for_region = []
    for n, i in enumerate (value):
        if i < -0.005:
            list_for_region.append(i)
            if start_pos is None:
                start_pos = n
        if (i > 0.005 or n == (len(value) - 1)) and start_pos is not None:
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

    print ("Summary of '{0}':".format(key))

    # unfoldability, charge and phobic mean information about the input sequence
    hb_mean = hb_avg (sequence_dict[key], hb_data_dict)
    charge = net_charge (sequence_dict[key], pka_data_dict, ph, True, True)
    unfoldability = unfoldability_algorithm (hb_mean, charge/(len(sequence_dict[key])))
    print ('{3} residues, unfoldability {2:.3f} (Charge: {1:.3f}, Phobic: {0:.3f})'
        .format(hb_mean, charge/(len(sequence_dict[key])), unfoldability, len(sequence_dict[key])))

    # information about the disordered regions in the input sequence
    print ('Number of Disordered Regions: {0}\nLongest Disordered Region: {1}\nNumber of Disordered Residues: {2}'
        .format(len(disorder_dict), longest_region, number_disordered_res))
    for key1, value1 in disorder_dict.items():
        mean_ = np.mean(value1)
        std_dev = np.std(value1)
        print ('Predicted disorder segment: {0}-{1} length: {2} score: {3:.3f} ± {4:.2f}'
            .format(key1[0], key1[1], len(value1), mean_, std_dev))

# generating unfoldability figures for each sequence
# todo: plot the sequence with green/red highlights of the letters;
# check C and N terminus calculations; decimal module
# exporting of numerical information to csv etc. (use pandas); move it to a function?
phobicity = True
charges = True

fig_counter = 1
for key, value in unfold_plot_dict.items():
    # removing first and last k/2 values
    y_axis = value
    for i in range(0, int(k/2 - 1)):
        y_axis[i] = np.nan
    for i in range(int(len(value) - k/2), len(value)):
        y_axis[i] = np.nan
    # the initial setup
    x_axis = np.arange(1, len(value)+1)
    fig, ax = plt.subplots()
    plt.figure(fig_counter)
    plt.plot(x_axis, value, color='k', linewidth=2, zorder=2)
    # plotting phobicity and charge
    for i in range(0, int(k/2 - 1)):
        hb_plot_dict[key][i] = np.nan
    for i in range(int(len(value) - k/2), len(value)):
        hb_plot_dict[key][i] = np.nan
    if phobicity is True:
        for i in range(0, int(k/2 - 1)):
            hb_plot_dict[key][i] = np.nan
        for i in range(int(len(value) - k/2), len(value)):
            hb_plot_dict[key][i] = np.nan
        plt.plot(x_axis, hb_plot_dict[key], color='blue', linewidth=1, zorder=3)
    if charges is True:
        for i in range(0, int(k/2 - 1)):
            charge_plot_dict[key][i] = np.nan
        for i in range(int(len(value) - k/2), len(value)):
            charge_plot_dict[key][i] = np.nan
        plt.plot(x_axis, charge_plot_dict[key], color='pink', linewidth=1, zorder=5)
    # labels and grid
    plt.title ('{0}'.format(key))
    plt.ylabel ("Unfoldability")
    plt.xlabel ("Residue Number")
    plt.grid(zorder=10)
    # filling neg/pos values with green/red colour
    x_axis_new = x_axis[int(k/2 - 1):int(len(value) - k/2)]
    y_axis_new = y_axis[int(k/2 - 1):int(len(value) - k/2)]
    plt.fill_between(x_axis_new, 0, y_axis_new , where=y_axis_new >0, interpolate=True, color='g', zorder=4)
    plt.fill_between(x_axis_new, 0, y_axis_new , where=y_axis_new <0, interpolate=True, color='r', zorder=4)
    # legend
    red_patch = mpatches.Patch(color='r', label='unfolded')
    green_patch = mpatches.Patch(color='g', label='folded')
    plt.legend(handles=[red_patch, green_patch])
    # x axis tick
    locs = ax.xaxis.get_ticklocs()
    ax.set_xticks(np.append(locs[1:], len(value)))
    fig_counter += 1
    plt.show()
    # plt.savefig('fold_predict_{0}.png'.format(key.split('|')[0].replace('>', '') + str(fig_counter)), format='png', dpi=1000, figsize=(8,4))
