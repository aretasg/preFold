import sys
import numpy as np
import matplotlib.pyplot as plt

# to do: modify fold index function to accept charge as arg
# fasta file parsing
# sliding window
# matplotlib

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

# a function to find a net charge of a sequence
# ph value min 0.1; max 14; terminus groups?
# algorithm from innovagen
def net_charge (sequence, pka_table, ph, terminus):

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

    if terminus is True:
        # N terminus
        if pka_table[sequence[0]][1] > ph:
            pos_sum += pos_charge_calc (ph, pka_table[sequence[0]][1])
        elif pka_table[sequence[0]][1] < ph:
            pass

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
# 1. remove stop codon; 2. fix the creation of an empty key 3. add additional file format checks
sequence_dict = {} # here all the peptide sequence are stored
try:
    with open (sys.argv[1], 'r') as fasta_file:
        file_2_string = fasta_file.read().split('>')
        for sequence in file_2_string:
            sequence_dict[sequence.split('\n')[0]] = ''.join(sequence.split('\n')[1:]).upper()
        del sequence_dict['']
        print ('Found {0} peptide sequence(s) in the FASTA file.'.format(len(sequence_dict)))
except IOError:
    print('Cannot open the file! The format must be in FASTA!')

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

# finding hydrophobicity mean and the mean net charge for the whole sequence
ph = 7
for key, value in sequence_dict.items():
    hb_mean = hb_avg (value, hb_data_dict)
    charge = net_charge (value, pka_data_dict, ph, True)
    unfoldability = unfoldability_algorithm (hb_mean, charge/(len(value)))
    print (hb_mean, charge/(len(value)), unfoldability, len(value))

# calculating the foldability using a sliding window
# calculating hb and charge means uses two loops; can solve it with one loop
# to add n and c terminal charges to the calculations
unfold_plot_dict = {} # tag : np array

for key, value in sequence_dict.items():

    win_start = 0
    window_size = 20
    step = 1
    win_end = win_start + window_size
    unfold_array = np.zeros (len(value)- int(window_size/2))

    while win_end <= len(value):
        hb_mean = hb_avg (value[win_start:win_end], hb_data_dict)
        charge_mean = net_charge (value[win_start:win_end], pka_data_dict, ph, False) / len(value[win_start:win_end])
        unfold_i = unfoldability_algorithm (hb_mean, charge_mean)
        if win_start == 0:
            for i in range (int(window_size/2)):
                unfold_array[i] = unfold_i
        else:
            unfold_array[win_start + int(window_size/2 - 1)] = unfold_i
        # print (win_start + 1, value[win_start], unfold_i)
        win_start += step
        win_end += step

    unfold_plot_dict[key] = unfold_array
    # print (unfold_array)

# some info about the unfolded regions
# matplotlib and writting output
# add hydrophobicity and net charge curves; remove first residue redundancy
for key, value in unfold_plot_dict.items():
    x_axis = [i for i in range (1, len(value) + 1)]
    plt.figure(1)
    plt.plot(x_axis, value, color = 'k', linewidth=2)
    plt.ylabel ("Unfoldability")
    plt.xlabel ("Residue Number")
    plt.grid(zorder=10)
    plt.fill_between(x_axis, 0, value, where=value>0, interpolate=True, color='g', zorder=5)
    plt.fill_between(x_axis, 0, value, where=value<0, interpolate=True, color='r', zorder=5)
    plt.show()
