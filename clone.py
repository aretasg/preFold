import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

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

# finding hydrophobicity mean and the mean net charge for the whole sequence
ph = 7
phobicity = False
charges = False
for key, value in sequence_dict.items():
    hb_mean = hb_avg (value, hb_data_dict)
    charge = net_charge (value, pka_data_dict, ph, True, True)
    unfoldability = unfoldability_algorithm (hb_mean, charge/(len(value)))
    print ('{3} residues, unfoldability {2:.3f} (Charge: {1:.3f}, Phobic: {0:.3f})'
        .format(hb_mean, charge/(len(value)), unfoldability, len(value)))

# calculating the foldability using a sliding window
# move it to a function?
unfold_plot_dict = {} # tag : np array
hb_plot_dict = {}
charge_plot_dict = {}

window_size = 20
step = 1

for key, value in sequence_dict.items():

    win_start = 0
    win_end = win_start + window_size

    unfold_array = np.zeros (len(value)- int(window_size/2))
    hb_array = np.zeros (len(value)- int(window_size/2))
    charge_array = np.zeros (len(value)- int(window_size/2))

    while win_end <= len(value):
        hb_mean = hb_avg (value[win_start:win_end], hb_data_dict)
        charge_mean = net_charge (value[win_start:win_end], pka_data_dict, ph,
            False, False) / len(value[win_start:win_end])
        unfold_i = unfoldability_algorithm (hb_mean, charge_mean)
        if win_start == 0:
            for i in range (int(window_size/2)): # populates the first k/2 arrays with the first window result
                unfold_array[i] = unfold_i
                hb_array[i] = hb_mean
                if i == 0: # adds N terminus charge into the calculation
                    charge_array[i] = net_charge (value[win_start:win_end],
                        pka_data_dict, ph, False, True) / len(value[win_start:win_end])
                else:
                    charge_array[i] = charge_mean
        else:
            unfold_array[win_start + int(window_size/2 - 1)] = unfold_i
            hb_array[win_start + int(window_size/2 - 1)] = hb_mean
            if win_end == len(value) - 1: # adds C ter charge into the calculation
                charge_array[win_start + int(window_size/2 - 1)] = net_charge (
                    value[win_start:win_end],pka_data_dict, ph, True, False
                    ) / len(value[win_start:win_end])
            else:
                charge_array[win_start + int(window_size/2 - 1)] = charge_mean
        win_start += step
        win_end += step

    unfold_plot_dict[key] = unfold_array
    hb_plot_dict[key] = hb_array
    charge_plot_dict[key] = charge_array

# some info about the unfolded regions
# add SD, fill in the missing end of the sequence
disorder_dict = {}
for key, value in unfold_plot_dict.items():
    summ = 0
    counter = 0
    start_pos = None
    end_pos = 0
    for n, i in enumerate (value):
        if i < -0.005:
            summ += i
            counter += 1
            if start_pos is None:
                start_pos = n
        if i > -0.005 or n == (len(value) - 1):
            end_pos = n
            if summ != 0 and end_pos - (start_pos + 1) > 4:
                region_mean = summ/counter
                #disorder_dict[key] = dickt[(start_pos + 1, end_pos)] = region_mean
                disorder_dict[(start_pos + 1, end_pos)] = region_mean
            start_pos = None
            summ = 0
            counter = 0
print (disorder_dict)

# matplotlib and writting output
# remove first k/2 residues plotting
fig_num = 1
for key, value in unfold_plot_dict.items():
    x_axis = np.arange(1, len(value)+1)
    plt.figure(fig_num)
    plt.plot(x_axis, value, color='k', linewidth=2, zorder=2)
    if phobicity is True:
        plt.plot(x_axis, hb_plot_dict[key], color='blue', linewidth=1, zorder=3)
    if charges is True:
        plt.plot(x_axis, charge_plot_dict[key], color='pink', linewidth=1, zorder=5)
    plt.title ('{}'.format(key))
    plt.ylabel ("Unfoldability")
    plt.xlabel ("Residue Number")
    plt.grid(zorder=10)
    plt.fill_between(x_axis, 0, value, where=value>0, interpolate=True, color='g', zorder=4)
    plt.fill_between(x_axis, 0, value, where=value<0, interpolate=True, color='r', zorder=4)
    red_patch = mpatches.Patch(color='r', label='unfolded')
    green_patch = mpatches.Patch(color='g', label='folded')
    plt.legend(handles=[red_patch, green_patch])
    # plt.show()
    fig_num += 1
    # plt.savefig('fold_predict_{0}.png'.format(key.split('|')[0].replace('>', '')))
