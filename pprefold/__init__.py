import os
import sys

# pka and hydrophobicity data import
def data_import (pka_file='PKA_DATA_VOET.DAT'):
    # retrieving data file paths
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
        sys.exit('Could not open the pKa table file! Please make sure {0} is in the right format.'
            .format(pka_file))

    return [hb_data_dictionary, pka_data_dictionary]

hb_data_dict = data_import()[0]
pka_data_dict = data_import()[1]

# an algorithm to describe peptide charge and hydrophobicity relantionship (Uversky et al, 2000)
def uversky_algorithm (hb_mean):
    r = 2.785 * hb_mean - 1.151
    return abs(r)

# a function to calculate the foldability of a peptide sequence (Prilsuky and Felder et al, 2005)
def unfoldability_algorithm (hb_mean, q_mean):
    index = 2.785 * hb_mean - abs(q_mean) - 1.151
    return index

# a function to calculate hydrophobicity average of a peptide sequence
def hb_avg (sequence, hb_dict=hb_data_dict):
    hb_sum = 0
    for i in sequence:
        hb_sum += hb_dict[i]
    return hb_sum / len(sequence)

# a function to find a net charge of a sequence (Moore, 1985)
# ph value min 0.1 max 14;
def net_charge (sequence, pka_table=pka_data_dict, ph=7.4, ter_c=True, ter_n=True):

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
def unfoldty_sliding_win (seq, ph_level=7.4, win_size=50, step=1, n_ter=True,
    c_ter=True, pka_table=pka_data_dict, hb_dict=hb_data_dict):

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

    dict_all = {}
    dict_all['hydrophobicity'] = hb_array
    dict_all['charge'] = charge_array
    dict_all['unfoldability'] = unfold_array

    dt = np.round(pd.DataFrame(dict_all), decimals=6)

    return dt
