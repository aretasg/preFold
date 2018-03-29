from __future__ import print_function
import sys
import pprefold
import pytest

fasta = open('foo.fasta')
tag = fasta.readline()
sequence = fasta.read().replace('\n', '')


def import_module(seq):
    mean_value_of_charge = pprefold.net_charge(seq) / len(seq)
    return round(float(mean_value_of_charge), 3)

def test_module_import_and_net_charge_function():
    assert import_module(sequence) == -0.003

def test_data_import():
    assert pprefold.data_import() != None

def test_sliding_window_is_empty():
    assert len(pprefold.unfoldty_sliding_win(sequence).index) != 0

def test_fasta_parser():
    assert pprefold.fasta_file_parsing('foo.fasta', 1) != {}

def test_multiple_sequences_in_fasta_parsing():
    assert len(pprefold.fasta_file_parsing('foo2.fasta')) == 4

def test_file_parsing_nuceotide_fasta(capsys):
    with pytest.raises(SystemExit):
        pprefold.fasta_file_parsing('foo_nucleotide.fasta.txt')
    out, err = capsys.readouterr()
    # assert out == "The input '{0}' sequence(s) is not a peptide sequence."
        # .format('1AGJ:A|PDBID|CHAIN|SEQUENCE')
    print(out, err)

fasta.close()
