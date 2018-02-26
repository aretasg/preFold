This is an alpha Fold Index clone more info is coming soon!

A CLI tool to predict peptide sequence foldability. The program accepts a FASTA file with multiple unlimited number of peptide sequences as input. Requires Python 3.6 to run (possibly no issues using it with Python 2.7).
This tool is a Python clone of a FoldIndex web app (https://fold.weizmann.ac.il) and surpasses the original in accuracy and availability of options:
allows the user to modify pH value;
use of residue pKa value tables of one's preference;
inclusion of N and C terminal charges into the calculation;
use of floating point values for increased precion in calculations through out all calculation (i.e the values are only rounded when messages are printed to the user in the CLI);
the boundry line ± 0.005 for calling disordered peptide regions is included in this tool in contrast to the FoldIndex (the authors in the original publication falsley claimed it to be included). In addion, this boundry line can be modified by the user.

example useage: 'python clone.py 1tsr.fasta.txt'

Note that the current version is still not particulary user friendly. This will change in the upcoming releases
