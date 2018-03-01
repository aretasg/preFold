# preFold
## TL;DR
A CLI tool to predict foldability of a peptide sequence. The program accepts a FASTA file with multiple number of peptide sequences as input.

Example useage: 'python preFold.py -i foo.fasta'. This will print out some information about the sequence, disordered regions in CLI and generate an .png file showing the predicted foldability of the peptide sequence.

**The script was written in and intended to be used with Python 3.6 use it with caution when executing with 2.7 interpreter. You will need a SciPy stack to be installed to run the program.**

## Usage and arguments

* -i <file name> a flag to specify the filename and/or the directory of a file (note that this is the only required option);
* -ph <int> a flag to specify a pH value to be used in the calculation (default = 7.4);
* -pka <file name> a flag to specify a file containing pKa values of amino acids. Please use the default file as an example for the input file format. (default = PKA_DATA_VOET.DAT);
* -s <int> a flag to specify a step for the sliding window to be used in the calculation (default = 1);
* -k <int> a flag to specify the window size (default = 20);
* -csv a flag to generate numerical output in .csv format;
* -z Specify the flag for charge of the sequence to be plotted on the figure;
* -hb Specify the flag for hydrophobicity of the sequence to be plotted on the figure.
* -ter a flag to exclude N and C terminal charges from the calculation;

Please use --help flag for more information on all options and parameters.

## Features
This tool is a Python clone of a [FoldIndex](https://fold.weizmann.ac.il) web app. preFold surpasses the original in accuracy, availability of options and fixed mistakes. Here are listed a few of them:
*Allows the user to modify a pH value (this feature is crucial for crystallographers);
*Allows the use of a residue pKa value table of one's preference;
*Allows inclusion of N and C terminal charges into the calculation;
*Uses floating point values for increased precion throughout all the calculation (i.e the values are only rounded when messages are printed to the user in the CLI);
*The boundry line ± 0.005 for calling disordered peptide regions is included in this tool in contrast to the FoldIndex (It seems the FoldIndex authors in the publication falsley claimed it to be included). In addion, this boundry line can be modified by the user.
*Improved charge plotting (negative values are plotted)

### Limitations
The tool does not account for the pertrubed pKa values of residues by the neighbouring residue groups and does not assume any disuphide bridges.

### Upcoming features:
*Coloured sequence highlights of the disordered regions
*Increases in speed/cython compiled functions
*Web app

### Disclaimer
preFold is an original work and does not copy any elements or principles of FoldIndex source code. It was built/cloned only from using the FoldIndex (Prilsuky and Felder et al, 2005) and algorithm (Uversky et al, 2000) publications in mind, preFold was further improved using these ideas as a foundation.
The peptide charge calculation method is used as described by Moore, 1985.
