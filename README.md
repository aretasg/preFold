# preFold

A CLI tool to predict foldability of a peptide sequence. The program accepts a FASTA file with a multiple number of peptide sequences as input.

## Getting Started

### Dependecies
**The script was written in and intended to be used with Python 3.6 use it with caution when executing with 2.7 interpreter. You will need a SciPy stack and colorama installed to run the program with a full functionality (see [requirements.txt](requirements.txt) file).**

To install required packages use pip (recommended together with virtualenv):
```
pip install -r requirements.txt
```

### Example usage:
```
python preFold.py -i foo.fasta
```
This will print out some information about the sequence, disordered regions in CLI and generate .png file showing the predicted foldability of the peptide sequence.

### Figure example using foo.fasta as input
![alt text](https://raw.githubusercontent.com/aretas2/preFold/master/example/foo.png)

### CLI message example using foo.fasta
![alt text](https://raw.githubusercontent.com/aretas2/preFold/master/example/CLI_message_example.png)

## Usage and arguments
* -i <file name> a flag to specify the filename and/or the directory of a file (note that this is the only required option);
* -ph <int> a flag to specify a pH value to be used in the calculation (default = 7.4);
* -pka <file name> a flag to specify a file containing pKa values of amino acids. Please use the default file as an example for the input file format. (default = PKA_DATA_VOET.DAT);
* -s <int> a flag to specify a step for the sliding window to be used in the calculation (default = 1);
* -k <int> a flag to specify the window size (default = 50);
* -csv a flag to generate numerical output in .csv format;
* -z Specify the flag for charge of the sequence to be plotted on the figure;
* -hb Specify the flag for hydrophobicity of the sequence to be plotted on the figure.
* -ter a flag to exclude N and C terminal charges from the calculation;
* -b Specify the boundry for calling disordered regions of peptide sequence (default=± 0.005);
* -f Specify the dpi (resolution) of a figure (default = 200).

Please use --help flag for more information on all options and parameters.

## Features
This tool is a Python clone of a [FoldIndex](https://fold.weizmann.ac.il) web app. preFold surpasses the original in accuracy, availability of options and fixed mistakes. Here are listed a few of them:
* Allows the user to modify pH value (this feature is crucial for crystallographers);
* The use of a residue pKa value table of one's preference (default is as provided in Voet & Voet);
* Optional inclusion of N and C terminal charges into the calculation;
* Uses floating point values for increased precision throughout the calculation (i.e the values are only rounded when messages are printed to the user in the CLI);
* The boundry line ± 0.005 for calling disordered peptide regions is included in this tool in contrast to the FoldIndex (it seems the FoldIndex authors in the publication mistakenly claimed it to be included). In additon, this boundry line can be modified by the user;
* Generation of publication quality figures with an option to regulate the resolution of the output figure;
* .pyx file with static typing for compilation using Cython for performence increase;
* Improved peptide charge plotting (negative values are plotted instead of absolute values).

### Upcoming features:
* Web app

### Limitations
* The tool does not account for the pertrubed pKa values of residues by the neighbouring residue groups and does not assume any disuphide bridges.
* You might also notice differences in results produced by FoldIndex and preFold this is attributed to the set boundry condition (-b), pH and most importantly the pKa values used for the calculations.
* Coloured schema of disordered regions is properly rendered only in the terminal window.

## Authors
* **Aretas Gaspariunas**

## License
This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

### Acknowledgments & Disclaimer
preFold is an original work and does not copy any code elements or principles of FoldIndex source code. It was built/cloned only using ideas presented in the FoldIndex (Prilsuky and Felder et al, 2005) and the foldability prediction algorithm (Uversky et al, 2000) publications, preFold was further improved using these ideas as a foundation.
The peptide charge calculation method is used as described by Moore, 1985.
