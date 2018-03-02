# preFold

A CLI tool to predict foldability of a peptide sequence. The program accepts a FASTA file with multiple number of peptide sequences as input.

## Getting Started

### Dependecies
**The script was written in and intended to be used with Python 3.6 use it with caution when executing with 2.7 interpreter. You will need a SciPy stack and colorama installed to run the program with full functionality.**

### Example useage:
```
python preFold.py -i foo.fasta
```
This will print out some information about the sequence, disordered regions in CLI and generate .png file showing the predicted foldability of the peptide sequence.

### CLI message example using foo.fasta
```
Found 1 peptide sequence(s) in the FASTA file.
Summary of '1AGJ:A|PDBID|CHAIN|SEQUENCE':
242 residues, unfoldability 0.028 (Charge: -0.003, Phobic: 0.425)
Number of Disordered Regions: 6
Longest Disordered Region: 29
Number of Disordered Residues: 100
Predicted disorder segment: 1-22 length: 22 score: -0.237 ± 0.08
Predicted disorder segment: 27-48 length: 22 score: -0.206 ± 0.10
Predicted disorder segment: 77-81 length: 5 score: -0.118 ± 0.07
Predicted disorder segment: 83-111 length: 29 score: -0.282 ± 0.15
Predicted disorder segment: 154-164 length: 11 score: -0.112 ± 0.05
Predicted disorder segment: 168-178 length: 11 score: -0.069 ± 0.03

   1 EVSAEEIKKH EEKWNKYYGV NAFNLPKELF SKVDEKDRQK YPYNTIGNVF
  51 VKGQTSATGV LIGKNTVLTN RHIAKFANGD PSKVSFRPSI NTDDNGNTET
 101 PYGEYEVKEI LQEPFGAGVD LALIRLKPDQ NGVSLGDKIS PAKIGTSNDL
 151 KDGDKLELIG YPFDHKVNQM HRSEIELTTL SRGLRYYGFT VPGNSGSGIF
 201 NSNGELVGIH SSKVSHLDRE HQINYGVGIG NYVKRIINEK NE

(Predicted disordered segment)
```

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
* -b Specify the boundry for calling disordered regions of peptide sequence (default=± 0.005);
* -f Specify the dpi of a figure (default = 750).

Please use --help flag for more information on all options and parameters.

## Features
This tool is a Python clone of a [FoldIndex](https://fold.weizmann.ac.il) web app. preFold surpasses the original in accuracy, availability of options and fixed mistakes. Here are listed a few of them:
* Allows the user to modify pH value (this feature is crucial for crystallographers);
* The use of a residue pKa value table of one's preference;
* Inclusion of N and C terminal charges into the calculation;
* Uses floating point values for increased precion throughout the calculation (i.e the values are only rounded when messages are printed to the user in the CLI);
* The boundry line ± 0.005 for calling disordered peptide regions is included in this tool in contrast to the FoldIndex (it seems the FoldIndex authors in the publication falsley claimed it to be included). In addion, this boundry line can be modified by the user;
* Publication quality figure generation
* Improved peptide charge plotting (negative values are plotted).

### Upcoming features:
* .pyx file with static typing for Cython compilation for performence increase
* Web app

### Limitations
The tool does not account for the pertrubed pKa values of residues by the neighbouring residue groups and does not assume any disuphide bridges.

## Authors
* **Aretas Gaspariunas**

## License
This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

### Acknowledgments & Disclaimer
preFold is an original work and does not copy any code elements or principles of FoldIndex source code. It was built/cloned only using ideas presented in the FoldIndex (Prilsuky and Felder et al, 2005) and algorithm (Uversky et al, 2000) publications, preFold was further improved using these ideas as a foundation.
The peptide charge calculation method is used as described by Moore, 1985.
