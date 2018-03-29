#!/usr/bin/env bash

if [ $# -eq 0 ]
  then
    echo "Please specify a Python interpreter to be used for installation!"
    exit 1
    # echo "Python interpreter was not specified for the installation. Using default: python"
    # ${1:-python}
fi

{
    # sudo interpreter -m pip install -r requirements.txt &&
    chmod +x setup.py &&
    $1 ./setup.py install &&
    echo "
***
CLI tool preFold for peptide sequence foldability prediction installed with $1 distribution.
Example usage: 'prefold -i <fasta_file>'
***"
} || {
    echo "
***
Failed to install the package! Please check if pip and the Python distribution are installed on the machine. Alternatively, try to install using virtualenv/venv/conda etc!
***"
}
