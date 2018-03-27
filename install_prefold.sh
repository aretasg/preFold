#!/bin/bash

sudo $1 -m pip install -r requirements.txt
chmod +x setup.py
$1 ./setup.py install
echo "***"
echo "CLI tool prefold installed with $1 distribution.
    Example usage: 'prefold -i <fasta_file>'"
echo "***"
