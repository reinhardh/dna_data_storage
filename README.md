Summary
=======

This readme file gives instructions on encoding and decoding information on DNA segments. The code was written by Reinhard Heckel for the publication:

**``Robust chemical preservation of digital information on DNA in silica with error-correcting codes''**, by R. N. Grass, R. Heckel, M. Puddu, D. Paunescu, and W. J. Stark, Angewandte Chemie International Edition, Vol. 54, Nr. 8, pp. 2552–2555, February 2015. 



Example
=======

The following describes an example that translates a part of Archimedes Book, encodes it to DNA segments, and back. 

To start, change to the folder simulate and execute the following commands:

- compile texttodna: (first adopt Makefile to make sure that the correct path to BOOST_LIB is provided) 
	
	`make texttodna' 

The program texttodna can be used to translate text to DNA segments and vice versa as follows.

1. The following command takes the text in the file *archimedes_short.txt*, encodes it on segments of DNA, and stores the segments in *arch_short_dna.txt*:

	`./texttodna --encode --input=../data/archimedes_short.txt --output=../data/arch_short_dna.txt`

2. The following command takes random lines (DNA segments) from *arch_short_dna.txt*, introduces errors, and saves the resulting file to *arch_short_drawnseg.txt*:

	`./texttodna --disturb --input=../data/arch_short_dna.txt --output=../data/arch_short_drawnseg.txt`

3. The following command decodes based on the perturbed data in *arch_short_dna.txt*, and writes the result to *arch_short_rec.txt*; this should recover the original text:

	`./texttodna --decode --numblocks=3 --input=../data/arch_short_drawnseg.txt --output=../data/arch_short_rec.txt`



Installation
============

The code is written in C++, and compillation requires installation of the boost libary. 


Installation of required software on Linux (not tested)
-------------------------------------------------------
	sudo apt-get install gcc
	sudo apt-get install make
	sudo apt-get install libboost-all-dev

Licence
==========

All files are provided under the terms of the Apache License, Version 2.0, see the included file "apache_licence_20" for details.