#!/bin/bash
#
#########################################################################
#########################################################################
#########################################################################
#########################################################################
#
#	THIS SCRIPT IS THE MAIN CALL OF THE CODE:
#		
#		RUN sh scipt.sh
#
#
#########################################################################
#########################################################################
#########################################################################
#########################################################################
#
#    comments
# make dir for diagns
rm -rf ADiag
mkdir ADiag
#
rm *.o
rm *.dat
make
./main