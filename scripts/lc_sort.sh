#!/bin/bash

if [ $# -ne 1 ]
then
    echo "Error. Invalid number of arguments"
fi

filename=$(basename ${1})   # get the file name without the path
extension=${filename##*.}   # extract the file extension
filename=${filename%.*}     # extract the file name

# if the file has no extension then filename == extension
if [ $filename == $extension ]
then
    extension=""
fi
#echo $filename
#echo $extension


# all extensions will be converted to txt

# convert file contents to lower case
tr A-Z a-z < $1 > ${filename}_lc.txt

# sort file contents
sort ${filename}_lc.txt > ${filename}_lc_sorted.txt

# remove temporary file
rm -f ${filename}_lc.txt

