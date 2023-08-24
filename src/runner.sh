#!/bin/bash

# Check if both path1 and path2 are provided as arguments
if [ $# -ne 2 ]; then
  echo "Usage: $0 path1 path2; where path1 is the absolute path to the KmerCo folder (containing KmerCo.c), and path2 is the absolute path to the folder containing the extracted sequences from the datasets in .txt format."
  exit 1
fi

# Store the arguments in variables
path1=$1
path2=$2

mkdir ${path1}/Results

# Compile KmerCo.c
gcc "${path1}/KmerCo.c" -o "${path1}/KmerCo.out" -lm

# Check if the compilation was successful
if [ $? -ne 0 ]; then
  echo "Compilation failed"
  exit 1
fi

# Loop through each filename.txt in path2
for file in "${path2}"/*.txt; do
  # Extract the filename without the extension
  filenamebase=$(basename "$file" .txt)
  # Run KmerCo with the current filename.txt
  ${path1}/KmerCo.out ${path2}/$filenamebase.txt
	
# Define an array of filenames
filenames=("linecount" "Distinct" "Erroneous" "Trustworthy" "Result" "Found" "NotFound")

# Loop through each filename in the array
for filename in "${filenames[@]}"; do

  # Rename the output file using the current filename
  fnappend=`echo $filenamebase | cut -d_ -f2`
  mv "${path1}/${filename}.txt" "${path1}/Results/${filename}_${fnappend}.txt"


done
done



