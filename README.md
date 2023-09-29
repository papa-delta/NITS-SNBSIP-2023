# A Custom Error Correction Algorithm for DNA Assembly

### Abstract

This repository is built to house our code for project work involving the design and construction of a custom error correction algorithm for DNA sequences. Our algorithm integrates both [KmerCo](https://github.com/patgiri/KmerCo-Main) and [RobustBF](https://github.com/patgiri/robustBF) within itself.

We have investigated the performance of our custom error correction algorithm by conducting rigorous experiments using DNA sequences from four organisms - the African
forest elephant (_Loxodonta cyclotis_), the Sunda flying lemur (_Galeopterus variegatus_), the
gray mouse lemur (_Microcebus murinus_), and the common minke whale (_Balaenoptera acutorostrata_).

We were able to experimentally show that our algorithm was able to reduce the erroneous rate and increase the trustworthy rate. However, the algorithm wasn’t efficient
and upto current standards of other error correcting algorithms ([Lighter](https://github.com/mourisl/Lighter), for example). Our project uncovers a lot of future avenues that can be pursued in this regard, which will help improve the correctness and efficiency of the overall process.

This project is undertaken by us as part of the Satyendranath Bose Summer Internship Programme 2023 at National Institue of Technology, Silchar.

### Steps to run the error correction code:

1. Download FASTQ files for the four datasets from [here](https://drive.google.com/drive/folders/1zM8VxICK28C0U_05r0hAktD_mW4kfIPE). Alternatively, one can also download it from [here](https://github.com/patgiri/KmerCo-Main/blob/master/README.md#dataset).

2. Retrieve the sequences from the FASTQ file into a text file. The following `AWK` and `sed` commands can be used for the same:

    ```$ awk '{if(NR%4==2)print $0}' dataset.fastq > sequence.txt```

    ```$ sed -i ':a; N; s/\\n/ /; ta' sequence_dataset.txt```

    Extract the sequences for all the datasets similarly, and keep them in the same folder. The next step takes care of running the algorithm on all the different extracted sequences automatically.

3. Run the `runner.sh` script, with the inputs `path1` and `path2`; where path1 is the absolute path to the repository folder (containing `main.c`), and path2 is the absolute path to the folder containing the extracted sequences from the datasets in .txt format.

    ```./runner.sh path1 path2```

    Wait for the script to complete running. It will create a separate results folder inside the directory pointed to by `path1`, and put all the obtained outputs from the various datasets in that folder.

Alternatively, one can manually compile and run `main.c` with a text file containing the extracted sequences from a FASTQ file as input.

### References

List of references, resources, and research materials we used for the project are given below:

- Nayak, Sabuzima & Patgiri, Ripon. (2023). KmerCo: A lightweight K-mer counting technique with a tiny memory footprint. [arXiv:2305.07545](https://doi.org/10.48550/arXiv.2305.07545)

- Sabuzima Nayak and Ripon Patgiri. 2021. [robustBF: A High Accuracy and Memory Efficient 2D Bloom Filter](https://doi.org/10.48550/arXiv.2106.04365).

- Sabuzima Nayak and Ripon Patgiri. 2021. countBF: A general-purpose high accuracy and space efficient counting bloom filter. [2021 17th International Conference on Network and Service Management (CNSM). IEEE, Izmir, Turkey, 355–359](https://ieeexplore.ieee.org/document/9615556 ).

- Sabuzima Nayak and Ripon Patgiri. 2019. A Review on Role of Bloom Filter on DNA Assembly. [IEEE Access 7 (2019), 66939–66954](https://doi.org/10.1109/ACCESS.2019.2910180).

- Song, L., Florea, L. & Langmead, B. Lighter: fast and memory-efficient sequencing error correction without counting. [Genome Biol 15, 509 (2014)](https://doi.org/10.1186/s13059-014-0509-9).

- [Github repository of Lighter](https://github.com/mourisl/Lighter), a fast and memory-efficient sequencing error corrector
  
- [Github repository of RobustBF](https://github.com/patgiri/robustBF), a memory-efficient and highly accurate 2D Bloom Filter

- [Github repository of KmerCo](https://github.com/patgiri/KmerCo-Main),  a lightweight k-mer counting technique with a tiny memory footprint
- [Github repository of CountBF](https://github.com/patgiri/KmerCo-Main),  a lightweight k-mer counting technique with a tiny memory footprint


### Acknowledgements: 

We thank and acknowledge our guide, Dr. Ripon Patgiri, and supervisor, Ms. Sabuzima Nayak, for supporting and guiding us during the internship. We also acknowledge and thank everyone associated with the Satyendranath Bose Summer Internship at the Department of CSE, NIT Silchar, including the relevant authorities
and faculties.

### Authors:

- [Preetodeep Dev](https://github.com/papa-delta/), B.Tech (CSE), Assam University, Silchar.
  
- [Swayampakula Kedharnath](https://github.com/SWAYAMPAKULA-KEDHARNATH), B.Tech (CSE), National Institute of Technology, Silchar.
