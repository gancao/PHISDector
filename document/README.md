# PHISDetector: a great tool to detect and systematically study diverse in silico phage-host interaction signals
[![standard-readme compliant](https://img.shields.io/badge/readme%20style-standard-brightgreen.svg?style=flat-square)](https://github.com/iwillnever/PHISDector)

## Table of Contents

- [Background](#background)
- [Requirements](requirements)
- [Install](#install)
- [Usage](#usage)
- [Contributors](#contributors)
- [License](#license)

## Background

Phage-host interactions are appealing systems to study bacterial adaptive evolution and are increasingly recognized as playing an important role in human health and diseases, which may contribute to novel therapeutic agents, such as phage therapy to combat multi-drug resistant infections.PHISDetector is the first comprehensive tool to detect and systematically study diverse in silico phage-host interaction signals, including analyses for oligonucleotide profile/sequence composition, CRISPR-targeting, prophages,alignment-based similarity, protein-protein interactions, special gene check and co-occurrence/co-abundance patterns.Further more,PHISDector provides fancy visualizations for users in [http://www.microbiome-bigdata.com/PHISDetector/index/](http://www.microbiome-bigdata.com/PHISDetector/index/)


## Requirements ##
The source code is written by python3,c++ and wrapped by PyInstaller. Thus it requires python3 and c++ compiler. It works under Linux environment. The following python packages 
are required.
Requirements: 

1. numpy
 
2. Biopython
 
3. joblib

## Development & Funding ##
The development of PHISDector was supported as a collaboration of School of Life Science and Technology, Harbin Institute of Technology.

## Install ##

- step1:Download the whole packages and partial profiles from [https://github.com/iwillnever/PHISDector](https://github.com/iwillnever/PHISDector)

- step2:Download the rest of large profiles and databases from [http://www.microbiome-bigdata.com/PHISDetector/index/download](http://www.microbiome-bigdata.com/PHISDetector/index/download)


- step3:Unpack the corresponding files downloaded from step2 and put these large profiles and databases in the directory named database in step1.


- step4:When you organize the whole files well,the corresponding directory structure are displayed as shown below. 
![](https://github.com/iwillnever/PHISDector/blob/master/standalone_directory_structure.png)


## Usage

The standalone version of PHISDector can predict the hosts of query phages in fasta ,multi-fasta or GenBank format,combining 5 different signals including sequence composition, CRISPR-targeting, prophages, regions of genetic homology and protein-protein or domain-domain interactions.The program of prediction of hosts was divided into two steps.First,detect all the candidate hosts in terms of CRISPR-targeting, prophages and regions of genetic homology,satisfying the **loose** criteria referred in the paper [http://www.microbiome-bigdata.com/PHISDetector/index/](http://www.microbiome-bigdata.com/PHISDetector/index/).Second,the candidate hosts satisfying the strict criteria will be directly output and the rest of candidate hosts will be sent to the machine learning models including RandomForest(RF), Decision Trees (DT), Logistic Regression (LR), and Support Vector Machines (SVM) with RBF kernel and linear kernel, Gaussian Naive Bayes, and Bernoulli Naive Bayes.Finally,the consensus results will be written to a text file and intermediate results of each signal will be saved in corresponding directory.

- step1:Add the PHIS directory named bin into your environment PATH or you can use the absolute path to run the program.
- step2:To run a prediction, you should proceed like the following instructions.The required input is input file path in fasta or multi-fasta or Genbank format and output directory.The optional input is signal name(all,crispr,prophage,protein_protein_interaction,blast),the default option is all.If you want to see more parameters,please excute `PHIS` in command line.The visualization html page is generated if you run all signals.
- step3:Predict the hosts of your query phages using all the signals with a simple command:<br>
`PHIS --input <file path> --output <folder name>`
- step4:Predict the hosts of your query phages using single signal with a simple command:
`PHIS --input <file path> --output <folder name> --model all/crispr/prophage/protein_protein_interaction/blast`
- step5:Predict the hosts of your query phages using the criterias you want with a simple example command,x is integer:
`PHIS --input <file path> --output <folder name> --min_mis_crispr <x> --min_cov_crispr <x>`

## Contributors

This project exists thanks to all the people who contribute.

## License

[MIT](LICENSE) © Richard Littauer
