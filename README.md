# TFM_code
Enhancing Alignment Accuracy through Core Alignment Refinement: A Comparative Tool-Based Study 

---

This repository contains the code used to carry out the master's thesis in bioinfomatics of Jordi Sevilla Fortuny at the University of Valencia.

To run the code first you must compile the [alnlen.c](software/alnlen.c) script with:
```
cd software
gcc -O2 -o alnlen alnlen.c
cd ..
```

Then add the [software](software/) folder to your PATH:
```
soft_dir=$(realpath software/)
echo "export PATH=\$soft_fir:\$PATH" >> $HOME/.bash_profile
source $HOME/.bash_profile
```
Make sure that you install the defined conda environment for each section prior to run the code.

---

The [figures](figures/) folder includes the data and the code necessary to generate the master's thesis figures.
