# Photosystem II (PSII) network analysis
Python code to cumpute excitation energy transfer between chlorophylls in PSII supercomplex

# Author
Heetae Kim, Data Science Institute, Universidad del Desarrollo, Chile
Eunchul Kim, National Institute for Basic Biology, Japan

# License
This work is licensed under the GNU General Public License v3.0. For details please read LICENSE.txt.

# General notes
This repository hosts python code to compute network effect of the excitation energy transfer of PSII SC.
The analysis is based on the Förster resonance energy transfer rate between chlorophylls.

# Requirements
numpy and pandas

# Usage
A user should have all seven data files (`FRET_all_a.csv`, `FRET_all_b.csv`, `FRET_major_a.csv`, `FRET_major_b.csv`, `FRET_minor_a.csv`,`FRET_minor_b.csv`, `FRET_natural.csv`) before running `main.py`.
Excuting `main.py` by python in console will automatically load the data file and generate the time series excitation energy trace as Pandas dataframe.
```bash
python main.py
```
The result `.pkl` files are in pyhton pickle format.
