# Photosystem II (PSII) network analysis
Python code to cumpute excitation energy transfer between chlorophylls in PSII supercomplex

# Author
- Heetae Kim, Data Science Institute, Universidad del Desarrollo, Chile
- Eunchul Kim, National Institute for Basic Biology, Japan

# License
This work is licensed under the GNU General Public License v3.0. For details please read LICENSE.txt.

# General notes
This repository hosts python code to compute FRET rate and the network effect of the excitation energy transfer of PSII SC.
The python code `fret.py` analyzes the Förster resonance energy transfer rate between chlorophylls. The output files are used for the main analysis code `main.py`.
In this repository, all data and output files already generated.

# Requirements
numpy and pandas

# Usage
A user should have all seven raw data files (`database_all_a.csv`, `database_all_b.csv`, `database_major_a.csv`, `database_major_b.csv`, `database_minor_a.csv`,`database_minor_b.csv`, `database_natural.csv`) before running `fret.py`.
The resulting files (`FRET_all_a.csv`, `FRET_all_b.csv`, `FRET_major_a.csv`, `FRET_major_b.csv`, `FRET_minor_a.csv`,`FRET_minor_b.csv`, `FRET_natural.csv`) need to be generated before running `main.py`.
Excuting `fret.py` or `main.py` by python in console will automatically load the necessary files and generate the result:
```bash
python fret.py
python main.py
```
The result `.pkl` files are in pyhton pickle format and it includes the time series excitation energy trace as Pandas dataframe.
