# Photosystem II (PSII) network analysis
Python code to cumpute excitation energy transfer between chlorophylls in PSII supercomplex

# Author
Heetae Kim

# License
This work is licensed under the GNU General Public License v3.0. For details please read LICENSE.txt.

# General notes
This repository hosts python code to compute network effect of the excitation energy transfer of PSII SC.
The analysis is based on the Förster resonance energy transfer rate between chlorophylls.

# Requirements
numpy and pandas

# Usage
A user should have all seven data files prepared before running `main.py`.
Excuting `main.py` by python in console will automatically load the data file and generate the time series excitation energy trace as Pandas dataframe.
```bash
python main.py
```
