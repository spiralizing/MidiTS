# MidiTS
<a href="https://zenodo.org/badge/latestdoi/191198911"><img src="https://zenodo.org/badge/191198911.svg" alt="DOI"></a>
MIDI - CSV - Time Series, this repository contains functions for constructing time series of music scores from a given CSV file.
### READ: It only works with Julia v0.6, will not work with Julia 1.0 yet.

### Usage
First you need to convert the .mid file to a .csv with midi to csv software freely available here: https://www.fourmilab.ch/webtools/midicsv/

Then just by reading the .csv file for Julia v0.6:

```
f = readcsv("File.csv")

```
If you want to get only the pitch sequences for the piece:

```
ps = get_pitch_seq(f)

```
Or the actual time series (with subdivision):

```
ts = get_time_series(f)

```

Both outputs would be a vector containing the n voices of the piece. 
