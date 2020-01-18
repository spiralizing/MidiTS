# MidiTS
<a href="https://zenodo.org/badge/latestdoi/191198911"><img src="https://zenodo.org/badge/191198911.svg" alt="DOI"></a>

MIDI - CSV - Time Series, this repository contains functions for constructing time series of music scores from a given CSV file.
### READ: Now it works with Julia 1.x versions

### Usage
First you need to convert the .mid file to a .csv with midi to csv software freely available here: https://www.fourmilab.ch/webtools/midicsv/

Then open the .csv file (for Julia v0.6):

```
f = readdlm("File.csv", ',')

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
