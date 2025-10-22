# wqsreg
### A Stata command for Weighted Quantile Sum regression

- Current version: `0.0.1` 
- Release date: `01Nov2025`

---

### Description

`wqsreg` allows estimating Weighted Quantile Sum regression for continuous, binary, and count outcomes

Further details can be found in Discacciati et al. (2018) and in the help file.

### Installation

- To install the current version of `wqsreg` directly from GitHub, run:
```Stata
net install med4way, from("https://raw.githubusercontent.com/PonzanoMarta/wqsreg/master/") replace
``` 
from within a web-aware Stata (version 13+).

- For older versions of Stata, download and extract the [zip file](https://github.com/PonzanoMarta/wqsreg/archive/master.zip) and then run:
```Stata
net install med4way, from(mydir) replace 
```
from within Stata, where *mydir* is the directory that containes the extracted files.

- After installation, see the help file:
```Stata
help wqsreg
```
- To download in the current working directory the datasets needed to run the example code in the help file, type:
```Stata
net get wqsreg, from("https://raw.githubusercontent.com/PonzanoMarta/wqsreg/master/")
```

### Authors

Marta Ponzano (1), Stefano Renzetti (2), Andrea Bellavia (3,4)

*(1) University of Genoa, Italy (2) University of Parma, Italy (3) Harvard T.H. Chan School of Public Health, Boston, MA, USA (4) Harvard Medical School, Boston, MA, USA*

### References

Carrico C, Gennings C, Wheeler DC, Factor-Litvak P. Characterization of weighted quantile sum regression for highly correlated data in a risk analysis setting. Journal of agricultural, biological, and environmental statistics. 2015 Mar;20(1):100-20.
