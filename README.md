# wqsreg
### A Stata command for Weighted Quantile Sum regression

- Current version: `V.1.0.0` 
- Release date: `12dec2025`

---

### Description

`wqsreg` allows estimating Weighted Quantile Sum regression for continuous, binary, and count outcomes

Further details can be found in the help file.

### Installation

- To install the current version of `wqsreg` directly from GitHub, run:
```Stata
net install wqsreg, from("https://raw.githubusercontent.com/PonzanoMarta/wqsreg/master/") replace
``` 
from within a web-aware Stata (version 13+).

- For older versions of Stata, download and extract the [zip file](https://github.com/PonzanoMarta/wqsreg/archive/master.zip) and then run:
```Stata
net install wqsreg, from(mydir) replace 
```
from within Stata, where *mydir* is the directory that containes the extracted files.

- After installation, see the help file:
```Stata
help wqsreg
```


### Authors

Marta Ponzano (1, 2), Stefano Renzetti (3), Andrea Bellavia (4,5)

*(1) University of Genoa, Genoa, Italy (2) Link Campus University, Rome, Italy (3) University of Parma, Parma, Italy (4) Harvard T.H. Chan School of Public Health, Boston, MA, USA (5) Harvard Medical School, Boston, MA, USA*

### References

Carrico C, Gennings C, Wheeler DC, Factor-Litvak P. Characterization of weighted quantile sum regression for highly correlated data in a risk analysis setting. Journal of agricultural, biological, and environmental statistics. 2015 Mar;20(1):100-20.

Tanner EM, Bornehag CG, Gennings C. Repeated holdout validation for weighted quantile sum regression. MethodsX. 2019 Nov 22;6:2855-2860. doi: 10.1016/j.mex.2019.11.008. PMID: 31871919; PMCID: PMC6911906.
