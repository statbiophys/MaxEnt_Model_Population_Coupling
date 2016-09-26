# Maximum Entropy Models for population coupling
Models are presented in the article [A tractable method for describing complex couplings between neurons and population rate, Gardella, Marre and Mora, 2016, eNeuro](http://eneuro.org/content/3/4/ENEURO.0160-15.2016)

This repository allows you to learn the minimal, linear-coupling and complete coupling models in Matlab. It is then possible to compute the predictions used in the article. Examples are provided in script `EXAMPLE.m`

## Warning
The code uses .mex functions, allowing Matlab to run C code. The .mex files must be compiled before using them on a new computer, running script `COMPILE_mex_files.m.`
