# PSPManalysis 0.3.4

* Corrected all "Format overflow" and "may be uninitialized" error messages in the code that showed up using gcc (from Rtools) on Windows

# PSPManalysis 0.3.3

* Added a C as well as an R implementation of another example model (Salmon.h and Salmon.h), used in the publication on PSPManalysis

* Added more informative output messages of computations in case no output was generated

* removed the compiler directive "-Wno-format-overflow" from buildSO.R and PSPMecodyn.R as this one is unknown to the clang compiler

* csbread() now returns a list of states when the state argument is of the form "State-4.0" and the value 4.0 occurs more than once in the CSB file

# PSPManalysis 0.3.2

* Corrected return statements without () in `PSPMecodyn.R` and `PSPMevodyn.R` as reported on the CRAN results pages for fedora-clang and fedora-gcc.

* Adapted the building process of the compiled module in `buildSO.R` and `PSPMecodyn.R` to deal with situations in which the current working directory is not the directory in which the model files are stored	

# PSPManalysis 0.3.1

* Updated to version 0.3.1 - Changed compilation command to work with R 4.0.0 under Windows

* Added the manual in bookdown format as a ZIP file 

# PSPManalysis 0.3.0

* Changed the manual to the bookdown version that can be opened in the Rstudio viewer	





