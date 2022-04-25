# Configurable Genetic Algorithm


This project provides a configurable framework of genetic algorithms.

If you are using the tool for the first time, please download or clone this branch, and run `cmake .`; `make` at the directory where the project locates.
* If you plan to install the package, please run `make install`.
* If you need to set up the install directory, please run `cmake -DCMAKE_INSTALL_PREFIX=your/path .` before installation.
* An exectuable file `main` will be generated.

<!-- After installation, you can compile your project as follow (with linking configGA library):
```
g++ $CMPL_FLAGS -o PROJECT_NAME PROJECT_NAME.cpp -configGA
```
If the library and header files are not installed as the default settings of your compiler, you may need to add additional include and libarary directory as 
```
g++ $CMPL_FLAGS -o PROJECT_NAME PROJECT_NAME.cpp -I/INCLUDE_INSTALL_PATH/ -L/LIBRARY_INSTALL_PATH -configGA
``` -->
