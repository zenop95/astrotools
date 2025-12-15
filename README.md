
# astrotools

Propagators:
* AIDA: a numerical propagator for Earth-orbiting object
* SADA: a semi-analytical propagator

## Getting started

### Prerequisites
This is a list of libraries that are required for compiling, linking and runnig the software, follow Page 31 of Manual to install them easily:
* [JSONcpp](https://github.com/open-source-parsers/jsoncpp)
* [DLib](http://dlib.net/)
* [Eigen](https://eigen.tuxfamily.org/dox/index.html)
* [CSPICE](https://naif.jpl.nasa.gov/naif/toolkit.html)\
  :open_file_folder: We suggest to install the CSPICE includes in ```/usr/local/include/cspice```, the binaries in ```/usr/local/bin```, and the libraries ```cspice.a``` and ```csupport.a``` in ```/usr/local/lib```.
* [DACE](https://github.com/dacelib/dace.git)\
  :exclamation: The Differential Algebra Computational toolbox shall be compiled with option ```-DWITH_ALGEBRAICMATRIX=ON``` to enable the algebraic matrices support.

To generate unit tests (optional) this additional dependency is required:
* [CPPunit](https://www.freedesktop.org/wiki/Software/cppunit/)

### Building and Installing
0. Download or clone the repository:
```bash
git clone "https://github.com/zenop95/astrotools.git"
```
1. Navigate to the repository folder:
```bash
cd astrotools
```
2. Run cmake and create a folder where to build the repository, e.g.:
```bash
cmake -S . -B ./_build
```
This will build the repository libraries and tools with the predefined build type (```RelWithDebInfo```) and with unit tests enabled. To select another build flavour, change compiler, or to disable unit tests you can play around with the following options:
```bash
cmake -S . -B ./_build -DCMAKE_BUILD_TYPE=Debug -DBUILD_TESTING=NO -DCMAKE_CXX_COMPILER=clang
```

3. Go to the build directory and run
```bash
cd _build
make
```
4. Install the library
```bash
sudo make install
```
Use ```sudo``` to get the required permissions to install into your system directories (usually ```/usr/local```), drop it when installing in your $HOME folder.

## Running the tests
Unit tests are located in the ```tests``` folder. To run test run from the ```_build``` folder:
```bash
make test
```
or use directly the ctest executable (which provides more options)
```bash
ctest
```

## Deployment

## Built with
* [CMake](https://cmake.org/)
