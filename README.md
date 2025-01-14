# Custom Mie Scattering Process for Geant4

## Getting source files
The easiest way is to create a folder for this project:
```
mkdir Norma && cd Norma
```
Clone this repository:
```
git clone -b custom_mie https://github.com/molnarmatyas/mie_with_geant
```


## Patch and build Geant4
Install the necessary packages (Ubuntu example):
```
sudo apt-get -y install bc gfortran libpcre3-dev doxygen \
xlibmesa-glu-dev libglew-dev libftgl-dev \
libmysqlclient-dev libfftw3-dev libcfitsio-dev \
graphviz-dev libavahi-compat-libdnssd-dev \
libldap2-dev python3-dev libxml2-dev libkrb5-dev \
libgsl0-dev qtwebengine5-dev dpkg-dev cmake libx11-dev libxpm-dev \
libxft-dev libxext-dev python3 libssl-dev binutils gcc g++ libxmu-dev
```

This patch is based on `geant4-11.2.0`, but may work with newer versions too.
```
wget https://github.com/Geant4/geant4/archive/refs/tags/v11.2.0.tar.gz
tar -xzf v11.2.0.tar.gz
cd geant4-11.2.0
```

Apply the [patch](patch/numeric_mie.patch):
```
git apply ../mie_with_geant/patch/numeric_mie.patch
```
or, if in this directory (make sure geant4-* is in the .gitignore):
```
git apply ../patch/numeric_mie.patch
```
Or, rather this one works REALLY, after copying it in the geant4-11.2.0 directory (and *maybe* performing `cp src/ThreadSafeWriter.cc geant4-11.2.0/source/processes/optical/src/`, `cp include/ThreadSafeWriter.hh geant4-11.2.0/source/processes/optical/include/` and `cp patch/numeric_mie_v2.patch geant4-11.2.0/`, in `cd geant4-11.2.0`):
```
patch -p1 --fuzz=3 --ignore-whitespace --merge --verbose < numeric_mie.patch
```



Now everything is ready to build Geant4:
```
cd ../
mkdir GeantBuild && cd GeantBuild
```
or:
```
mkdir geant4-install
mkdir geant4-build && cd geant4-build
```

Run CMake:
```
cmake -DCMAKE_INSTALL_PREFIX=/opt/geant4/ -DGEANT4_USE_GDML=ON -DXERCESC_ROOT_DIR=/opt/xerces-c/ -DGEANT4_INSTALL_DATA=ON -DGEANT4_USE_QT=ON -DGEANT4_USE_OPENGL_X11=ON -DGEANT4_USE_RAYTRACER_X11=ON -DGEANT4_USE_SYSTEM_CLHEP=ON -DCLHEP_INCLUDE_DIR=/opt/CLHEP/include/ -DCLHEP_LIBRARY=/opt/CLHEP/lib/ -DGEANT4_INSTALL_EXAMPLES=ON -DCLHEP_ROOT_DIR=/opt/CLHEP/ -DGEANT4_BUILD_MULTITHREADED=ON $geant_extra_args ../geant4-11.2.0
```
or:
```
cmake -DCMAKE_INSTALL_PREFIX=../geant4-install -DGEANT4_USE_GDML=ON -DGEANT4_INSTALL_DATA=ON -DGEANT4_USE_QT=ON -DGEANT4_USE_OPENGL_X11=ON -DGEANT4_USE_RAYTRACER_X11=ON -DGEANT4_USE_SYSTEM_CLHEP=ON -DCLHEP_INCLUDE_DIR=../../../GEANT/CLHEP/install/include/ -DCLHEP_LIBRARY=../../../GEANT/CLHEP/install/include/lib -DGEANT4_INSTALL_EXAMPLES=ON -DCLHEP_ROOT_DIR=../../../GEANT/CLHEP/install/ -DGEANT4_BUILD_MULTITHREADED=ON $geant_extra_args ../geant4-11.2.0
```

Use `make -jN` command. Replace `N` with the number of jobs you prefer. For the best performance use your physical core count. For example on an 8 core CPU use:
```
make -j8
```

Install Geant4:
```
sudo make install
```

If you prefer, add Geant4 source script to your `bash.bashrc`. If you skip this part, you must call this source script everytime when you want to use the modded version.
```
sudo grep -q "source /opt/geant4/bin/geant4.sh" /etc/bash.bashrc
```
or just call from `mie_with_geant` folder:
```
source geant4-install/bin/geant4.sh
```

If you make any changes in Geant4 source, repeat the `make -jN` and `sudo make install` steps. Build your simulation too.

## Build simulation
Change directory to `mie_with_geant`
```
cd mie_with_geant
mkdir build && cd build
cmake ../
make -j8
```

## Using the simulation
After a successful build of Geant4 and simulation you must need the numeric simulation's result file. Set `NUMERIC_MIE_FPATH` environment variable to the file path which you want to use, set `CELL_RADIUS_UM` environment variable to set the radius for the numerical simulation selection, then run the software:
```
NUMERIC_MIE_FPATH=diff_cross_0-180_10000pt.txt CELL_RADIUS_UM=3.0 ./Norma
```

If `NUMERIC_MIE_FPATH` is not set, it will use the original `OPMieHG` process.

### Command line arguments:
- `-m [MACRO_PATH]` - Specify a macro to run
- `-r [SEED]` - Sets a new seed
- `-b [BUBBLE_RADIUS]` - Sets the bubble radius (nm)

### Output
Currently it is using [ThreadSafeWriter](include/ThreadSafeWriter.hh) to write the calculated theta and generated theta to a file in the following format:
```
[calculated_theta], [generated_theta]
...
[calculated_theta], [generated_theta]
```
