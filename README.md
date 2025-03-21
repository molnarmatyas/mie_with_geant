# Custom Mie Scattering Process for Geant4

## Getting source files
The easiest way is to create a folder for this project:
```
mkdir Norma && cd Norma
```
Clone this repository:
```
git clone -b 3D_complete_custom_mie https://github.com/molnarmatyas/mie_with_geant
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

Apply the [patch](patch/numeric_mie_v4.patch):
Before applying, one needs to copy the patch to the `geant4-11.2.0` directory.
You should be on the level where you see both `geant4-11.2.0` directory as well as `mie_with_geant`.
```
cp mie_with_geant/patch/numeric_mie_v4.patch geant4-11.2.0/
cd geant4-11.2.0
patch -p1 --fuzz=3 --ignore-whitespace --merge --verbose < numeric_mie_v4.patch
```


Now everything is ready to build Geant4. Be one level above the `geant4-11.2.0` directory:
```
mkdir GeantBuild && mkdir GeantInstall && cd GeantBuild
```

Run CMake:
```
cmake -DCMAKE_INSTALL_PREFIX=../geant4-install -DGEANT4_USE_GDML=OFF -DGEANT4_INSTALL_DATA=ON -DGEANT4_USE_QT=ON -DGEANT4_USE_OPENGL_X11=ON -DGEANT4_USE_RAYTRACER_X11=ON -DGEANT4_USE_SYSTEM_CLHEP=OFF -DGEANT4_INSTALL_EXAMPLES=ON -DGEANT4_BUILD_MULTITHREADED=ON $geant_extra_args ../geant4-11.2.0
```
alternatively, if you have CLHEP in /opt, for example in Docker or other VM:
```
cmake -DCMAKE_INSTALL_PREFIX=/opt/geant4/ -DGEANT4_USE_GDML=ON -DXERCESC_ROOT_DIR=/opt/xerces-c/ -DGEANT4_INSTALL_DATA=ON -DGEANT4_USE_QT=ON -DGEANT4_USE_OPENGL_X11=ON -DGEANT4_USE_RAYTRACER_X11=ON -DGEANT4_USE_SYSTEM_CLHEP=ON -DCLHEP_INCLUDE_DIR=/opt/CLHEP/include/ -DCLHEP_LIBRARY=/opt/CLHEP/lib/ -DGEANT4_INSTALL_EXAMPLES=ON -DCLHEP_ROOT_DIR=/opt/CLHEP/ -DGEANT4_BUILD_MULTITHREADED=ON $geant_extra_args ../geant4-11.2.0
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
sudo grep -q "source /PATH/TO/FOLDER/GeantInstall/bin/geant4.sh" /etc/bash.bashrc
```
or just call from `mie_with_geant` folder:
```
source /PATH/TO/FOLDER/GeantInstall/bin/geant4.sh
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
After a successful build of Geant4 and simulation you must need the numeric simulation's result file. Set `NUMERIC_MIE_FPATH` environment variable to the file path which you want to use. Set then the  `CELL_RADIUS_UM` environment variable to the radius of the cell for the numerical simulation selection, then say:
```
NUMERIC_MIE_FPATH=diff_cross_0-180_10000pt.txt CELL_RADIUS_UM=3.0 ./Norma -m run2.mac
```

If `NUMERIC_MIE_FPATH` is not set, it will use the original `OPMieHG` process.

### Command line arguments:
- `-m [MACRO_PATH]` - Specify a macro to run
- `-r [SEED]` - Sets a new seed
- `-b [BUBBLE_RADIUS]` - Sets the bubble radius (um)
- `-t [THREADS]` - Sets number of threads to use

### Output
Currently it is using [ThreadSafeWriter](include/ThreadSafeWriter.hh) to write the calculated theta and generated theta to a file in the following format:
```
[calculated_theta], [generated_theta], [distance on 2D plane], [detector hit position X],  [detector hit position Y],  [detector hit position Z], [momentum direction at the cell - alpha angle (not used anymore)], [manually calculated alpha angle (not used anymore)], [phi angle], [detector number: 0 - CCD, 1 - LA, 2 - HA] 
...
[calculated_theta], [generated_theta], [distance on 2D plane], [detector hit position X],  [detector hit position Y],  [detector hit position Z], [momentum direction at the cell - alpha angle (not used anymore)], [manually calculated alpha angle (not used anymore)], [phi angle], [detector number: 0 - CCD, 1 - LA, 2 - HA] 
```
