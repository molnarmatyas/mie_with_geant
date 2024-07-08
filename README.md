# mie_with_geant
ELTE-Norma együttműködés, GEANT csoport

Most csinálok ilyet először, szóval lesz itt ám minden...


Compile:

`mkdir build`

`cd build`

`cmake -DGeant4_DIR=/opt/Geant4 ... # a Geant4 install directory helye`

`make -j8 # elerheto CPU magok szama`

Futtat:

`./exampleB1 run2.mac`
