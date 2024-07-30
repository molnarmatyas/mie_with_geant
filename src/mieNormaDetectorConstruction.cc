//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: mieNormaDetectorConstruction.cc 75117 2013-10-28 09:38:37Z gcosmo $
//
/// \file mieNormaDetectorConstruction.cc
/// \brief Implementation of the mieNormaDetectorConstruction class

#include "mieNormaDetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4SystemOfUnits.hh"
#include "G4Element.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

mieNormaDetectorConstruction::mieNormaDetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(0)
{
  DefineCommands();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

mieNormaDetectorConstruction::~mieNormaDetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* mieNormaDetectorConstruction::Construct()
{  
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
  
  // Envelope parameters
  //
  G4double env_sizeXY = 20*cm, env_sizeZ = 30*cm;
  G4Material* env_mat = nist->FindOrBuildMaterial("G4_Galactic");
   
  // Option to switch on/off checking of volumes overlaps
  //
  checkOverlaps = true;

  // ------------ Generate & Add Material Properties Table ------------
  G4double photonEnergy[] ={0.01 * eV, 0.1 * eV, 1. * eV, 10. * eV};
  const G4int nEntries = sizeof (photonEnergy) / sizeof (G4double);
  G4double absorptionLength[nEntries] ={1.*m, 1.*m, 1.*m, 1.*m};

  // ------------- Materials -------------
  G4double a, z, density;
  G4int nelements;
  // Air
  auto N = new G4Element("Nitrogen", "N", z = 7, a = 14.01 * g / mole);
  auto O = new G4Element("Oxygen", "O", z = 8, a = 16.00 * g / mole);
  auto air = new G4Material("Air", density = 1.29 * mg / cm3, nelements = 2);
  air->AddElement(N, 70. * perCent);
  air->AddElement(O, 30. * perCent);
  //
  // Water
  auto H = new G4Element("Hydrogen", "H", z = 1, a = 1.01 * g / mole);
  auto water = new G4Material("Water", density = 1.0 * g / cm3, nelements = 2);
  water->AddElement(H, 2);
  water->AddElement(O, 1);


  // ------------ Specific Properties for This Material ---------------
  G4double refractiveIndex0[nEntries] ={1.331, 1.331, 1.331, 1.331};
  G4MaterialPropertiesTable* myMPT0 = new G4MaterialPropertiesTable();
  myMPT0->AddProperty("RINDEX", photonEnergy, refractiveIndex0, nEntries); //->SetSpline(true);
  myMPT0->AddProperty("ABSLENGTH", photonEnergy, absorptionLength, nEntries); //->SetSpline(true);

  // World
  G4double world_sizeXY = 1.2*env_sizeXY;
  G4double world_sizeZ  = 1.2*env_sizeZ;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_Galactic");
  world_mat->SetMaterialPropertiesTable(myMPT0);
  
  G4Box* solidWorld =    
    new G4Box("World",                       //its name
       0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);     //its size
      
  G4LogicalVolume* logicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        "World");            //its name
                                   
  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking
                     
  //     
  // Envelope
  //  
  solidEnv =    
    new G4Box("Envelope",                    //its name
        0.5*env_sizeXY, 0.5*env_sizeXY, 0.5*env_sizeZ); //its size
      
  logicEnv =                         
    new G4LogicalVolume(solidEnv,            //its solid
                        env_mat,             //its material
                        "Envelope");         //its name

   physEnv = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicEnv,              //its logical volume
                      "Envelope",            //its name
                      logicWorld,            //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking

  // ------------ Specific Properties for This Material ---------------
  G4double refractiveIndex1[nEntries] ={1.352, 1.352, 1.352, 1.352};
  G4MaterialPropertiesTable* myMPT1 = new G4MaterialPropertiesTable();
  myMPT1->AddProperty("RINDEX", photonEnergy, refractiveIndex1, nEntries); //->SetSpline(true);
  myMPT1->AddProperty("ABSLENGTH", photonEnergy, absorptionLength, nEntries); //->SetSpline(true);

  //energies for MIE scattering 
  G4double energy_water[] = {
     1.56962*eV, 1.58974*eV, 1.61039*eV, 1.63157*eV,
     1.65333*eV, 1.67567*eV, 1.69863*eV, 1.72222*eV,
     1.74647*eV, 1.77142*eV, 1.7971 *eV, 1.82352*eV,
     1.85074*eV, 1.87878*eV, 1.90769*eV, 1.93749*eV,
     1.96825*eV, 1.99999*eV, 2.03278*eV, 2.06666*eV,
     2.10169*eV, 2.13793*eV, 2.17543*eV, 2.21428*eV,
     2.25454*eV, 2.29629*eV, 2.33962*eV, 2.38461*eV,
     2.43137*eV, 2.47999*eV, 2.53061*eV, 2.58333*eV,
     2.63829*eV, 2.69565*eV, 2.75555*eV, 2.81817*eV,
     2.88371*eV, 2.95237*eV, 3.02438*eV, 3.09999*eV,
     3.17948*eV, 3.26315*eV, 3.35134*eV, 3.44444*eV,
     3.54285*eV, 3.64705*eV, 3.75757*eV, 3.87499*eV,
     3.99999*eV, 4.13332*eV, 4.27585*eV, 4.42856*eV,
     4.59258*eV, 4.76922*eV, 4.95999*eV, 5.16665*eV,
     5.39129*eV, 5.63635*eV, 5.90475*eV, 6.19998*eV
  };
  const G4int numentries_water = sizeof(energy_water)/sizeof(G4double);

  //assume 100 times larger than the rayleigh scattering for now.
  G4double mie_water[] = {
     167024.4*nm, 158726.7*nm, 150742  *nm,
     143062.5*nm, 135680.2*nm, 128587.4*nm,
     121776.3*nm, 115239.5*nm, 108969.5*nm,
     102958.8*nm, 97200.35*nm, 91686.86*nm,
     86411.33*nm, 81366.79*nm, 76546.42*nm,
     71943.46*nm, 67551.29*nm, 63363.36*nm,
     59373.25*nm, 55574.61*nm, 51961.24*nm,
     48527.00*nm, 45265.87*nm, 42171.94*nm,
     39239.39*nm, 36462.50*nm, 33835.68*nm,
     31353.41*nm, 29010.30*nm, 26801.03*nm,
     24720.42*nm, 22763.36*nm, 20924.88*nm,
     19200.07*nm, 17584.16*nm, 16072.45*nm,
     14660.38*nm, 13343.46*nm, 12117.33*nm,
     10977.70*nm, 9920.416*nm, 8941.407*nm,
     8036.711*nm, 7202.470*nm, 6434.927*nm,
     5730.429*nm, 5085.425*nm, 4496.467*nm,
     3960.210*nm, 3473.413*nm, 3032.937*nm,
     2635.746*nm, 2278.907*nm, 1959.588*nm,
     1675.064*nm, 1422.710*nm, 1200.004*nm,
     1004.528*nm, 833.9666*nm, 686.1063*nm
  };
  assert(sizeof(mie_water) == sizeof(energy_water));

  // gforward, gbackward, forward backward ratio
  G4double mie_water_const[3]={0.99,0.99,1.0};
  //G4double mie_photon_energy[2] = {1.8 * eV, 2 * eV}; 
  //G4double mieScatteringLength[2] = {10.0 * um, 10.0 * um};
  //const G4int numentries_mie_water = sizeof(mie_photon_energy)/sizeof(G4double);
 
  //energy_water is an array containing the photon energy values.
  //mie_water is an array containing the corresponding Mie scattering phase function values.
  //num_entries is the number of entries in the arrays.
  myMPT1->AddProperty("MIEHG",energy_water, mie_water, numentries_water, false, true);

  myMPT1->AddConstProperty("MIEHG_FORWARD",mie_water_const[0]);
  myMPT1->AddConstProperty("MIEHG_BACKWARD",mie_water_const[1]);
  myMPT1->AddConstProperty("MIEHG_FORWARD_RATIO",mie_water_const[2]);

  // Shape 1
  //G4Material* shape1_mat = nist->FindOrBuildMaterial("G4_A-150_TISSUE");
  G4Material* shape1_mat = nist->FindOrBuildMaterial("G4_WATER");
  shape1_mat->SetMaterialPropertiesTable(myMPT1);
  water->SetMaterialPropertiesTable(myMPT1);
  //G4ThreeVector pos1 = G4ThreeVector(0.1*micrometer, 0.1*micrometer, 0.1*micrometer);
  G4ThreeVector pos1 = G4ThreeVector(0, 0, 5*cm);

  //G4double shape1_radius = 0.5*micrometer;
  solidShape1 = new G4Orb("Shape1", fRadius);

  logicShape1 =
      new G4LogicalVolume(solidShape1,         //its solid
          shape1_mat,          //its material
          "Shape1");           //its name
  
  physShape1 =
      new G4PVPlacement(0,                     //no rotation
          pos1,                  //at position
          logicShape1,           //its logical volume
          "Shape1",              //its name
          logicEnv,              //its mother  volume
          false,                 //no boolean operation
          0,                     //copy number
          checkOverlaps);        //overlaps checking
  
  Shape1Wrap = new G4OpticalSurface("Shape1Wrap");
  new G4LogicalBorderSurface("Shape1Wrap", physShape1, physEnv, Shape1Wrap);

  // Set Shape1 as scoring volume
  //
  fScoringVolume = logicShape1;

  //
  //always return the physical World
  //
  return physWorld;
}

void mieNormaDetectorConstruction::UpdateSphere()
{
    if (logicShape1) {
        // Remove the old physical volume of the sphere
        delete physShape1;
        delete logicShape1;
        delete solidShape1;

        // Recreate the sphere with the new radius
        solidShape1 = new G4Orb("Shape1", fRadius);

        logicShape1 =
            new G4LogicalVolume(solidShape1,         //its solid
                shape1_mat,          //its material
                "Shape1");           //its name

        physShape1 =
            new G4PVPlacement(0,                     //no rotation
                pos1,                  //at position
                logicShape1,           //its logical volume
                "Shape1",              //its name
                logicEnv,              //its mother  volume
                false,                 //no boolean operation
                0,                     //copy number
                checkOverlaps);        //overlaps checking
        
        new G4LogicalBorderSurface("Shape1Wrap", physShape1, physEnv, Shape1Wrap);
    }
}

void mieNormaDetectorConstruction::SetRadius(G4double value)
{
    fRadius = value;
    UpdateSphere();
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}


void mieNormaDetectorConstruction::DefineCommands()
{
    // Define /B5/detector command directory using generic messenger class
    fMessenger = new G4GenericMessenger(this,
        "/sphere/",
        "Radius of the sphere.");

    // radius command
    auto& radiusCmd
        = fMessenger->DeclareMethodWithUnit("radius", "micrometer",
            &B1DetectorConstruction::SetRadius,
            "Set radius of the sphere.");
    radiusCmd.SetParameterName("radius", true);
    radiusCmd.SetRange("radius>=0.");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
