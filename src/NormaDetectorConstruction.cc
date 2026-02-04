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
/// \file Norma/src/NormaDetectorConstruction.cc
/// \brief Implementation of the NormaDetectorConstruction class
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "NormaDetectorConstruction.hh"
#include "NormaDetectorMessenger.hh"

#include "G4Box.hh"
#include "G4Element.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4OpticalSurface.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include <G4Orb.hh>


#include "CADMesh.hh"
#include "G4VisAttributes.hh"
#include "G4Tubs.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
NormaDetectorConstruction::NormaDetectorConstruction() : G4VUserDetectorConstruction()
{
  // create a messenger for this class
  fDetectorMessenger = new NormaDetectorMessenger(this);
  fBubble_additional_offset = G4ThreeVector(0, 0, 0);
}

void NormaDetectorConstruction::SetCellOffset(G4double x, G4double y, G4double z) 
{
  fBubble_additional_offset.set(x, y, z);
}

void NormaDetectorConstruction::UpdateGeometry() 
{
  G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
NormaDetectorConstruction::~NormaDetectorConstruction()
{
  delete fDetectorMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VPhysicalVolume *NormaDetectorConstruction::Construct()
{
  G4bool checkOverlaps = true;
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

  // Get nist material manager
  G4NistManager *nist = G4NistManager::Instance();
  G4Material *sil = nist->FindOrBuildMaterial("G4_Si");
  G4Material *copper = nist->FindOrBuildMaterial("G4_Cu");
  G4Material *alu_oxide = nist->FindOrBuildMaterial("G4_ALUMINUM_OXIDE");

  // ------------ Generate & Add Material Properties Table ------------
  //
  std::vector<G4double> photonEnergy ={0.01 * eV, 0.1 * eV, 1. * eV, 10. * eV};
  const G4int nEntries = 4;
  std::vector<G4double> absorptionLength ={1.*m, 1.*m, 1.*m, 1.*m};
  std::vector<G4double> refractiveIndexWorld ={1.000277, 1.000277, 1.000277, 1.000277}; // instead of 1.331

  // Water
  // Material properties can be added as arrays. However, in this case it is
  // up to the user to make sure both arrays have the same number of elements.
  auto myMPT1 = new G4MaterialPropertiesTable();

  std::vector<G4double> photonEnergyMie = {1.8 *eV, 1.884*eV, 1.9*eV, 2.0*eV};
  std::vector<G4double> refractiveIndexCell = {1.592, 1.592, 1.592, 1.592};
  myMPT1->AddProperty("RINDEX", photonEnergyMie, refractiveIndexCell, nEntries); //->SetSpline(true);
  myMPT1->AddProperty("ABSLENGTH", photonEnergyMie, absorptionLength, false, true);

  // Values can be added to the material property table individually.
  // With this method, spline interpolation cannot be set. Arguments
  // createNewKey and spline both take their default values of false.
  // Need to specify the number of entries (1) in the arrays, as an argument
  // to AddProperty.

  // Check that group velocity is calculated from RINDEX
  if (myMPT1->GetProperty("RINDEX")->GetVectorLength() != myMPT1->GetProperty("GROUPVEL")->GetVectorLength())
  {
    G4ExceptionDescription ed;
    ed << "Error calculating group velocities. Incorrect number of entries "
      "in group velocity material property vector.";
    G4Exception("Norma::NormaDetectorConstruction", "Norma001", FatalException, ed);
  }

  std::vector<G4double> energy_water = {
    1.56962 * eV, 1.58974 * eV, 1.61039 * eV, 1.63157 * eV, 1.65333 * eV, 1.67567 * eV, 1.69863 * eV, 1.72222 * eV,
    1.74647 * eV, 1.77142 * eV, 1.7971 * eV,  1.82352 * eV, 1.85074 * eV, 1.87878 * eV, 1.90769 * eV, 1.93749 * eV,
    1.96825 * eV, 1.99999 * eV, 2.03278 * eV, 2.06666 * eV, 2.10169 * eV, 2.13793 * eV, 2.17543 * eV, 2.21428 * eV,
    2.25454 * eV, 2.29629 * eV, 2.33962 * eV, 2.38461 * eV, 2.43137 * eV, 2.47999 * eV, 2.53061 * eV, 2.58333 * eV,
    2.63829 * eV, 2.69565 * eV, 2.75555 * eV, 2.81817 * eV, 2.88371 * eV, 2.95237 * eV, 3.02438 * eV, 3.09999 * eV,
    3.17948 * eV, 3.26315 * eV, 3.35134 * eV, 3.44444 * eV, 3.54285 * eV, 3.64705 * eV, 3.75757 * eV, 3.87499 * eV,
    3.99999 * eV, 4.13332 * eV, 4.27585 * eV, 4.42856 * eV, 4.59258 * eV, 4.76922 * eV, 4.95999 * eV, 5.16665 * eV,
    5.39129 * eV, 5.63635 * eV, 5.90475 * eV, 6.19998 * eV};

  // Rayleigh scattering length is calculated by G4OpRayleigh

  // Mie: assume 100 times larger than the rayleigh scattering
  std::vector<G4double> mie_water = {
    167024.4 * 0.001 * nm, 158726.7 * 0.001 * nm, 150742 * 0.001 * nm,	 143062.5 * 0.001 * nm,
    135680.2 * 0.001 * nm, 128587.4 * 0.001 * nm, 121776.3 * 0.001 * nm, 115239.5 * 0.001 * nm,
    108969.5 * 0.001 * nm, 102958.8 * 0.001 * nm, 97200.35 * 0.001 * nm, 91686.86 * 0.001 * nm,
    86411.33 * 0.001 * nm, 81366.79 * 0.001 * nm, 76546.42 * 0.001 * nm, 71943.46 * 0.001 * nm,
    67551.29 * 0.001 * nm, 63363.36 * 0.001 * nm, 59373.25 * 0.001 * nm, 55574.61 * 0.001 * nm,
    51961.24 * 0.001 * nm, 48527.00 * 0.001 * nm, 45265.87 * 0.001 * nm, 42171.94 * 0.001 * nm,
    39239.39 * 0.001 * nm, 36462.50 * 0.001 * nm, 33835.68 * 0.001 * nm, 31353.41 * 0.001 * nm,
    29010.30 * 0.001 * nm, 26801.03 * 0.001 * nm, 24720.42 * 0.001 * nm, 22763.36 * 0.001 * nm,
    20924.88 * 0.001 * nm, 19200.07 * 0.001 * nm, 17584.16 * 0.001 * nm, 16072.45 * 0.001 * nm,
    14660.38 * 0.001 * nm, 13343.46 * 0.001 * nm, 12117.33 * 0.001 * nm, 10977.70 * 0.001 * nm,
    9920.416 * 0.001 * nm, 8941.407 * 0.001 * nm, 8036.711 * 0.001 * nm, 7202.470 * 0.001 * nm,
    6434.927 * 0.001 * nm, 5730.429 * 0.001 * nm, 5085.425 * 0.001 * nm, 4496.467 * 0.001 * nm,
    3960.210 * 0.001 * nm, 3473.413 * 0.001 * nm, 3032.937 * 0.001 * nm, 2635.746 * 0.001 * nm,
    2278.907 * 0.001 * nm, 1959.588 * 0.001 * nm, 1675.064 * 0.001 * nm, 1422.710 * 0.001 * nm,
    1200.004 * 0.001 * nm, 1004.528 * 0.001 * nm, 833.9666 * 0.001 * nm, 686.1063 * 0.001 * nm};


  // Mie: gforward, gbackward, forward backward ratio
  G4double mie_water_const[3] = {mieFg, 0.99, 0.99};

  myMPT1->AddProperty("MIEHG", energy_water, mie_water, false, true);
  myMPT1->AddConstProperty("MIEHG_FORWARD", mie_water_const[0]);
  myMPT1->AddConstProperty("MIEHG_BACKWARD", mie_water_const[1]);
  myMPT1->AddConstProperty("MIEHG_FORWARD_RATIO", mie_water_const[2]);

  G4cout << "Water G4MaterialPropertiesTable:" << G4endl;

  water->SetMaterialPropertiesTable(myMPT1);

  // Set the Birks Constant for the Water scintillator
  water->GetIonisation()->SetBirksConstant(0.126 * mm / MeV);

  // Air

  auto myMPT0 = new G4MaterialPropertiesTable();
  myMPT0->AddProperty("RINDEX", photonEnergy, refractiveIndexWorld);
  myMPT0->AddProperty("ABSLENGTH", photonEnergy, absorptionLength, false, true);

  /*
     myMPT0->AddProperty("MIEHG", energy_water, mie_water, false, true);
     myMPT0->AddConstProperty("MIEHG_FORWARD", mie_water_const[0]);
     myMPT0->AddConstProperty("MIEHG_BACKWARD", mie_water_const[1]);
     myMPT0->AddConstProperty("MIEHG_FORWARD_RATIO", mie_water_const[2]);
   */

  G4cout << "Air G4MaterialPropertiesTable:" << G4endl;

  air->SetMaterialPropertiesTable(myMPT0);


  // Mirror material
  auto mirrorMaterial = nist->FindOrBuildMaterial("G4_Al");

  G4MaterialPropertiesTable* myMPT3 = new G4MaterialPropertiesTable();
  std::vector<G4double> photonEnergyMirror = {1.8*eV, 1.85*eV, 1.9*eV, 1.95*eV };
  mirrorMaterial->SetMaterialPropertiesTable(myMPT3);

  // Steel
  G4Material* steel = G4NistManager::Instance()->FindOrBuildMaterial("G4_STAINLESS-STEEL");
  std::vector<G4double> steelreflectivity = {0.7, 0.7, 0.7 , 0.7}; // Approximate steel reflectivity
  std::vector<G4double> steelabsorption = {0.1 * mm, 0.1 * mm, 0.1 * mm, 0.1 * mm}; // Very high absorption

  G4MaterialPropertiesTable* steelMPT = new G4MaterialPropertiesTable();
  //steelMPT->AddProperty("REFLECTIVITY", photonEnergyMirror, steelreflectivity, nEntries);
  //steelMPT->AddProperty("ABSLENGTH", photonEnergyMirror, steelabsorption, nEntries);
  steel->SetMaterialPropertiesTable(steelMPT);


  // Shielding
  auto shieldMaterial = new G4Material("BlackPlastic", 0.94*g/cm3, 1);
  shieldMaterial->AddMaterial(nist->FindOrBuildMaterial("G4_POLYETHYLENE"), 1.0);
  G4MaterialPropertiesTable* myMPT4 = new G4MaterialPropertiesTable();


  std::vector<G4double> absLengthShield = { 0.0000001*mm, 0.0000001*mm, 0.0000001*mm, 0.0000001*mm };
  myMPT4->AddProperty("ABSLENGTH", photonEnergyMirror, absLengthShield, false, false);


  shieldMaterial->SetMaterialPropertiesTable(myMPT4);


  // Aspherical lens - B270 material, https://www.schott.com/en-us/products/b-270-p1000313/technical-details
  // https://refractiveindex.info/?shelf=specs&book=SCHOTT-misc&page=B270
  G4Material* lensMaterial = new G4Material("B270", 2.56*g/cm3, 1);
  lensMaterial->AddMaterial(nist->FindOrBuildMaterial("G4_Pyrex_Glass"), 1.0); // should do it for now
  G4MaterialPropertiesTable* myMPT5 = new G4MaterialPropertiesTable();

  std::vector<G4double> rindexLens = {1.5203, 1.5203, 1.5203, 1.5203};
  myMPT5->AddProperty("RINDEX", photonEnergyMirror, rindexLens, nEntries);

  std::vector<G4double> absLengthLens = {100.0*mm, 100.0*mm, 100.0*mm, 100.0*mm};
  myMPT5->AddProperty("ABSLENGTH", photonEnergyMirror, absLengthLens, false, false);

  /*
     std::vector<G4double> transmittanceB270 = {0.99,  0.99,  0.99,  0.99};//{0.91320,  0.91320,  0.91320,  0.91320};
     myMPT5->AddProperty("TRANSMITTANCE", photonEnergyMirror, transmittanceB270, nEntries);
   */

  lensMaterial->SetMaterialPropertiesTable(myMPT5);


  // LB1258-A lens - N-BK7 material
  // https://www.schott.com/shop/medias/schott-datasheet-n-bk7-eng.pdf?context=bWFzdGVyfHJvb3R8NjkxODAwfGFwcGxpY2F0aW9uL3BkZnxoZTUvaDM4Lzg4MTAzMTYxMDM3MTAucGRmfGJjNmI4ZjFmY2Q1NjMxMTE0MjkzMTUwOGRmMTUzOTg2NWJjZTgzMjA0OTc2NTNiMThjN2RhMjI4NGZmMWM4MWU
  G4Material* lensMaterialNBK7 = new G4Material("NBK7", 2.51*g/cm3, 1);
  lensMaterialNBK7->AddMaterial(nist->FindOrBuildMaterial("G4_Pyrex_Glass"), 1.0); // should do it for now
  G4MaterialPropertiesTable* MPT_nbk7 = new G4MaterialPropertiesTable();

  std::vector<G4double> rindexLensNBK7 = {1.51680, 1.51680, 1.51680, 1.51680};
  MPT_nbk7->AddProperty("RINDEX", photonEnergyMirror, rindexLensNBK7, nEntries);

  std::vector<G4double> absLengthLensNBK7 = {100.0*mm, 100.0*mm, 100.0*mm, 100.0*mm};
  MPT_nbk7->AddProperty("ABSLENGTH", photonEnergyMirror, absLengthLensNBK7, false, false);

  /*
     std::vector<G4double> transmittanceB270 = {0.99,  0.99,  0.99,  0.99};//{0.91320,  0.91320,  0.91320,  0.91320};
     myMPT5->AddProperty("TRANSMITTANCE", photonEnergyMirror, transmittanceB270, nEntries);
   */

  lensMaterialNBK7->SetMaterialPropertiesTable(MPT_nbk7);

  //OD 0.5 ND filter - NG4 material
  //https://www.sydor.com/wp-content/uploads/SCHOTT-NG4-Neutral-Density-Filter.pdf
  G4Material* lensMaterialNG4 = new G4Material("NG4", 2.43*g/cm3, 1);
  lensMaterialNG4->AddMaterial(nist->FindOrBuildMaterial("G4_Pyrex_Glass"), 1.0); // should do it for now
  G4MaterialPropertiesTable* MPT_ng4 = new G4MaterialPropertiesTable();

  std::vector<G4double> rindexLensNG4 = {1.510, 1.510, 1.510, 1.510};
  MPT_ng4->AddProperty("RINDEX", photonEnergyMirror, rindexLensNG4, nEntries);

  std::vector<G4double> absLengthLensNG4 = {100.0*mm, 100.0*mm, 100.0*mm, 100.0*mm};
  MPT_ng4->AddProperty("ABSLENGTH", photonEnergyMirror, absLengthLensNG4, false, false);

  lensMaterialNG4->SetMaterialPropertiesTable(MPT_ng4);

  // Helma flow cell - made of quartz glass
  G4Material* flowcellMaterial = nist->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
  G4MaterialPropertiesTable* myMPT6 = new G4MaterialPropertiesTable();

  std::vector<G4double> rindexQuartz = {1.4585, 1.4585, 1.4585, 1.4585};
  myMPT6->AddProperty("RINDEX", photonEnergyMirror, rindexQuartz, nEntries);

  std::vector<G4double> absLengthQuartz = {100.0*mm, 100.0*mm, 100.0*mm, 100.0*mm};
  myMPT6->AddProperty("ABSLENGTH", photonEnergyMirror, absLengthQuartz, false, false);

  flowcellMaterial->SetMaterialPropertiesTable(myMPT6);

  // Physiological saline solution
  G4Material* saltwater = nist->FindOrBuildMaterial("G4_WATER");
  G4MaterialPropertiesTable* saline_MPT = new G4MaterialPropertiesTable();

  //std::vector<G4double> rindexSaline = {1.34, 1.34, 1.34, 1.34}; // FIXME more precise!
  std::vector<G4double> rindexSaline = {1.3309, 1.3309, 1.3309, 1.3309}; //new values, more precise
  saline_MPT->AddProperty("RINDEX", photonEnergyMirror, rindexSaline, nEntries);

  std::vector<G4double> absLengthSaline = {100.0*mm, 100.0*mm, 100.0*mm, 100.0*mm};
  saline_MPT->AddProperty("ABSLENGTH", photonEnergyMirror, absLengthSaline, false, false);

  saltwater->SetMaterialPropertiesTable(saline_MPT);

  // OD 0.5 ND filter; partially based on https://www.3doptix.com/catalog/optics/filter/edmund-optics/88-277
  G4Material* NDglass = nist->FindOrBuildMaterial("G4_Pyrex_Glass");
  G4MaterialPropertiesTable* ND_MPT = new G4MaterialPropertiesTable();

  std::vector<G4double> rindexND = {1.458448271212196, 1.458448271212196, 1.458448271212196, 1.458448271212196};
  ND_MPT->AddProperty("RINDEX", photonEnergyMirror, rindexND, nEntries);

  std::vector<G4double> absLengthND = {5.0*m, 5.0*m, 5.0*m, 5.0*m};
  ND_MPT->AddProperty("ABSLENGTH", photonEnergyMirror, absLengthND, false, false);

  std::vector<G4double> transmittanceND = {100.0*mm, 100.0*mm, 100.0*mm, 100.0*mm};
  ND_MPT->AddProperty("TRANSMITTANCE", photonEnergyMirror, transmittanceND, false, false);

  NDglass->SetMaterialPropertiesTable(ND_MPT);

  // PMMA for injector, catcher tube and flowcella
  G4Material* PMMA = nist->FindOrBuildMaterial("G4_PLEXIGLASS");

  // ruby for capillary based on https://en.wikipedia.org/wiki/Corundum and https://en.wikipedia.org/wiki/Ruby
  auto ruby = new G4Material("Ruby", density = 4.02 * g / cm3, nelements = 2);
  auto Al = new G4Element("Aluminium", "Al", z = 13, a = 26.98 * g / mole);
  ruby -> AddElement(Al, 2);
  ruby -> AddElement(O, 3);

  G4MaterialPropertiesTable* ruby_MPT = new G4MaterialPropertiesTable();

  std::vector<G4double> rindexRUBY = {1.768, 1.768, 1.768, 1.768};
  ruby_MPT->AddProperty("RINDEX", photonEnergyMirror, rindexRUBY, nEntries);

  ruby->SetMaterialPropertiesTable(ruby_MPT);

  // salty water for flowcell-front- and backsheat
  G4Material* saltwater_frontbacksheat = nist->FindOrBuildMaterial("G4_WATER");
  G4MaterialPropertiesTable* saline_frontbacksheat_MPT = new G4MaterialPropertiesTable();

  std::vector<G4double> rindexSaline_frontbacksheat = {1.332, 1.332, 1.332, 1.332};
  saline_MPT->AddProperty("RINDEX", photonEnergyMirror, rindexSaline_frontbacksheat, nEntries);

  saline_MPT->AddProperty("ABSLENGTH", photonEnergyMirror, absLengthSaline, false, false);

  saltwater->SetMaterialPropertiesTable(saline_frontbacksheat_MPT);





  // ------------- Volumes --------------
  //
  // The world
  auto world_box = new G4Box("World", fWorld_x, fWorld_y, fWorld_z);
  auto world_log = new G4LogicalVolume(world_box, air, "World");
  G4VPhysicalVolume *world_phys =
    new G4PVPlacement(nullptr, G4ThreeVector(), world_log, "world", nullptr, false, 0, checkOverlaps);

  /*
  // The Bubble

  G4double Pz[11] = {-fBubble_r * 1.0, -fBubble_r * 0.8, -fBubble_r * 0.6, -fBubble_r * 0.4, -fBubble_r * 0.2, 0,
  fBubble_r * 0.2,	 fBubble_r * 0.4,  fBubble_r * 0.6,	 fBubble_r * 0.8,  fBubble_r * 1.0};
  G4double Prin[11] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
  G4double Prout[11] = {fBubble_r * 0.01, fBubble_r * 0.25, fBubble_r * 0.55, fBubble_r * 0.8,
  fBubble_r * 0.9,	fBubble_r * 0.99, fBubble_r * 0.92, fBubble_r * 0.8,
  fBubble_r * 0.5,	fBubble_r * 0.3,  fBubble_r * 0.05};
  std::cout << fBubble_r / CLHEP::um << std::endl;
  for (int i = 0; i < 11; i++)
  {
  std::cout << Pz[i] / CLHEP::um << " " << Prin[i] / CLHEP::um << " " << Prout[i] / CLHEP::um << std::endl;
  }
  auto bubbleWP = new G4Polycone("Bubble", 0., 360. * deg, 11, Pz, Prin, Prout);
  auto bubbleWP_log = new G4LogicalVolume(bubbleWP, water, "Bubble");

  auto bubbleW = new G4Orb("Bubble", fBubble_r);
  auto bubbleW_log = new G4LogicalVolume(bubbleW, water, "Bubble");

  if (isPolycone)
  {
  bubble_phys = new G4PVPlacement(nullptr, G4ThreeVector(14.09 * mm, 96.2425 * mm, -137.51 * mm), bubbleWP_log,
  "Bubble_dis_bnd_proc", world_log, false, 0);
  }
  else
  {
  bubble_phys = new G4PVPlacement(nullptr, G4ThreeVector(14.09 * mm, 96.2425 * mm, -137.51 * mm), bubbleW_log,
  "Bubble_dis_bnd_proc", world_log, false, 0);
  }
   */

  // screen
  double box_position = 395;
  auto screen_box = new G4Box("Screen", fScreen_x, fScreen_y, fScreen_z);
  auto screen_log = new G4LogicalVolume(screen_box, sil, "Screen");
  G4VPhysicalVolume *screen_phys = new G4PVPlacement(nullptr, G4ThreeVector(box_position * mm, 0.0 * mm, 0.0 * mm),
      screen_log, "Screen", world_log, false, 0);

  auto screenB_box = new G4Box("ScreenB", fScreen_x, fScreen_y, fScreen_z);
  auto screenB_log = new G4LogicalVolume(screenB_box, sil, "ScreenB");
  G4VPhysicalVolume *screenB_phys = new G4PVPlacement(nullptr, G4ThreeVector(-box_position * mm, 0.0 * mm, 0.0 * mm),
      screenB_log, "ScreenB", world_log, false, 0);

  auto screenL_box = new G4Box("ScreenL", fScreen_y, fScreen_y, fScreen_x);
  auto screenL_log = new G4LogicalVolume(screenL_box, sil, "ScreenL");
  G4VPhysicalVolume *screenL_phys = new G4PVPlacement(nullptr, G4ThreeVector(0.0 * mm, 0.0 * mm, box_position * mm),
      screenL_log, "ScreenL", world_log, false, 0);

  auto screenR_box = new G4Box("ScreenR", fScreen_y, fScreen_y, fScreen_x);
  auto screenR_log = new G4LogicalVolume(screenR_box, sil, "ScreenR");
  G4VPhysicalVolume *screenR_phys = new G4PVPlacement(nullptr, G4ThreeVector(0.0 * mm, 0.0 * mm, -box_position * mm),
      screenR_log, "ScreenR", world_log, false, 0);

  auto screenU_box = new G4Box("ScreenU", fScreen_y, fScreen_x, fScreen_y);
  auto screenU_log = new G4LogicalVolume(screenU_box, sil, "ScreenU");
  G4VPhysicalVolume *screenU_phys = new G4PVPlacement(nullptr, G4ThreeVector(0.0 * mm, box_position * mm, 0.0 * mm),
      screenU_log, "ScreenU", world_log, false, 0);

  auto screenD_box = new G4Box("ScreenD", fScreen_y, fScreen_x, fScreen_y);
  auto screenD_log = new G4LogicalVolume(screenD_box, sil, "ScreenD");
  G4VPhysicalVolume *screenD_phys = new G4PVPlacement(nullptr, G4ThreeVector(0.0 * mm, -box_position * mm, 0.0 * mm),
      screenD_log, "ScreenD", world_log, false, 0);


  //3D Modell load
  auto mesh = CADMesh::TessellatedMesh::FromOBJ("./I5R_GEANT_20260119_1.obj");
  G4cout << " MESH NAME: " << mesh->GetFileName() << G4endl;
  mesh->SetScale(1.0);
  std::vector<G4VSolid*> solids; // = mesh->GetSolids();

  //Simplified
  solids.push_back(mesh->GetSolid("beam_splitter"));
  solids.push_back(mesh->GetSolid("mirror"));
  solids.push_back(mesh->GetSolid("COHERENT_MINI-701L-660S"));
  solids.push_back(mesh->GetSolid("Flowcell"));
  solids.push_back(mesh->GetSolid("ACL12708U"));
  solids.push_back(mesh->GetSolid("GS3-U3-23S6M-C_sensor"));
  solids.push_back(mesh->GetSolid("BST04_BeamSplitter"));
  solids.push_back(mesh->GetSolid("Direct_beam_stop_2"));
  solids.push_back(mesh->GetSolid("vbpw34s_1"));
  solids.push_back(mesh->GetSolid("vbpw34s_2"));
  solids.push_back(mesh->GetSolid("LB1258-A"));
  solids.push_back(mesh->GetSolid("LA_mirror"));
  solids.push_back(mesh->GetSolid("LA_HA_separator")); 
  //solids.push_back(mesh->GetSolid("Direct_beam_stop")); //currently disabled
  solids.push_back(mesh->GetSolid("Direct_beam_stop_0.75"));
  solids.push_back(mesh->GetSolid("HA_mirror"));
  //solids.push_back(mesh->GetSolid("GS3-U3-23S6M-C_sensor_housing_PRIM"));
  //solids.push_back(mesh->GetSolid("shield"));
  solids.push_back(mesh->GetSolid("Saltywater"));

  //Complete
  solids.push_back(mesh->GetSolid("lense_outer_housing"));
  solids.push_back(mesh->GetSolid("filter_adapter"));
  solids.push_back(mesh->GetSolid("ND_filter_housing"));
  solids.push_back(mesh->GetSolid("BST04_BeamSplitter_housing"));
  solids.push_back(mesh->GetSolid("LA_HA_housing"));
  solids.push_back(mesh->GetSolid("LA_HA_holder"));
  solids.push_back(mesh->GetSolid("ND_filter"));
  //solids.push_back(mesh->GetSolid("LA_HA_separator_underpart"));
  //solids.push_back(mesh->GetSolid("HA_mirror_underpart"));
  //solids.push_back(mesh->GetSolid("LA_mirror_underpart"));
  //solids.push_back(mesh->GetSolid("BST04_BeamSplitter_underpart"));
  solids.push_back(mesh->GetSolid("GS3-U3-23S6M-C_sensor_housing"));
  solids.push_back(mesh->GetSolid("camera_outer_window"));
  solids.push_back(mesh->GetSolid("sensor_package_window"));
  solids.push_back(mesh->GetSolid("vbpw34s_1_sensor"));
  solids.push_back(mesh->GetSolid("vbpw34s_2_sensor"));
  solids.push_back(mesh->GetSolid("injector"));
  solids.push_back(mesh->GetSolid("catcher_tube"));
  solids.push_back(mesh->GetSolid("capillary"));
  solids.push_back(mesh->GetSolid("Flowcellwater-backsheath"));
  solids.push_back(mesh->GetSolid("Flowcellwater-frontsheath"));
  //solids.push_back(mesh->GetSolid("BeamProfiler")); // FIXME only for beam profile testing

  std::vector<G4LogicalVolume*> argosz_log(solids.size());
  std::vector<G4VPhysicalVolume*> argosz_phys(solids.size());

  std::vector<G4Material*> argosz_mat(solids.size());

  /* 0 -> 16
     Simplified 3D model
     1  solid name: beam_splitter
     2  solid name: mirror
     3  solid name: COHERENT_MINI-701L-660S
     4  solid name: Flowcell
     5  solid name: ACL12708U
     6  solid name: GS3-U3-23S6M-C_sensor
     7  solid name: BST04_BeamSplitter
     8  solid name: Direct_beam_stop_2
     9  solid name: vbpw34s_1
     10 solid name: vbpw34s_2
     11 solid name: LB1258-A
     12 solid name: LA_Mirror
     13 solid name: LA_HA_mirror
     14 solid name: Direct_beam_stop
     15 solid name: HA_mirror
     16 solid name: Saltywater
   */
  /*
     Complete 3D model
     0  solid name: ACL12708U
     1  solid name: lense_outer_housing
     2  solid name: LB1258-A
     3  solid name: COHERENT_MINI-701L-660S
     4  solid name: sensor_shield_1
     5  solid name: half_shield
     6  solid name: BST04_BeamSplitter_housing
     7  solid name: LA_HA_housing
     8  solid name: LA_HA_holder
     9  solid name: GS3-U3-23S6M-C_sensor
     10  solid name: LA_HA_mirror
     11  solid name: half_shield_2
     12  solid name: BST04_BeamSplitter
     13  solid name: Direct_beam_stop
     14  solid name: beam_splitter_1
     15  solid name: mirror_1
     16  solid name: Flowcell
     17  solid name: HA_mirror
     18  solid name: LA_mirror
     19  solid name: vbpw34s_1_sensor
     20  solid name: vbpw34s_2_sensor
     21  solid name: Direct_beam_stop_2
     22  solid name: vbpw34s_1
     23  solid name: vbpw34s_2
     24  solid name: GS3-U3-23S6M-C_sensor_housing
     25  solid name: injector
     26  solid name: catcher_tube
     27  solid name: capillary
     28  solid name: Flowcellwater-backsheath
     29  solid name: Flowcellwater-frontsheath
     30  solid name: Saltywater
   */


  // beam_splitter
  argosz_mat[0] = mirrorMaterial;
  // mirror
  argosz_mat[1] = mirrorMaterial;
  // COHERENT_MINI-701L-660S
  argosz_mat[2] = mirrorMaterial;
  // Flowcell
  argosz_mat[3] = PMMA; 
  // ACL12708U
  argosz_mat[4] = lensMaterial;
  // GS3-U3-23S6M-C_sensor
  argosz_mat[5] = sil;
  // BST04_BeamSplitter
  argosz_mat[6] = lensMaterial;
  // Direct_beam_stop_2
  argosz_mat[7] = shieldMaterial;
  // vbpw34s_1
  argosz_mat[8] = lensMaterial;
  // vbpw34s_2
  argosz_mat[9] = lensMaterial;
  // LB1258-A
  argosz_mat[10] = lensMaterial;
  // LA_mirror
  argosz_mat[11] = steel; //mirrorMaterial;
  // LA_HA_separator_
  argosz_mat[12] = shieldMaterial;
  // Direct_beam_stop_0.75
  argosz_mat[13] = shieldMaterial;
  // HA_mirror
  argosz_mat[14] = steel; //mirrorMaterial;
  // Saltywater
  argosz_mat[15] = saltwater;
  // lense_outer_housing
  argosz_mat[16] = lensMaterial;
  // filter_adapter
  argosz_mat[17] = shieldMaterial;
  // ND_filter_housing
  argosz_mat[18] = shieldMaterial;
  // BST04_BeamSplitter_housing
  argosz_mat[19] = shieldMaterial;
  // LA_HA_housing
  argosz_mat[20] = shieldMaterial;
  // LA_HA_holder
  argosz_mat[21] = shieldMaterial;
  // ND_filter
  argosz_mat[22] = NDglass;
  // BST04_BeamSplitter_underpart
  //argosz_mat[23] = shieldMaterial;
  // GS3-U3-23S6M-C_sensor_housing
  argosz_mat[23] = shieldMaterial;
  // camera_outer_window
  argosz_mat[24] = lensMaterial;
  // sensor_package_window
  argosz_mat[25] = lensMaterial;
  // vbpw34s_1_sensor001
  argosz_mat[26] = sil;
  // vbpw34s_2_sensor001
  argosz_mat[27] = sil;
  // BeamProfiler (CCD also) FIXME only for beam profile testing
  //argosz_mat[32] = air;
  // injector
  argosz_mat[28] = PMMA;
  // catcher_tube
  argosz_mat[29] = PMMA;
  // capillary
  argosz_mat[30] = ruby;
  // Flowcellwater-backsheath
  argosz_mat[31] = saltwater_frontbacksheat;
  // Flowcellwater-frontsheath
  argosz_mat[32] = saltwater_frontbacksheat;

  int isolid = 0;
  G4double dbshift = -0.050 * mm; // shift direct beam stop to make CCD image symmetrical
  for (auto solid : solids)
  {
    /*
       if(isolid == 14 || isolid == 16){
       isolid++
       continue;
       } 
     */
    G4cout << "solid name: " << solid->GetName() << G4endl;
    argosz_log[isolid]  = new G4LogicalVolume( solid
        , argosz_mat[isolid]
        , solid->GetName()//"logical"
        , 0, 0, 0
        );
    switch(isolid) {
    case 3:
    case 15:
      argosz_phys[isolid] = new G4PVPlacement( 0
          , G4ThreeVector(-0.40, 0, 0)
          , argosz_log[isolid]
          , solid->GetName()
          , world_log
          , false, 0
          );
      break;
    case 13: // Direct_beam_stop_0.75
      argosz_phys[isolid] = new G4PVPlacement( 0
          , G4ThreeVector(0.8660254038*dbshift, 0, 0.5*dbshift)
          , argosz_log[isolid]
          , solid->GetName()
          , world_log
          , false, 0
          );
      break;
    case 30:
      argosz_phys[isolid] = new G4PVPlacement( 0
          , G4ThreeVector(0.0, 0, 0)
          , argosz_log[isolid]
          , solid->GetName()
          , world_log 
          , false, 0 
          );
      break;
    case 31:
      argosz_phys[isolid] = new G4PVPlacement( 0
          , G4ThreeVector(0.0, 0, 0)
          , argosz_log[isolid]
          , solid->GetName()
          , world_log
          , false, 0
          );
      break;
    default:
      argosz_phys[isolid] = new G4PVPlacement( 0
          , G4ThreeVector(0, 0, 0)
          , argosz_log[isolid]
          , solid->GetName()
          , world_log
          , false, 0
          );
      break;
    }

    isolid++;
    //if(isolid == 16) break;
  }


  // The Bubble IN SALINE SOLUTION

  G4double Pz[11] = {-fBubble_r * 1.0, -fBubble_r * 0.8, -fBubble_r * 0.6, -fBubble_r * 0.4, -fBubble_r * 0.2, 0,
    fBubble_r * 0.2,	 fBubble_r * 0.4,  fBubble_r * 0.6,	 fBubble_r * 0.8,  fBubble_r * 1.0};
  G4double Prin[11] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
  G4double Prout[11] = {fBubble_r * 0.01, fBubble_r * 0.25, fBubble_r * 0.55, fBubble_r * 0.8,
    fBubble_r * 0.9,	fBubble_r * 0.99, fBubble_r * 0.92, fBubble_r * 0.8,
    fBubble_r * 0.5,	fBubble_r * 0.3,  fBubble_r * 0.05};
  std::cout << fBubble_r / CLHEP::um << std::endl;
  for (int i = 0; i < 11; i++)
  {
    std::cout << Pz[i] / CLHEP::um << " " << Prin[i] / CLHEP::um << " " << Prout[i] / CLHEP::um << std::endl;
  }
  auto bubbleWP = new G4Polycone("Bubble", 0., 360. * deg, 11, Pz, Prin, Prout);
  auto bubbleWP_log = new G4LogicalVolume(bubbleWP, water, "Bubble");

  auto bubbleW = new G4Orb("Bubble", fBubble_r);
  auto bubbleW_log = new G4LogicalVolume(bubbleW, water, "Bubble");

  G4double shift = 0.0 * mm;//0.035 * mm; // to make resulting CCD image symmetrical
  if (isPolycone)
  {
    // old position 14.09 * mm, 96.2425 * mm, -137.51 * mm 
    // 14.4999505 * mm, 96.250088 * mm, -137.4700015 * mm
    fBubble_def_pos = G4ThreeVector(14.0999505 * mm, 96.250088 * mm, -137.4700015 * mm + shift);
    bubble_phys = new G4PVPlacement(nullptr, fBubble_def_pos + fBubble_additional_offset, bubbleWP_log,
        "Bubble_dis_bnd_proc", argosz_log[15], false, 0);
  }
  else
  {
    //old model
    //fBubble_def_pos = G4ThreeVector(14.4999505 * mm, 96.250088 * mm, -137.4700015 * mm + shift);
    fBubble_def_pos = G4ThreeVector(-33.96995 * mm, 11.0412 * mm, -1.99995 * mm + shift);
    bubble_phys = new G4PVPlacement(nullptr, fBubble_def_pos + fBubble_additional_offset, bubbleW_log,
        "Bubble_dis_bnd_proc", argosz_log[15], false, 0);
  }

  if (IsVerbose())
  {
    std::cout << "Cell default position: " << fBubble_def_pos << std::endl;
    std::cout << "Cell offset: " << fBubble_additional_offset << std::endl;
    std::cout << "Cell current position: " << fBubble_additional_offset + fBubble_def_pos << std::endl;
  }

  /*
  //   -----  DETECTOR ARRAY  -----
  //

  // SENSITIVE DETECTOR MATERIAL PROPERTIES
  G4Material* detectorMaterial = nist->FindOrBuildMaterial("G4_Si");
  G4MaterialPropertiesTable* mptDetector = new G4MaterialPropertiesTable();

  // Define absorption length for optical photons
  const G4int detNumEntries = 2;
  G4double detPhotonEnergy[detNumEntries] = {1.5*eV, 2.0*eV}; // Photon energies
  G4double detAbsorptionLength[detNumEntries] = {.00001*mm, .00001*mm}; // Absorption length at those energies

  // Assign the material properties table to the detector material
  detectorMaterial->SetMaterialPropertiesTable(mptDetector);
  //mptDetector->AddProperty("ABSLENGTH", detPhotonEnergy, absorptionLength, detNumEntries);

  // SENSITIVE DETECTOR CONSTRUCTION
  G4double fullDetectorWidth = 1*mm;   // Width of the full detector
  G4double fullDetectorHeight = 1*mm;  // Height of the full detector
  G4double fullDetectorThickness = 0.1*mm; // Thickness of the detector
                                           // Number of pixels in X and Y directions
                                           G4int nPixelsZ = 10;
                                           G4int nPixelsY = 10;
  // Therefore, the pixel dimensions:
  G4double pixelWidth = fullDetectorWidth / nPixelsZ;
  G4double pixelHeight = fullDetectorHeight / nPixelsY;
  G4double pixelThickness = fullDetectorThickness;  // Same thickness for each pixe

  // Single pixel definition
  G4Box *solidDetector = new G4Box("solidDetector", fullDetectorWidth/(2*nPixelsZ), fullDetectorHeight/(2*nPixelsY), fullDetectorThickness/2);
  logicDetector = new G4LogicalVolume(solidDetector, detectorMaterial, "logicDetector");
  // Creating individual physical pixel instances
  for(G4int iz=0; iz<nPixelsZ; iz++)
  {
  for(G4int iy=0; iy<nPixelsY; iy++)
  {
  // X, Y, and Z position of the pixel
  G4double posX = 5*CLHEP::mm; // so far the whole plane at the same Z pos.
  G4double posY = -fullDetectorHeight/2 + pixelWidth*(iy +0.5);
  G4double posZ = -fullDetectorWidth/2 + pixelWidth*(iz +0.5);
  // Pixel placement
  G4VPhysicalVolume * physDetector = new G4PVPlacement(0, G4ThreeVector(posX, posY, posZ),
  logicDetector, "physDetector", world_log, false, iz+iy*nPixelsZ, true); // or logicWorld???
  }
  }
   */
  // ------------- Surfaces --------------

  // Water Tank
  auto opWaterSurface = new G4OpticalSurface("WaterSurface");
  opWaterSurface = new G4OpticalSurface("WaterSurface", glisur, polished, x_ray);

  //auto waterSurface = new G4LogicalBorderSurface("WaterSurface", world_phys, bubble_phys, opWaterSurface);
  auto waterSurface = new G4LogicalBorderSurface("WaterSurface", argosz_phys[17], bubble_phys, opWaterSurface); // in SALINE SOLUTION

  auto opticalSurface =
    dynamic_cast<G4OpticalSurface *>(waterSurface->GetSurface(argosz_phys[17], bubble_phys)->GetSurfaceProperty()); // in SALINE SOLUTION
  if (opticalSurface)
    opticalSurface->DumpInfo();

  // Mirror

  // Create optical surface
  G4OpticalSurface* opticalSurfaceMirror = new G4OpticalSurface("MirrorSurface");
  opticalSurfaceMirror = new G4OpticalSurface("MirrorSurface", unified, polishedfrontpainted, dielectric_dielectric);

  // Define reflection and transmission properties
  std::vector<G4double> reflectivity = { 0.91320,  0.91320,  0.91320,  0.91320};
  //std::vector<G4double> transmittance = {0.00, 0.00, 0.00, 0.00};
  std::vector<G4double> refractiveIndexMirror = {1.0972, 1.0972, 1.0972, 1.0972};
  G4MaterialPropertiesTable* SMPT = new G4MaterialPropertiesTable();

  SMPT->AddProperty("RINDEX", photonEnergyMirror, refractiveIndexMirror, nEntries);
  SMPT->AddProperty("REFLECTIVITY", photonEnergyMirror, reflectivity, nEntries);
  //SMPT->AddProperty("TRANSMITTANCE", photonEnergyMirror, transmittance, nEntries);

  opticalSurfaceMirror->SetMaterialPropertiesTable(SMPT);

  G4LogicalBorderSurface* mirrorSurface_1 = new G4LogicalBorderSurface("MirrorBorderSurface_1", world_phys, argosz_phys[0], opticalSurfaceMirror);
  G4LogicalBorderSurface* mirrorSurface_2 = new G4LogicalBorderSurface("MirrorBorderSurface_2", world_phys, argosz_phys[1], opticalSurfaceMirror);
  G4LogicalBorderSurface* mirrorSurface_11 = new G4LogicalBorderSurface("MirrorBorderSurface_11", world_phys, argosz_phys[11], opticalSurfaceMirror);
  G4LogicalBorderSurface* mirrorSurface_14 = new G4LogicalBorderSurface("MirrorBorderSurface_14", world_phys, argosz_phys[14], opticalSurfaceMirror);

  // Lens
  // Create optical surface
  G4OpticalSurface* opticalSurfaceLens = new G4OpticalSurface("LensSurface");
  opticalSurfaceLens = new G4OpticalSurface("LensSurface", unified, polished, dielectric_dielectric);

  //G4LogicalBorderSurface* flowcellSurface1 = new G4LogicalBorderSurface("LensBorderSurface1", world_phys, argosz_phys[3], opticalSurfaceLens);
  // ACL12708U
  G4LogicalBorderSurface* lensSurface1_in = new G4LogicalBorderSurface("LensBorderSurface2", world_phys, argosz_phys[4], opticalSurfaceLens);
  G4LogicalBorderSurface* lensSurface1_out = new G4LogicalBorderSurface("LensBorderSurface2_1", argosz_phys[4], world_phys, opticalSurfaceLens);
  // LB1258-A
  G4LogicalBorderSurface* lensSurface2_in = new G4LogicalBorderSurface("LensBorderSurface2", world_phys, argosz_phys[10], opticalSurfaceLens);
  G4LogicalBorderSurface* lensSurface2_out = new G4LogicalBorderSurface("LensBorderSurface2_1", argosz_phys[10], world_phys, opticalSurfaceLens);

  // Beam splitter
  // BST04_BeamSplitter
  G4OpticalSurface* splitterSurface_front = new G4OpticalSurface("SplitterSurface", unified, polished, dielectric_dielectric);
  G4OpticalSurface* splitterSurface_back = new G4OpticalSurface("SplitterSurface", unified, polished, dielectric_dielectric);

  G4MaterialPropertiesTable* surfaceMPT_front = new G4MaterialPropertiesTable();
  G4MaterialPropertiesTable* surfaceMPT_back = new G4MaterialPropertiesTable();
  std::vector<G4double> reflectivitySplitter_front = {0.5, 0.5, 0.5, 0.5};  // 50% reflection
  std::vector<G4double> reflectivitySplitter_back = {0., 0., 0., 0.};  // 0% reflection
  std::vector<G4double> transmittanceSplitter_front = {0.5, 0.5, 0.5, 0.5}; // 50% transmission
  std::vector<G4double> transmittanceSplitter_back = {1.0, 1.0, 1.0, 1.}; // 100% transmission
                                                                          //surfaceMPT_front->AddProperty("REFLECTIVITY", photonEnergyMirror, reflectivitySplitter_front);
                                                                          //surfaceMPT_back->AddProperty("REFLECTIVITY", photonEnergyMirror, reflectivitySplitter_back);
  surfaceMPT_front->AddProperty("TRANSMITTANCE", photonEnergyMirror, transmittanceSplitter_front);
  surfaceMPT_back->AddProperty("TRANSMITTANCE", photonEnergyMirror, transmittanceSplitter_back);

  splitterSurface_front->SetMaterialPropertiesTable(surfaceMPT_front);
  splitterSurface_back->SetMaterialPropertiesTable(surfaceMPT_back);
  //this is not used now in the simulation
  //G4LogicalBorderSurface* splitterSurface1 = new G4LogicalBorderSurface("splitterBorderSurface1", world_phys, argosz_phys[0], opticalSurfaceLens);
  G4LogicalBorderSurface* splitterBorderSurface_front = new G4LogicalBorderSurface("splitterBorderSurface_front", world_phys, argosz_phys[6], splitterSurface_front);
  G4LogicalBorderSurface* splitterBorderSurface_back = new G4LogicalBorderSurface("splitterBorderSurface_back", argosz_phys[6], world_phys, splitterSurface_back);

  //flowcell - seems to be acceptable?
  G4OpticalSurface* flowcellSurface = new G4OpticalSurface("flowcellSurface", unified, polished, dielectric_dielectric);
  G4LogicalBorderSurface* flowcellBorderSurface_in_out = new G4LogicalBorderSurface("flowcellBorderSurface_in_out", argosz_phys[3], world_phys, splitterSurface_back);
  G4LogicalBorderSurface* flowcellBorderSurface_out_in = new G4LogicalBorderSurface("flowcellBorderSurface_out_in", world_phys, argosz_phys[3], splitterSurface_back);
  
  // ND filter
  G4OpticalSurface* NDFilterSurface = new G4OpticalSurface("NDFilterSurface");
  NDFilterSurface->SetType(dielectric_dielectric);
  NDFilterSurface->SetFinish(polished); // or ground if it’s a diffusing ND filter
  NDFilterSurface->SetModel(unified);

  G4LogicalBorderSurface* ND_in = new G4LogicalBorderSurface("NDFilterBorder", world_phys, argosz_phys[22], NDFilterSurface);
  G4LogicalBorderSurface* ND_out = new G4LogicalBorderSurface("NDFilterBorder", argosz_phys[22], world_phys, NDFilterSurface);




  ////////////////////////////////////////////////////////////////////////////
  // test user-defined properties
  G4String ed;
  if (myMPT1->GetProperty("USERDEFINED") != nullptr)
  {
    ed = "USERDEFINED != nullptr";
    PrintError(ed);
  }
  myMPT1->AddProperty("USERDEFINED", energy_water, mie_water, true, true);
  if (myMPT1->GetProperty("USERDEFINED") == nullptr)
  {
    ed = "USERDEFINED == nullptr";
    PrintError(ed);
  }
  [[maybe_unused]] G4int index_userdefined = -1;
  if (myMPT1->GetProperty("USERDEFINED") != nullptr)
  {
    index_userdefined = myMPT1->GetPropertyIndex("USERDEFINED");
  }
  if (!(index_userdefined >= 0 && index_userdefined < (G4int)myMPT1->GetMaterialPropertyNames().size()))
  {
    ed = "USERDEFINED index out of range";
    PrintError(ed);
  }
  myMPT1->RemoveProperty("USERDEFINED");
  if (myMPT1->GetProperty("USERDEFINED") != nullptr)
  {
    ed = "USERDEFINED != nullptr at end";
    PrintError(ed);
  }

  if (myMPT1->ConstPropertyExists("USERDEFINEDCONST") == true)
  {
    ed = "ConstProperty USERDEFINEDCONST already exists.";
    PrintError(ed);
  }
  myMPT1->AddConstProperty("USERDEFINEDCONST", 3.14, true);
  if (myMPT1->ConstPropertyExists("USERDEFINEDCONST") == false)
  {
    ed = "ConstProperty USERDEFINEDCONST doesn't exist.";
    PrintError(ed);
  }
  [[maybe_unused]] G4int index_pi = -1;
  if (myMPT1->ConstPropertyExists("USERDEFINEDCONST") == true)
  {
    index_pi = myMPT1->GetConstPropertyIndex("USERDEFINEDCONST");
  }
  if (!(index_pi >= 0 && index_pi < (G4int)myMPT1->GetMaterialConstPropertyNames().size()))
  {
    ed = "ConstProperty USERDEFINEDCONST index out of range.";
    PrintError(ed);
  }
  myMPT1->RemoveConstProperty("USERDEFINEDCONST");
  if (myMPT1->ConstPropertyExists("USERDEFINEDCONST") == true)
  {
    ed = "ConstProperty USERDEFINEDCONST still exists.";
    PrintError(ed);
  }
  // done testing user-defined properties
  ////////////////////////////////////////////////////////////////////////////

  return world_phys;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void NormaDetectorConstruction::SetVerbose(G4bool val)
{
  fVerbose = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool NormaDetectorConstruction::IsVerbose() const
{
  return fVerbose;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void NormaDetectorConstruction::PrintError(G4String ed)
{
  G4Exception("NormaDetectorConstruction:MaterialProperty test", "op001", FatalException, ed);
}

void NormaDetectorConstruction::SetParameters(Parameters p)
{
  fBubble_r = p.bubbleRadius;
  std::cout << "Bubble radius set to " << fBubble_r << std::endl;
  mieFg = p.g;
  std::cout << "mie forward " << mieFg << std::endl;
  isPolycone = p.p;
  std::cout << "is Polycone " << mieFg << std::endl;
}

/*
   void NormaDetectorConstruction::ConstructSDandField()
   {
   G4cout << "ConstructSDandField called, setting sensitive detector..." << G4endl;
   PhotoSensitiveDetector *sensDet = new PhotoSensitiveDetector("SensitiveDetector");
   logicDetector->SetSensitiveDetector(sensDet);
   }
 */
