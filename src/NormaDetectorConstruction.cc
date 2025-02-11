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
	myMPT1->DumpTable();

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
	myMPT0->DumpTable();

	air->SetMaterialPropertiesTable(myMPT0);


  // Mirror material
  auto mirrorMaterial = nist->FindOrBuildMaterial("G4_Al");
  
  G4MaterialPropertiesTable* myMPT3 = new G4MaterialPropertiesTable();
  std::vector<G4double> photonEnergyMirror = {1.8*eV, 1.85*eV, 1.9*eV, 1.95*eV };
  mirrorMaterial->SetMaterialPropertiesTable(myMPT3);

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




	// ------------- Volumes --------------
	//
	// The world
	auto world_box = new G4Box("World", fWorld_x, fWorld_y, fWorld_z);
	auto world_log = new G4LogicalVolume(world_box, air, "World");
	G4VPhysicalVolume *world_phys =
		new G4PVPlacement(nullptr, G4ThreeVector(), world_log, "world", nullptr, false, 0, checkOverlaps);

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
  auto mesh = CADMesh::TessellatedMesh::FromOBJ("./Argosz_sensorupdate_with_shield.obj");
  G4cout << " MESH NAME: " << mesh->GetFileName() << G4endl;;
  mesh->SetScale(1.0);
  std::vector<G4VSolid*> solids = mesh->GetSolids();

  std::vector<G4LogicalVolume*> argosz_log(solids.size());
  std::vector<G4VPhysicalVolume*> argosz_phys(solids.size());

	std::vector<G4Material*> argosz_mat(solids.size());

  /* 1 -> 16
  0 - solid name: Beamstop_flowcell
  1 - solid name: beam_splitter
  2 - solid name: mirror
  3 - solid name: COHERENT_MINI-701L-660S
  4 - solid name: Hellma_flowcell_131-814-40
  5 - solid name: ACL12708U
  6 - solid name: GS3-U3-23S6M-C_sensor
  7 - solid name: BST04_BeamSplitter
  8 - solid name: Direct_beam_stop_2
  9 - solid name: vbpw34s_1
  10 - solid name: vbpw34s_2
  11 - solid name: LB1258-A
  12 - solid name: LA_Mirror
  13 - solid name: LA_HA_mirror
  14 - solid name: Direct_beam_stop
  15 - solid name: HA_mirror
  16 - solid name: GS3-U3-23S6M-C_sensor_housing
  */
 argosz_mat[1] = mirrorMaterial;
 argosz_mat[2] = mirrorMaterial;
 argosz_mat[3] = mirrorMaterial;
 argosz_mat[4] = lensMaterial; // let us try this, larger abs. length
 argosz_mat[5] = lensMaterial;
 argosz_mat[6] = sil;
 argosz_mat[7] = lensMaterial;
 argosz_mat[8] = shieldMaterial;
 argosz_mat[9] = lensMaterial;
 argosz_mat[10] = lensMaterial;
 argosz_mat[11] = lensMaterial;
 argosz_mat[12] = mirrorMaterial;
 argosz_mat[13] = mirrorMaterial;
 argosz_mat[14] = shieldMaterial;
 argosz_mat[15] = mirrorMaterial;
 argosz_mat[16] = lensMaterial;
 argosz_mat[0] = shieldMaterial;

  int isolid = 0;
  for (auto solid : mesh->GetSolids())
  {
    if(isolid == 14 || isolid == 16) continue;
    G4cout << "solid name: " << solid->GetName() << G4endl;
    argosz_log[isolid]  = new G4LogicalVolume( solid
                                        , argosz_mat[isolid]
                                        , solid->GetName()//"logical"
                                        , 0, 0, 0
    );
      argosz_phys[isolid] = new G4PVPlacement( 0
                        , G4ThreeVector(0, 0, 0)
                        , argosz_log[isolid]
                        , solid->GetName()
                        , world_log
                        , false, 0
      );
    
    isolid++;
    //if(isolid == 16) break;
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

	auto waterSurface = new G4LogicalBorderSurface("WaterSurface", world_phys, bubble_phys, opWaterSurface);

	auto opticalSurface =
		dynamic_cast<G4OpticalSurface *>(waterSurface->GetSurface(world_phys, bubble_phys)->GetSurfaceProperty());
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

	G4LogicalBorderSurface* mirrorSurface_1 = new G4LogicalBorderSurface("MirrorBorderSurface_1", world_phys, argosz_phys[1], opticalSurfaceMirror);
	G4LogicalBorderSurface* mirrorSurface_2 = new G4LogicalBorderSurface("MirrorBorderSurface_2", world_phys, argosz_phys[2], opticalSurfaceMirror);
	G4LogicalBorderSurface* mirrorSurface_11 = new G4LogicalBorderSurface("MirrorBorderSurface_11", world_phys, argosz_phys[12], opticalSurfaceMirror);
	G4LogicalBorderSurface* mirrorSurface_12 = new G4LogicalBorderSurface("MirrorBorderSurface_12", world_phys, argosz_phys[13], opticalSurfaceMirror);
	G4LogicalBorderSurface* mirrorSurface_14 = new G4LogicalBorderSurface("MirrorBorderSurface_14", world_phys, argosz_phys[15], opticalSurfaceMirror);

  // Lens
	// Create optical surface
  G4OpticalSurface* opticalSurfaceLens = new G4OpticalSurface("LensSurface");
  opticalSurfaceLens = new G4OpticalSurface("LensSurface", unified, polished, dielectric_dielectric);

  // Define reflection and transmission properties
  std::vector<G4double> transmittanceLens = {0.99,  0.99,  0.99,  0.99};//{0.91320,  0.91320,  0.91320,  0.91320};
  std::vector<G4double> reflectivityLens = {0.01, 0.01, 0.01, 0.01};//{0.00, 0.00, 0.00, 0.00};
	std::vector<G4double> refractiveIndexLens = {1.52, 1.52, 1.52, 1.52};//{1.0972, 1.0972, 1.0972, 1.0972};
  G4MaterialPropertiesTable* SMPTlens = new G4MaterialPropertiesTable();
  

	G4LogicalBorderSurface* lensSurface1 = new G4LogicalBorderSurface("LensBorderSurface1", world_phys, argosz_phys[4], opticalSurfaceLens);
  
	G4LogicalBorderSurface* lensSurface2 = new G4LogicalBorderSurface("LensBorderSurface2", world_phys, argosz_phys[5], opticalSurfaceLens);
	G4LogicalBorderSurface* lensSurface2_1 = new G4LogicalBorderSurface("LensBorderSurface2_1", argosz_phys[5], world_phys, opticalSurfaceLens);
//  opticalSurfaceLens->SetMaterialPropertiesTable(SMPTlens);//myMPT5);

  // Beam splitter
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
	G4LogicalBorderSurface* splitterBorderSurface_front = new G4LogicalBorderSurface("splitterBorderSurface_front", world_phys, argosz_phys[7], splitterSurface_front);
	G4LogicalBorderSurface* splitterBorderSurface_back = new G4LogicalBorderSurface("splitterBorderSurface_back", argosz_phys[7], world_phys, splitterSurface_back);

  //flowcell - not working atm FIXME
  G4OpticalSurface* flowcellSurface = new G4OpticalSurface("flowcellSurface", unified, polished, dielectric_dielectric);
	G4LogicalBorderSurface* flowcellBorderSurface_in_out = new G4LogicalBorderSurface("flowcellBorderSurface_in_out", argosz_phys[4], world_phys, splitterSurface_back);
	G4LogicalBorderSurface* flowcellBorderSurface_out_in = new G4LogicalBorderSurface("flowcellBorderSurface_out_in", world_phys, argosz_phys[4], splitterSurface_back);



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
