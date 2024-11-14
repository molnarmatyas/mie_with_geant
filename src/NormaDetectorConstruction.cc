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
	std::vector<G4double> photonEnergy = {
		2.034 * eV, 2.068 * eV, 2.103 * eV, 2.139 * eV, 2.177 * eV, 2.216 * eV, 2.256 * eV, 2.298 * eV,
		2.341 * eV, 2.386 * eV, 2.433 * eV, 2.481 * eV, 2.532 * eV, 2.585 * eV, 2.640 * eV, 2.697 * eV,
		2.757 * eV, 2.820 * eV, 2.885 * eV, 2.954 * eV, 3.026 * eV, 3.102 * eV, 3.181 * eV, 3.265 * eV,
		3.353 * eV, 3.446 * eV, 3.545 * eV, 3.649 * eV, 3.760 * eV, 3.877 * eV, 4.002 * eV, 4.136 * eV};

	// Water
	std::vector<G4double> refractiveIndex1 = {1.3435, 1.344,  1.3445, 1.345,  1.3455, 1.346,  1.3465, 1.347,
											  1.3475, 1.348,  1.3485, 1.3492, 1.35,	  1.3505, 1.351,  1.3518,
											  1.3522, 1.3530, 1.3535, 1.354,  1.3545, 1.355,  1.3555, 1.356,
											  1.3568, 1.3572, 1.358,  1.3585, 1.359,  1.3595, 1.36,	  1.3608};
	std::vector<G4double> absorption = {
		3.448 * m,	4.082 * m,	6.329 * m,	9.174 * m,	12.346 * m, 13.889 * m, 15.152 * m, 17.241 * m,
		18.868 * m, 20.000 * m, 26.316 * m, 35.714 * m, 45.455 * m, 47.619 * m, 52.632 * m, 52.632 * m,
		55.556 * m, 52.632 * m, 52.632 * m, 47.619 * m, 45.455 * m, 41.667 * m, 37.037 * m, 33.333 * m,
		30.000 * m, 28.500 * m, 27.000 * m, 24.500 * m, 22.000 * m, 19.500 * m, 17.500 * m, 14.500 * m};

	// Material properties can be added as arrays. However, in this case it is
	// up to the user to make sure both arrays have the same number of elements.
	G4double scintilFastArray[]{1.0, 1.0};
	G4double energyArray[]{2.034 * eV, 4.136 * eV};
	G4int lenArray = 2;

	std::vector<G4double> scintilSlow = {0.01, 1.00, 2.00, 3.00, 4.00, 5.00, 6.00, 7.00, 8.00, 9.00, 8.00,
										 7.00, 6.00, 4.00, 3.00, 2.00, 1.00, 0.01, 1.00, 2.00, 3.00, 4.00,
										 5.00, 6.00, 7.00, 8.00, 9.00, 8.00, 7.00, 6.00, 5.00, 4.00};

	auto myMPT1 = new G4MaterialPropertiesTable();

	// Values can be added to the material property table individually.
	// With this method, spline interpolation cannot be set. Arguments
	// createNewKey and spline both take their default values of false.
	// Need to specify the number of entries (1) in the arrays, as an argument
	// to AddProperty.
	G4int numEntries = 1;
	myMPT1->AddProperty("RINDEX", &photonEnergy[0], &refractiveIndex1[0], numEntries);

	for (size_t i = 1; i < photonEnergy.size(); ++i)
	{
		myMPT1->AddEntry("RINDEX", photonEnergy[i], refractiveIndex1[i]);
	}

	// Check that group velocity is calculated from RINDEX
	if (myMPT1->GetProperty("RINDEX")->GetVectorLength() != myMPT1->GetProperty("GROUPVEL")->GetVectorLength())
	{
		G4ExceptionDescription ed;
		ed << "Error calculating group velocities. Incorrect number of entries "
			  "in group velocity material property vector.";
		G4Exception("Norma::NormaDetectorConstruction", "Norma001", FatalException, ed);
	}

	// Adding a property from two std::vectors. Argument createNewKey is false
	// and spline is true.
	myMPT1->AddProperty("ABSLENGTH", photonEnergy, absorption, false, true);

	// Adding a property using a C-style array.
	// Spline interpolation isn't used for scintillation.
	// Arguments spline and createNewKey both take default value false.
	myMPT1->AddProperty("SCINTILLATIONCOMPONENT1", energyArray, scintilFastArray, lenArray);

	myMPT1->AddProperty("SCINTILLATIONCOMPONENT2", photonEnergy, scintilSlow, false, true);
	myMPT1->AddConstProperty("SCINTILLATIONYIELD", 50. / MeV);
	myMPT1->AddConstProperty("RESOLUTIONSCALE", 1.0);
	myMPT1->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 1. * ns);
	myMPT1->AddConstProperty("SCINTILLATIONTIMECONSTANT2", 10. * ns);
	myMPT1->AddConstProperty("SCINTILLATIONYIELD1", 0.8);
	myMPT1->AddConstProperty("SCINTILLATIONYIELD2", 0.2);
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

	// std::vector<G4double> mie_water = {
	// 	1670.244 * nm, 1587.267 * nm, 1507.42 * nm, 1430.625 * nm, 1356.802 * nm,
	// 	1285.874 * nm, 1217.763 * nm, 1152.395 * nm, 1089.69 * nm, 1029.588 * nm,
	// 	972.0035 * nm, 916.8686 * nm, 864.1133 * nm, 813.6679 * nm, 765.4642 * nm,
	// 	71943.46 * nm, 67551.29 * nm, 63363.36 * nm, 59373.25 * nm, 55574.61 * nm,
	// 	51961.24 * nm, 48527.00 * nm, 45265.87 * nm, 42171.94 * nm, 39239.39 * nm,
	// 	36462.50 * nm, 33835.68 * nm, 31353.41 * nm, 29010.30 * nm, 26801.03 * nm,
	// 	24720.42 * nm, 22763.36 * nm, 20924.88 * nm, 19200.07 * nm, 17584.16 * nm,
	// 	16072.45 * nm, 14660.38 * nm, 13343.46 * nm, 12117.33 * nm, 10977.70 * nm,
	// 	9920.416 * nm, 8941.407 * nm, 8036.711 * nm, 7202.470 * nm, 6434.927 * nm,
	// 	5730.429 * nm, 5085.425 * nm, 4496.467 * nm, 3960.210 * nm, 3473.413 * nm,
	// 	3032.937 * nm, 2635.746 * nm, 2278.907 * nm, 1959.588 * nm, 1675.064 * nm,
	// 	1422.710 * nm, 1200.004 * nm, 1004.528 * nm, 833.9666 * nm, 686.1063 * nm};

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
	std::vector<G4double> refractiveIndex2 = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
											  1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
											  1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

	auto myMPT2 = new G4MaterialPropertiesTable();
	myMPT2->AddProperty("RINDEX", photonEnergy, refractiveIndex2);

	/*
	myMPT2->AddProperty("MIEHG", energy_water, mie_water, false, true);
	myMPT2->AddConstProperty("MIEHG_FORWARD", mie_water_const[0]);
	myMPT2->AddConstProperty("MIEHG_BACKWARD", mie_water_const[1]);
	myMPT2->AddConstProperty("MIEHG_FORWARD_RATIO", mie_water_const[2]);
	*/

	G4cout << "Air G4MaterialPropertiesTable:" << G4endl;
	myMPT2->DumpTable();

	air->SetMaterialPropertiesTable(myMPT2);

	// ------------- Volumes --------------
	//
	// The world
	auto world_box = new G4Box("World", fWorld_x, fWorld_y, fWorld_z);
	auto world_log = new G4LogicalVolume(world_box, air, "World");
	G4VPhysicalVolume *world_phys =
		new G4PVPlacement(nullptr, G4ThreeVector(), world_log, "world", nullptr, false, 0, checkOverlaps);

	// The experimental Hall
	/*
	auto expHall_box = new G4Box("expHall", fExpHall_x, fExpHall_y, fExpHall_z);
	auto expHall_log = new G4LogicalVolume(expHall_box, air, "expHall");
	G4VPhysicalVolume *expHall_phys = new G4PVPlacement(
	nullptr, G4ThreeVector(), expHall_log, "expHall", world_log, false, 0);
	*/

	// The Water Tank
	/*
	auto waterTank_box = new G4Box("Tank", fTank_x, fTank_y, fTank_z);
	auto waterTank_log = new G4LogicalVolume(waterTank_box, water, "Tank");
	G4VPhysicalVolume *waterTank_phys = new G4PVPlacement(
	nullptr, G4ThreeVector(), waterTank_log, "Tank", expHall_log,
	false, 0);
	*/

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
		bubble_phys = new G4PVPlacement(nullptr, G4ThreeVector(0.0 * cm, 0.0 * cm, 0.0 * cm), bubbleWP_log,
										"Bubble_dis_bnd_proc", world_log, false, 0);
	}
	else
	{
		bubble_phys = new G4PVPlacement(nullptr, G4ThreeVector(0.0 * cm, 0.0 * cm, 0.0 * cm), bubbleW_log,
										"Bubble_dis_bnd_proc", world_log, false, 0);
	}

	// screen
	auto screen_box = new G4Box("Screen", fScreen_x, fScreen_y, fScreen_z);
	auto screen_log = new G4LogicalVolume(screen_box, sil, "Screen");
	G4VPhysicalVolume *screen_phys = new G4PVPlacement(nullptr, G4ThreeVector(19.0 * mm, 0.0 * mm, 0.0 * mm),
													   screen_log, "Screen", world_log, false, 0);

	auto screenB_box = new G4Box("ScreenB", fScreen_x, fScreen_y, fScreen_z);
	auto screenB_log = new G4LogicalVolume(screenB_box, sil, "ScreenB");
	G4VPhysicalVolume *screenB_phys = new G4PVPlacement(nullptr, G4ThreeVector(-19.0 * mm, 0.0 * mm, 0.0 * mm),
														screenB_log, "ScreenB", world_log, false, 0);

	auto screenL_box = new G4Box("ScreenL", fScreen_y, fScreen_y, fScreen_x);
	auto screenL_log = new G4LogicalVolume(screenL_box, sil, "ScreenL");
	G4VPhysicalVolume *screenL_phys = new G4PVPlacement(nullptr, G4ThreeVector(0.0 * mm, 0.0 * mm, 19.0 * mm),
														screenL_log, "ScreenL", world_log, false, 0);

	auto screenR_box = new G4Box("ScreenR", fScreen_y, fScreen_y, fScreen_x);
	auto screenR_log = new G4LogicalVolume(screenR_box, sil, "ScreenR");
	G4VPhysicalVolume *screenR_phys = new G4PVPlacement(nullptr, G4ThreeVector(0.0 * mm, 0.0 * mm, -19.0 * mm),
														screenR_log, "ScreenR", world_log, false, 0);

	auto screenU_box = new G4Box("ScreenU", fScreen_y, fScreen_x, fScreen_y);
	auto screenU_log = new G4LogicalVolume(screenU_box, sil, "ScreenU");
	G4VPhysicalVolume *screenU_phys = new G4PVPlacement(nullptr, G4ThreeVector(0.0 * mm, 19.0 * mm, 0.0 * mm),
														screenU_log, "ScreenU", world_log, false, 0);

	auto screenD_box = new G4Box("ScreenD", fScreen_y, fScreen_x, fScreen_y);
	auto screenD_log = new G4LogicalVolume(screenD_box, sil, "ScreenD");
	G4VPhysicalVolume *screenD_phys = new G4PVPlacement(nullptr, G4ThreeVector(0.0 * mm, -19.0 * mm, 0.0 * mm),
														screenD_log, "ScreenD", world_log, false, 0);

	// ------------- Surfaces --------------

	// Water Tank
	auto opWaterSurface = new G4OpticalSurface("WaterSurface");
	opWaterSurface = new G4OpticalSurface("WaterSurface", glisur, polished, x_ray);
	// opWaterSurface->SetType(dielectric_LUTDAVIS);
	// opWaterSurface->SetFinish(Rough_LUT);
	// opWaterSurface->SetModel(DAVIS);

	auto waterSurface = new G4LogicalBorderSurface("WaterSurface", world_phys, bubble_phys, opWaterSurface);

	auto opticalSurface =
		dynamic_cast<G4OpticalSurface *>(waterSurface->GetSurface(world_phys, bubble_phys)->GetSurfaceProperty());
	if (opticalSurface)
		opticalSurface->DumpInfo();

	// Air Bubble
	/*
	auto opAirSurface = new G4OpticalSurface("AirSurface");
	opAirSurface->SetType(dielectric_dielectric);
	opAirSurface->SetFinish(polished);
	opAirSurface->SetModel(glisur);

	auto airSurface =
	new G4LogicalSkinSurface("AirSurface", bubbleAir_log, opAirSurface);

	opticalSurface = dynamic_cast<G4OpticalSurface *>(
	airSurface->GetSurface(bubbleAir_log)->GetSurfaceProperty());
	if (opticalSurface)
	opticalSurface->DumpInfo();
	*/
	// Generate & Add Material Properties Table attached to the optical surfaces
	//
	std::vector<G4double> ephoton = {2.034 * eV, 4.136 * eV};

	// OpticalAirSurface
	std::vector<G4double> reflectivity = {0.3, 0.5};
	std::vector<G4double> efficiency = {0.8, 1.0};

	/*
	auto myST2 = new G4MaterialPropertiesTable();

	myST2->AddProperty("REFLECTIVITY", ephoton, reflectivity);
	myST2->AddProperty("EFFICIENCY", ephoton, efficiency);
	if (fVerbose)
	{
	G4cout << "Air Surface G4MaterialPropertiesTable:" << G4endl;
	myST2->DumpTable();
}
opAirSurface->SetMaterialPropertiesTable(myST2);
*/

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
