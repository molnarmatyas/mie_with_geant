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
/// \file Norma/src/NormaPrimaryGeneratorAction.cc
/// \brief Implementation of the NormaPrimaryGeneratorAction class
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "NormaPrimaryGeneratorAction.hh"
#include "G4Event.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "NormaPrimaryGeneratorMessenger.hh"
#include "Randomize.hh"
#include <numeric>  // for std::accumulate
#include <fstream>
#include <sstream>


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
NormaPrimaryGeneratorAction::NormaPrimaryGeneratorAction() : G4VUserPrimaryGeneratorAction(), fParticleGun(nullptr)
{
	G4int n_particle = 1;
	fParticleGun = new G4ParticleGun(n_particle);
	// create a messenger for this class
	fGunMessenger = new NormaPrimaryGeneratorMessenger(this);
	// default kinematic
	//
	G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
	G4ParticleDefinition *particle = particleTable->FindParticle("opticalphoton");

	fParticleGun->SetParticleDefinition(particle);
	fParticleGun->SetParticleTime(0.0 * ns);
	// G4double y = (G4UniformRand()-0.5) * LEDSize;
	G4double y = (0.0);// - (1.5 * CLHEP::um));
	// std::cout<<"startpos photon "<<y/CLHEP::um<<std::endl;
	//fParticleGun->SetParticlePosition(G4ThreeVector(-17 * CLHEP::um, y, 0*mm));
  //fParticleGun->SetParticlePosition(G4ThreeVector(-29.35 * CLHEP::mm, 96.076 * CLHEP::mm, -157.841 * CLHEP::mm)); //laser
  G4double shift = .0;//0.035 * mm; // to make resulting CCD image symmetrical
	fParticleGun->SetParticlePosition(G4ThreeVector(4.4999505 * mm, 96.250088 * mm, -137.4700015 * mm + shift)); //14.49 //center of cell
		
  /*
	G4double dy = (G4UniformRand() - 0.5) * 0.1;
	G4double dz = (G4UniformRand() - 0.5) * 0.1;
	// fParticleGun->SetParticleMomentumDirection(G4ThreeVector(1., dy, dz));
  */
	fParticleGun->SetParticleMomentumDirection(G4ThreeVector(1., 0.0, 0.)); // -1 backwards
	fParticleGun->SetParticleEnergy(1.884 * eV);
  //G4cout << "Gun from " << fParticleGun->GetParticlePosition().x() << ", " << fParticleGun->GetParticlePosition().y() << ", " << fParticleGun->GetParticlePosition().z() << G4endl; 

  //fParticleGun->SetParticlePosition(G4ThreeVector(0,0,0));
  //fParticleGun->GeneratePrimaryVertex(anEvent);

  // Initialize intensity profile
  std::string intensityFile = "N958260004_10.csv"; //"RayCi8_180mm.csv"; // You can pass this via a macro later
  G4double leftrightZ = 0.0*mm;
  //old model
  //profileCenterWorld = G4ThreeVector(14.4999505 * mm + -0.4*mm, 96.250088 * mm - 0.23*mm, -137.4700015 * mm - 0.9*mm + leftrightZ);
  //new model
  profileCenterWorld = G4ThreeVector(-33.96995 * mm -10*mm, 11.0412  * mm, -1.99995 * mm + 0.57*mm + leftrightZ); // 0.57mm shift needed to put profile in the middle; let it start 10 mm from cell, before the flowcell

  InitializeIntensityProfile(intensityFile, pixelSize);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
NormaPrimaryGeneratorAction::~NormaPrimaryGeneratorAction()
{
	delete fParticleGun;
	delete fGunMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void NormaPrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent)
{
	auto [i, j] = SamplePixel();
	int idx = i * numCols + j;
	double localZ = (j - numCols / 2.0) * pixelSize;
	double localY = (i - numRows / 2.0) * pixelSize; //z-y plane!
  
	G4ThreeVector emissionPoint = profileCenterWorld + G4ThreeVector(0, localY, localZ);
	G4ThreeVector emissionDirection = GetDirectionForPixel(idx);
	
	fParticleGun->SetParticlePosition(emissionPoint);
	fParticleGun->SetParticleMomentumDirection(emissionDirection);
    G4cout << "Emitting photon from pixel (" << i << ", " << j << ") at " 
           << emissionPoint.x() << ", " << emissionPoint.y() << ", " << emissionPoint.z() 
           << " with direction " << emissionDirection.x() << ", " << emissionDirection.y() << ", " << emissionDirection.z() 
           << G4endl;
 
    fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void NormaPrimaryGeneratorAction::SetOptPhotonPolar()
{
	G4double angle = G4UniformRand() * 360.0 * deg;
	SetOptPhotonPolar(angle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void NormaPrimaryGeneratorAction::SetOptPhotonPolar(G4double angle)
{
	if (fParticleGun->GetParticleDefinition()->GetParticleName() != "opticalphoton")
	{
		G4ExceptionDescription ed;
		ed << "Warning: the particleGun is not an opticalphoton";
		G4Exception("NormaPrimaryGeneratorAction::SetOptPhotonPolar()", "Norma_010", JustWarning, ed);
		return;
	}

	G4ThreeVector normal(1., 0., 0.);
	G4ThreeVector kphoton = fParticleGun->GetParticleMomentumDirection();
	G4ThreeVector product = normal.cross(kphoton);
	G4double modul2 = product * product;

	G4ThreeVector e_perpend(0., 0., 1.);
	if (modul2 > 0.)
		e_perpend = (1. / std::sqrt(modul2)) * product;
	G4ThreeVector e_paralle = e_perpend.cross(kphoton);

	G4ThreeVector polar = std::cos(angle) * e_paralle + std::sin(angle) * e_perpend;
	fParticleGun->SetParticlePolarization(polar);
}

void NormaPrimaryGeneratorAction::SetParameters(Parameters p)
{
	LEDSize = p.bubbleRadius * 2;
	std::cout << "LEDSzize " << LEDSize << std::endl;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// For csv intensity profile
void NormaPrimaryGeneratorAction::LoadIntensityCSV(const std::string& filename) {
    std::ifstream file(filename);
    std::string line;

    // Skip the first (header) line
    std::getline(file, line);

    std::vector<double> tempData;
    numRows = 0;

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string cell;
        int colCount = 0;

        while (std::getline(ss, cell, ',')) {
            double val = std::stod(cell);
            tempData.push_back(val);
            ++colCount;
        }

        if (colCount > 0) {
            ++numRows;
            numCols = colCount; // assume consistent width
        }
    }

    intensityMap = std::move(tempData);
}


void NormaPrimaryGeneratorAction::LoadDirectionsCSV(const std::string& filename) {
    std::ifstream file(filename);
    
    if (!file.is_open()) {
        G4cout << "Warning: Could not open photon directions file: " << filename << G4endl;
        G4cout << "Falling back to default direction (1., 0., 0.)" << G4endl;
        directionsLoaded = false;
        return;
    }
    
    std::string line;

    // Skip the first (header) line
    std::getline(file, line);

    std::vector<double> tempData;

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string cell;

        while (std::getline(ss, cell, ',')) {
            double val = std::stod(cell);
            tempData.push_back(val);
        }
    }

    file.close();
    directionMap = std::move(tempData);
    directionsLoaded = true;
}


void NormaPrimaryGeneratorAction::NormalizeAndBuildCDF() {
    // First, zero out values below threshold (e.g., 0.01)
    double threshold = 0.05;
    for (auto& val : intensityMap) {
        if (val < threshold) val = 0.0;
    }
    // Normalize
    // Find the maximum value for normalization
    double maxValue = *std::max_element(intensityMap.begin(), intensityMap.end());
    
    // Renormalize only if max value is greater than zero
    if (maxValue > 0.0) {
        for (auto& value : intensityMap) {
            value /= maxValue;
        }
    }
    /*
    double sum = std::accumulate(intensityMap.begin(), intensityMap.end(), 0.0);
    if (sum == 0.0) {
        G4Exception("NormaPrimaryGeneratorAction::NormalizeAndBuildCDF",
                    "ZeroIntensityMap", FatalException,
                    "All intensity values are zero after thresholding.");
    }

    for (auto& val : intensityMap) val /= sum;
    */
}


G4ThreeVector NormaPrimaryGeneratorAction::ComputeProfileCenter() {
	// Find center of mass of the intensity profile
    double xSum = 0, ySum = 0, total = 0;
    for (int i = 0; i < numRows; ++i) {
        for (int j = 0; j < numCols; ++j) {
            double w = intensityMap[i * numCols + j];
            xSum += j * w;
            ySum += i * w;
            total += w;
        }
    }

    if (total == 0.0) {
        G4Exception("NormaPrimaryGeneratorAction::ComputeProfileCenter",
                    "ZeroIntensityMap", FatalException,
                    "Cannot compute center of mass: all intensity values are zero.");
    }

    double cx = xSum / total;
    double cy = ySum / total;

    // Convert to world coordinates
    double worldX = 0.0; // X is fixed
    double worldZ = (cx - numCols / 2.0) * pixelSize; // in z-y plane!
    double worldY = (cy - numRows / 2.0) * pixelSize;
    /*
    G4cout << "Coordinates of profile center in profile's coordinates: " << G4endl;
    G4cout << "Y: " << worldY << "\tZ: " << worldZ << G4endl;
    */
    
    worldZ = cx * pixelSize;
    worldY = cy * pixelSize;

    return G4ThreeVector(worldX, worldY, worldZ) + profileCenterWorld;
}

std::pair<int, int> NormaPrimaryGeneratorAction::SamplePixel() {
	// Sampling from the 2D intensity map
    //this directly returns the idx based on the intensity map
    int idx = fDist_disc(fGen);
    
    // logging
    /*
    ++map_test[idx];
    for (const auto& [num, count] : map_test)
        std::cout << num << " generated " << std::setw(4) << count << " times\n";
    
    G4cout << "idx: " << idx << G4endl;
    */
    
    // to coordinates
    int i = idx / numCols;
    int j = idx % numCols;
    return {i, j};
}

G4ThreeVector NormaPrimaryGeneratorAction::GetDirectionForPixel(int idx) {
    // If directions were not loaded, return default direction
    if (!directionsLoaded) {
        return G4ThreeVector(1., 0., 0.);
    }
    
    // Each pixel has 3 direction components (x, y, z)
    if (idx < 0 || (idx * 3 + 2) >= static_cast<int>(directionMap.size())) {
        G4Exception("NormaPrimaryGeneratorAction::GetDirectionForPixel",
                    "InvalidPixelIndex", FatalException,
                    "Pixel index out of range for direction map.");
    }

    double dx = directionMap[idx * 3];
    double dy = directionMap[idx * 3 + 1];
    double dz = directionMap[idx * 3 + 2];

    return G4ThreeVector(dx, dy, dz);
}

void NormaPrimaryGeneratorAction::InitializeIntensityProfile(const std::string& filename, double pixelSizeMM) {
    //pixelSize = pixelSizeMM * CLHEP::mm;
    //G4cout << "pixelSize normal: " << pixelSize << G4endl;


    LoadIntensityCSV(filename);
    NormalizeAndBuildCDF();
    G4ThreeVector beamCenter = ComputeProfileCenter();

	fDist_disc = std::discrete_distribution<>(intensityMap.begin(), intensityMap.end());
	//fDist_disc = std::piecewise_linear_distribution<>(intensityMap.begin(), intensityMap.end(), intensityMap.begin());
	fGen.seed(time(0)); // if you want different results from different runs

    // Load directions from the corresponding CSV file
    LoadDirectionsCSV("photon_directions_full_matrix_sequential.csv");

   // G4cout << "→ Laser profile loaded from: " << filename << G4endl;
   // G4cout << "→ Profile dimensions: " << numCols << " × " << numRows << G4endl;
   // G4cout << "→ Computed center of mass: " << beamCenter / mm << " mm" << G4endl;
}
