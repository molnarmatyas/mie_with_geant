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
	//fParticleGun->SetParticlePosition(G4ThreeVector(4.4999505 * mm, 96.250088 * mm, -137.4700015 * mm + shift)); //14.49 //center of cell
		
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
  std::string intensityFile = "RayCi8_180mm.csv"; // You can pass this via a macro later
  G4ThreeVector profileCenterWorld = G4ThreeVector(14.4999505 * mm, 96.250088 * mm, -137.4700015 * mm);

  InitializeIntensityProfile(intensityFile, profileCenterWorld, pixelSize);
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
	/*
	G4double t = 2*CLHEP::pi  * G4UniformRand();
	G4double r = std::sqrt(0.1/CLHEP::pi) * (G4UniformRand());
	G4double x0 = 14.09 * CLHEP::mm;
	G4double y0 = r * std::cos(t) + 96.2425 * CLHEP::mm;
	G4double z0 = r * std::sin(t) + -137.51 * CLHEP::mm;
	//G4cerr << "Gun from " << x0 << ", " << y0 << ", " << z0 << G4endl; 
	fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0)*CLHEP::mm); // so then 0.1 mm^2 is the source
		G4double dy = 0.0;//(G4UniformRand() - 0.5) * 0.1;
		G4double dz = 0.0;//(G4UniformRand() - 0.5) * 0.1;
		fParticleGun->SetParticleMomentumDirection(G4ThreeVector(1., dy, dz));
	  */
	
	auto [i, j] = SamplePixel();
	double localZ = (j - numCols / 2.0) * pixelSize;
	double localY = (i - numRows / 2.0) * pixelSize; //z-y plane!
  G4cout << "selected index: " << i << "\t" << j << G4endl;
  G4cout << "selected pixel: " << localZ << "\t" << localY << G4endl;
	G4ThreeVector emissionPoint = G4ThreeVector(0, localY, localZ) + profileCenterWorld;
	
	fParticleGun->SetParticlePosition(emissionPoint);
	// Set direction etc.
		
 
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


void NormaPrimaryGeneratorAction::NormalizeAndBuildCDF() {
    // First, zero out values below threshold (e.g., 0.01)
    double threshold = 2e-05;
    for (auto& val : intensityMap) {
        if (val < threshold) val = 0.0;
    }
    // Normalize
    double sum = std::accumulate(intensityMap.begin(), intensityMap.end(), 0.0);
    if (sum == 0.0) {
        G4Exception("NormaPrimaryGeneratorAction::NormalizeAndBuildCDF",
                    "ZeroIntensityMap", FatalException,
                    "All intensity values are zero after thresholding.");
    }

    for (auto& val : intensityMap) val /= sum;

    cumulativeProb.resize(intensityMap.size());
    cumulativeProb[0] = intensityMap[0];
    for (size_t i = 1; i < intensityMap.size(); ++i) {
        cumulativeProb[i] = cumulativeProb[i - 1] + intensityMap[i];
    }
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
    double worldX = profileCenterWorld.x(); // X is fixed
	double worldZ = (cx - numCols / 2.0) * pixelSize; // in z-y plane!
    double worldY = (cy - numRows / 2.0) * pixelSize;

    return G4ThreeVector(worldX, worldY, worldZ) + profileCenterWorld;
}

std::pair<int, int> NormaPrimaryGeneratorAction::SamplePixel() {
	// Sampling from the 2D intensity map
		double r = discgenerate(cumulativeProb);
    //to see how many of each were generated
    ++map_test[r];
    for (const auto& [num, count] : map_test)
        std::cout << num << " generated " << std::setw(4) << count << " times\n";

    //double r = G4UniformRand();
    auto it = std::lower_bound(cumulativeProb.begin(), cumulativeProb.end(), r);
    int idx = std::distance(cumulativeProb.begin(), it);
    G4cout << "r: " << r << "\tidx: " << idx << G4endl;
    int i = idx / numCols;
    int j = idx % numCols;
    return {i, j};
}

void NormaPrimaryGeneratorAction::InitializeIntensityProfile(const std::string& filename, const G4ThreeVector& centerWorld, double pixelSizeMM) {
    profileCenterWorld = centerWorld;
    //pixelSize = pixelSizeMM * CLHEP::mm;
    G4cout << "pixelSize normal: " << pixelSize << G4endl;


    LoadIntensityCSV(filename);
    NormalizeAndBuildCDF();
    G4ThreeVector beamCenter = ComputeProfileCenter();
		fDist_disc = std::discrete_distribution<>(intensityMap.begin(), intensityMap.end());
		fGen.seed(time(0)); // if you want different results from different runs

    G4cout << "→ Laser profile loaded from: " << filename << G4endl;
    G4cout << "→ Profile dimensions: " << numCols << " × " << numRows << G4endl;
    G4cout << "→ Computed center of mass: " << beamCenter / mm << " mm" << G4endl;
}

double NormaPrimaryGeneratorAction::discgenerate(std::vector<double> &values) // discrete
{
  if (values.empty()) {
    G4cerr << "Error: Values vector is empty!" << G4endl;
    return -1;
  }

  int index = fDist_disc(fGen);
  if (index < 0 || index >= values.size()) {
    G4cerr << "Error: Generated index out of bounds!" << G4endl;
    return -1;
  }
	return values[index]; // Draw a number from the values vector
}
