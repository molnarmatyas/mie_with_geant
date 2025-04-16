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
/// \file Norma/include/NormaPrimaryGeneratorAction.hh
/// \brief Definition of the NormaPrimaryGeneratorAction class
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef NormaPrimaryGeneratorAction_h
#define NormaPrimaryGeneratorAction_h 1

#include "globals.hh"
#include "G4ParticleGun.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "ParametersStruct.hh"
#include <random>

class G4Event;
class NormaPrimaryGeneratorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class NormaPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  NormaPrimaryGeneratorAction();
  ~NormaPrimaryGeneratorAction() override;

  void GeneratePrimaries(G4Event*) override;

  void SetOptPhotonPolar();
  void SetOptPhotonPolar(G4double);

  G4ParticleGun* GetParticleGun() { return fParticleGun; }

  void SetParameters(Parameters p);

private:
  G4ParticleGun* fParticleGun = nullptr;
  NormaPrimaryGeneratorMessenger* fGunMessenger = nullptr;
  double LEDSize = 5 * CLHEP::um;
  // For intensity profile
  void InitializeIntensityProfile(const std::string& filename, const G4ThreeVector& centerWorld, double pixelSizeMM);
  void LoadIntensityCSV(const std::string& filename);
  void NormalizeAndBuildCDF();
  G4ThreeVector ComputeProfileCenter();
  std::pair<int, int> SamplePixel();

  std::vector<double> intensityMap;
  std::vector<double> cumulativeProb;
	std::mt19937 fGen;
	std::discrete_distribution<int> fDist_disc;
  std::map<double, int> map_test;
  int numRows = 0;
  int numCols = 0;
  double pixelSize = 4.5 * CLHEP::um; // CinCam CMOS-1203
  G4ThreeVector profileCenterWorld = G4ThreeVector(); // Where the center should be in Geant4
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif /*NormaPrimaryGeneratorAction_h*/
