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
/// \file Norma/include/NormaRun.hh
/// \brief Definition of the NormaRun class
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef NormaRun_h
#define NormaRun_h 1

#include "G4Run.hh"

class G4ParticleDefinition;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class NormaRun : public G4Run
{
 public:
  NormaRun() = default;
  ~NormaRun() override = default;

  void SetPrimary(G4ParticleDefinition* particle, G4double energy);

  void AddCerenkov(G4double n)
  {
    fCerenkovCounter += n;
    fCerenkov2 += n * n;
  };
  void AddScintillation(G4double n)
  {
    fScintillationCounter += n;
    fScintillation2 += n * n;
  };
  void AddRayleigh(G4double n)
  {
    fRayleighCounter += n;
    fRayleigh2 += n * n;
  };
  void AddAbsorption(G4double n)
  {
    fAbsorptionCounter += n;
    fAbsorption2 += n * n;
  };
  void AddMie(G4double n)
  {
    fMieCounter += n;
    fMie2 += n * n;
  };
  void AddBoundary(G4double n)
  {
    fBoundaryCounter += n;
    fBoundary2 += n * n;
  };

  void Merge(const G4Run*) override;
  void EndOfRun();

 private:
  G4ParticleDefinition* fParticle = nullptr;

  G4double fCerenkovCounter = 0.;
  G4double fCerenkov2 = 0.;
  G4double fScintillationCounter = 0.;
  G4double fScintillation2 = 0.;
  G4double fRayleighCounter = 0.;
  G4double fRayleigh2 = 0.;
  G4double fAbsorptionCounter = 0.;
  G4double fAbsorption2 = 0.;
  G4double fMieCounter = 0.;
  G4double fMie2 = 0.;
  G4double fBoundaryCounter = 0.;
  G4double fBoundary2 = 0.;
  G4double fEnergy = -1.;
};
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#endif
