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
/// \file optical/Norma/include/NormaEventAction.hh
/// \brief Definition of the NormaEventAction class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef NormaEventAction_h
#define NormaEventAction_h 1

#include "globals.hh"
#include "G4UserEventAction.hh"
#include "ThreadSafeWriter.hh"
#include <sstream>      // std::stringstream

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class NormaEventAction : public G4UserEventAction
{
public:
  NormaEventAction() = default;
  ~NormaEventAction() override = default;

  void BeginOfEventAction(const G4Event*) override;
  void EndOfEventAction(const G4Event*) override;

  void AddRayleigh() { ++fRayleigh; }
  void AddAbsorption() { ++fAbsorption; }
  void AddMie() { ++fMie; }
  void AddBoundary() { ++fBoundary; }
  void SetMIE() {isMIE=1;}
  G4int GetMIE() {return isMIE;}

  void SetBeforeTheta(G4double theta) {beforeTheta = theta;}
  void SetAfterTheta(G4double theta) {afterTheta = theta;}
  void SetBeforePhi(G4double phi) {beforePhi = phi;}
  void SetAfterPhi(G4double phi) {afterPhi = phi;}
  void SetWriter(ThreadSafeWriter *writer);
  G4int isFilled = 0;


private:
  G4int fRayleigh = 0;
  G4int fAbsorption = 0;
  G4int fMie = 0;
  G4int fBoundary = 0;
  G4int isMIE = 0;
  G4double beforeTheta = 0;
  G4double afterTheta = 0;
  G4double beforePhi = 0;
  G4double afterPhi = 0;
  ThreadSafeWriter *writer;
  std::stringstream ss;


};
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#endif
