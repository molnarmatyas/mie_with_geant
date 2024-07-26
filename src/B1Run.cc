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
// $Id: B1Run.cc 66536 2012-12-19 14:32:36Z ihrivnac $
//
/// \file B1Run.cc
/// \brief Implementation of the B1Run class

#include "B1Run.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1Run::B1Run() : G4Run() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1Run::~B1Run() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1Run::Merge(const G4Run *run) {
  const B1Run *localRun = static_cast<const B1Run *>(run);
  magvector.insert(magvector.end(), localRun->magvector.begin(),
                   localRun->magvector.end());
  phivector.insert(phivector.end(), localRun->phivector.begin(),
                   localRun->phivector.end());
  thetavector.insert(thetavector.end(), localRun->thetavector.begin(),
                     localRun->thetavector.end());

  fParticle = localRun->fParticle;
  fEnergy = localRun->fEnergy;
  fRayleighCounter += localRun->fRayleighCounter;
  fRayleigh2 += localRun->fRayleigh2;
  fAbsorptionCounter += localRun->fAbsorptionCounter;
  fAbsorption2 += localRun->fAbsorption2;
  fMieCounter += localRun->fMieCounter;
  fMie2 += localRun->fMie2;
  fBoundaryCounter += localRun->fBoundaryCounter;
  fBoundary2 += localRun->fBoundary2;

  G4Run::Merge(run);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1Run::SaveAnglesToVector(G4double thismag, G4double thisphi,
                               G4double thistheta) {
  magvector.push_back(thismag);
  phivector.push_back(thisphi);
  thetavector.push_back(thistheta);
}

void B1Run::PrintAngles() const {
  auto man = G4AnalysisManager::Instance();
  G4cout << "Printing angles list" << G4endl;
  for (unsigned int ientry = 0; ientry < phivector.size(); ientry++) {
    /*G4cout << "mag,phi,theta:\t" << magvector.at(ientry) << "\t"
           << phivector.at(ientry) << "\t" << thetavector.at(ientry) << G4endl;*/
    man->FillNtupleIColumn(0, 1);
    man->FillNtupleDColumn(1, magvector.at(ientry));
    man->FillNtupleDColumn(2, phivector.at(ientry));
    man->FillNtupleDColumn(3, thetavector.at(ientry));
    man->AddNtupleRow(0);
  }
}
void B1Run::EndOfRun() {
  if (numberOfEvent == 0)
    return;
  auto TotNbofEvents = G4double(numberOfEvent);
  G4double totalRayleigh = fRayleighCounter;
  G4double totalAbs = fAbsorptionCounter;
  G4double totalMie = fMieCounter;
  G4double totalBoundary = fBoundaryCounter;

  fRayleighCounter /= TotNbofEvents;
  fRayleigh2 /= TotNbofEvents;
  G4double rmsRayleigh = fRayleigh2 - fRayleighCounter * fRayleighCounter;
  if (rmsRayleigh > 0.)
    rmsRayleigh = std::sqrt(rmsRayleigh);
  else
    rmsRayleigh = 0.;

  fAbsorptionCounter /= TotNbofEvents;
  fAbsorption2 /= TotNbofEvents;
  G4double rmsAbsorption =
      fAbsorption2 - fAbsorptionCounter * fAbsorptionCounter;
  if (rmsAbsorption > 0.)
    rmsAbsorption = std::sqrt(rmsAbsorption);
  else
    rmsAbsorption = 0.;
  fMieCounter /= TotNbofEvents;
  fMie2 /= TotNbofEvents;
  G4double rmsMie = fMie2 - fMieCounter * fMieCounter;
  if (rmsMie > 0.)
    rmsMie = std::sqrt(rmsMie);
  else
    rmsMie = 0.;

  fBoundaryCounter /= TotNbofEvents;
  fBoundary2 /= TotNbofEvents;
  G4double rmsBoundary = fBoundary2 - fBoundaryCounter * fBoundaryCounter;
  if (rmsBoundary > 0.)
    rmsBoundary = std::sqrt(rmsBoundary);
  else
    rmsBoundary = 0.;

  G4int prec = G4cout.precision(3);
  G4cout << "\n ======================== run summary ======================\n";

  //G4cout << "Primary particle was: " << fParticle->GetParticleName()
  //       << " with energy " << G4BestUnit(fEnergy, "Energy") << "." << G4endl;
  G4cout << "Number of events: " << numberOfEvent << G4endl;
  G4cout << "Total number of optical Rayleigh interactions: " << totalRayleigh << "\tper event: "
         << fRayleighCounter << " +- " << rmsRayleigh << G4endl;
  G4cout << "Total number of optical absorption interactions: " << totalAbs << "\tper event: "
         << fAbsorptionCounter << " +- " << rmsAbsorption << G4endl;
  G4cout << "Total number of optical Mie interactions: " << totalMie <<"\tper event: "
         << fMieCounter << " +- " << rmsMie << G4endl;
  G4cout << "Total number of optical boundary interactions: " << totalBoundary << "\tper event: "
         << fBoundaryCounter << " +- " << rmsBoundary << G4endl;
  G4cout << G4endl;
  G4cout.precision(prec);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
