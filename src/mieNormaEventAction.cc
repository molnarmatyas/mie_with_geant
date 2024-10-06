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
// $Id: mieNormaEventAction.cc 75117 2013-10-28 09:38:37Z gcosmo $
//
/// \file mieNormaEventAction.cc
/// \brief Implementation of the mieNormaEventAction class

#include "mieNormaEventAction.hh"
#include "mieNormaRun.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

mieNormaEventAction::mieNormaEventAction()
: G4UserEventAction(),
  fMag(0.),
  fPhi(0.),
  fTheta(0.)
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

mieNormaEventAction::~mieNormaEventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void mieNormaEventAction::BeginOfEventAction(const G4Event*)
{    
  fMag = 0.;
  fPhi = 0.;
  fTheta = 0.;
  fRayleigh   = 0;
  fAbsorption = 0;
  fMie        = 0;
  fBoundary   = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void mieNormaEventAction::EndOfEventAction(const G4Event*)
{   
  // accumulate statistics in mieNormaRun
  mieNormaRun* run = static_cast<mieNormaRun*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  run->SaveAnglesToVector(fMag,fPhi,fTheta);
  run->AddRayleigh(fRayleigh);
  run->AddAbsorption(fAbsorption);
  run->AddMie(fMie);
  run->AddBoundary(fBoundary);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
