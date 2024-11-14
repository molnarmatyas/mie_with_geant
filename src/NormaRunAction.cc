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
/// \file Norma/src/NormaRunAction.cc
/// \brief Implementation of the NormaRunAction class
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#include "NormaRunAction.hh"
#include "G4ParticleDefinition.hh"
#include "G4Run.hh"
#include "NormaPrimaryGeneratorAction.hh"
#include "NormaRun.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
NormaRunAction::NormaRunAction(NormaPrimaryGeneratorAction *prim) : G4UserRunAction(), fRun(nullptr), fPrimary(prim)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4Run *NormaRunAction::GenerateRun()
{
	fRun = new NormaRun();
	return fRun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void NormaRunAction::BeginOfRunAction(const G4Run *)
{
	if (fPrimary)
	{
		G4ParticleDefinition *particle = fPrimary->GetParticleGun()->GetParticleDefinition();
		G4double energy = fPrimary->GetParticleGun()->GetParticleEnergy();
		fRun->SetPrimary(particle, energy);
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void NormaRunAction::EndOfRunAction(const G4Run *)
{
	if (isMaster)
		fRun->EndOfRun();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
