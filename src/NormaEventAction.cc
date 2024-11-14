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
/// \file optical/Norma/src/NormaEventAction.cc
/// \brief Implementation of the NormaEventAction class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "NormaEventAction.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "NormaRun.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void NormaEventAction::BeginOfEventAction(const G4Event *)
{
	fRayleigh = 0;
	fAbsorption = 0;
	fMie = 0;
	fBoundary = 0;
	isMIE = 0;
	beforeTheta = 0;
	afterTheta = 0;
	beforePhi = 0;
	afterPhi = 0;
	isFilled = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void NormaEventAction::EndOfEventAction(const G4Event *event)
{
	auto run = static_cast<NormaRun *>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
	run->AddRayleigh(fRayleigh);
	run->AddAbsorption(fAbsorption);
	run->AddMie(fMie);
	run->AddBoundary(fBoundary);
	/*
	if (isMIE){
	  ss.str(""); //ss.clear();
	  ss<< "MIEEA "<<event->GetEventID()<<" "<<  beforeTheta<<" "<<beforePhi<<" "<<afterTheta<<" "<<afterPhi;
	  if (writer) writer -> write(ss.str());
	  else std::cout << "EA File writer is null" << std::endl;
	}

	ss.str(""); //ss.clear();
	ss<< "ALLEA "<<event->GetEventID()<<" "<<  beforeTheta<<" "<<beforePhi<<" "<<afterTheta<<" "<<afterPhi;
	if (writer) writer -> write(ss.str());
	else std::cout << "EA File writer is null" << std::endl;
  */
}

void NormaEventAction::SetWriter(ThreadSafeWriter *writer)
{
	this->writer = writer;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
