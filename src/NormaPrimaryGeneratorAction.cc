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
	fParticleGun->SetParticlePosition(G4ThreeVector(-29.35 * CLHEP::mm, 96 * CLHEP::mm, -158.98 * CLHEP::mm));
	G4double dy = (G4UniformRand() - 0.5) * 0.1;
	G4double dz = (G4UniformRand() - 0.5) * 0.1;
	// fParticleGun->SetParticleMomentumDirection(G4ThreeVector(1., dy, dz));
	fParticleGun->SetParticleMomentumDirection(G4ThreeVector(-1., 0, 0));
	fParticleGun->SetParticleEnergy(1.884 * eV);
  //G4cout << "Gun from " << fParticleGun->GetParticlePosition().x() << ", " << fParticleGun->GetParticlePosition().y() << ", " << fParticleGun->GetParticlePosition().z() << G4endl; 

  //fParticleGun->SetParticlePosition(G4ThreeVector(0,0,0));
  //fParticleGun->GeneratePrimaryVertex(anEvent);
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
  G4double x0 = -17 * CLHEP::um;
  G4double y0 = r * std::cos(t);
  G4double z0 = r * std::sin(t);
  //G4cerr << "Gun from " << x0 << ", " << y0 << ", " << z0 << G4endl; 
  fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0)*CLHEP::mm); // so then 0.1 mm^2 is the source
  */
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
