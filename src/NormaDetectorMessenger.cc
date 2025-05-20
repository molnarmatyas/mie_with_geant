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

#include "NormaDetectorMessenger.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIdirectory.hh"
#include "NormaDetectorConstruction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

NormaDetectorMessenger::NormaDetectorMessenger(G4VUserDetectorConstruction *detcon)
	: G4UImessenger(), fNormaDetCon(detcon)
{
	fDetConDir = new G4UIdirectory("/Norma/DetectorConstruction/");
	fDetConDir->SetGuidance("Configuring Detector Construction");

	fVerboseCmd = new G4UIcmdWithABool("/Norma/DetectorConstruction/enableVerbose", this);
	fVerboseCmd->SetGuidance("Set flag for enabling verbose diagnostic printout");
	fVerboseCmd->SetDefaultValue(false);
	fVerboseCmd->AvailableForStates(G4State_PreInit);

	setOffsetCmd = new G4UIcmdWith3VectorAndUnit("/Norma/DetectorConstruction/Cell/setOffset", this);
    setOffsetCmd->SetGuidance("Set X/Y/Z offsets for the Cell.");
    setOffsetCmd->SetParameterName("x", "y", "z", false);
    setOffsetCmd->SetUnitCategory("Length");
    setOffsetCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

NormaDetectorMessenger::~NormaDetectorMessenger()
{
	delete fDetConDir;
	delete fVerboseCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void NormaDetectorMessenger::SetNewValue(G4UIcommand *command, G4String newValue)
{
	auto dc1 = dynamic_cast<NormaDetectorConstruction *>(fNormaDetCon);
	if (dc1 != nullptr)
	{
		if (command == fVerboseCmd)
		{
			dc1->SetVerbose(fVerboseCmd->GetNewBoolValue(newValue));
		}
		if (command == setOffsetCmd) 
		{
			G4ThreeVector offset = setOffsetCmd->GetNew3VectorValue(newValue);
			dc1->SetCellOffset(offset.x(), offset.y(), offset.z());
			dc1->UpdateGeometry();  // Rebuild geometry
		}
	}
	
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
