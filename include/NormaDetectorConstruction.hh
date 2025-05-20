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
/// \file Norma/include/NormaDetectorConstruction.hh
/// \brief Definition of the NormaDetectorConstruction class
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef NormaDetectorConstruction_h
#define NormaDetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include <CLHEP/Units/SystemOfUnits.h>
#include "ParametersStruct.hh"
#include "G4NistManager.hh"
#include "G4Polycone.hh"
#include "G4Ellipsoid.hh"
#include "G4RunManager.hh"

//#include "NormaSensor.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class NormaDetectorMessenger;

class NormaDetectorConstruction : public G4VUserDetectorConstruction
{
public:
  NormaDetectorConstruction();
  ~NormaDetectorConstruction() override;

  G4VPhysicalVolume *Construct() override;

  void SetVerbose(G4bool verbose);
  G4bool IsVerbose() const;

  void SetParameters(Parameters p);
  void SetCellOffset(G4double x, G4double y, G4double z);
  void UpdateGeometry();

private:
  void PrintError(G4String);

  NormaDetectorMessenger *fDetectorMessenger = nullptr;

  G4double fWorld_x = 40 * CLHEP::cm;
  G4double fWorld_y = 40 * CLHEP::cm;
  G4double fWorld_z = 40 * CLHEP::cm;

  G4double fScreen_x = 0.5 * CLHEP::cm;
  G4double fScreen_y = 40 * CLHEP::cm;
  G4double fScreen_z = 40 * CLHEP::cm;

  G4VPhysicalVolume *bubble_phys;

  /*
  G4double fExpHall_x = 10. * CLHEP::m;
  G4double fExpHall_y = 10. * CLHEP::m;
  G4double fExpHall_z = 10. * CLHEP::m;

  G4double fTank_x = 5. * CLHEP::m;
  G4double fTank_y = 5. * CLHEP::m;
  G4double fTank_z = 5. * CLHEP::m;
  */

  G4double fBubble_r = 3 * CLHEP::um;
  G4ThreeVector fBubble_def_pos; //Cell default position
  G4ThreeVector fBubble_additional_offset;
  G4double mieFg = 0.99;
  G4int isPolycone = 0;
  G4bool fVerbose = false;

  G4LogicalVolume *logicDetector;     // For sensitive detector array and field (EM etc.)
  //virtual void ConstructSDandField();
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif /*NormaDetectorConstruction_h*/
