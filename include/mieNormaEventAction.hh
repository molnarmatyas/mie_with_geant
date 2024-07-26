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
// $Id: mieNormaEventAction.hh 75216 2013-10-29 16:08:11Z gcosmo $
//
/// \file mieNormaEventAction.hh
/// \brief Definition of the mieNormaEventAction class

#ifndef mieNormaEventAction_h
#define mieNormaEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include <vector>

/// Event action class
///

class mieNormaEventAction : public G4UserEventAction
{
  public:
    mieNormaEventAction();
    virtual ~mieNormaEventAction();
    
    virtual void BeginOfEventAction(const G4Event* event);
    virtual void EndOfEventAction(const G4Event* event);

    void SaveAngles(G4double thismag, G4double thisphi, G4double thistheta) { fMag = thismag; fPhi = thisphi; fTheta = thistheta; }
    void AddRayleigh() { ++fRayleigh; }
    void AddAbsorption() { ++fAbsorption; }
    void AddMie() { ++fMie; }
    void AddBoundary() { ++fBoundary; }
  
  private:
    G4int fRayleigh = 0;
    G4int fAbsorption = 0;
    G4int fMie = 0;
    G4int fBoundary = 0;
    G4double fMag;
    G4double fPhi;
    G4double fTheta;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
