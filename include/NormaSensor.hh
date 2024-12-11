#ifndef MIENORMASENSOR_HH
#define MIENORMASENSOR_HH

#include "G4VSensitiveDetector.hh"
#include "G4RunManager.hh"
#include "G4AnalysisManager.hh"

class PhotoSensitiveDetector : public G4VSensitiveDetector
{
public:
    PhotoSensitiveDetector(G4String);
    ~PhotoSensitiveDetector();

private:
    virtual G4bool ProcessHits(G4Step *, G4TouchableHistory *);
};

#endif
