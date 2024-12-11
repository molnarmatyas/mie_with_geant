#include "NormaSensor.hh"
//#include "G4AutoLock.hh"
//namespace { G4Mutex myMutex = G4MUTEX_INITIALIZER; }

PhotoSensitiveDetector::PhotoSensitiveDetector(G4String name) : G4VSensitiveDetector(name)
{
  //G4cout << "Sensitive detector created." << G4endl;
  //G4cout << "Detector active: " << isActive() << G4endl;
}

PhotoSensitiveDetector::~PhotoSensitiveDetector()
{}

G4bool PhotoSensitiveDetector::ProcessHits(G4Step *aStep, G4TouchableHistory *ROhist)
{
    //G4AutoLock lock(&myMutex);
    G4Track *track = aStep->GetTrack();

    G4StepPoint *preStepPoint = aStep->GetPreStepPoint();
    G4StepPoint *postStepPoint = aStep->GetPostStepPoint();

    G4ThreeVector posPhoton = preStepPoint->GetPosition();

    //G4cout << "Photon position in detector: " << posPhoton << G4endl;
    
    const G4VTouchable *touchable = aStep->GetPreStepPoint()->GetTouchable();
    G4int copyNo = touchable->GetCopyNumber();
    //G4cout << "Copy number of detector pos.: " << copyNo << G4endl;
    
    G4VPhysicalVolume *physVol = touchable->GetVolume();
    G4ThreeVector posDetector = physVol->GetTranslation();
    G4cout << "Detector position of scattered photon: " << posDetector << G4endl;

     
    G4int evt = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
    
    //G4AnalysisManager *man = G4AnalysisManager::Instance();
    auto man = G4AnalysisManager::Instance();
    //G4cout << "DEBUG fEvt: " << evt << G4endl;
    //man->SetNtupleMerging(true); // TESTING, does not resolve the segfault
    /* This part causes segfault idk why
    man->FillNtupleDColumn(1, 0, evt);
    G4cout << "Did I get here?" << G4endl; // nope
    man->FillNtupleDColumn(1, 1, posDetector.x());
    man->FillNtupleDColumn(1, 2, posDetector.y());
    man->FillNtupleDColumn(1, 3, posDetector.z());
    man->AddNtupleRow(1);
    */
   
    
    // Return just because yes
    return true; 
}
