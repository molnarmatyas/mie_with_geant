//
// $Id: mieNormaNumericReader.cc  $
//
/// \file mieNormaNumericReader.cc
/// \brief Implementation of the mieNormaNumericReader class

#include "mieNormaNumericReader.hh"
#include "mieNormaDetectorConstruction.hh"

#include "G4Event.hh"
#include "G4LogicalVolume.hh"
#include "G4RunManager.hh"
#include "G4Step.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

mieNormaNumericReader::mieNormaNumericReader()
// fScoringVolume(0)
{
  const mieNormaDetectorConstruction* detectorConstruction
        = static_cast<const mieNormaDetectorConstruction*>
          (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  setXSect("diff_cross_0-180_10000pt.txt");
  G4cout << detectorConstruction->GetRadius()*1000 << G4endl;
	int colNum = (int)(detectorConstruction->GetRadius()*10000) / 5 - 5;

  for (int i = 0; i < 10000; i++) {
    fWeights.push_back(fMieXSect[colNum + (11) * i]);
    //G4cout << fMieXSect[colNum + (11) * i] << G4endl;
  }
	//std::random_device rd;
  fGen = std::mt19937(time(0));
  fDist = std::discrete_distribution<>(fWeights.begin(), fWeights.end());
  fGen.seed(time(0)); // if you want different results from different runs
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

mieNormaNumericReader::~mieNormaNumericReader() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void mieNormaNumericReader::setXSect(const char *filename) {
  int column_number = 12;
  // Create arrays to hold the data for each column
  double theta;
  double values; // Assuming column_number includes the angle column

  std::ifstream file(filename);
  std::string line;
  while (getline(file, line)) {
    std::stringstream ss(line);

    // Read the first column as the angle
    ss >> theta;
    fTheta.push_back(theta);

    // Read the rest of the columns as values
    for (int i = 0; i < column_number - 1; i++) {
      ss >> values;
      fMieXSect.push_back(values);
    }
  }
}

double mieNormaNumericReader::generate(std::vector<double> &values) {
  return values[fDist(fGen)]; // Draw a number from the values vector
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
