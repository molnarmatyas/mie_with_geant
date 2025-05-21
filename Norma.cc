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
/// \file Norma/Norma.cc
/// \brief Main program of the Norma example
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
// Description: Test of Continuous Process G4Cerenkov
//              and RestDiscrete Process G4Scintillation
//              -- Generation Cerenkov Photons --
//              -- Generation Scintillation Photons --
//              -- Transport of optical Photons --
// Version:     5.0
// Created:     1996-04-30
// Author:      Juliet Armstrong
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "FTFP_BERT.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4OpticalPhysics.hh"
#include "G4RunManagerFactory.hh"
#include "G4Types.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "NormaActionInitialization.hh"
#include "NormaDetectorConstruction.hh"
#include "ThreadSafeWriter.hh"
#include <execinfo.h>
#include <filesystem>
#include <signal.h>
#include "G4AutoDelete.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
namespace
{
void PrintUsage()
{
	G4cerr << " Usage: " << G4endl;
	G4cerr << " Norma  [-m macro ] [-u UIsession] [-t nThreads] [-r seed] " << G4endl;
	G4cerr << "   note: -t option is available only for multi-threaded mode." << G4endl;
	G4cerr << "   advanced note: CELL_RADIUS_UM and -b must be the same number." << G4endl;
}
} // namespace

std::string getAvailableFilename(const std::string &filename)
{
	// Split filename into name and extension
	std::string name = filename;
	std::string extension;

	std::size_t pos = filename.find_last_of(".");
	if (pos != std::string::npos)
	{
		name = filename.substr(0, pos);
		extension = filename.substr(pos);
	}

	// Check if file exists and find the next available filename
	std::string newFilename = filename;
	int iterator = 0;
	while (std::filesystem::exists(newFilename))
	{
		newFilename = name + "_" + std::to_string(iterator++) + extension;
	}

	return newFilename;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void handler(int sig)
{
	void *array[50];
	size_t size;

	// get void*'s for all entries on the stack
	size = backtrace(array, 50);

	// print out all the frames to stderr
	fprintf(stderr, "Error: signal %d:\n", sig);
	backtrace_symbols_fd(array, size, STDERR_FILENO);
	exit(1);
}

int main(int argc, char **argv)
{
	// Evaluate arguments
	//
	signal(SIGSEGV, handler);
	if (argc > 11)
	{
		PrintUsage();
		return 1;
	}
	G4String macro;
	G4String session;
#ifdef G4MULTITHREADED
	G4int nThreads = 0;
#endif

	G4long myseed = 345354;

	/* A Parameters struktúrába adhatsz hozzá új változókat
	 * A DetectorConstruction SetParameters függvényében el tudod menteni
	 * egy belső változóba.
	 */
	Parameters p;

	p.bubbleRadius = 3 * CLHEP::um;

	p.g = 0.99;

	p.p = 0; // by default it's orb if it's 1 it is a ploycone

	for (G4int i = 1; i < argc; i++)
	{
		std::cout << argv[i] << std::endl;
	}

	for (G4int i = 1; i < argc; i = i + 2)
	{
		if (G4String(argv[i]) == "-m")
			macro = argv[i + 1];
		else if (G4String(argv[i]) == "-u")
			session = argv[i + 1];
		else if (G4String(argv[i]) == "-r")
			myseed = atoi(argv[i + 1]);
		// parancssori argumentumból lehet változtatni a radiust pl itt és így:
		else if (G4String(argv[i]) == "-b")
			p.bubbleRadius = atof(argv[i + 1]) * CLHEP::um;
		else if (G4String(argv[i]) == "-g")
			p.g = atof(argv[i + 1]) / 1000;
		else if (G4String(argv[i]) == "-p")
			p.p = atoi(argv[i + 1]);
#ifdef G4MULTITHREADED
		else if (G4String(argv[i]) == "-t")
		{
			nThreads = G4UIcommand::ConvertToInt(argv[i + 1]);
		}
#endif
		else
		{
			PrintUsage();
			return 1;
		}
	}
  /*
  if(p.bubbleRadius != std::stof(std::getenv("CELL_RADIUS_UM")) * CLHEP::um) {
    G4cerr << "Radius variables must match" << G4endl;
		PrintUsage();
    return 1;
  }
  */

	const char* numeric_simulation_file_path = nullptr; 
  std::string filename;
	if ((numeric_simulation_file_path =  std::getenv("NUMERIC_MIE_FPATH"))) {
    std::string numeric_path = numeric_simulation_file_path;
	  filename = "output" + numeric_path + std::to_string((int)(p.bubbleRadius * 1000000)) + "_" +
	  					   std::to_string((int)(p.g * 1000)) + "_" + std::to_string((int)(p.p)) + ".txt";
  }
  else {
	  filename = "output" + std::to_string((int)(p.bubbleRadius * 1000000)) + "_" +
	  					   std::to_string((int)(p.g * 1000)) + "_" + std::to_string((int)(p.p)) + ".txt";
  }
	std::string availableFilename = getAvailableFilename(filename);

	// ez létrehoz egy output_x.txt-t és abba ír, thread safe, át van passzolva a stepping actionnek
	// ott van is egy példa
	// std::stringstream -el lehet sztringet csinálni és beleírni
//  ThreadSafeWriter *writer = new ThreadSafeWriter(availableFilename);
	ThreadSafeWriter writer(availableFilename);
  	G4AutoDelete::Register(&writer);

	// Instantiate G4UIExecutive if interactive mode
	G4UIExecutive *ui = nullptr;
	if (macro.size() == 0)
	{
		ui = new G4UIExecutive(argc, argv);
	}

	// Construct the default run manager
	auto runManager = G4RunManagerFactory::CreateRunManager(/*G4RunManagerType::SerialOnly*/);
#ifdef G4MULTITHREADED
	if (nThreads > 0)
		runManager->SetNumberOfThreads(nThreads);
#endif

	// Seed the random number generator manually
	G4Random::setTheSeed(myseed);

	// Set mandatory initialization classes
	//
	// Detector construction
	NormaDetectorConstruction *detConst = new NormaDetectorConstruction();
	detConst->SetParameters(p);

	runManager->SetUserInitialization(detConst);
	// Physics list
	G4VModularPhysicsList *physicsList = new FTFP_BERT;
	physicsList->ReplacePhysics(new G4EmStandardPhysics_option4());
	auto opticalPhysics = new G4OpticalPhysics();
	physicsList->RegisterPhysics(opticalPhysics);
	runManager->SetUserInitialization(physicsList);

	NormaActionInitialization *actionInit = new NormaActionInitialization();

	actionInit->SetWriter(&writer);

	runManager->SetUserInitialization(actionInit);

	G4VisManager *visManager = new G4VisExecutive("Quiet");
	visManager->Initialize();

	G4UImanager *UImanager = G4UImanager::GetUIpointer();

	if (macro.size())
	{
		G4String command = "/control/execute ";
		UImanager->ApplyCommand(command + macro);
	}
	else // Define UI session for interactive mode
	{
		UImanager->ApplyCommand("/control/execute simplevis.mac");
		ui->SessionStart();
		delete ui;
	}

	runManager = G4RunManager::GetRunManager();
	if (visManager) delete visManager;
	if (runManager) delete runManager;

	return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
