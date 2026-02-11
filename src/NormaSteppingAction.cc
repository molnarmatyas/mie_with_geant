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
//
/// \file NormaSteppingAction.cc
/// \brief Implementation of the NormaSteppingAction class

#include "NormaSteppingAction.hh"
#include "G4Event.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4OpMieHG.hh"
#include "G4OpticalPhoton.hh"
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "NormaRun.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

NormaSteppingAction::NormaSteppingAction(NormaEventAction *event) : G4UserSteppingAction(), fEventAction(event)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void NormaSteppingAction::UserSteppingAction(const G4Step *step)
{
	static G4ParticleDefinition *opticalphoton = G4OpticalPhoton::OpticalPhotonDefinition();

	const G4ParticleDefinition *particleDef = step->GetTrack()->GetDynamicParticle()->GetParticleDefinition();

	G4Track *track = step->GetTrack();
	G4Material *material = track->GetMaterial();
	G4MaterialPropertiesTable *mpt = material->GetMaterialPropertiesTable();
	try
	{
		if (mpt->ConstPropertyExists(kMIEHG_FORWARD_RATIO) && false)
		{

			std::cout << "kMIEHG_FORWARD_RATIO const property for Material:  " << material->GetName() << std::endl;
			std::cout << "forwardRatio " << mpt->GetConstProperty(kMIEHG_FORWARD_RATIO) << std::endl;
		}
		else if (false)
		{
			std::cout << "kMIEHG_FORWARD_RATIO const property not found for Material: " << material->GetName()
					  << std::endl;
			std::cout << "-----------------------------------------------DONE--------"
						 "---------------------------------------"
					  << std::endl;
		}
	}
	catch (const std::exception &e)
	{
		std::cout << e.what() << '\n';
	}

	G4OpMieHG *mieProcess;
	G4ProcessManager *processManager = track->GetDefinition()->GetProcessManager();
	if (processManager)
	{
		G4int numProcesses = processManager->GetProcessListLength();
		G4ProcessVector *processes = processManager->GetProcessList();
		for (G4int i = 0; i < numProcesses; ++i)
		{
			G4VProcess *process = (*processes)[i];
			if (process->GetProcessName() == "OpMieHG")
			{
				mieProcess = dynamic_cast<G4OpMieHG *>(process);
        /*
				if (mieProcess)
				{
					G4cout << "Mie process found" << G4endl;
				}
        */
			}
		}
	}

	G4StepPoint *endPoint = step->GetPostStepPoint();
	const G4VProcess *pds = endPoint->GetProcessDefinedStep();
	G4String procname = pds->GetProcessName();
	G4StepPoint *startPoint = step->GetPreStepPoint();
	const G4VProcess *pds2 = startPoint->GetProcessDefinedStep();

	G4String procname2 = pds->GetProcessName();
	// trackID
	G4int trackID = track->GetTrackID();
	// pre step porsition
	G4double preX = step->GetPreStepPoint()->GetPosition().x();
	G4double preY = step->GetPreStepPoint()->GetPosition().y();
	G4double preZ = step->GetPreStepPoint()->GetPosition().z();

	// Fixed preStep position
	//preX = 0.0;
	//preY = 0.0;
	//preZ = 0.0;

	// post step porsition
  G4ThreeVector postPos = step->GetPostStepPoint()->GetPosition();
	G4double postX = postPos.x();
	G4double postY = postPos.y();
	G4double postZ = postPos.z();

	G4VPhysicalVolume *prevolume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
	G4VPhysicalVolume *postvolume = step->GetPostStepPoint()->GetTouchableHandle()->GetVolume();

	// seems to be useless for theta and phi, they are something different
	G4ThreeVector momentumdirection = track->GetMomentumDirection();
	G4double theta = momentumdirection.theta();
	G4double phi = momentumdirection.phi();

  G4double gunPosX = 0.0;
  G4double gunPosY = 0.0;
  G4double gunPosZ = 0.0;

  //Threevectors related to gun
  G4ThreeVector gunPosition = G4RunManager::GetRunManager()->GetCurrentEvent()->GetPrimaryVertex()->GetPosition();
  G4ThreeVector momdir = G4ThreeVector(1, 0., 0.); 
  gunPosX = gunPosition.x();
  gunPosY = gunPosition.y();
  gunPosZ = gunPosition.z();

  G4double alpha = momdir.angle(momentumdirection);

//  G4cout << "PrimVertex from " << gunPosX << ", " << gunPosY << ", " << gunPosZ << G4endl; 

	//G4cout << "MomentumDirTheta Phi: " << theta << " | " << phi << G4endl;

	G4double phi2 = 0.0;
	G4double theta2;

	if (abs(postX - preX) < 0.00001)
	{
		theta2 = 0.0;
	}
	else
	{
		if ((postX - preX) > 0.0 && (postY - preY) > 0.0)
		{
			theta2 = atan((postY - preY) / (postX - preX)) / M_PI * 360;
		}
		if ((postX - preX) < 0.0 && (postY - preY) > 0.0)
		{
			theta2 = 180 - atan((postY - preY) / -(postX - preX)) / M_PI * 360;
		}
		if ((postX - preX) < 0.0 && (postY - preY) < 0.0)
		{
			theta2 = 180 + atan(-(postY - preY) / -(postX - preX)) / M_PI * 360;
		}
		if ((postX - preX) > 0.0 && (postY - preY) < 0.0)
		{
			theta2 = atan(-(postY - preY) / (postX - preX)) / M_PI * 360;
		}
	}
	if (abs(postY - preY) < 0.00001)
	{
		phi2 == 0.0;
	}
	else
	{
		phi2 = atan((postZ - preZ) / (postY - preY)) / M_PI * 360;
	}

	if (procname == "OpMieHG")
	{
		genTheta = ((G4OpMieHG *)pds)->radThetaGen;
	}
  if(prevolume->GetName() == "Bubble_dis_bnd_proc" && postvolume->GetName() == "Saltywater") {
	  //G4cout << "Zero point based preZ, postZ, preY, postY: " << preZ << " | " << postZ << " | " << preY << " | " << postY << G4endl;
		phi3 = atan2((postY - preY),(postZ - preZ)) / M_PI * 180.0;
    G4cout << "genTheta: " << genTheta * 180.0 / M_PI;
    G4cout << "phi3: " << phi3 << G4endl;
  }
  /*
	std::cout << G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID() << " " << prevolume->GetName() << " "
			  << postvolume->GetName() << " momdir: " << momentumdirection.x() << " " << momentumdirection.y() << " "
			  << momentumdirection.z() << " "
			  << " " << theta2 << " " << phi2 << " Pre " << preX / CLHEP::um << " " << preY / CLHEP::um << " "
			  << preZ / CLHEP::um << " Post " << postX / CLHEP::um << " " << postY / CLHEP::um << " "
			  << postZ / CLHEP::um << " " << procname << " " << material->GetName() << std::endl; // 0 world Bubble Air
  */

	if (fEventAction->isFilled == 0 && prevolume->GetName() == "world" && postvolume->GetName() == "Bubble")
	{
		fEventAction->SetBeforeTheta(theta2);
		fEventAction->SetBeforePhi(phi2);
		fEventAction->isFilled = 1;
	}

	if (postvolume->GetName() == "world")
	{
		if (mieProcess)
		{
			mieProcess->ResetMieCounter();
			//G4cout << "Mie counter explicitly reset to zero" << G4endl;
		}
	}
  G4ThreeVector Vn = G4ThreeVector(1, 1, 1);

  //calculating angle after the track exits the bubble, should be equal to theta...
	if (prevolume->GetName() == "Bubble_dis_bnd_proc" && postvolume->GetName() == "world")
	{
    G4ThreeVector posAfterBubble = postPos;
    G4ThreeVector momAfterBubble = momentumdirection;
    /*
	  std::cout << "Getting direction for alpha" << G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID() << " " << prevolume->GetName() << " "
			  << postvolume->GetName() << " posdir: " << momAfterBubble.x() << " " << momAfterBubble.y() << " "
			  << momAfterBubble.z() << std::endl;
    */
    angleMan = acos((momdir).dot( (momAfterBubble / momAfterBubble.mag() ) ) );
    G4ThreeVector cross = momdir.cross(momAfterBubble);
    if (Vn.dot(cross) < 0) { // Or > 0
      angleMan = -angleMan;
    }
    //G4cout << "momdir angle: " << angleMan << " | cross: " << cross << G4endl; 
  }

  int det_num = 0;

	if (
      ( prevolume->GetName() == "world" && postvolume->GetName() == "GS3-U3-23S6M-C_sensor") || 
      ( prevolume->GetName() == "vbpw34s_1" && postvolume->GetName() == "vbpw34s_1_sensor") ||
      ( prevolume->GetName() == "vbpw34s_2" && postvolume->GetName() == "vbpw34s_2_sensor") ||
      ( prevolume->GetName() == "Saltywater" && postvolume->GetName() == "Screen-at-cell") // FIXME beam profiler testing
     ) 
	{
		ss.str("");
		double theta3 = asin(sqrt(postY * postY + postZ * postZ) / sqrt(postX * postX + postY * postY + postZ * postZ));
    double r = sqrt(postY * postY + postZ * postZ);
    G4cout << "Zero point based theta and phi3: " << genTheta*180.0 / M_PI << " | " << phi3 << G4endl;
		//G4cout << "Zero point based x, y, z: " << postX << " | " << postY << " | " << postZ << G4endl;
		//G4cout << "Prim vertex based x, y, z: " << gunPosX << " | " << gunPosY << " | " << gunPosZ << G4endl;
    //G4cout << "Angle between prim vertex and mie: " << gunPosition.angle(posAfterBubble) * 180 / CLHEP::pi <<G4endl;
    //G4cout << "Angle between momdir prim vertex and mie: " << momdir.angle(momentumdirection) * 180 / CLHEP::pi <<G4endl;
    //G4cout << "Angle between prim vertex and mie manual 2: " << angleMan <<G4endl;
    //G4cout << "Angle Theta: " << theta3 * 180 / CLHEP::pi <<G4endl;

    if(postvolume->GetName() == "GS3-U3-23S6M-C_sensor") det_num = 0;
    if(postvolume->GetName() == "vbpw34s_1_sensor") det_num = 1;
    if(postvolume->GetName() == "vbpw34s_2_sensor") det_num = 2;
    if(postvolume->GetName() == "Screen-at-cell") det_num = 3;

		if (((G4OpMieHG *)pds)->generated)
		{
			ss << theta3 << ", " << genTheta << ", " << r << ", " << postX << ", " << postY << ", " << postZ << ", " << alpha << ", " << angleMan << ", " << phi3 << ", " << det_num;
			if (writer)
				writer->write(ss.str());
		}

		// fEventAction->SetAfterTheta(theta2);fEventAction->SetAfterPhi(phi2);
	}

	if (particleDef == opticalphoton)
	{

		if (procname == "OpRayleigh")
			fEventAction->AddRayleigh();
		else if (procname == "OpAbsorption")
			fEventAction->AddAbsorption();
		else if (procname == "OpMieHG")
		{
			fEventAction->SetMIE();
			// std::cout << "OpMie process detected" << "\t|\t";
			// std::cout << "Material type: " << material->GetName();
			double myMieAMX = ((G4OpMieHG *)pds)->myAMX;
			double myMieAMY = ((G4OpMieHG *)pds)->myAMY;
			double myMieAMZ = ((G4OpMieHG *)pds)->myAMZ;
			double myMieBMX = ((G4OpMieHG *)pds)->myBMX;
			double myMieBMY = ((G4OpMieHG *)pds)->myBMY;
			double myMieBMZ = ((G4OpMieHG *)pds)->myBMZ;
			double myCosTheta = ((G4OpMieHG *)pds)->cosTheta;
			double myG = ((G4OpMieHG *)pds)->theG;

			double scAngle = acos((myMieAMX * myMieBMX + myMieAMY * myMieBMY + myMieAMZ * myMieBMZ) /
								  (sqrt(myMieAMX * myMieAMX + myMieAMY * myMieAMY + myMieAMZ * myMieAMZ) *
								   sqrt(myMieBMX * myMieBMX + myMieBMY * myMieBMY + myMieBMZ * myMieBMZ)));
			double scAngleDeg = 360 * scAngle / (2.0 * CLHEP::pi);
			// a writer thread safe szóval írhatsz vele korlátlanul fájlba
			if (scAngleDeg < 175)
			{
				// ss.str(""); //ss.clear();
				//  ss <<"SA "<< G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID()
				//  <<" " << trackID <<" newmom " << myMieAMX <<" " << myMieAMY <<" " << myMieAMZ
				//  <<" oldmom "	<< myMieBMX<<" " << myMieBMY <<" " << myMieBMZ <<" "
				//  <<" scAngle "<<scAngle
				//  <<" Pre "<<preX/CLHEP::um<<" "<<preY/CLHEP::um<<" "<<preZ/CLHEP::um
				//  <<" Post "<<postX/CLHEP::um<<" "<<postY/CLHEP::um<<" "<<postZ/CLHEP::um
				//  <<" " << procname2 <<" "<< material->GetName()<<std::endl;
				// ss <<"MIEXS "<<G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID()
				//<<" "<<scAngleDeg<<" "<<(1- myG * myG)/pow(1 + myG * myG - 2 * myG * myCosTheta,1.5);//<<std::endl;
				//  !!!!!!!!!!!!!!!!!OpMIE variable how to write!!!!!!!!!!!!!!!
				//  ss << "myMieVar = " << myMieVar << std::endl;
				// if (writer) writer -> write(ss.str());
				// else std::cout << "SE File writer is null" << std::endl;
			}
			fEventAction->AddMie();
		}

		// for boundary scattering, process name in 'transportation'.
		// Need to check differently:
		if (endPoint->GetStepStatus() == fGeomBoundary)
		{
			G4OpBoundaryProcessStatus theStatus = Undefined;
			G4ProcessManager *opManager = opticalphoton->GetProcessManager();
			G4int n_proc = opManager->GetPostStepProcessVector(typeDoIt)->entries();
			G4ProcessVector *postStepDoItVector = opManager->GetPostStepProcessVector(typeDoIt);
			for (G4int i = 0; i < n_proc; ++i)
			{
				G4VProcess *currentProcess = (*postStepDoItVector)[i];

				auto opProc = dynamic_cast<G4OpBoundaryProcess *>(currentProcess);
				if (opProc)
					theStatus = opProc->GetStatus();
			}
			if (theStatus != Undefined && theStatus != NotAtBoundary && theStatus != StepTooSmall)
			{
				fEventAction->AddBoundary();
			}
		}
	}
}

void NormaSteppingAction::SetWriter(ThreadSafeWriter *writer)
{
	this->writer = writer;
}

G4ThreeVector NormaSteppingAction::GetGunPosition() {
    return G4ThreeVector(0,0,0);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
