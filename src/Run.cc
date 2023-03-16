
#include "Run.hh"

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4MTRunManager.hh"
#include "G4SDManager.hh"
#include "G4VPrimitiveScorer.hh"

#include "PenOutputPrinter.hh"

#include <map>
#include <vector>


//==============================================================================
Run::Run(PenOutputPrinter* outPrinter)
:G4Run()
{
  G4cout << "Creating Run object" << G4endl;
  fOutputPrinter = outPrinter;
  // 2d-dose-map
  fGridRZDoseRun = 0;
  fGridRZDoseRun2 = 0;
  // depth-dose
  fDepthDoseRun = 0;
  fDepthDoseRun2 = 0;

  // spc-enddet-##
  fSpcEnddetRun = 0;
  fSpcEnddetRun2 = 0;
  fSpcEnddetBinDEVec = 0;

  // emerging particles
  // polar-angle
  fPolarAngleElectRun = 0;
  fPolarAngleElectRun2 = 0;
  fPolarAngleGammaRun = 0;
  fPolarAngleGammaRun2 = 0;
  fPolarAnglePositRun = 0;
  fPolarAnglePositRun2 = 0;
  // energy-up
  fEmergEneUpElectRun = 0;
  fEmergEneUpElectRun2 = 0;
  fEmergEneUpGammaRun = 0;
  fEmergEneUpGammaRun2 = 0;
  fEmergEneUpPositRun = 0;
  fEmergEneUpPositRun2 = 0;
  // energy-down
  fEmergEneDownElectRun = 0;
  fEmergEneDownElectRun2 = 0;
  fEmergEneDownGammaRun = 0;
  fEmergEneDownGammaRun2 = 0;
  fEmergEneDownPositRun = 0;
  fEmergEneDownPositRun2 = 0;
  

  // 2d-dose-map & depth-dose info
  if (fOutputPrinter->GetRZDoseMapFlag()) {
    G4int gridRZnoBins = fOutputPrinter->GetGridRZnoBins();
    fGridRZDoseRun = new std::vector<G4double>;
    fGridRZDoseRun->reserve(gridRZnoBins);
    for (G4int ii = 0; ii < gridRZnoBins; ii++) {
      fGridRZDoseRun->push_back(0.0);
    }
    fGridRZDoseRun2 = new std::vector<G4double>;
    fGridRZDoseRun2->reserve(gridRZnoBins);
    for (G4int ii = 0; ii < gridRZnoBins; ii++) {
      fGridRZDoseRun2->push_back(0.0);
    }
    G4int gridZnoBins = fOutputPrinter->GetGridZnoBins();
    fDepthDoseRun = new std::vector<G4double>;
    fDepthDoseRun->reserve(gridZnoBins);
    for (G4int ii = 0; ii < gridZnoBins; ii++) {
      fDepthDoseRun->push_back(0.0);
    }
    fDepthDoseRun2 = new std::vector<G4double>;
    fDepthDoseRun2->reserve(gridZnoBins);
    for (G4int ii = 0; ii < gridZnoBins; ii++) {
      fDepthDoseRun2->push_back(0.0);
    }
  }

  // spc-enddet-## info
  if (fOutputPrinter->GetSpcEnddetFlag()) {
    fSpcEnddetRun = new std::vector< std::vector<G4double>* >;
    fSpcEnddetRun2 = new std::vector< std::vector<G4double>* >;
    fSpcEnddetBinDEVec = new std::vector<G4double>;

    size_t nEnddet = fOutputPrinter->GetNoSpcEnddet();

    for (size_t ii = 0; ii < nEnddet; ii++) {
      G4int nBins = fOutputPrinter->GetSpcEnddetNoBins(ii);
      std::vector<G4double>* aSpcEnddetRun = new std::vector<G4double>;
      std::vector<G4double>* aSpcEnddetRun2 = new std::vector<G4double>;
      aSpcEnddetRun->reserve(nBins);
      aSpcEnddetRun2->reserve(nBins);
      for (G4int jj = 0; jj < nBins; jj++) {
	aSpcEnddetRun->push_back(0.0);
	aSpcEnddetRun2->push_back(0.0);
      }
      fSpcEnddetRun->push_back( aSpcEnddetRun );
      fSpcEnddetRun2->push_back( aSpcEnddetRun2 );

      // Once vectors are created, we calculate and store binDE
      G4double emax = fOutputPrinter->GetSpcEnddetMax(ii);
      G4double emin = fOutputPrinter->GetSpcEnddetMin(ii);
      fSpcEnddetBinDEVec->push_back( (emax-emin)/nBins );
    }
  }

  // emerging particles info
  // -----------------------
  if (fOutputPrinter->GetEmergingFlag()) {
    // polar-angle info
    G4int nBinsPolarAngle = fOutputPrinter->GetPolarAngleNoBins();
    fPolarAngleElectRun = new std::vector<G4double>;
    fPolarAngleGammaRun = new std::vector<G4double>;
    fPolarAnglePositRun = new std::vector<G4double>;
    fPolarAngleElectRun->reserve(nBinsPolarAngle);
    fPolarAngleGammaRun->reserve(nBinsPolarAngle);
    fPolarAnglePositRun->reserve(nBinsPolarAngle);
    for (G4int ii = 0; ii < nBinsPolarAngle; ii++) {
      fPolarAngleElectRun->push_back(0.0);
      fPolarAngleGammaRun->push_back(0.0);
      fPolarAnglePositRun->push_back(0.0);
    }
    fPolarAngleElectRun2 = new std::vector<G4double>;
    fPolarAngleGammaRun2 = new std::vector<G4double>;
    fPolarAnglePositRun2 = new std::vector<G4double>;
    fPolarAngleElectRun2->reserve(nBinsPolarAngle);
    fPolarAngleGammaRun2->reserve(nBinsPolarAngle);
    fPolarAnglePositRun2->reserve(nBinsPolarAngle);
    for (G4int ii = 0; ii < nBinsPolarAngle; ii++) {
      fPolarAngleElectRun2->push_back(0.0);
      fPolarAngleGammaRun2->push_back(0.0);
      fPolarAnglePositRun2->push_back(0.0);
    }
    G4int nBinsEmergEne = fOutputPrinter->GetEmergEneNoBins();
    // energy-up & energy-down info
    fEmergEneUpElectRun = new std::vector<G4double>;
    fEmergEneUpGammaRun = new std::vector<G4double>;
    fEmergEneUpPositRun = new std::vector<G4double>;
    fEmergEneUpElectRun->reserve(nBinsEmergEne);
    fEmergEneUpGammaRun->reserve(nBinsEmergEne);
    fEmergEneUpPositRun->reserve(nBinsEmergEne);
    fEmergEneDownElectRun = new std::vector<G4double>;
    fEmergEneDownGammaRun = new std::vector<G4double>;
    fEmergEneDownPositRun = new std::vector<G4double>;
    fEmergEneDownElectRun->reserve(nBinsEmergEne);
    fEmergEneDownGammaRun->reserve(nBinsEmergEne);
    fEmergEneDownPositRun->reserve(nBinsEmergEne);
    for (G4int ii = 0; ii < nBinsEmergEne; ii++) {
      fEmergEneUpElectRun->push_back(0.0);
      fEmergEneUpGammaRun->push_back(0.0);
      fEmergEneUpPositRun->push_back(0.0);
      fEmergEneDownElectRun->push_back(0.0);
      fEmergEneDownGammaRun->push_back(0.0);
      fEmergEneDownPositRun->push_back(0.0);
    }
    fEmergEneUpElectRun2 = new std::vector<G4double>;
    fEmergEneUpGammaRun2 = new std::vector<G4double>;
    fEmergEneUpPositRun2 = new std::vector<G4double>;
    fEmergEneUpElectRun2->reserve(nBinsEmergEne);
    fEmergEneUpGammaRun2->reserve(nBinsEmergEne);
    fEmergEneUpPositRun2->reserve(nBinsEmergEne);
    fEmergEneDownElectRun2 = new std::vector<G4double>;
    fEmergEneDownGammaRun2 = new std::vector<G4double>;
    fEmergEneDownPositRun2 = new std::vector<G4double>;
    fEmergEneDownElectRun2->reserve(nBinsEmergEne);
    fEmergEneDownGammaRun2->reserve(nBinsEmergEne);
    fEmergEneDownPositRun2->reserve(nBinsEmergEne);
    for (G4int ii = 0; ii < nBinsEmergEne; ii++) {
      fEmergEneUpElectRun2->push_back(0.0);
      fEmergEneUpGammaRun2->push_back(0.0);
      fEmergEneUpPositRun2->push_back(0.0);
      fEmergEneDownElectRun2->push_back(0.0);
      fEmergEneDownGammaRun2->push_back(0.0);
      fEmergEneDownPositRun2->push_back(0.0);
    }
  }
}



//==============================================================================
Run::Run()
:G4Run()
{
  G4cout << "Run void constructor not meant to be used!" << G4endl;
}


//==============================================================================
// Destructor
//    clear all data members.
Run::~Run()
{
  G4cout << "Destroying Run object" << G4endl;
  if (fGridRZDoseRun) {
    fGridRZDoseRun->clear();
    delete fGridRZDoseRun;
  }
  if (fGridRZDoseRun2) {
    fGridRZDoseRun2->clear();
    delete fGridRZDoseRun2;
  }

  if (fDepthDoseRun) {
    fDepthDoseRun->clear();
    delete fDepthDoseRun;
  }
  if (fDepthDoseRun2) {
    fDepthDoseRun2->clear();
    delete fDepthDoseRun2;
  }

  if (fSpcEnddetRun) {
    size_t nDet = fSpcEnddetRun->size();
    for (size_t ii = 0; ii < nDet; ii++) {
      (*fSpcEnddetRun)[ii]->clear();
      delete (*fSpcEnddetRun)[ii];
    }
    fSpcEnddetRun->clear();
    delete fSpcEnddetRun;
  }
  if (fSpcEnddetRun2) {
    size_t nDet = fSpcEnddetRun2->size();
    for (size_t ii = 0; ii < nDet; ii++) {
      (*fSpcEnddetRun2)[ii]->clear();
      delete (*fSpcEnddetRun2)[ii];
    }
    fSpcEnddetRun2->clear();
    delete fSpcEnddetRun2;
  }
  if (fSpcEnddetBinDEVec) {
    fSpcEnddetBinDEVec->clear();
    delete fSpcEnddetBinDEVec;
  }

  if (fPolarAngleElectRun) {
    fPolarAngleElectRun->clear();
    delete fPolarAngleElectRun;
  }
  if (fPolarAngleElectRun2) {
    fPolarAngleElectRun2->clear();
    delete fPolarAngleElectRun2;
  }
  if (fPolarAngleGammaRun) {
    fPolarAngleGammaRun->clear();
    delete fPolarAngleGammaRun;
  }
  if (fPolarAngleGammaRun2) {
    fPolarAngleGammaRun2->clear();
    delete fPolarAngleGammaRun2;
  }
  if (fPolarAnglePositRun) {
    fPolarAnglePositRun->clear();
    delete fPolarAnglePositRun;
  }
  if (fPolarAnglePositRun2) {
    fPolarAnglePositRun2->clear();
    delete fPolarAnglePositRun2;
  }
  if (fEmergEneUpElectRun) {
    fEmergEneUpElectRun->clear();
    delete fEmergEneUpElectRun;
  }
  if (fEmergEneUpElectRun2) {
    fEmergEneUpElectRun2->clear();
    delete fEmergEneUpElectRun2;
  }
  if (fEmergEneUpGammaRun) {
    fEmergEneUpGammaRun->clear();
    delete fEmergEneUpGammaRun;
  }
  if (fEmergEneUpGammaRun2) {
    fEmergEneUpGammaRun2->clear();
    delete fEmergEneUpGammaRun2;
  }
  if (fEmergEneUpPositRun) {
    fEmergEneUpPositRun->clear();
    delete fEmergEneUpPositRun;
  }
  if (fEmergEneUpPositRun2) {
    fEmergEneUpPositRun2->clear();
    delete fEmergEneUpPositRun2;
  }
  if (fEmergEneDownElectRun) {
    fEmergEneDownElectRun->clear();
    delete fEmergEneDownElectRun;
  }
  if (fEmergEneDownElectRun2) {
    fEmergEneDownElectRun2->clear();
    delete fEmergEneDownElectRun2;
  }
  if (fEmergEneDownGammaRun) {
    fEmergEneDownGammaRun->clear();
    delete fEmergEneDownGammaRun;
  }
  if (fEmergEneDownGammaRun2) {
    fEmergEneDownGammaRun2->clear();
    delete fEmergEneDownGammaRun2;
  }
  if (fEmergEneDownPositRun) {
    fEmergEneDownPositRun->clear();
    delete fEmergEneDownPositRun;
  }
  if (fEmergEneDownPositRun2) {
    fEmergEneDownPositRun2->clear();
    delete fEmergEneDownPositRun2;
  }

}



//==============================================================================
//  RecordEvent() is a method called at the end of each event.
//  For scoring purpose, the resultant quantity in a event,
//  is accumulated during a Run.
void Run::RecordEvent(const G4Event* aEvent)
{
  G4Run::RecordEvent(aEvent);   // Mandatory to increment 'numberOfEvent'
  //  G4cout << "numberOfEvent = " << numberOfEvent << G4endl;

  // Dump accum info to vectors form evt-persistent to run-persistent objects
  // ------------------------------------------------------------------------

  // 2d-dose-map & depth-dose info
  if (fOutputPrinter->GetRZDoseMapFlag()) {
    std::vector<G4double>* gridRZDoseEvt = fOutputPrinter->GetGridRZDoseEvt();
    G4int gridRZnoBins = fOutputPrinter->GetGridRZnoBins();
    for (G4int ii = 0; ii < gridRZnoBins; ii++) {
      G4double evtDose = (*gridRZDoseEvt)[ii];
      (*fGridRZDoseRun)[ii] += evtDose;
      (*fGridRZDoseRun2)[ii] += (evtDose*evtDose);
    }
    std::vector<G4double>* depthDoseEvt = fOutputPrinter->GetDepthDoseEvt();
    G4int gridZnoBins = fOutputPrinter->GetGridZnoBins();
    for (G4int ii = 0; ii < gridZnoBins; ii++) {
      G4double evtDose = (*depthDoseEvt)[ii];
      (*fDepthDoseRun)[ii] += evtDose;
      (*fDepthDoseRun2)[ii] += (evtDose*evtDose);
    }
  }

  // spc-enddet-## data
  if (fOutputPrinter->GetSpcEnddetFlag()) {
    auto nSpcEnddet = fSpcEnddetRun->size();
    for (size_t ii = 0; ii < nSpcEnddet; ii++) {
      G4double evtEdep = fOutputPrinter->GetSpcEnddetEvt(ii);
      G4double edMin = fOutputPrinter->GetSpcEnddetMin(ii);
      G4double edMax = fOutputPrinter->GetSpcEnddetMax(ii);
      //DEBUG!!
      // G4cout << "\tevtEdep[eV] = " << evtEdep/eV << G4endl;
      // G4cout << "\tedMin[eV] = " << edMin/eV;
      // G4cout << "\tedMax[eV] = " << edMax/eV << G4endl;
    
      if (evtEdep > edMin && evtEdep < edMax) {
	G4int edBin =
	  static_cast<G4int>( (evtEdep-edMin)/((*fSpcEnddetBinDEVec)[ii]) );
	(*(*fSpcEnddetRun)[ii])[edBin] += 1.;
	(*(*fSpcEnddetRun2)[ii])[edBin] += 1.;
	//DEBUG!!
	// G4cout <<"edBin= " << edBin << G4endl;
	// G4cout << "(*(*fSpcEnddetRun)[ii])[edBin] = "
	//           << (*(*fSpcEnddetRun)[ii])[edBin]
	// 	     << G4endl;
      }
    }
  }
  
  // emerging particles
  // ----------------
  if (fOutputPrinter->GetEmergingFlag()) {
    // polar-angle data
    G4int nBinsT = fOutputPrinter->GetPolarAngleNoBins();
    std::vector<G4double>* evtPolarAngleE
      = fOutputPrinter->GetPolarAngleElectEvt();
    std::vector<G4double>* evtPolarAngleG
      = fOutputPrinter->GetPolarAngleGammaEvt();
    std::vector<G4double>* evtPolarAngleP
      = fOutputPrinter->GetPolarAnglePositEvt();
    for (G4int ii = 0; ii < nBinsT; ii++) {
      G4double wE = (*evtPolarAngleE)[ii];
      G4double wG = (*evtPolarAngleG)[ii];
      G4double wP = (*evtPolarAngleP)[ii];
      (*fPolarAngleElectRun)[ii] += wE;
      (*fPolarAngleElectRun2)[ii] += (wE*wE);
      (*fPolarAngleGammaRun)[ii] += wG;
      (*fPolarAngleGammaRun2)[ii] += (wG*wG);
      (*fPolarAnglePositRun)[ii] += wP;
      (*fPolarAnglePositRun2)[ii] += (wP*wP);
    }

    // energy-up & energy-down data
    G4int nBinsE = fOutputPrinter->GetEmergEneNoBins();
    std::vector<G4double>* evtEmergEneUpE
      = fOutputPrinter->GetEmergEneUpElectEvt();
    std::vector<G4double>* evtEmergEneUpG
      = fOutputPrinter->GetEmergEneUpGammaEvt();
    std::vector<G4double>* evtEmergEneUpP
      = fOutputPrinter->GetEmergEneUpPositEvt();
    std::vector<G4double>* evtEmergEneDownE
      = fOutputPrinter->GetEmergEneDownElectEvt();
    std::vector<G4double>* evtEmergEneDownG
      = fOutputPrinter->GetEmergEneDownGammaEvt();
    std::vector<G4double>* evtEmergEneDownP
      = fOutputPrinter->GetEmergEneDownPositEvt();
    for (G4int ii = 0; ii < nBinsE; ii++) {
      G4double wuE = (*evtEmergEneUpE)[ii];
      G4double wuG = (*evtEmergEneUpG)[ii];
      G4double wuP = (*evtEmergEneUpP)[ii];
      G4double wdE = (*evtEmergEneDownE)[ii];
      G4double wdG = (*evtEmergEneDownG)[ii];
      G4double wdP = (*evtEmergEneDownP)[ii];
      (*fEmergEneUpElectRun)[ii] += wuE;
      (*fEmergEneUpElectRun2)[ii] += (wuE*wuE);
      (*fEmergEneUpGammaRun)[ii] += wuG;
      (*fEmergEneUpGammaRun2)[ii] += (wuG*wuG);
      (*fEmergEneUpPositRun)[ii] += wuP;
      (*fEmergEneUpPositRun2)[ii] += (wuP*wuP);
      (*fEmergEneDownElectRun)[ii] += wdE;
      (*fEmergEneDownElectRun2)[ii] += (wdE*wdE);
      (*fEmergEneDownGammaRun)[ii] += wdG;
      (*fEmergEneDownGammaRun2)[ii] += (wdG*wdG);
      (*fEmergEneDownPositRun)[ii] += wdP;
      (*fEmergEneDownPositRun2)[ii] += (wdP*wdP);
    }

  }
}




//============================================================================
// Method to be overwritten by the user for merging local Run object to the
// global Run object.
void Run::Merge(const G4Run* aRun) 
{
  G4cout << "Run::Merge() started" << G4endl;

  const Run* localRun = static_cast<const Run*>(aRun);

  // Merge 2d-dose-map & depth-dose info
  if (fOutputPrinter->GetRZDoseMapFlag()) {
    G4int gridRZnoBins = fOutputPrinter->GetGridRZnoBins();
    for (G4int ii = 0; ii < gridRZnoBins; ii++) {
      (*fGridRZDoseRun)[ii] += (*(localRun->fGridRZDoseRun))[ii];
      (*fGridRZDoseRun2)[ii] += (*(localRun->fGridRZDoseRun2))[ii];
    }
    G4int gridZnoBins = fOutputPrinter->GetGridZnoBins();
    for (G4int ii = 0; ii < gridZnoBins; ii++) {
      (*fDepthDoseRun)[ii] += (*(localRun->fDepthDoseRun))[ii];
      (*fDepthDoseRun2)[ii] += (*(localRun->fDepthDoseRun2))[ii];
    }
  }

  // Merge spc-enddet-## info
  if (fOutputPrinter->GetSpcEnddetFlag()) {
    size_t nDet = fOutputPrinter->GetNoSpcEnddet();
    // G4cout << "nDet= " << nDet << G4endl;
    // G4cout << "fSpcEnddetRun->size()= " << fSpcEnddetRun->size() << G4endl;
    // G4cout << "localRun->fSpcEnddetRun->size()= "
    // 	   << localRun->fSpcEnddetRun->size() << G4endl;
    for (size_t jj = 0; jj < nDet; jj++) {
      size_t edNbin = fOutputPrinter->GetSpcEnddetNoBins(jj);
      // G4cout << "edNbin spc-enddet # " << jj+1 << " = " << edNbin << G4endl;
      // G4cout << "(*fSpcEnddetRun)[jj]->size()= "
      // 	     << (*fSpcEnddetRun)[jj]->size() << G4endl;
      // G4cout << "(*(localRun->fSpcEnddetRun))[jj]->size()= "
      // 	     << (*(localRun->fSpcEnddetRun))[jj]->size() << G4endl;
      for (size_t ii = 0; ii < edNbin; ii++) {
	(*(*fSpcEnddetRun)[jj])[ii]
	  += (*(*(localRun->fSpcEnddetRun))[jj])[ii];
	(*(*fSpcEnddetRun2)[jj])[ii]
	  += (*(*(localRun->fSpcEnddetRun2))[jj])[ii];
      }
    }
  }
	  
  // Merge emerging particles info
  if (fOutputPrinter->GetEmergingFlag()) {
    // Merge polar-angle info
    G4int nBinsT = fOutputPrinter->GetPolarAngleNoBins();
    for (G4int ii = 0; ii < nBinsT; ii++) {
      (*fPolarAngleElectRun)[ii] += (*(localRun->fPolarAngleElectRun))[ii];
      (*fPolarAngleElectRun2)[ii] += (*(localRun->fPolarAngleElectRun2))[ii];
      (*fPolarAngleGammaRun)[ii] += (*(localRun->fPolarAngleGammaRun))[ii];
      (*fPolarAngleGammaRun2)[ii] += (*(localRun->fPolarAngleGammaRun2))[ii];
      (*fPolarAnglePositRun)[ii] += (*(localRun->fPolarAnglePositRun))[ii];
      (*fPolarAnglePositRun2)[ii] += (*(localRun->fPolarAnglePositRun2))[ii];
    }

    G4int nBinsE = fOutputPrinter->GetEmergEneNoBins();
    for (G4int ii = 0; ii < nBinsE; ii++) {
      // Merge energy-up info
      (*fEmergEneUpElectRun)[ii] += (*(localRun->fEmergEneUpElectRun))[ii];
      (*fEmergEneUpElectRun2)[ii] += (*(localRun->fEmergEneUpElectRun2))[ii];
      (*fEmergEneUpGammaRun)[ii] += (*(localRun->fEmergEneUpGammaRun))[ii];
      (*fEmergEneUpGammaRun2)[ii] += (*(localRun->fEmergEneUpGammaRun2))[ii];
      (*fEmergEneUpPositRun)[ii] += (*(localRun->fEmergEneUpPositRun))[ii];
      (*fEmergEneUpPositRun2)[ii] += (*(localRun->fEmergEneUpPositRun2))[ii];
      // Merge energy-down info
      (*fEmergEneDownElectRun)[ii]+=(*(localRun->fEmergEneDownElectRun))[ii];
      (*fEmergEneDownElectRun2)[ii]+=(*(localRun->fEmergEneDownElectRun2))[ii];
      (*fEmergEneDownGammaRun)[ii]+=(*(localRun->fEmergEneDownGammaRun))[ii];
      (*fEmergEneDownGammaRun2)[ii]+=(*(localRun->fEmergEneDownGammaRun2))[ii];
      (*fEmergEneDownPositRun)[ii]+=(*(localRun->fEmergEneDownPositRun))[ii];
      (*fEmergEneDownPositRun2)[ii]+=(*(localRun->fEmergEneDownPositRun2))[ii];
    }
  }

  G4Run::Merge(aRun);
}
