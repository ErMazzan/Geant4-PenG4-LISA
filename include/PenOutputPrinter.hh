
#ifndef PenOutputPrinter_hh
#define PenOutputPrinter_hh 1

#include "globals.hh"
#include <fstream>
#include <vector>

class G4Step;
class G4Track;
class DetectorConstruction;
class PenOutputMessenger;
class Run;

class PenOutputPrinter
{
private:
  PenOutputPrinter();

public:
  PenOutputPrinter(const DetectorConstruction* );
  ~PenOutputPrinter();

public:
  void CreateFiles();  // Only Master
  void BookForMaster();    // Initialize DM for Master
  void BookForWorker();    // Initialize DM for Worker
  void PrintPreTrackInfo(const G4Track* );
  void PrintStepInfo(const G4Step* );
  void StoreStepInfo(const G4Step* );  // Only Worker
  void ResetEventInfo();  // Only Worker
  void DumpAll(const Run* );    // Only Master
  void CloseAllFiles();  // Only Master

private:
  void CreateRZDoseMapFile();
  void CreateSpcEnddetFile();
  void CreatePolarAngleFile();
  void CreateEmergEneFiles();
  void CreateTrkFiles();
  void DumpRZDoseMap(const Run* );
  void DumpDepthDoseMap(const Run* );
  void DumpSpcEnddet(const Run* );
  void DumpPolarAngle(const Run* );
  void DumpEmergEne(const Run* );

public: // Get or Set methods
  inline void SetRZDoseMapFlag(G4bool choice) { fRZDoseMapFlag = choice; }
  inline void SetGridZmin(G4double val)       { fGridZmin = val; }
  inline void SetGridZmax(G4double val)       { fGridZmax = (1.+1.e-8)*val; }
  inline void SetGridZnoBins(G4int val)       { fGridZnoBins = val; }
  inline void SetGridRmax(G4double val)       { fGridRmax = (1.+1.e-8)*val; }
  inline void SetGridRnoBins(G4int val)       { fGridRnoBins = val; }

  inline void SetSpcEnddetFlag(G4bool choice) { fSpcEnddetFlag = choice; }
  void AddSpcEnddetPVName(G4String name);
  inline void AddSpcEnddetMin(G4double val)
  { fSpcEnddetMinVec->push_back(val); }
  inline void AddSpcEnddetMax(G4double val)
  { fSpcEnddetMaxVec->push_back( (1.+1.e-8)*val ); }
  inline void AddSpcEnddetNoBins(G4int val)
  { fSpcEnddetNoBinsVec->push_back(val);}

  inline void SetEmergingFlag(G4bool choice)  { fEmergingFlag = choice; }
  inline void SetPolarAngleNoBins(G4int val)  { fPolarAngleNoBins = val; }
  inline void SetEmergEneMin(G4double val)    { fEmergEneMin = val; }
  inline void SetEmergEneMax(G4double val)    { fEmergEneMax = (1.+1.e-8)*val;}
  inline void SetEmergEneNoBins(G4int val)    { fEmergEneNoBins = val; }

  inline void SetTrkFlag(G4bool choice)       { fTrkFlag = choice; }


  inline G4bool GetRZDoseMapFlag()      { return fRZDoseMapFlag; }
  inline G4int GetGridRnoBins()         { return fGridRnoBins; }
  inline G4int GetGridZnoBins()         { return fGridZnoBins; }
  inline G4int GetGridRZnoBins()        { return (fGridRnoBins*fGridZnoBins); }
  inline std::vector<G4double>* GetGridRZDoseEvt()  { return fGridRZDoseEvt; }
  inline std::vector<G4double>* GetDepthDoseEvt()   { return fDepthDoseEvt; }

  inline G4bool GetSpcEnddetFlag()    { return fSpcEnddetFlag; }
  inline size_t GetNoSpcEnddet()      { return fNoSpcEnddet; }
  inline G4double GetSpcEnddetMin(size_t ii) { return (*fSpcEnddetMinVec)[ii]; }
  inline G4double GetSpcEnddetMax(size_t ii) { return (*fSpcEnddetMaxVec)[ii]; }
  inline G4int GetSpcEnddetNoBins(size_t ii)
  { return (*fSpcEnddetNoBinsVec)[ii]; }
  inline G4double GetSpcEnddetEvt(size_t ii) { return (*fSpcEnddetEvtVec)[ii]; }

  inline G4bool GetEmergingFlag()      { return fEmergingFlag; }
  inline G4int GetPolarAngleNoBins()   { return fPolarAngleNoBins; }
  inline std::vector<G4double>* GetPolarAngleElectEvt()
  { return fPolarAngleElectEvt; }
  inline std::vector<G4double>* GetPolarAngleGammaEvt()
  { return fPolarAngleGammaEvt; }
  inline std::vector<G4double>* GetPolarAnglePositEvt()
  { return fPolarAnglePositEvt; }
  inline G4double GetEmergEneMin()  { return fEmergEneMin; }
  inline G4double GetEmergEneMax()  { return fEmergEneMax; }
  inline G4int GetEmergEneNoBins()  { return fEmergEneNoBins; }
  inline std::vector<G4double>* GetEmergEneUpElectEvt()
  { return fEmergEneUpElectEvt; }
  inline std::vector<G4double>* GetEmergEneUpGammaEvt()
  { return fEmergEneUpGammaEvt; }
  inline std::vector<G4double>* GetEmergEneUpPositEvt()
  { return fEmergEneUpPositEvt; }
  inline std::vector<G4double>* GetEmergEneDownElectEvt()
  { return fEmergEneDownElectEvt; }
  inline std::vector<G4double>* GetEmergEneDownGammaEvt()
  { return fEmergEneDownGammaEvt; }
  inline std::vector<G4double>* GetEmergEneDownPositEvt()
  { return fEmergEneDownPositEvt; }


private:

  const DetectorConstruction* fDetConst;
  PenOutputMessenger*            fMessenger;

  // 2d-dose-map.dat, z-dose.dat & depth-dose.dat
  G4double fGridZmin;
  G4double fGridZmax;
  G4int    fGridZnoBins;
  G4double fGridRmax;
  G4int    fGridRnoBins;
  G4double fGridBinDZ;      // to save operations at step level
  G4double fGridBinDR;      // to save operations at step level
  std::vector<G4double>* fGridRZDoseEvt;
  std::vector<G4double>* fDepthDoseEvt;
  std::ofstream fRZDoseMapFile;
  std::ofstream fDoseZFile;
  std::ofstream fDepthDoseFile;
  G4bool        fRZDoseMapFlag;

  // spc-enddet-##.dat
  size_t                 fNoSpcEnddet;            // To save time at step level
  std::vector<G4String>* fSpcEnddetNameVec;
  std::vector<G4double>* fSpcEnddetEvtVec;
  std::vector<std::ofstream*>* fSpcEnddetFileVec;
  std::vector<G4double>* fSpcEnddetMinVec;
  std::vector<G4double>* fSpcEnddetMaxVec;
  std::vector<G4int>*    fSpcEnddetNoBinsVec;
  G4bool fSpcEnddetFlag;

  // emerging particles
  // ------------------
  G4bool fEmergingFlag;
  // polar-angle.dat
  G4int    fPolarAngleNoBins;
  G4double fPolarAngleBinDT;
  std::vector<G4double>* fPolarAngleElectEvt;
  std::vector<G4double>* fPolarAngleGammaEvt;
  std::vector<G4double>* fPolarAnglePositEvt;
  std::ofstream fPolarAngleFile;
  // energy-up.dat & energy-down.dat
  G4double fEmergEneMin;
  G4double fEmergEneMax;
  G4int    fEmergEneNoBins;
  G4double fEmergEneBinDE;  // Must be defined from variables above
  std::vector<G4double>* fEmergEneUpElectEvt;
  std::vector<G4double>* fEmergEneUpGammaEvt;
  std::vector<G4double>* fEmergEneUpPositEvt;
  std::vector<G4double>* fEmergEneDownElectEvt;
  std::vector<G4double>* fEmergEneDownGammaEvt;
  std::vector<G4double>* fEmergEneDownPositEvt;
  std::ofstream fEmergEneUpFile;
  std::ofstream fEmergEneDownFile;

  // tracking files
  // --------------
  std::ofstream fElectronTrkFile;
  std::ofstream fGammaTrkFile;
  std::ofstream fPositronTrkFile;
  G4bool        fTrkFlag;

};

#endif
