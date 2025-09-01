// ============================================================================
// GlobalPolarizationAnalysis_FilePrep.h
//
// Header file for Global Polarization Analysis (LHCb Pb+Pb 2024)
//
// Defines data structures used in the analysis: EventPlane, Event, Lambda,
// and Daughter classes, each corresponding to a set of variables stored in
// ROOT trees. These classes facilitate I/O and cut-based selection for
// polarization studies.
//
// Author: Maria Stefaniak
// Affiliation: The Ohio State University
// Year: 2025
// ============================================================================

#ifndef GlobalPolarizationAnalysis_FilePrep_h
#define GlobalPolarizationAnalysis_FilePrep_h

#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH2.h>
#include <iostream>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TList.h>
#include <TString.h>
#include <vector>
#include <string>
#include <TRegexp.h>

// ============================================================================
// Class storing event plane information matched by (RUNNUMBER, EVENTNUMBER)
// ============================================================================
class EventPlane {
public:
    ULong64_t EVENTNUMBER;   // Unique event ID
    UInt_t    RUNNUMBER;     // Run ID
    Double_t  Psi1Full;      // Full Psi1 event plane angle
    Double_t  Psi2Full;      // Full Psi2 event plane angle
    Double_t  PsiBack[2];    // Psi1 and Psi2 from backward side
    Double_t  PsiFor[2];     // Psi1 and Psi2 from forward side
    Double_t  r1;            // Resolution for Psi1
    Double_t  r2;            // Resolution for Psi2

    EventPlane() : EVENTNUMBER(0), RUNNUMBER(0) {}
};

// ============================================================================
// Class storing global event-level quantities for a Lambda event
// ============================================================================
class Event {
public:
    ULong64_t EVENTNUMBER;
    UInt_t    RUNNUMBER;
    Double_t  Psi1Full, Psi2Full;
    Double_t  Psi1back, Psi2back;
    Double_t  Psi1for, Psi2for;
    Double_t  r1, r2;
    Float_t   PVX, PVY, PVZ;               // Primary vertex position
    Int_t     nBackTracks, nVeloTracks;    // Event multiplicity info
    Int_t     nEcalClusters;               // Not used in cuts (but stored)

    Event() : EVENTNUMBER(0), RUNNUMBER(0) {}
};

// ============================================================================
// Lambda candidate reconstructed from proton and pion daughters
// ============================================================================
class Lambda {
public:
    Int_t     ID;
    Float_t   ETA, PHI;
    Double_t  MASS;
    Float_t   PT, PX, PY, PZ;
    Float_t   BPVIPCHI2;     // Impact parameter chi2 w.r.t. PV
    Float_t   BPVFDCHI2;     // Flight distance chi2 w.r.t. PV
    Float_t   B_PV_Z, B_PV_X, B_PV_Y; // PV position
    Float_t   BPVDIRA;       // Direction angle cosine

    Lambda() : ID(0) {}
};

// ============================================================================
// Daughter particle (proton or pion) from Lambda decay
// ============================================================================
class Daughter {
public:
    Int_t     ID;
    Float_t   ETA, PHI;
    Double_t  MASS;
    Float_t   PT, PX, PY, PZ;
    Double_t  BPVIPCHI2;
    Double_t  GHOSTPROB;

    Daughter() : ID(0) {}
};

/*
// OPTIONAL: Function to read ROOT files from directory (not used currently)
std::vector<std::string> GetRootFilesInDirectory(const std::string& dirPath, int file_nr) {
    std::vector<std::string> fileNames;

    TSystemDirectory dir("dir", dirPath.c_str());
    TList* files = dir.GetListOfFiles();
    if (!files) {
        std::cerr << "Could not open directory: " << dirPath << std::endl;
        return fileNames;
    }

    TIter next(files);
    TSystemFile* file;
    while ((file = (TSystemFile*)next())) {
        TString fname = file->GetName();
        if (!file->IsDirectory() && fname.EndsWith(".root")) {
            fileNames.push_back(dirPath + "/" + fname.Data());
        }
    }

    delete files;
    return fileNames;
}
*/

#endif // GlobalPolarizationAnalysis_FilePrep_h
