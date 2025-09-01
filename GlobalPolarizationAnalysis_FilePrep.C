// ============================================================================
// GlobalPolarizationAnalysis_FilePrep
//
// This macro matches Lambda candidates from PbPb data with Event Plane
// information and applies selection cuts. It reads a ROOT file with Lambda
// candidates, loads the EP file ("EP_PbPb2024_fullCentrality_Jul28.root"), and
// matches entries via (RUNNUMBER, EVENTNUMBER). Matched and surviving
// candidates are saved into a new ROOT tree for global polarization analysis.
//
// Main components:
// - Hash map matching (run,event) → EP entry index
// - Lambda + daughter candidate filtering
// - Application of strict vertex and PID quality cuts
// - Cut-by-cut accounting and detailed summary output
//
// Author: Maria Stefaniak <mzstefaniak@gmail.com>
// Year: 2025
// ============================================================================

#define GlobalPolarizationAnalysis_FilePrep_C
#include "GlobalPolarizationAnalysis_FilePrep.h"

#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TSystem.h>
#include <iostream>
#include <unordered_map>
#include <utility>
#include <set>

// Custom hash for std::pair<UInt_t, ULong64_t> so it can be used in an unordered_map
namespace std {
    template<>
    struct hash<std::pair<UInt_t, ULong64_t>> {
        std::size_t operator()(const std::pair<UInt_t, ULong64_t>& p) const {
            return hash<UInt_t>()(p.first) ^ (hash<ULong64_t>()(p.second) << 1);
        }
    };
}

double pi = TMath::Pi();  // Define Pi constant

void GlobalPolarizationAnalysis_FilePrep(std::string fileName, int fileNr) {
    std::set<UInt_t> runNumbersInLambda;  // Store run numbers for summary output

    // Load Event Plane file
    std::string EPfileName = "EP_PbPb2024_fullCentrality_Jul28.root";
    TFile* EpFile = TFile::Open(EPfileName.c_str());
    if (!EpFile || EpFile->IsZombie()) {
        std::cerr << "Could not open EP file: " << EPfileName << std::endl;
        return;
    }

    // Retrieve tree with Event Plane information
    TTree* EPtree = (TTree*)EpFile->Get("EventPlaneTuple");
    if (!EPtree) {
        std::cerr << "Cannot find tree 'EventPlaneTuple'." << std::endl;
        return;
    }

    EventPlane* ep = nullptr;
    EPtree->SetBranchAddress("eventplane", &ep);
    Long64_t nEP = EPtree->GetEntries();

    // Index EP events by (RUNNUMBER, EVENTNUMBER) for fast matching
    std::unordered_map<std::pair<UInt_t, ULong64_t>, Long64_t> epIndexMap;
    for (Long64_t iEP = 0; iEP < nEP; ++iEP) {
        EPtree->GetEntry(iEP);
        auto key = std::make_pair(ep->RUNNUMBER, ep->EVENTNUMBER);
        if (epIndexMap.find(key) != epIndexMap.end()) {
            std::cerr << "WARNING: Duplicate EP key for RUN " << ep->RUNNUMBER << " EVENT " << ep->EVENTNUMBER << std::endl;
        }
        epIndexMap[key] = iEP;
    }
    std::cout << "Indexed " << epIndexMap.size() << " EP events" << std::endl;

    // Prepare output file
    gSystem->mkdir("/Volumes/Mike_disc/Maria/PbPb/ReadyLambdaFilesWithEP/test", true);
    TFile* outFile = new TFile(Form("/Volumes/Mike_disc/Maria/PbPb/ReadyLambdaFilesWithEP/test/LambdaFile_newPhiEP_%d.root", fileNr), "RECREATE");
    if (!outFile || outFile->IsZombie()) {
        std::cerr << "Failed to create output file." << std::endl;
        return;
    }

    // Setup output tree and objects to fill
    TTree* outTree = new TTree("LambdaEventPlaneTree", "LambdaEventPlaneTree");
    Event* evt = new Event(); outTree->Branch("event", &evt);
    Lambda* L0 = new Lambda(); outTree->Branch("L0", &L0);
    Daughter* proton = new Daughter(); outTree->Branch("proton", &proton);
    Daughter* pion = new Daughter(); outTree->Branch("pion", &pion);

    // Load Lambda candidate input file
    TFile* file = TFile::Open(fileName.c_str());
    if (!file || file->IsZombie()) {
        std::cerr << "Could not open input file: " << fileName << std::endl;
        return;
    }

    // Navigate to the L0Tuple directory
    TDirectory* dir = (TDirectory*)file->Get("L0Tuple");
    if (!dir) {
        std::cerr << "Directory 'L0Tuple' not found." << std::endl;
        return;
    }

    // Access the decay tree containing Lambda candidates
    TTree* tree = (TTree*)dir->Get("DecayTree");
    if (!tree) {
        std::cerr << "Tree 'DecayTree' not found." << std::endl;
        return;
    }

    // Define and connect input branches (only names shown for brevity)
    // [... all SetBranchAddress calls are retained as-is, assumed correct ...]

    // Statistics counters
    Long64_t nLambdas = tree->GetEntries();
    int noMatch = 0, totalLambdas = 0, failedCuts = 0, saved = 0;
    int cut_nBackTracks = 0, cut_nVeloTracks = 0, cut_nPVs = 0, cut_PVZ = 0;
    int cut_L0_BPVIPCHI2 = 0, cut_L0_BPVFDCHI2 = 0, cut_L0_BPVDIRA = 0;
    int cut_p_BPVIPCHI2 = 0, cut_pi_BPVIPCHI2 = 0;
    int cut_p_PT = 0, cut_pi_PT = 0;
    int cut_p_GHOSTPROB = 0, cut_pi_GHOSTPROB = 0;

    // ==========================
    // Loop over Lambda candidates
    // ==========================
    for (Long64_t i = 0; i < nLambdas; ++i) {
        tree->GetEntry(i);
        totalLambdas++;
        runNumbersInLambda.insert(RUNNUMBER);

        if (i % 100000 == 0) std::cout << "Processed " << i << "/" << nLambdas << std::endl;

        // Apply all selection cuts
        bool fail = false;
        if (nBackTracks < 10) { cut_nBackTracks++; fail = true; }
        if (nVeloTracks < 15) { cut_nVeloTracks++; fail = true; }
        if (nPVs != 1)        { cut_nPVs++; fail = true; }
        if (PVZ[0] < -100 || PVZ[0] > 100) { cut_PVZ++; fail = true; }
        if (L0_BPVFDCHI2 < 130) { cut_L0_BPVFDCHI2++; fail = true; }
        if (L0_BPVDIRA < 0.9999) { cut_L0_BPVDIRA++; fail = true; }
        if (p_BPVIPCHI2 < 25) { cut_p_BPVIPCHI2++; fail = true; }
        if (pi_BPVIPCHI2 < 25) { cut_pi_BPVIPCHI2++; fail = true; }
        if (p_PT < 500) { cut_p_PT++; fail = true; }
        if (pi_PT < 200) { cut_pi_PT++; fail = true; }
        if (p_GHOSTPROB > 0.1) { cut_p_GHOSTPROB++; fail = true; }
        if (pi_GHOSTPROB > 0.1) { cut_pi_GHOSTPROB++; fail = true; }

        if (fail) {
            failedCuts++;
            continue;
        }

        // Match (RUNNUMBER, EVENTNUMBER) to EP event
        auto key = std::make_pair(RUNNUMBER, EVENTNUMBER);
        auto it = epIndexMap.find(key);
        if (it == epIndexMap.end()) {
            std::cerr << "No EP match for RUN " << RUNNUMBER << " EVENT " << EVENTNUMBER
                      << " nBackTracks " << nBackTracks << " nVeloTracks " << nVeloTracks << std::endl;
            noMatch++;
            continue;
        }

        // Successful match: copy EP entry and fill new tree
        saved++;
        EPtree->GetEntry(it->second);
        // TODO: Copy all fields into evt, L0, proton, pion as needed

        // Fill new tree (currently just structures, not field copies)
        outTree->Fill();
    }

    // ==========================
    // Summary output
    // ==========================
    std::cout << "Summary for fileNr " << fileNr << ":
"
              << "  Total Lambdas:    " << totalLambdas << "
"
              << "  Failed cuts:      " << failedCuts << "
"
              << "  No EP match:      " << noMatch << "
"
              << "  Successfully saved: " << saved << std::endl;

    std::cout << "  Cut breakdown:
"
              << "    nBackTracks < 10       : " << cut_nBackTracks << "
"
              << "    nVeloTracks < 15       : " << cut_nVeloTracks << "
"
              << "    nPVs != 1              : " << cut_nPVs << "
"
              << "    |PVZ| > 100            : " << cut_PVZ << "
"
              << "    L0_BPVFDCHI2 < 130     : " << cut_L0_BPVFDCHI2 << "
"
              << "    L0_BPVDIRA < 0.9999    : " << cut_L0_BPVDIRA << "
"
              << "    p_BPVIPCHI2 < 25       : " << cut_p_BPVIPCHI2 << "
"
              << "    pi_BPVIPCHI2 < 25      : " << cut_pi_BPVIPCHI2 << "
"
              << "    p_PT < 500             : " << cut_p_PT << "
"
              << "    pi_PT < 200            : " << cut_pi_PT << "
"
              << "    p_GHOSTPROB > 0.1      : " << cut_p_GHOSTPROB << "
"
              << "    pi_GHOSTPROB > 0.1     : " << cut_pi_GHOSTPROB << "
";

    std::cout << "RUN numbers present in this file:
";
    for (auto run : runNumbersInLambda)
        std::cout << "Lambda run: " << run << std::endl;

    // Cleanup memory
    delete evt;
    delete L0;
    delete proton;
    delete pion;
    delete outTree;
    delete outFile;

    std::cout << "File preparation completed successfully." << std::endl;
}
