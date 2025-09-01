// ============================================================================
// This is the original ROOT macro used for calculating event-plane Q-vectors
// using VELO tracks from PbPb collisions collected by LHCb in 2024.
//
// The code reads ROOT files from Analysis Production, applies event-level
// and track-level quality cuts, computes Q-vectors (1st and 2nd harmonics)
// in backward and forward η regions (with and without η-weighting), and stores
// the results into a structured TTree for further analysis.
//
// Author: Maria Stefaniak <mzstefaniak@gmail.com>
// Year: 2025
//
// Main structural components:
//
// - Header includes and Event class definition
// - Function: GetFilteredRootFiles()
//     - Scans input directory and selects ROOT files matching AP filename pattern
// - Function: EventPlaneAnalysis()
//     - Loads and loops over VELO AP files
//     - Applies event- and track-level cuts
//     - Calculates Q-vectors in 3 forward η bins and 1 backward bin
//     - Stores Q-vectors (with/without η weight), multiplicities, and event info
//       into a TTree saved in a new ROOT file
// ============================================================================

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
#include <TMath.h>

// ==========================
// Event class definition
// ==========================
// Holds event-level information and Q-vectors that will be saved into the output tree.
class Event {
    public:
        // General event information
        ULong64_t       outGPSTIME;
        ULong64_t       outEVENTNUMBER;
        Float_t         outPVX;   
        Float_t         outPVY;  
        Float_t         outPVZ;  
        UInt_t          outRUNNUMBER;

        // Global multiplicities
        Int_t           outnBackTracks;
        Int_t           outnVeloClusters;
        Int_t           outnVeloTracks;
        Int_t           outnEcalClusters;
        Int_t           outECalETot;
        Int_t           outnLongTracks;
        Int_t           outnVPClusters;

        // Q-vectors (harmonic orders 1 and 2)
        Double_t        outQx_back[2];
        Double_t        outQy_back[2];
        Double_t        outQx_for[2][4];   // 4 eta bins
        Double_t        outQy_for[2][4];

        // Q-vectors with eta weight
        Double_t        outQx_back_wEta[2];
        Double_t        outQy_back_wEta[2];
        Double_t        outQx_for_wEta[2][4];
        Double_t        outQy_for_wEta[2][4];

        // Multiplicity of forward/backward tracks in each eta bin
        Int_t           out_Qmulti[4];
    
        Event() : outGPSTIME(0), outEVENTNUMBER(0), outRUNNUMBER(0){}
};

// ==========================
// GetFilteredRootFiles
// ==========================
// Scans a directory and selects ROOT files matching a specific regex pattern.
// Only files of the form: 00274156_00000###_1.tuple_pbpb2024.root are accepted.
std::vector<std::string> GetFilteredRootFiles(const std::string& dirPath) {
    std::vector<std::string> fileNames;

    TSystemDirectory dir("dir", dirPath.c_str());
    TList* files = dir.GetListOfFiles();
    if (!files) {
        std::cerr << "Could not open or read directory: " << dirPath << std::endl;
        return fileNames;
    }

    TIter next(files);
    TSystemFile* file;

    // Regex pattern for PbPb 2024 AP file naming scheme
    TRegexp pattern("^00274156_00000[0-4][0-9][0-9]_1\\.tuple_pbpb2024\\.root$");
    while ((file = (TSystemFile*)next())) {
        TString fname = file->GetName();
        if (!file->IsDirectory() && fname.EndsWith(".root") && fname.Contains(pattern)) {
            std::string fullpath = dirPath + "/" + fname.Data();
            std::cout << "Adding file: " << fullpath << std::endl;
            fileNames.push_back(fullpath);
        }
    }

    delete files;
    return fileNames;
}

// ==========================
// EventPlaneAnalysis
// ==========================
// Main function that:
// 1. Reads AP files
// 2. Applies event and track cuts
// 3. Calculates Q-vectors for different eta ranges
// 4. Stores event information into an output ROOT tree
void EventPlaneAnalysis(){

    const int maxTracks = 10000;  // max number of tracks per event
    std::string eosDir = "/eos/lhcb/grid/prod/lhcb/anaprod/lhcb/LHCb/Lead24/TUPLE_PBPB2024.ROOT/00274156/0000";

    // Collect input ROOT files from EOS
    std::vector<std::string> fileNames = GetFilteredRootFiles(eosDir);

    // Output ROOT file + tree
    TFile* outFile = new TFile("centTests/weights_event_plane_pbpb_localtest.root", "RECREATE");
    TTree* outTree = new TTree("EventPlaneTuple", "Event Plane");

    // Event object linked to tree
    Event* evt = new Event();
    outTree->Branch("event", &evt);

    // ==========================
    // Loop over input files
    // ==========================
    for (const auto& fileName : fileNames) { 

        TFile *file = TFile::Open(fileName.c_str());
        if (!file || file->IsZombie()) {
            std::cerr << "Could not open file: " << fileName << std::endl;
            continue;
        }

        // Navigate to EventTuplePV directory
        TDirectory* dir = (TDirectory*)file->Get("EventTuplePV");
        if (!dir) {
            std::cerr << "Directory 'EventTuplePV' not found in file: " << fileName << std::endl;
            file->Close();
            continue;
        }

        // Retrieve tree from directory
        TTree* tree = (TTree*)dir->Get("EventTuplePV");
        if (!tree) {
            std::cerr << "Tree not found in directory 'EventTuplePV' in file: " << fileName << std::endl;
            file->Close();
            continue;
        }
        std::cerr << "file: " << fileName << " opened" << std::endl;

        // ==========================
        // Declare branch variables
        // ==========================
        ULong64_t       GPSTIME;
        ULong64_t       EVENTNUMBER;
        Float_t         PVX[maxTracks];   
        Float_t         PVY[maxTracks];   
        Float_t         PVZ[maxTracks];   
        UInt_t          RUNNUMBER;
        Float_t         VELOTRACK_BIPCHI2[maxTracks];   
        Float_t         VELOTRACK_ETA[maxTracks];  
        Float_t         VELOTRACK_ISBACKWARD[maxTracks];  
        Float_t         VELOTRACK_NVPHITS[maxTracks];   
        Float_t         VELOTRACK_PHI[maxTracks];   
        Float_t         VELOTRACK_PX[maxTracks];   
        Float_t         VELOTRACK_PY[maxTracks];  
        Float_t         VELOTRACK_PZ[maxTracks];  
        Float_t         VELOTRACK_RHO[maxTracks];   
        Int_t           nBackTracks;
        Int_t           nPVs;
        Int_t           nVeloClusters;
        Int_t           nVeloTracks;
        Int_t           nEcalClusters;
        Int_t           ECalETot;
        Int_t           nLongTracks;
        Int_t           nVPClusters;

        // ==========================
        // Connect branches
        // ==========================
        tree->SetBranchAddress("GPSTIME", &GPSTIME);
        tree->SetBranchAddress("EVENTNUMBER", &EVENTNUMBER);
        tree->SetBranchAddress("PVX", &PVX);
        tree->SetBranchAddress("PVY", &PVY);
        tree->SetBranchAddress("PVZ", &PVZ);
        tree->SetBranchAddress("RUNNUMBER", &RUNNUMBER);
        tree->SetBranchAddress("VELOTRACK_BIPCHI2", &VELOTRACK_BIPCHI2);
        tree->SetBranchAddress("VELOTRACK_ETA", &VELOTRACK_ETA);
        tree->SetBranchAddress("VELOTRACK_ISBACKWARD", &VELOTRACK_ISBACKWARD);
        tree->SetBranchAddress("VELOTRACK_NVPHITS", &VELOTRACK_NVPHITS);
        tree->SetBranchAddress("VELOTRACK_PHI", &VELOTRACK_PHI);
        tree->SetBranchAddress("VELOTRACK_PX", &VELOTRACK_PX);
        tree->SetBranchAddress("VELOTRACK_PY", &VELOTRACK_PY);
        tree->SetBranchAddress("VELOTRACK_PZ", &VELOTRACK_PZ);
        tree->SetBranchAddress("VELOTRACK_RHO", &VELOTRACK_RHO);
        tree->SetBranchAddress("nBackTracks", &nBackTracks);
        tree->SetBranchAddress("nPVs", &nPVs);
        tree->SetBranchAddress("nVeloClusters", &nVeloClusters);
        tree->SetBranchAddress("nVeloTracks", &nVeloTracks);
        tree->SetBranchAddress("nEcalClusters", &nEcalClusters);
        tree->SetBranchAddress("ECalETot", &ECalETot);
        tree->SetBranchAddress("nLongTracks", &nLongTracks);
        tree->SetBranchAddress("nVPClusters", &nVPClusters);

        // ==========================
        // Loop over events in tree
        // ==========================
        Long64_t nEvents = tree->GetEntries();
        for (Long64_t iEvent = 0; iEvent < nEvents; ++iEvent) {
            tree->GetEntry(iEvent);

            // ------------------
            // Event selection
            // ------------------
            if(nPVs!=1) continue;                // Only single collision events
            if(nBackTracks < 10) continue;       // Ensure PbPb collision, not background
            if(PVZ[0]<-100 || PVZ[0]>100) continue;  // Vertex within fiducial region
            if(nVeloTracks<15) continue;         // Require enough tracks

            // Track multiplicities
            int nPrimaryForTracks = 0; 
            int nPrimaryBackTracks = 0; 
            int nPrimaryForTracks_1 = 0; 
            int nPrimaryForTracks_2 = 0; 
            int nPrimaryForTracks_3 = 0;

            // Initialize Q-vectors (two harmonics: n=1,2)
            double Qx_back[2], Qy_back[2], Qx_for[2][4], Qy_for[2][4];
            double Qx_back_wEta[2], Qy_back_wEta[2], Qx_for_wEta[2][4], Qy_for_wEta[2][4];
            for(int i = 0; i < 2; i++){
                Qx_back[i] = Qy_back[i] = 0;
                Qx_back_wEta[i] = Qy_back_wEta[i] = 0;
                for(int j = 0; j < 4; j++){
                    Qx_for[i][j] = Qy_for[i][j] = 0;
                    Qx_for_wEta[i][j] = Qy_for_wEta[i][j] = 0;
                }
            }

            // ------------------
            // Loop over tracks
            // ------------------
            for(int iTrack = 0; iTrack < nVeloTracks; iTrack++){
                
                // Basic track quality cut
                if(VELOTRACK_BIPCHI2[iTrack]>1.5) continue;

                // Handle backward tracks: flip eta and phi
                if(VELOTRACK_ISBACKWARD[iTrack]==1) {
                    VELOTRACK_ETA[iTrack] *= -1;
                    VELOTRACK_PHI[iTrack] += TMath::Pi();
                }

                // Count forward/backward tracks
                if(VELOTRACK_ISBACKWARD[iTrack]!=1) nPrimaryForTracks++;
                if(VELOTRACK_ISBACKWARD[iTrack]==1) nPrimaryBackTracks++;

                // Weights
                double w1=1;              // default weight = 1
                double w2=1;              // for n=2 harmonic
                double wEta1=VELOTRACK_ETA[iTrack]; // eta-weighted option

                // Backward eta region
                if(VELOTRACK_ETA[iTrack] < -0.5){
                    // Harmonic n=1,2
                    Qx_back[0] += w1 * cos(1 * VELOTRACK_PHI[iTrack]);
                    Qx_back[1] += w2 * cos(2 * VELOTRACK_PHI[iTrack]);
                    Qy_back[0] += w1 * sin(1 * VELOTRACK_PHI[iTrack]);
                    Qy_back[1] += w2 * sin(2 * VELOTRACK_PHI[iTrack]);

                    // Eta-weighted
                    Qx_back_wEta[0] += wEta1 * cos(1 * VELOTRACK_PHI[iTrack]);
                    Qy_back_wEta[0] += wEta1 * sin(1 * VELOTRACK_PHI[iTrack]);
                }

                // Forward eta bin 1: 0.5–2.5
                if(VELOTRACK_ETA[iTrack] > 0.5 && VELOTRACK_ETA[iTrack] <= 2.5 ){
                    Qx_for[0][0] += w1 * cos(1 * VELOTRACK_PHI[iTrack]);
                    Qx_for[1][0] += w2 * cos(2 * VELOTRACK_PHI[iTrack]);
                    Qy_for[0][0] += w1 * sin(1 * VELOTRACK_PHI[iTrack]);
                    Qy_for[1][0] += w2 * sin(2 * VELOTRACK_PHI[iTrack]);

                    Qx_for_wEta[0][0] += wEta1 * cos(1 * VELOTRACK_PHI[iTrack]);
                    Qy_for_wEta[0][0] += wEta1 * sin(1 * VELOTRACK_PHI[iTrack]);

                    nPrimaryForTracks_1++;
                }

                // Forward eta bin 2: 2.5–4.0
                if(VELOTRACK_ETA[iTrack] > 2.5 && VELOTRACK_ETA[iTrack] <= 4.0 ){
                    Qx_for[0][1] += w1 * cos(1 * VELOTRACK_PHI[iTrack]);
                    Qx_for[1][1] += w2 * cos(2 * VELOTRACK_PHI[iTrack]);
                    Qy_for[0][1] += w1 * sin(1 * VELOTRACK_PHI[iTrack]);
                    Qy_for[1][1] += w2 * sin(2 * VELOTRACK_PHI[iTrack]);

                    Qx_for_wEta[0][1] += wEta1 * cos(1 * VELOTRACK_PHI[iTrack]);
                    Qy_for_wEta[0][1] += wEta1 * sin(1 * VELOTRACK_PHI[iTrack]);

                    nPrimaryForTracks_2++;
                }

                // Forward eta bin 3: 4.0–6.0
                if(VELOTRACK_ETA[iTrack] > 4.0 && VELOTRACK_ETA[iTrack] <= 6.0 ){
                    Qx_for[0][2] += w1 * cos(1 * VELOTRACK_PHI[iTrack]);
                    Qx_for[1][2] += w2 * cos(2 * VELOTRACK_PHI[iTrack]);
                    Qy_for[0][2] += w1 * sin(1 * VELOTRACK_PHI[iTrack]);
                    Qy_for[1][2] += w2 * sin(2 * VELOTRACK_PHI[iTrack]);

                    Qx_for_wEta[0][2] += wEta1 * cos(1 * VELOTRACK_PHI[iTrack]);
                    Qy_for_wEta[0][2] += wEta1 * sin(1 * VELOTRACK_PHI[iTrack]);

                    nPrimaryForTracks_3++;
                }

                // Forward inclusive eta: 0.5–6.0
                if(VELOTRACK_ETA[iTrack] > 0.5 && VELOTRACK_ETA[iTrack] <= 6.0 ){
                    Qx_for[0][3] += w1 * cos(1 * VELOTRACK_PHI[iTrack]);
                    Qx_for[1][3] += w2 * cos(2 * VELOTRACK_PHI[iTrack]);
                    Qy_for[0][3] += w1 * sin(1 * VELOTRACK_PHI[iTrack]);
                    Qy_for[1][3] += w2 * sin(2 * VELOTRACK_PHI[iTrack]);

                    Qx_for_wEta[0][3] += wEta1 * cos(1 * VELOTRACK_PHI[iTrack]);
                    Qy_for_wEta[0][3] += wEta1 * sin(1 * VELOTRACK_PHI[iTrack]);

                    nPrimaryForTracks++; 
                }
            }// end of track loop
            
            // ------------------
            // Require enough tracks in all subevents
            // ------------------
            if(nPrimaryForTracks<5)     continue;
            if(nPrimaryBackTracks<5)    continue;
            if(nPrimaryForTracks_1<5)   continue;
            if(nPrimaryForTracks_2<5)   continue;
            if(nPrimaryForTracks_3<5)   continue;

            // ------------------
            // Fill event object
            // ------------------
            evt->outGPSTIME       = GPSTIME;
            evt->outEVENTNUMBER   = EVENTNUMBER;
            evt->outPVX           = PVX[0];
            evt->outPVY           = PVY[0]; 
            evt->outPVZ           = PVZ[0];  
            evt->outRUNNUMBER     = RUNNUMBER;
            evt->outnBackTracks   = nBackTracks;
            evt->outnVeloClusters = nVeloClusters;
            evt->outnVeloTracks   = nVeloTracks;
            evt->outnEcalClusters = nEcalClusters;
            evt->outECalETot      = ECalETot;
            evt->outnLongTracks   = nLongTracks;
            evt->outnVPClusters   = nVPClusters;

            evt->out_Qmulti[0] = nPrimaryForTracks_1;
            evt->out_Qmulti[1] = nPrimaryForTracks_2;
            evt->out_Qmulti[2] = nPrimaryForTracks_3;
            evt->out_Qmulti[3] = nPrimaryBackTracks;

            for(int i = 0; i < 2; i++){
                evt->outQx_back[i] = Qx_back[i]; 
                evt->outQy_back[i] = Qy_back[i]; 
                evt->outQx_back_wEta[i] = Qx_back_wEta[i]; 
                evt->outQy_back_wEta[i] = Qy_back_wEta[i]; 

                for(int j = 0; j < 4; j++){
                    evt->outQx_for[i][j] = Qx_for[i][j]; 
                    evt->outQy_for[i][j] = Qy_for[i][j];
                    evt->outQx_for_wEta[i][j] = Qx_for_wEta[i][j]; 
                    evt->outQy_for_wEta[i][j] = Qy_for_wEta[i][j];
                }
            }

            // Fill the event tree
            outTree->Fill(); 
        } // End of event loop

        // Close and clean up per-file
        file->Close();
        delete file;
    } // End of input file loop

    // ==========================
    // Finalize output
    // ==========================
    outFile->cd();           // Go to output file directory
    outTree->Write();        // Write the tree to file
    outFile->Close();        // Close the ROOT file

    std::cout << "Done. Saved to event_plane_pbpb.root" << std::endl;

    // Clean up the event object
    delete evt; 
}


