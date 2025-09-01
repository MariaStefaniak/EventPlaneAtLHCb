
# Event Plane Determination in PbPb Collisions at LHCb (2024)

The goal of this project is to determine the Event Plane from the PbPb collision data collected by LHCb in 2024, using VELO tracks. The input dataset consists of multiple terabytes of VELO and Long tracks stored in the Analysis Production (AP): https://lhcb-productions.web.cern.ch/ana-prod/productions/?wg=ift&analysis=velotracks_pbpb2024. The primary code that runs over all AP output ROOT files is `EventPlaneAnalysis.cxx`. You can also find this code in the public CERN work folder at:
/afs/cern.ch/work/m/mastefan/public/EventPlanePbPb

## Step 1: Q-vector Extraction from VELO Tracks

To run the analysis:
> root -l  
> .L EventPlaneAnalysis.cxx++  
> EventPlaneAnalysis()

It will take a while to process all files. You can optionally run the code in parallel jobs or split the dataset and analyze it in parts. To do that, modify the `GetFilteredRootFiles()` function and adjust the file-matching pattern:
TRegexp pattern("^00274156_00000[0-4][0-9][0-9]_1\.tuple_pbpb2024\.root$");
Modify the [0-4][0-9][0-9] part to define which files are read in.

To run locally, define the list of ROOT files and output file explicitly:
std::vector<std::string> fileNames = {...};  
TFile* outFile = new TFile("output.root", "RECREATE");

This macro calculates Q-vectors for each event using VELO tracks. Applied event cuts:
- nPVs == 1
- nBackTracks > 10
- -100 < PVZ < 100 mm
- nVeloTracks > 15
- nPrimaryVELOtracks (in forward and backward) > 5

Backward track correction: if IS_BACKWARD == 1 then
  eta = -eta  
  phi = phi + π

Two sets of Q1 vectors are calculated: one with weight w=1, one with w=eta. Be cautious:
- For w=1 → sum of Q-vectors does not need normalization.
- For w=eta → normalization IS required. (Default is w=1, so normalization is skipped.)

Q-vectors are calculated for:
- all backward tracks
- forward tracks in 3 eta bins:  
  Bin 1: 0.5–2.5  
  Bin 2: 2.5–4.0  
  Bin 3: 4.0–6.0

Each event also stores:
GPSTIME, EVENTNUMBER, PVX, PVY, PVZ, RUNNUMBER, nBackTracks, nVeloClusters, nVeloTracks, nEcalClusters, ECalETot, nLongTracks, nVPClusters, nPrimaryForTracks_1, nPrimaryForTracks_2, nPrimaryForTracks_3, nPrimaryBackTracks

Example output:
/afs/cern.ch/work/m/mastefan/public/EventPlanePbPb/July28/EventPlanes_PbPb2024.root

## Step 2: Calculate Event Plane

This step runs locally on the Q-vector output. Files:
- calculateEventPlane.cpp
- ExecuteEPcalculations.cpp

Adjust file paths in `calculateEventPlane.cpp` and then run:
> root -l ExecuteEPcalculations.cpp

Edit the following parameters:
1. Centrality bins by nVeloTracks:
   int CentralityBins[nrCentBins+1] = {14, 126, 270, 2000};
2. Eta bin used for EP:  
   int iEta = 1; // 1 = midEta (0.5–2.5)
3. Input file:
   TFile* file = TFile::Open("Qvectors_input.root");
4. Output EP file:  
   TFile* outFile = new TFile("EP_output.root", "RECREATE");
5. Set weights file name consistently in 2 places:
   - TFile *fWeights = new TFile("weights.root", "READ");
   - TFile *weightsFile = new TFile("weights.root", "RECREATE");

Three steps:
- Step 1: Create Qx/Qy histograms, store in weights file
- Step 2: Center Q-vectors and prepare Ψ-shift histograms (save in weights file)
- Step 3: Shift Ψ and determine final EP angles and resolution

Output tree stores:
  EVENTNUMBER, RUNNUMBER, Psi1Full, Psi2Full, r1, r2, PsiBack[0/1], PsiFor[0/1]

Note on flipping sign:  
int a = -1; // forward  
int b = 1;  // backward  
→ our Q-vectors have weight 1, but v1 is negative in forward η, positive in backward η

## Step 3: Matching EP with Lambda Candidates

Matching is based on RUNNUMBER and EVENTNUMBER. Files:
- ExecuteGlobalPolarizationAnalysisFilePrep.C
- GlobalPolarizationAnalysis_FilePrep.C
- GlobalPolarizationAnalysis_FilePrep.h

Set path in ExecuteGlobalPolarizationAnalysisFilePrep.C:  
std::string fileName = Form("your_path/...");
Adjust loop over files and run:  
> root -l ExecuteGlobalPolarizationAnalysisFilePrep.C

Mapping with std::map is used for fast matching. Objects created: EventPlane, Event, Lambda, Daughter.

Debug prints are included — currently ~20% match rate, probably because VELO AP does not contain all triggered events. Consider relaxing Lambda cuts or verifying event coverage.

## Contact

For any questions, please contact:  
Maria Zofia Stefaniak – mzstefaniak@gmail.com
