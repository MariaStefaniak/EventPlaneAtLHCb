// ============================================================================
// ExecuteGlobalPolarizationAnalysisFilePrep
// 
// This macro is responsible for executing the Global Polarization analysis
// file preparation function (`GlobalPolarizationAnalysis_FilePrep`) on a set 
// of ROOT files containing Lambda candidates from PbPb collisions.
//
// It dynamically loads the implementation (`GlobalPolarizationAnalysis_FilePrep.C`) 
// and then loops over numbered input files, invoking the file preparation 
// macro on each.
//
// Adjust the path and loop index range as needed.
//
// Author: Maria Stefaniak <mzstefaniak@gmail.com>
// Year: 2025
// ============================================================================

void ExecuteGlobalPolarizationAnalysisFilePrep() {
    // Load the GlobalPolarizationAnalysis_FilePrep.C file with ACLiC compilation
    gROOT->ProcessLine(".L GlobalPolarizationAnalysis_FilePrep.C+");

    // Loop over file indices (currently set to only run for i = 0)
    for (int i = 0; i <= 0; ++i) {
        // Construct the input filename
        std::string fileName = Form("/Volumes/Mike_disc/Maria/PbPb/pbpb_%d.root", i);

        // Construct the ROOT command to call the analysis macro with filename and file index
        std::string command = Form("GlobalPolarizationAnalysis_FilePrep(\"%s\", %d)", fileName.c_str(), i);

        // Execute the command via ROOT interpreter
        gROOT->ProcessLine(command.c_str());
    }
}