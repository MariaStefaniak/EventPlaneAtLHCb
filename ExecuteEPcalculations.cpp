//////////////////////////////////////////////////////////
// ExecuteEPcalculations Macro
// LHCb: Pb+Pb 2024 Event Plane Analysis
// Executes the calculateEventPlane() function over selected input files
// Author: Maria Stefaniak (The Ohio State University)
//////////////////////////////////////////////////////////

void ExecuteEPcalculations() {
    // Loop over the desired range of input file indices
    // In this example, it only processes file index 1 (ii = 1)
    for (int ii = 1; ii <= 1; ii++) {

        // Dynamically load and compile the calculateEventPlane.cpp script
        gROOT->ProcessLine(".L calculateEventPlane.cpp+");

        // Construct the command to call calculateEventPlane with the given index
        std::string command = Form("calculateEventPlane(%d)", ii);

        // Execute the command in the ROOT interpreter
        gROOT->ProcessLine(command.c_str());
    }
}