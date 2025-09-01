//////////////////////////////////////////////////////////////
// Event Plane Calculation Macro
// Pb+Pb 2024 at LHCb - using Q-vectors from VELO tracks
// Author: Maria Stefaniak, The Ohio State University
// Description: Applies centering and shifting corrections to 
//              calculate final Event Plane angles and resolution.
//////////////////////////////////////////////////////////////

#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TVector2.h>
#include <TProfile2D.h>
#include <TCanvas.h>
#include <iostream>


double pi = TMath::Pi();


class Event {
    public:
        ULong64_t       outGPSTIME;
        ULong64_t       outEVENTNUMBER;
        Float_t         outPVX;   
        Float_t         outPVY;  
        Float_t         outPVZ;  
        UInt_t          outRUNNUMBER;
        Int_t           outnBackTracks;
        Int_t           outnVeloClusters;
        Int_t           outnVeloTracks;
        Int_t           outnEcalClusters;

        Int_t           outECalETot;
        Int_t           outnLongTracks;
        Int_t           outnVPClusters;


        Double_t        outQx_back[2];
        Double_t        outQy_back[2];
        Double_t        outQx_for[2][4];
        Double_t        outQy_for[2][4];

        Double_t        outQx_back_wEta[2];
        Double_t        outQy_back_wEta[2];
        Double_t        outQx_for_wEta[2][4];
        Double_t        outQy_for_wEta[2][4];


        Int_t           out_Qmulti[4];
    
        Event() : outGPSTIME(0), outEVENTNUMBER(0), outRUNNUMBER(0){}
    };

class EventPlane {
    public:
        ULong64_t       EVENTNUMBER;
        UInt_t          RUNNUMBER;
        Double_t        Psi1Full;
        Double_t        Psi2Full;
        Double_t        PsiBack[2];
        Double_t        PsiFor[2];
        Double_t        r1;
        Double_t        r2;          

        EventPlane() :  EVENTNUMBER(0), RUNNUMBER(0){}
    };

         
double makeShift(double psi, TH2D *hShiftSin, TH2D *hShiftCos, int iep, int n){
     
    double PsiShifted=psi;
    for(int j = 1; j < 9; j++)  PsiShifted += (2.0/(float)j) * (- hShiftSin->GetBinContent(iep+1, j)*cos(j*n*psi)+ hShiftCos->GetBinContent(iep+1, j)*sin(j*n*psi)) / n;
    return PsiShifted;
}
     
double keepPsiInPi(double psi){
   
    double psiNew = psi;
    if(abs(psi) < pi) psiNew = psi;
    else{
       if(psi < 0) psiNew = 2*pi - psi;
       if(psi > 0) psiNew = psi - 2*pi;
    }
    return psiNew;
}
double keepPsiInHalfPi(double psi){

    double psiNew = psi;
    if(abs(psi) < 0.5*pi) psiNew = psi;
    else{
        if(psi < 0) psiNew = pi + psi;
        if(psi > 0) psiNew = psi - pi;
    }
    return psiNew;
}

void calculateEventPlane(int EP_correction=1){

    cout << "EP_correction "<< EP_correction << endl;
    const int nrCentBins = 3;
// Define centrality bins based on number of VELO tracks (nVeloTracks)
    int CentralityBins[nrCentBins+1] = {14  , 126 ,   270 ,  2000};//350};
   // int CentralityBins[nrCentBins+1] = {400 ,  500,  600, 700};
   // int CentralityBins[nrCentBins+1] = {14, 200, 10000};

// Input ROOT file containing Q vectors from VELO tracks
    TFile* file = TFile::Open("/Users/stefaniak.9/OneDriveOSU/LHCb_Maria/EventPlane/PbPb2024_VELO/output/event_plane_pbpb_fulleta_weq1.root"); //centTests/weights_event_plane_pbpb_localtest.root");//
    if (!file || file->IsZombie()) {
        std::cerr << "Cannot open file." << std::endl;
        return;
    }
    // Get the tree
    TTree* tree = (TTree*)file->Get("EventPlaneTuple");
    if (!tree) {
        std::cerr << "Cannot find tree 'EventPlaneTuple'." << std::endl;
        return;
    }
      // Set up pointer to the Event object in the tree
    Event* evt = nullptr;
    tree->SetBranchAddress("event", &evt);

    // Create output file and tree
// Create output ROOT file for final event plane results
    TFile* outFile = new TFile("EP_PbPb2024_calculated_midEtaBin_test.root", "RECREATE");
    TTree* outTree = new TTree("EventPlaneTuple", "Event Plane");
    EventPlane *ep = nullptr;
    outTree->Branch("eventplane", &ep);


// Select which eta bin to use for event plane determination (1, 2, or 3)
    int iEta = 1; // mid Eta bin! 
    // some test histos:
    TH1D *hQdotQ[3][nrCentBins];
    TH1D *hQdotQback[3][nrCentBins];
    for(int iCent = 0; iCent < nrCentBins; iCent++){
        for(int ii = 0; ii < 3; ii ++){
            hQdotQ[ii][iCent] = new TH1D(Form("QdotQ_%d_cent%d", ii, iCent),Form("QdotQ_%d_cent%d", ii, iCent), 300, -1.5, 1.5);
            hQdotQback[ii][iCent] = new TH1D(Form("QdotQback_%d_cent%d", ii, iCent),Form("QdotQback_%d_cent%d", ii, iCent), 300, -1.5, 1.5);
        }
    }

    TH1D *hCentrality =    new TH1D("hCentrality", "Centrality; Centrality bin; Events", 500, 0, 5000);
    TH2D * hVeloClusters_EcalClusters = new TH2D("hVeloClusters_EcalClusters", "Velo clusters vs Ecal clusters; Velo Tracks; Ecal clusters", 300, 0, 5000, 300, 0, 800);
    TH2D * hVPClusters_EcalClusters = new TH2D("hnVPClusters_EcalClusters", "nVPClusters vs Ecal clusters; VP clusters; Ecal clusters", 300, 0, 5000, 300, 0, 800);
    TH2D * hnVeloTracks_outECalETot = new TH2D("hnVeloTracks_EcalClusters", "nVeloTracks vs Ecal clusters; Velo Tracks; outECalETot", 300, 0, 2000, 300, 0, 12000);



    // QxQy correction:
    TH2D *hQxQy_back[2][nrCentBins], *hQxQy_for[2][nrCentBins]; // iCent: 0-nbins-1 centrality bins and nbis = min bias // 0-2 eta bins, 3 - full forward eta
    TH2D *hQxQy_full[2][nrCentBins];
    int qbins = 20; double qmin = -10; double qmax = 10;
    for(int iCent = 0; iCent < nrCentBins; iCent++){
        for(int in = 0; in <2; in++){
            hQxQy_back[in][iCent] = new TH2D(Form("hQxQy_back_n%d_cent%d",in, iCent), Form("hQxQy_back_n%d_cent%d; Qx; Qy;",in,iCent), qbins, qmin, qmax, qbins,qmin, qmax );
            hQxQy_for[in][iCent]  = new TH2D(Form("hQxQy_for_n%d_cent%d_Eta%d",in, iCent, iEta), Form("hQxQy_for_n%d_cent%d_eta%d; Qx; Qy;",in,iCent, iEta), qbins, qmin, qmax, qbins,qmin, qmax );
            hQxQy_full[in][iCent] = new TH2D(Form("hQxQy_full_n%d_cent%d_Eta%d",in ,iCent, iEta), Form("hQxQy_full_n%d_cent%d_eta%d; Qx; Qy;",in,iCent, iEta), qbins, qmin, qmax, qbins,qmin, qmax );  
        }
    }
    // Psi calculations:
    TH1D *hPsi_back[2][nrCentBins], *hPsi_for[2][nrCentBins], *hPsi_full[2][nrCentBins];
    for(int in = 0; in < 2; in++){
        for(int iCent = 0; iCent < nrCentBins; iCent++){
            hPsi_back[in][iCent] = new TH1D(Form("hPsi_back_n%d_cent_%d", in, iCent), Form("hPsi_back_n%d_cent_%d", in, iCent), 100, -2*TMath::Pi(), 2*TMath::Pi()+0.1 );
            hPsi_for[in][iCent]  = new TH1D(Form("hPsi_for_n%d_cent_%d_eta%d", in, iCent, iEta), Form("hPsi_for_n%d_cent_%d_eta%d", in, iCent, iEta), 100, -2*TMath::Pi(), 2*TMath::Pi()+0.1 );
            hPsi_full[in][iCent] = new TH1D(Form("hPsi_full_n%d_cent_%d_eta%d", in, iCent, iEta), Form("hPsi_full._n%d_cent_%d_eta%d", in, iCent, iEta), 100, -2*TMath::Pi(), 2*TMath::Pi()+0.1 );   
        }
    }

    // Test Psi correlations:
    TH2D *hPsi_back_for[2][nrCentBins];
    for(int in = 0; in < 2; in++){
        for(int iCent = 0; iCent < nrCentBins; iCent++){
                hPsi_back_for[in][iCent] = new TH2D(Form("hPsi_back_for_n%d_cent_%d_eta%d", in, iCent, iEta), Form("hPsi_back_for_n%d_cent_%d_eta%d", in, iCent, iEta), 50, -TMath::Pi(), TMath::Pi()+0.1, 50, -TMath::Pi(), TMath::Pi()+0.1 );
        }
    }


    //shifting the EP:
    TProfile2D  *hEPshift_sin[nrCentBins], *hEPshift_cos[nrCentBins];
    for(int iCent = 0; iCent < nrCentBins; iCent++){
        hEPshift_sin[iCent] = new TProfile2D(Form("hEPshift_sin_cent%d",iCent), "", 6.0,-0.5,5.5,  9,0.5,9.5,  -2.0,2.0,""); // 0 - psi1 back, 1 - psi1 for, 2-  psi 1 full, 3 - psi2 back, 4- psi2 for ,  5- psi2 full j- moments, 
        hEPshift_cos[iCent] = new TProfile2D(Form("hEPshift_cos_cent%d",iCent), "", 6.0,-0.5,5.5,  9,0.5,9.5,  -2.0,2.0,""); // forward/backward , j- moments,     
    }

    double Resolution1[3] = {0,0,0};
    double Resolution2[3] = {0,0,0};
    int nrR[4] = {0,0,0};


     // QxQy centering:
    double Qxmean_back[2][nrCentBins], Qxmean_for[2][nrCentBins], Qxmean_full[2][nrCentBins]; 
    double Qymean_back[2][nrCentBins], Qymean_for[2][nrCentBins], Qymean_full[2][nrCentBins]; 
    TH2D *hQxQy_back_corr[2][nrCentBins], *hQxQy_for_corr[2][nrCentBins], *hQxQy_full_corr[2][nrCentBins]; 
    TProfile2D *hEPshift_sinIN[nrCentBins], *hEPshift_cosIN[nrCentBins];

    if(EP_correction > 1){
// Open file with previously calculated Q-vector centering weights
        TFile *fWeights = new TFile("EP_PbPb2024_weights_test.root", "READ");

        for(int iCent = 0; iCent < nrCentBins; iCent++){
            for(int in = 0; in < 2; in++){
                hQxQy_back_corr[in][iCent]  = (TH2D*)fWeights->Get(Form("hQxQy_back_n%d_cent%d",in, iCent));
                Qxmean_back[in][iCent]      = hQxQy_back_corr[in][iCent]->GetMean(1);
                Qymean_back[in][iCent]      = hQxQy_back_corr[in][iCent]->GetMean(2);

                hQxQy_for_corr[in][iCent] = (TH2D*)fWeights->Get(Form("hQxQy_for_n%d_cent%d_Eta%d",in,iCent, iEta));
                Qxmean_for[in][iCent]           = hQxQy_for_corr[in][iCent]->GetMean(1);
                Qymean_for[in][iCent]           = hQxQy_for_corr[in][iCent]->GetMean(2);

                hQxQy_full_corr[in][iCent] = (TH2D*)fWeights->Get(Form("hQxQy_full_n%d_cent%d_Eta%d",in,iCent, iEta));
                Qxmean_full[in][iCent]           = hQxQy_full_corr[in][iCent]->GetMean(1);
                Qymean_full[in][iCent]           = hQxQy_full_corr[in][iCent]->GetMean(2);
            }
        }

            // shifting:
        
        for(int iCent = 0; iCent < nrCentBins; iCent++){
                hEPshift_sinIN[iCent] = (TProfile2D*)fWeights->Get(Form("hEPshift_sin_cent%d",iCent));
                hEPshift_cosIN[iCent] = (TProfile2D*)fWeights->Get(Form("hEPshift_cos_cent%d",iCent));
        }
        
    }
    // Fill histogram
    Long64_t nEntries = tree->GetEntries();
    cout << "nEntries " << nEntries << endl;
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        if (evt->outRUNNUMBER ==   310318 && evt->outEVENTNUMBER == 93971618) cout << "93971618 " <<  endl;
        if(evt->outnVeloTracks > 1000 && evt->outnEcalClusters < 480) continue;

        hCentrality->Fill(evt->outnVeloTracks);
        hVeloClusters_EcalClusters->Fill(evt->outnVeloTracks, evt->outnEcalClusters);
        hVPClusters_EcalClusters->Fill(evt->outnVPClusters, evt->outnEcalClusters);
        hnVeloTracks_outECalETot->Fill(evt->outnVeloTracks, evt->outECalETot);
// Eta sign flip: forward v1 is negative, backward is positive
        int a = -1; int b = 1;  // as the w for Q vectors are equal to 1, we can here modify if we want to same sign or opposite for weights
        int CentBin = 2;
        // Centrality:
        for(int iCent = 0; iCent < nrCentBins; iCent++){
            if(evt->outnVeloTracks > CentralityBins[iCent] && evt->outnVeloTracks <= CentralityBins[iCent+1]) 
                CentBin = iCent;
        }
        double Qx_back[2],  Qy_back[2], Qx_for[2], Qy_for[2]; // 0: n = 1, and 1: n=2, for forward the iEta will determine which bin we take
        double Qx_full[2],  Qy_full[2]; // for final Psi determination
        double Psi_back[2], Psi_for[2], Psi_full[2]; 
        //set everything to 0:
        for(int in = 0; in <2; in++){Qx_back[in] = 0; Qy_back[in] = 0; Qx_for[in] = 0; Qy_for[in] = 0;  Qx_full[in] = 0; Qy_full[in] = 0;}
        // n = 1:
        // w = eta
      //  Qx_back[0] =  b*evt->outQx_back_wEta[0];           Qy_back[0] =  b*evt->outQy_back_wEta[0]; 
      //  Qx_for[0]  =  a*evt->outQx_for_wEta[0][iEta];      Qy_for[0]  =  a*evt->outQy_for_wEta[0][iEta];

      //w = 1
        Qx_back[0] =  b*evt->outQx_back[0];           Qy_back[0] =  b*evt->outQy_back[0]; 
        Qx_for[0]  =  a*evt->outQx_for[0][iEta];      Qy_for[0]  =  a*evt->outQy_for[0][iEta];


        Qx_full[0] =  Qx_back[0] + Qx_for[0];         Qy_full[0] =  Qy_back[0] +  Qy_for[0];
        // n = 2:
        Qx_back[1] =  evt->outQx_back[1];             Qy_back[1] = evt->outQy_back[1];
        Qx_for[1]  =  evt->outQx_for[1][iEta];        Qy_for[1]  = evt->outQy_for[1][iEta];
        Qx_full[1] =  Qx_back[1] + Qx_for[1];         Qy_full[1] = Qy_back[1] + Qy_for[1];
        

        hQxQy_back[0][CentBin] -> Fill(Qx_back[0], Qy_back[0]);
        hQxQy_back[1][CentBin] -> Fill(Qx_back[1], Qy_back[1]);
 
        hQxQy_for[0][CentBin] -> Fill(Qx_for[0], Qy_for[0]);
        hQxQy_for[1][CentBin] -> Fill(Qx_for[1], Qy_for[1]);
 
        hQxQy_full[0][CentBin] -> Fill(Qx_full[0], Qy_full[0]);
        hQxQy_full[1][CentBin] -> Fill(Qx_full[1], Qy_full[1]);
 

        if(EP_correction<2) continue;
        for(int in = 0; in <2; in++){
            Qx_back[in] -= Qxmean_back[in][CentBin];
            Qy_back[in] -= Qymean_back[in][CentBin];

            Qx_for[in]  -= Qxmean_for[in][CentBin];
            Qy_for[in]  -= Qymean_for[in][CentBin];

            Qx_full[in] -= Qxmean_full[in][CentBin];
            Qy_full[in] -= Qymean_full[in][CentBin];
        }
       // =====================================


        //test:
        int order = 0; 

       /*TVector2 Q1(a*evt->outQx_for[order][0], a*evt->outQy_for[order][0]);
        TVector2 Q2(a*evt->outQx_for[order][1], a*evt->outQy_for[order][1]);
        TVector2 Q3(a*evt->outQx_for[order][2], a*evt->outQy_for[order][2]);
        TVector2 Qback(b*evt->outQx_back[order], b*evt->outQy_back[order]);*/

        TVector2 Q1(Qx_for[0], Qy_for[0]);
        TVector2 Q2(a*evt->outQx_for[order][1], a*evt->outQy_for[order][1]);
        TVector2 Q3(a*evt->outQx_for[order][2], a*evt->outQy_for[order][2]);
        TVector2 Qback(Qx_back[0], Qy_back[0]);

        Q1 = Q1.Unit();
        Q2 = Q2.Unit();
        Q3 = Q3.Unit();
        Qback = Qback.Unit();

        hQdotQ[0][CentBin] -> Fill(Q1*Q2);
        hQdotQ[1][CentBin] -> Fill(Q1*Q3);
        hQdotQ[2][CentBin] -> Fill(Q2*Q3);
        hQdotQback[0][CentBin] -> Fill(Q1*Qback);
        hQdotQback[1][CentBin] -> Fill(Q2*Qback);
        hQdotQback[2][CentBin] -> Fill(Q3*Qback);

        // =====================================
        // fliiping signes:
     //   Qx_back[0] *= (-1); /// backward needs to be positive
     //   Qy_back[0] *= (-1);

      //  Qx_for[0] *= (-1); /// forward negative
      //  Qy_for[0] *= (-1);

     //   Qx_full[0] *= (-1); /// if I flip both signes I can do the full also with * (-1)
    //    Qy_full[0] *= (-1);

        // end of flipping signes

        Psi_back[0] = atan2(Qy_back[0], Qx_back[0]);
        Psi_for[0]  = atan2(Qy_for[0],  Qx_for[0]) ;
        Psi_full[0] = atan2(Qy_full[0], Qx_full[0]);

        Psi_back[1] = atan2(Qy_back[1], Qx_back[1]) * 0.5;
        Psi_for[1]  = atan2(Qy_for[1],  Qx_for[1])  * 0.5;
        Psi_full[1] = atan2(Qy_full[1], Qx_full[1]) * 0.5;


        double FullPsi[6] = {Psi_back[0], Psi_for[0], Psi_full[0], Psi_back[1], Psi_for[1], Psi_full[1]};
        for(int j = 1; j < 9; j++){
            for(int iep = 0; iep < 3; iep++){
               hEPshift_sin[CentBin] -> Fill(iep, j, sin(j*FullPsi[iep]));
               hEPshift_cos[CentBin] -> Fill(iep, j, cos(j*FullPsi[iep]));
            }
            for(int iep = 3; iep < 6; iep++){
               hEPshift_sin[CentBin] -> Fill(iep, j, sin(j*2*FullPsi[iep]));
               hEPshift_cos[CentBin] -> Fill(iep, j, cos(j*2*FullPsi[iep]));
            }
         }
        if(EP_correction<3) continue; 
        // shift Psi:
        double PsiFullShifted[6] = {0,0,0,0,0,0};
        if(EP_correction==3){
            for(int iep = 0; iep < 3; iep++){
                PsiFullShifted[iep] = makeShift(FullPsi[iep], hEPshift_sinIN[CentBin], hEPshift_cosIN[CentBin], iep, 1 );
            }
            for(int iep = 3; iep < 6; iep++){
                PsiFullShifted[iep] = makeShift(FullPsi[iep], hEPshift_sinIN[CentBin], hEPshift_cosIN[CentBin], iep, 2 );
            }
        }
        PsiFullShifted[0] = keepPsiInPi(PsiFullShifted[0]); //backward psi 1
        PsiFullShifted[1] = keepPsiInPi(PsiFullShifted[1]); //forward psi 1
        PsiFullShifted[2] = keepPsiInPi(PsiFullShifted[2]); //full psi 1
        
        PsiFullShifted[3] = keepPsiInHalfPi(PsiFullShifted[3]); //backward psi 2
        PsiFullShifted[4] = keepPsiInHalfPi(PsiFullShifted[4]); //forward psi 2
        PsiFullShifted[5] = keepPsiInHalfPi(PsiFullShifted[5]); //full psi 2

        /*hPsi_back[0][CentBin]->Fill(FullPsi[0]);
        hPsi_back[1][CentBin]->Fill(FullPsi[3]);

        hPsi_for[0][CentBin]->Fill(FullPsi[1]);
        hPsi_for[1][CentBin]->Fill(FullPsi[4]);

        hPsi_full[0][CentBin]->Fill(FullPsi[2]);
        hPsi_full[1][CentBin]->Fill(FullPsi[5]);*/


        hPsi_back[0][CentBin]->Fill(PsiFullShifted[0]);
        hPsi_back[1][CentBin]->Fill(PsiFullShifted[3]);

        hPsi_for[0][CentBin]->Fill(PsiFullShifted[1]);
        hPsi_for[1][CentBin]->Fill(PsiFullShifted[4]);

        hPsi_full[0][CentBin]->Fill(PsiFullShifted[2]);
        hPsi_full[1][CentBin]->Fill(PsiFullShifted[5]);


        hPsi_back_for[0][CentBin]->Fill(PsiFullShifted[0], PsiFullShifted[1]);
        hPsi_back_for[1][CentBin]->Fill(PsiFullShifted[3], PsiFullShifted[4]);

        double r1 = cos(1*(PsiFullShifted[0]-PsiFullShifted[1]));
        double r2 = cos(2*(PsiFullShifted[3]-PsiFullShifted[4]));
        Resolution1[CentBin] += r1; 
        Resolution2[CentBin] += r2;  
        nrR[CentBin]++;


        
        //Save to tree:
// Save output variables to tree for each event
        ep->EVENTNUMBER     =   evt->outEVENTNUMBER;
        ep->RUNNUMBER       =   evt->outRUNNUMBER;
        ep->Psi1Full        =   PsiFullShifted[2];
        ep->Psi2Full        =   PsiFullShifted[5];
        ep->r1              =   r1;
        ep->r2              =   r2; 
        ep->PsiBack[0]      =   PsiFullShifted[0];
        ep->PsiBack[1]      =   PsiFullShifted[3];
        ep->PsiFor[0]       =   PsiFullShifted[1];
        ep->PsiFor[1]       =   PsiFullShifted[4];
        outTree->Fill(); 
        
    }//End of loop over Events
    
    for(int iCent = 0; iCent <nrCentBins; iCent++){
        Resolution1[iCent] = sqrt(2*Resolution1[iCent]/nrR[iCent]) *100;
        Resolution2[iCent] = sqrt(2*Resolution2[iCent]/nrR[iCent]) *100;
        cout << "Cent: " << iCent << endl;
        cout << "R1: " << Resolution1[iCent] << endl;
        cout << "R2: " << Resolution2[iCent] << endl;
    }

    if(EP_correction<3){
// Create file to store centering/shifting histograms for corrections
        TFile *weightsFile = new TFile("EP_PbPb2024_weights_test.root", "RECREATE");
        for(int iCent = 0; iCent <nrCentBins; iCent++){

            hEPshift_sin[iCent]->Write();
            hEPshift_cos[iCent]->Write();

            for(int in = 0; in < 2; in++){
                hQxQy_back[in][iCent]->Write();
                hQxQy_for[in][iCent]->Write();
                hQxQy_full[in][iCent]->Write();
            
        }}
        weightsFile->Close();
    }
    


    TCanvas *cTestQ = new TCanvas();
    cTestQ->Divide(3,2);
    cTestQ->cd(1);
    hQdotQ[0][0]->Draw();
    cTestQ->cd(2);
    hQdotQ[1][0]->Draw();
    cTestQ->cd(3);
    hQdotQ[2][0]->Draw();
    cTestQ->cd(4);
    hQdotQ[0][1]->Draw();
    cTestQ->cd(5);
    hQdotQ[1][1]->Draw();
    cTestQ->cd(6);
    hQdotQ[2][1]->Draw();

    TCanvas *cTestQback = new TCanvas();
    cTestQback->Divide(3,2);
    cTestQback->cd(1);
    hQdotQback[0][0]->Draw();
    cTestQback->cd(2);
    hQdotQback[1][0]->Draw();
    cTestQback->cd(3);
    hQdotQback[2][0]->Draw();
    cTestQback->cd(4);
    hQdotQback[0][2]->Draw();
    cTestQback->cd(5);
    hQdotQback[1][2]->Draw();
    cTestQback->cd(6);
    hQdotQback[2][2]->Draw();



    TCanvas *cCentrality = new TCanvas();
    cCentrality->SetLogy();
    hCentrality->Draw();

    TCanvas *cVeloClusters_EcalClusters = new TCanvas();
    cVeloClusters_EcalClusters->SetLogz();
    hVeloClusters_EcalClusters->Draw("colz");


    TCanvas *cVPClusters_EcalClusters = new TCanvas();
    cVPClusters_EcalClusters->SetLogz();
    hVPClusters_EcalClusters->Draw("colz");



    TCanvas *cnVeloTracks_outECalETot = new TCanvas();
    cnVeloTracks_outECalETot->SetLogz();
    hnVeloTracks_outECalETot->Draw("colz");


    int cc = 0;
    TCanvas *c1 = new TCanvas();
    c1->Divide(3,2);
    c1->cd(1);
    hQxQy_back[0][cc]->Draw("colz");
    c1->cd(2);
    hQxQy_for[0][cc]->Draw("colz");
    c1->cd(3);
    hQxQy_full[0][cc]->Draw("colz");
    c1->cd(4);
    hQxQy_back[1][cc]->Draw("colz");
    c1->cd(5);
    hQxQy_for[1][cc]->Draw("colz");
    c1->cd(6);
    hQxQy_full[1][cc]->Draw("colz");


    TCanvas *c2 = new TCanvas();
    c2->Divide(3,2);
    c2->cd(1);
    hPsi_back[0][cc]->Draw();
    c2->cd(2);
    hPsi_for[0][cc]->Draw();
    c2->cd(3);
    hPsi_full[0][cc]->Draw();
    c2->cd(4);
    hPsi_back[1][cc]->Draw();
    c2->cd(5);
    hPsi_for[1][cc]->Draw();
    c2->cd(6);
    hPsi_full[1][cc]->Draw();

    TCanvas *c3 = new TCanvas();
    c3->Divide(3,2);
    c3->cd(1);
    hPsi_back_for[0][0]->Draw("colz");
    c3->cd(2);
    hPsi_back_for[0][1]->Draw("colz");
    c3->cd(3);
    hPsi_back_for[0][2]->Draw("colz");
    c3->cd(4);
    hPsi_back_for[1][0]->Draw("colz");
    c3->cd(5);
    hPsi_back_for[1][1]->Draw("colz");
    c3->cd(6);
    hPsi_back_for[1][2]->Draw("colz");

    
    TCanvas *projectionTest = new TCanvas();
    hPsi_back_for[1][2]->ProjectionX()->Draw();


    // Save and close
    outFile->cd();
    outTree->Write();
    outFile->Close();
    std::cout << "Done. Saved to EventPlanes_forLowEtaBin_w1_FULL.root" << std::endl;
    delete ep;

}   
