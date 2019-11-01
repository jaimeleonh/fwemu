#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TPaveText.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TDatime.h"
#include "TAxis.h"
#include "TStyle.h"
#include "TH1.h"
#include "TH2.h"
#include "TSystem.h"
#include <string>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <map>
#include <fstream>
#include <cstring>
#include <string>
#include <sstream>
#include <time.h>
#include "TLegend.h"
#include "TEfficiency.h"
#ifdef __MAKECINT__
#pragma link C++ class vector<short>;
#pragma link C++ class vector<float>;
#endif
using namespace std;




/***************************************************************************************
 *
 * 				CONSTANTES INICIALES
 *
 ***************************************************************************************/

const bool debug = false; 

/***************************************************************************************
 *
 *			           SUBPROGRAMAS
 *
 * *************************************************************************************/

bool sameFit (std::vector <int> wires1, std::vector <int> wires2, std::vector <int> tdcs1, std::vector <int> tdcs2, std::vector <short> lateralities1, std::vector <short> lateralities2) {
  if ((wires1[0]!=wires2[0] || abs(tdcs1[0]-tdcs2[0])>1 || lateralities1[0]!=lateralities2[0])&&!(wires1[0]==-1 && wires2[0] == -1)) return false; 
  if ((wires1[1]!=wires2[1] || abs(tdcs1[1]-tdcs2[1])>1 || lateralities1[1]!=lateralities2[1])&&!(wires1[1]==-1 && wires2[1] == -1)) return false; 
  if ((wires1[2]!=wires2[2] || abs(tdcs1[2]-tdcs2[2])>1 || lateralities1[2]!=lateralities2[2])&&!(wires1[2]==-1 && wires2[2] == -1)) return false; 
  if ((wires1[3]!=wires2[3] || abs(tdcs1[3]-tdcs2[3])>1 || lateralities1[3]!=lateralities2[3])&&!(wires1[3]==-1 && wires2[3] == -1)) return false; 
  if ((wires1[4]!=wires2[4] || abs(tdcs1[4]-tdcs2[4])>1 || lateralities1[4]!=lateralities2[4])&&!(wires1[4]==-1 && wires2[4] == -1)) return false; 
  if ((wires1[5]!=wires2[5] || abs(tdcs1[5]-tdcs2[5])>1 || lateralities1[5]!=lateralities2[5])&&!(wires1[5]==-1 && wires2[5] == -1)) return false; 
  if ((wires1[6]!=wires2[6] || abs(tdcs1[6]-tdcs2[6])>1 || lateralities1[6]!=lateralities2[6])&&!(wires1[6]==-1 && wires2[6] == -1)) return false; 
  if ((wires1[7]!=wires2[7] || abs(tdcs1[7]-tdcs2[7])>1 || lateralities1[7]!=lateralities2[7])&&!(wires1[7]==-1 && wires2[7] == -1)) return false; 
  return true; 
}

int qualityGroup (int quality) {
  if (quality == 1 || quality == 2 || quality == 5) return 1; 
  if (quality == 3 || quality == 4 || quality == 7) return 2; 
  if (quality == 6) return 3; 
  if (quality == 8) return 4; 
  if (quality == 9) return 5;
  return -1;  
}



/**************************************************************************************
 * 				Declaramos histogramas
**************************************************************************************/


/***************************************************************************************
 *
 *				PROGRAMA PRINCIPAL
 *	
 * *************************************************************************************/

void comparator_fw(){



  gStyle->SetOptStat(1111111);

// Cargo la ntupla de las primitivas
TFile *f = new TFile("./primitives.root");
TTree *ntuple = (TTree*)f->Get("ntuple");

UInt_t nTrigsEm = 0;
UInt_t eventBX = 0; 
std::vector<float> *shift = 0;
std::vector<float> *positionEm= 0;
std::vector<float> *directionEm= 0;
std::vector<int> *bxTimeEm= 0;
std::vector<short> *qualityEm= 0;
std::vector<short> *wheelEm= 0;
std::vector<short> *sectorEm= 0;
std::vector<short> *stationEm= 0;
std::vector<std::vector<short>> *lateralitiesEm= 0;
std::vector<std::vector<int>> *wiresEm= 0;
std::vector<std::vector<int>> *tdcsEm= 0;

ntuple->SetBranchAddress("eventBX",  &eventBX);
ntuple->SetBranchAddress("numberOfTrigs",  &nTrigsEm);
ntuple->SetBranchAddress("Quality",&qualityEm);
ntuple->SetBranchAddress("Position",&positionEm);
ntuple->SetBranchAddress("Direction",&directionEm);
ntuple->SetBranchAddress("Time",&bxTimeEm);
ntuple->SetBranchAddress("Shift",&shift);
ntuple->SetBranchAddress("Wheel",&wheelEm);
ntuple->SetBranchAddress("Sector",&sectorEm);
ntuple->SetBranchAddress("Station",&stationEm);
ntuple->SetBranchAddress("Lateralities",&lateralitiesEm);
ntuple->SetBranchAddress("Cellnumber",&wiresEm);
ntuple->SetBranchAddress("tdcTime",&tdcsEm);

// Cargo la ntupla del firmware
TFile *f1 = new TFile("./outputfw.root");
TTree *ntuple1 = (TTree*)f1->Get("ntuple");

UInt_t nTrigsFw = 0;
std::vector<float> *positionFw= 0;
std::vector<float> *directionFw= 0;
std::vector<int> *bxTimeFw= 0;
std::vector<short> *qualityFw= 0;
std::vector<short> *wheelFw= 0;
std::vector<short> *sectorFw= 0;
std::vector<short> *stationFw= 0;
std::vector<std::vector<short>> *lateralitiesFw= 0;
std::vector<std::vector<int>> *wiresFw= 0;
std::vector<std::vector<int>> *tdcsFw= 0;

ntuple1->SetBranchAddress("numberOfTrigs",&nTrigsFw);
ntuple1->SetBranchAddress("Position",&positionFw);
ntuple1->SetBranchAddress("Direction",&directionFw);
ntuple1->SetBranchAddress("Time",&bxTimeFw);
ntuple1->SetBranchAddress("Quality",&qualityFw);
ntuple1->SetBranchAddress("Wheel",&wheelFw);
ntuple1->SetBranchAddress("Sector",&sectorFw);
ntuple1->SetBranchAddress("Station",&stationFw);
ntuple1->SetBranchAddress("Lateralities",&lateralitiesFw);
ntuple1->SetBranchAddress("Cellnumber",&wiresFw);
ntuple1->SetBranchAddress("tdcTime",&tdcsFw);

int n1 = ntuple->GetEntries(); 
int n2 = ntuple1->GetEntries(); 

int nEntries = -1; 
if (n1 < n2) nEntries = n1; 
else nEntries = n2;  

cout << "Emulator segments: " << ntuple->GetEntries() << " AB7 outputs " << ntuple1->GetEntries() << endl; 

/**************************************************************************************
 * 				Declaramos histogramas
**************************************************************************************/

std::map<std::string, TH1*> m_plots;
std::map<std::string, TH2*> m_plots2;

std::vector <std::string> qualityCategories = {"All", "sameQuality", "emulBetterQuality", "fwBetterQuality"};  
std::vector <std::string> qualityNumbers = {"Q1", "Q2", "Q3", "Q4", "Q6", "Q8", "Q9"};  

int nbinx = 21;    double xmin = 0.105;
int nbinTime = 21; double timemin = 10.5;
int nbinTan = 9;  double tanmin = 0.0225;
int nbinBX = 15; double bxmin = 7.5; 


for (auto & category : qualityCategories) {

  m_plots2["h_tanPhi2D"+category] = new TH2F (("h_tanPhi2D" + category).c_str(), "TanPhi firmware vs TanPhi emulator; TanPhi FW; TanPsi Emul", 400, -1,1, 400, -1, 1);
  m_plots2["h_pos2D"+category] = new TH2F (("h_pos2D" + category).c_str(), "Position firmware vs Position emulator; Position FW; Position Emul", 250, -250,250, 250, -250, 250);
  m_plots2["h_time2D"+category] = new TH2F (("h_time2D" + category).c_str(), "BxTime firmware vs BxTime emulator; BxTime FW; BxTime Emul", 1000, 0, 90000, 1000, 0, 90000);
  m_plots["h_time"+category] = new TH1F (("h_time" + category).c_str(), "BxTime firmware - BxTime emulator; #Delta BxTime (ns); Entries", nbinTime, -timemin, timemin);
  m_plots["h_time_whenGoodX"+category] = new TH1F (("h_time_whenGoodX" + category).c_str(), "BxTime firmware - BxTime emulator when goodX; #Delta BxTime (ns); Entries", nbinTime, -timemin, timemin);
  m_plots["h_BX"+category] = new TH1F (("h_BX" + category).c_str(), "BxTime firmware - BxTime emulator; #Delta BxTime (ns); Entries", 7, -87.5, 87.5);
  m_plots["h_tanPhi"+category] = new TH1F (("h_tanPhi" + category).c_str(), "TanPhi firmware - TanPhi emulator; #Delta TanPhi (adim); Entries", nbinTan, -tanmin, tanmin);
  m_plots["h_tanPhi_whenGoodX"+category] = new TH1F (("h_tanPhi_whenGoodX" + category).c_str(), "TanPhi firmware - TanPhi emulator when goodX; #Delta TanPhi (adim); Entries", nbinTan, -tanmin, tanmin);
  m_plots["h_pos"+category] = new TH1F (("h_pos" + category).c_str(), "Position firmware - Position emulator; #Delta Position (cm); Entries / 0.01 cm", nbinx, -xmin , xmin );
  m_plots["h_BXFW"+category] = new TH1F (("h_BXFW" + category).c_str(), "BX firmware - eventBX; #Delta BX (adim.); Entries / BX", nbinBX, -bxmin , bxmin );
  m_plots["h_BXEmul"+category] = new TH1F (("h_BXEmul" + category).c_str(), "BX emulator - eventBX; #Delta BX (adim.); Entries / BX", nbinBX, -bxmin , bxmin );

}

for (auto & category : qualityNumbers) {

  m_plots2["h_tanPhi2D"+category] = new TH2F (("h_tanPhi2D" + category).c_str(), "TanPhi firmware vs TanPhi emulator; TanPhi FW; TanPsi Emul", 400, -1,1, 400, -1, 1);
  m_plots2["h_pos2D"+category] = new TH2F (("h_pos2D" + category).c_str(), "Position firmware vs Position emulator; Position FW; Position Emul", 250, -250,250, 250, -250, 250);
  m_plots2["h_time2D"+category] = new TH2F (("h_time2D" + category).c_str(), "BxTime firmware vs BxTime emulator; BxTime FW; BxTime Emul", 1000, 0, 90000, 1000, 0, 90000);
  m_plots["h_time"+category] = new TH1F (("h_time" + category).c_str(), "BxTime firmware - BxTime emulator; #Delta BxTime (ns); Entries", nbinTime, -timemin, timemin);
  m_plots["h_time_whenGoodX"+category] = new TH1F (("h_time_whenGoodX" + category).c_str(), "BxTime firmware - BxTime emulator when goodX; #Delta BxTime (ns); Entries", nbinTime, -timemin, timemin);
  m_plots["h_BX"+category] = new TH1F (("h_BX" + category).c_str(), "BxTime firmware - BxTime emulator; #Delta BxTime (ns); Entries", 7, -87.5, 87.5);
  m_plots["h_tanPhi"+category] = new TH1F (("h_tanPhi" + category).c_str(), "TanPhi firmware - TanPhi emulator; #Delta TanPhi (adim); Entries", nbinTan, -tanmin, tanmin);
  m_plots["h_tanPhi_whenGoodX"+category] = new TH1F (("h_tanPhi_whenGoodX" + category).c_str(), "TanPhi firmware - TanPhi emulator when goodX; #Delta TanPhi (adim); Entries", nbinTan, -tanmin, tanmin);
  m_plots["h_pos"+category] = new TH1F (("h_pos" + category).c_str(), "Position firmware - Position emulator; #Delta Position (cm); Entries / 0.01 cm", nbinx, -xmin , xmin );
  m_plots["h_BXFW"+category] = new TH1F (("h_BXFW" + category).c_str(), "BX firmware - eventBX; #Delta BX (adim.); Entries / BX", nbinBX, -bxmin , bxmin );
  m_plots["h_BXEmul"+category] = new TH1F (("h_BXEmul" + category).c_str(), "BX emulator - eventBX; #Delta BX (adim.); Entries / BX", nbinBX, -bxmin , bxmin );
}

/**************************************************************************************
 * 					BUCLE
 *************************************************************************************/ 					

//for (Int_t i = 0; i < 1000; i++){
for (Int_t i = 0; i < nEntries; i++){
//cout << "************************************************************************" << endl;  
//  if (i!= 176) continue;
  ntuple->GetEntry(i);
  ntuple1->GetEntry(i);

  for (std::size_t iTrigEm = 0; iTrigEm < nTrigsEm; iTrigEm++) {
    int bestTrigFw = -1; 
    int bestTimeFw = 99999;     

    for (std::size_t iTrigFw = 0; iTrigFw < nTrigsFw; iTrigFw++) {

      if ( (wheelEm->at(iTrigEm) == wheelFw->at(iTrigFw)) && (stationEm->at(iTrigEm) == stationFw->at(iTrigFw)) && (sectorEm->at(iTrigEm) == sectorFw->at(iTrigFw)) ) {
        /*if (debug) {
          cout << "++++++++++++++++++++" << endl; 
	  for (int deb = 0; deb < 8; deb++){
            cout << wiresEm->at(iTrigEm)[deb] << " " << wiresFw->at(iTrigFw)[deb] << " " << tdcsEm->at(iTrigEm)[deb] << " "<< tdcsFw->at(iTrigFw)[deb] << " "<< lateralitiesEm->at(iTrigEm)[deb] << " "<< lateralitiesFw->at(iTrigFw)[deb] << endl;  
          } 
        }*/
        if (!sameFit(wiresEm->at(iTrigEm), wiresFw->at(iTrigFw), tdcsEm->at(iTrigEm), tdcsFw->at(iTrigFw), lateralitiesEm->at(iTrigEm), lateralitiesFw->at(iTrigFw))) continue;
        if (debug) {
          cout << "--------------------" << endl; 
	  cout << "Wh:" << wheelEm->at(iTrigEm) << " Se:" << sectorEm->at(iTrigEm) << " St:" << stationEm->at(iTrigEm) << endl; 
	  for (int deb = 0; deb < 8; deb++){
            cout << wiresEm->at(iTrigEm)[deb] << " " << wiresFw->at(iTrigFw)[deb] << " " << tdcsEm->at(iTrigEm)[deb] << " "<< tdcsFw->at(iTrigFw)[deb] << " "<< lateralitiesEm->at(iTrigEm)[deb] << " "<< lateralitiesFw->at(iTrigFw)[deb] << endl;  
	  } 
	    cout << "PosEmul "    << positionEm->at(iTrigEm) - shift->at(iTrigEm)  << " Pos FW "  << positionFw->at(iTrigFw) << endl;
	    cout << "TimeEmul "   << bxTimeEm->at(iTrigEm)    << " Time FW " << bxTimeFw->at(iTrigFw)                        << endl;
	    cout << "TanPsiEmul " << directionEm->at(iTrigEm) << " TanPsi FW "  << -directionFw->at(iTrigFw)                     << endl;
	}
        int deltaT0 = abs (bxTimeEm->at(iTrigEm) - bxTimeFw->at(iTrigFw));
        if (deltaT0 <= bestTimeFw) {
          bestTrigFw = iTrigFw; 
          bestTimeFw = deltaT0; 
        }
      } // If same chamber
    } // loop over fw

    // Filling plots
    std::vector <std::string> qualityCategories = {"All", "sameQuality", "emulBetterQuality", "fwBetterQuality"};  
    if (bestTrigFw != -1) {
      std::vector <std::string> outCategories = {"All"};
      if (qualityGroup(qualityEm->at(iTrigEm)) == qualityGroup(qualityFw->at(bestTrigFw))) 
      //if (qualityEm->at(iTrigEm)  == qualityFw->at(bestTrigFw)) 
        {
	  outCategories.push_back("sameQuality");
          if (qualityEm->at(iTrigEm) == 1) outCategories.push_back("Q1");
          else if (qualityEm->at(iTrigEm) == 2) outCategories.push_back("Q2");
          else if (qualityEm->at(iTrigEm) == 3) outCategories.push_back("Q3");
          else if (qualityEm->at(iTrigEm) == 4) outCategories.push_back("Q4");
          else if (qualityEm->at(iTrigEm) == 6) outCategories.push_back("Q6");
          else if (qualityEm->at(iTrigEm) == 8) outCategories.push_back("Q8");
          else if (qualityEm->at(iTrigEm) == 9) outCategories.push_back("Q9");
        } 	
      else if (qualityEm->at(iTrigEm) > qualityFw->at(bestTrigFw)) outCategories.push_back("emulBetterQuality"); 	
      else if (qualityEm->at(iTrigEm) < qualityFw->at(bestTrigFw)) outCategories.push_back("fwBetterQuality"); 

//  	if (fabs(positionFw->at(bestTrigFw) + shift->at(iTrigEm) - positionEm->at(iTrigEm)) > 0.5) { cout << i << endl;  }
  	//if (fabs(-directionFw->at(bestTrigFw) - directionEm->at(iTrigEm)) > 0.015) { cout << i << endl; }  
//        cout << positionFw->at(bestTrigFw) + shift->at(iTrigEm) - positionEm->at(iTrigEm) << endl; }
        if (debug) {
	  cout << "--------------------" << endl;
	  cout << "Wh:" << wheelEm->at(iTrigEm) << " Se:" << sectorEm->at(iTrigEm) << " St:" << stationEm->at(iTrigEm) << endl; 
	  for (int deb = 0; deb < 8; deb++){
            cout << wiresEm->at(iTrigEm)[deb] << " " << wiresFw->at(bestTrigFw)[deb] << " " << tdcsEm->at(iTrigEm)[deb] << " "<< tdcsFw->at(bestTrigFw)[deb] << " "<< lateralitiesEm->at(iTrigEm)[deb] << " "<< lateralitiesFw->at(bestTrigFw)[deb] << endl;  
	  } 
	    cout << "PosEmul "    << positionEm->at(iTrigEm) - shift->at(iTrigEm)  << " Pos FW "  << positionFw->at(bestTrigFw) << endl;
	    cout << "TimeEmul "   << bxTimeEm->at(iTrigEm)    << " Time FW " << bxTimeFw->at(bestTrigFw)                        << endl;
	    cout << "TanPsiEmul " << directionEm->at(iTrigEm) << " TanPsi FW "  << -directionFw->at(bestTrigFw)                     << endl;
	}
      for (auto & outCat : outCategories) {
        m_plots2["h_tanPhi2D"+outCat]->Fill(-directionFw->at(bestTrigFw), directionEm->at(iTrigEm)); 
        m_plots["h_tanPhi"+outCat]->Fill(-directionFw->at(bestTrigFw) - directionEm->at(iTrigEm)); 
        m_plots2["h_pos2D"+outCat]->Fill(positionFw->at(bestTrigFw)  + shift ->at(iTrigEm) , positionEm->at(iTrigEm)); 
        m_plots["h_pos"+outCat]->Fill(positionFw->at(bestTrigFw) + shift->at(iTrigEm) - positionEm->at(iTrigEm)); 
        m_plots2["h_time2D"+outCat]->Fill(bxTimeFw->at(bestTrigFw), bxTimeEm->at(iTrigEm)); 
        m_plots["h_time"+outCat]->Fill(bxTimeFw->at(bestTrigFw) - bxTimeEm->at(iTrigEm)); 
        m_plots["h_BX"+outCat]->Fill(bxTimeFw->at(bestTrigFw) - bxTimeEm->at(iTrigEm)); 
        m_plots["h_BXFW"+outCat]->Fill(round(bxTimeFw->at(bestTrigFw)/25.)  - eventBX);
        m_plots["h_BXEmul"+outCat]->Fill(round(bxTimeEm->at(iTrigEm)/25.)  - eventBX); 
        if ( fabs(positionFw->at(bestTrigFw) + shift->at(iTrigEm) - positionEm->at(iTrigEm) ) <= xmin / ((double) nbinx )){
	  m_plots["h_tanPhi_whenGoodX"+outCat]->Fill(-directionFw->at(bestTrigFw) - directionEm->at(iTrigEm));  
	  m_plots["h_time_whenGoodX"+outCat]->Fill( bxTimeFw->at(bestTrigFw) - bxTimeEm->at(iTrigEm) );  
        }
      } 
    } // if filling plots 

  } // loop over emul



} // loop over entries

TFile * file = new TFile ("outPlots.root", "RECREATE");


for (auto & category : qualityCategories) {

  m_plots2["h_tanPhi2D"+category] -> Write();
  m_plots2["h_pos2D"+category] -> Write();
  m_plots2["h_time2D"+category] -> Write();
  m_plots["h_time"+category] -> Write();
  m_plots["h_time_whenGoodX"+category] -> Write();
  m_plots["h_BX"+category] -> Write();
  m_plots["h_BXFW"+category] -> Write();
  m_plots["h_BXEmul"+category] -> Write();
  m_plots["h_tanPhi"+category] -> Write();
  m_plots["h_tanPhi_whenGoodX"+category] -> Write();
  m_plots["h_pos"+category] -> Write();
} 
for (auto & category : qualityNumbers) {

  m_plots2["h_tanPhi2D"+category] -> Write();
  m_plots2["h_pos2D"+category] -> Write();
  m_plots2["h_time2D"+category] -> Write();
  m_plots["h_time"+category] -> Write();
  m_plots["h_time_whenGoodX"+category] -> Write();
  m_plots["h_BX"+category] -> Write();
  m_plots["h_BXFW"+category] -> Write();
  m_plots["h_BXEmul"+category] -> Write();
  m_plots["h_tanPhi"+category] -> Write();
  m_plots["h_tanPhi_whenGoodX"+category] -> Write();
  m_plots["h_pos"+category] -> Write();
} 

delete file; 

/********************************************* PAINT AND SAVE HISTOS *********************************************************/


/*
TFile *file = new TFile (".//plots/histos.root","RECREATE");

gStyle->SetStatW(0.16);
gStyle->SetStatH(0.09);
gStyle->SetStatX(0.9);
gStyle->SetStatY(0.9);

//TFile *ff = new TFile ("./plots/2016/h_tanphi.root","RECREATE");
TCanvas *cTanphi = new TCanvas();
cTanphi->SetCanvasSize(1600, 1000);
cTanphi->SetWindowSize(1600, 1000);
h_tanphi->GetXaxis()->SetTitle("TanPhi firmware");
h_tanphi->GetYaxis()->SetTitle("TanPhi emulator");
h_tanphi->Draw("colz");
cTanphi->SaveAs("./plots/h_tanphi.png");
h_tanphi->Write();
//delete ff; 

//TFile *ffp = new TFile ("./plots/2016/h_pos.root","RECREATE");
TCanvas *cPos = new TCanvas();
cPos->SetCanvasSize(1600, 1000);
cPos->SetWindowSize(1600, 1000);
h_pos->GetXaxis()->SetTitle("Position firmware (cm)");
h_pos->GetYaxis()->SetTitle("Position emulator (cm)");
h_pos->Draw("colz");
cPos->SaveAs("./plots/h_pos.png");
h_pos->Write();
//delete ffp; 

//TFile *fft = new TFile ("./plots/2016/h_time.root","RECREATE");
TCanvas *cTime = new TCanvas();
cTime->SetCanvasSize(1600, 1000);
cTime->SetWindowSize(1600, 1000);
h_time->GetXaxis()->SetTitle("BxTime firmware (ns)");
h_time->GetYaxis()->SetTitle("BxTime emulator (ns)");
h_time->Draw("colz");
cTime->SaveAs("./plots/h_time.png");
h_time->Write();
//delete fft; 

//TFile *fft25 = new TFile ("./plots/2016/h_time25.root","RECREATE");
TCanvas *cTime25 = new TCanvas();
cTime25->SetCanvasSize(1600, 1000);
cTime25->SetWindowSize(1600, 1000);
h_time25->GetXaxis()->SetTitle("BxTime firmware (ns)");
h_time25->GetYaxis()->SetTitle("BxTime emulator (ns)");
h_time25->Draw("colztext");
cTime25->SaveAs("./plots/h_time25.png");
h_time25->Write();
//delete fft25; 

//TFile *fft1 = new TFile ("./plots/2016/h_t.root","RECREATE");
TCanvas *cT = new TCanvas();
cT->SetCanvasSize(1600, 1000);
cT->SetWindowSize(1600, 1000);
h_t->GetXaxis()->SetTitle("BxTime firmware - BxTime emulator (ns)");
h_t->GetYaxis()->SetTitle("Entries/ns");
h_t->Draw();
cT->SaveAs("./plots/h_t.png");
h_t->Write();
//delete fft1; 

TCanvas *cTP = new TCanvas();
cTP->SetCanvasSize(1600, 1000);
cTP->SetWindowSize(1600, 1000);
h_tP->GetXaxis()->SetTitle("tanPhi firmware - tanPhi emulator");
h_tP->GetYaxis()->SetTitle("Entries/0.01");
h_tP->Draw();
cTP->SaveAs("./plots/h_tP.png");
h_tP->Write();

TCanvas *cP = new TCanvas();
cP->SetCanvasSize(1600, 1000);
cP->SetWindowSize(1600, 1000);
h_p->GetXaxis()->SetTitle("Position firmware - Position emulator (cm)");
h_p->GetYaxis()->SetTitle("Entries");
h_p->Draw();
cP->SaveAs("./plots/h_p.png");
h_p->Write();

TCanvas *cq = new TCanvas();
cq->SetCanvasSize(1600, 1000);
cq->SetWindowSize(1600, 1000);
h_qFW->GetXaxis()->SetTitle("Qualities");
h_qFW->GetYaxis()->SetTitle("Entries");
h_qFW->Draw();
cq->SaveAs("./plots/h_qFW.png");
h_qFW->Write();

TCanvas *cqE = new TCanvas();
cqE->SetCanvasSize(1600, 1000);
cqE->SetWindowSize(1600, 1000);
h_qEmul->GetXaxis()->SetTitle("Qualities");
h_qEmul->GetYaxis()->SetTitle("Entries");
h_qEmul->Draw();
cqE->SaveAs("./plots/h_qEmul.png");
h_qEmul->Write();



gStyle->SetOptStat(111100);
gStyle->SetStatX(0.84);
gStyle->SetStatY(0.4);
//TFile *fft2 = new TFile ("./plots/2016/h_time1800.root","RECREATE");
TCanvas *cTime2 = new TCanvas();
cTime2->SetCanvasSize(1600, 1000);
cTime2->SetWindowSize(1600, 1000);
h_time2000->GetXaxis()->SetTitle("BxTime firmware mod 1800 (ns)");
h_time2000->GetYaxis()->SetTitle("BxTime emulator mod 1800 (ns)");
h_time2000->Draw("colztext");
cTime2->SaveAs("./plots/h_time1800.png");
h_time2000->Write();
//delete fft2; 

delete file; 


*/





}//end macro



