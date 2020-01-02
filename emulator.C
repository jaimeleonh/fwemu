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
#include <boost/iterator/zip_iterator.hpp>
#include <boost/range.hpp>
#ifdef __MAKECINT__
#pragma link C++ class vector<short>;
#pragma link C++ class vector<double>;
#endif
using namespace std;




/***************************************************************************************
 *
 * 				CONSTANTES INICIALES
 *
 ***************************************************************************************/

const bool debug = false; 
//const bool debug = true; 

const bool verLosMalos = false; 
//const bool verLosMalos = true; 

/***************************************************************************************
 *
 *			           SUBPROGRAMAS
 *
 * *************************************************************************************/

int numFromQuality (int quality){
  if (quality == 1) return 3; 
  if (quality == 2) return 3; 
  if (quality == 3) return 4; 
  if (quality == 4) return 4; 
  if (quality == 6) return 6; 
  if (quality == 8) return 7; 
  if (quality == 9) return 8;
  return -1;  
}

bool sameFit (std::vector <int> wires1, std::vector <int> wires2, std::vector <int> tdcs1, std::vector <int> tdcs2, std::vector <short> lateralities1, std::vector <short> lateralities2) {
  if ((wires1[0]!=wires2[0] || abs(tdcs1[0]-tdcs2[0])>0 || lateralities1[0]!=lateralities2[0]) && !(wires1[0]==-1 && wires2[0] == -1)) return false; 
  if ((wires1[1]!=wires2[1] || abs(tdcs1[1]-tdcs2[1])>0 || lateralities1[1]!=lateralities2[1]) && !(wires1[1]==-1 && wires2[1] == -1)) return false; 
  if ((wires1[2]!=wires2[2] || abs(tdcs1[2]-tdcs2[2])>0 || lateralities1[2]!=lateralities2[2]) && !(wires1[2]==-1 && wires2[2] == -1)) return false; 
  if ((wires1[3]!=wires2[3] || abs(tdcs1[3]-tdcs2[3])>0 || lateralities1[3]!=lateralities2[3]) && !(wires1[3]==-1 && wires2[3] == -1)) return false; 
  if ((wires1[4]!=wires2[4] || abs(tdcs1[4]-tdcs2[4])>0 || lateralities1[4]!=lateralities2[4]) && !(wires1[4]==-1 && wires2[4] == -1)) return false; 
  if ((wires1[5]!=wires2[5] || abs(tdcs1[5]-tdcs2[5])>0 || lateralities1[5]!=lateralities2[5]) && !(wires1[5]==-1 && wires2[5] == -1)) return false; 
  if ((wires1[6]!=wires2[6] || abs(tdcs1[6]-tdcs2[6])>0 || lateralities1[6]!=lateralities2[6]) && !(wires1[6]==-1 && wires2[6] == -1)) return false; 
  if ((wires1[7]!=wires2[7] || abs(tdcs1[7]-tdcs2[7])>0 || lateralities1[7]!=lateralities2[7]) && !(wires1[7]==-1 && wires2[7] == -1)) return false; 
  return true; 
}

int sharedHits (std::vector <int> wires1, std::vector <int> wires2, std::vector <int> tdcs1, std::vector <int> tdcs2, std::vector <short> lateralities1, std::vector <short> lateralities2) {
  int numberOfSharedHits = 0;   
  for (int i=0; i<8; i++){
    numberOfSharedHits += (wires1[i]==wires2[i] && tdcs1[i]==tdcs2[i] && lateralities1[i]==lateralities2[i] && wires1[i]!=-1 && wires2[i]!=-1);
//    cout << numberOfSharedHits << endl;  
  }
  return numberOfSharedHits; 
}

int qualityGroup (int quality) {
  if (quality == 1 || quality == 2 || quality == 5) return 1; 
  if (quality == 3 || quality == 4 || quality == 7) return 2; 
  if (quality == 6) return 3; 
  if (quality == 8) return 4; 
  if (quality == 9) return 5;
  return -1;  
}





/***************************************************************************************
 *										       *
 *		*****************PROGRAMA PRINCIPAL*******************		       *
 *										       *
 * *************************************************************************************/

void emulator(){



  gStyle->SetOptStat(1111111);

/**************************************************************************************
 * 			              NTUPLAS
**************************************************************************************/
// Cargo la ntupla de las primitivas
TFile *f = new TFile("./primitives_badchi2.root");
TTree *ntuple = (TTree*)f->Get("ntuple");

UInt_t nTrigsEm = 0;
UInt_t eventBX = 0; 
std::vector<double> *shift = 0;
std::vector<double> *positionEm= 0;
std::vector<double> *directionEm= 0;
std::vector<int> *bxTimeEm= 0;
std::vector<short> *qualityEm= 0;
std::vector<double> *chi2Em= 0;
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
ntuple->SetBranchAddress("chi2",&chi2Em);
ntuple->SetBranchAddress("Shift",&shift);
ntuple->SetBranchAddress("Wheel",&wheelEm);
ntuple->SetBranchAddress("Sector",&sectorEm);
ntuple->SetBranchAddress("Station",&stationEm);
ntuple->SetBranchAddress("Lateralities",&lateralitiesEm);
ntuple->SetBranchAddress("Cellnumber",&wiresEm);
ntuple->SetBranchAddress("tdcTime",&tdcsEm);

int n1 = ntuple->GetEntries();
cout << "Emulator segments: " << ntuple->GetEntries() << endl; 

/**************************************************************************************
 * 				Declaramos histogramas
**************************************************************************************/

std::map<std::string, TH1*> m_plots;
std::map<std::string, TH2*> m_plots2;

//std::vector <std::string> qualityCategories = {"All", "sameQuality", "emulBetterQuality", "fwBetterQuality","Q1", "Q2", "Q3", "Q4", "Q6", "Q8", "Q9"};  
std::vector <std::string> qualityCategories = {"Q1", "Q2", "Q3", "Q4", "Q6", "Q8", "Q9","Correlated"};  

//int nbinx = 21;    double xmin = 0.105;
int nbinx = 11;    double xmin = (1./160.)*(0.5+nbinx/2.);
int nbinTime = 21; double timemin = 10.5;
//int nbinTan = 21;  double tanmin = 0.0105;
int nbinTan = 11;  double tanmin = 2.5E-4 * (0.5 + nbinTan/2) ;
int nbinBX = 15; double bxmin = 7.5; 
int nbinChi2 = 15; double chi2min = 1E-5 * (0.5 + nbinChi2/2) ; 

//DISTRIBUTIONS
int nbinDistx = 1000; double xDistMin = 250;  
int nbinDistPsi = 1000; double psiDistMin = 2;  
int nbinDistTime = 300; double timeDistMin = 150;  
int nbinDistChi2 = 20000; double chi2Distmin = 1  ; 

//LABELS
std::vector <std::string> fitLabels = {"Same fit", "Different fit"}; 




for (auto & category : qualityCategories) {


  m_plots["h_chi2"+category] = new TH1F (("h_chi2" + category).c_str(), "Chi2 emulator; chi2 (cm2); Entries", nbinDistChi2,0,chi2Distmin);

}

/*
for (auto & category : qualityNumbers) {

  m_plots2["h_tanPhi2D"+category] = new TH2F (("h_tanPhi2D" + category).c_str(), "TanPhi firmware vs TanPhi emulator; TanPhi FW; TanPsi Emul", 400, -1,1, 400, -1, 1);
  m_plots2["h_pos2D"+category] = new TH2F (("h_pos2D" + category).c_str(), "Position firmware vs Position emulator; Position FW; Position Emul", 250, -250,250, 250, -250, 250);
  m_plots2["h_time2D"+category] = new TH2F (("h_time2D" + category).c_str(), "BxTime firmware vs BxTime emulator; BxTime FW; BxTime Emul", 1000, 0, 90000, 1000, 0, 90000);
  m_plots["h_time"+category] = new TH1F (("h_time" + category).c_str(), "BxTime firmware - BxTime emulator; #Delta BxTime (ns); Entries", nbinTime, -timemin, timemin);
  m_plots["h_time_whenGoodX"+category] = new TH1F (("h_time_whenGoodX" + category).c_str(), "BxTime firmware - BxTime emulator when goodX; #Delta BxTime (ns); Entries", nbinTime, -timemin, timemin);
  m_plots["h_BX"+category] = new TH1F (("h_BX" + category).c_str(), "BxTime firmware - BxTime emulator; #Delta BxTime (ns); Entries", 7, -87.5, 87.5);
  m_plots["h_tanPhi"+category] = new TH1F (("h_tanPhi" + category).c_str(), "TanPhi firmware - TanPhi emulator; #Delta TanPhi (adim); Entries", nbinTan, -tanmin, tanmin);
  m_plots["h_tanPhi_whenSameT0"+category] = new TH1F (("h_tanPhi_whenSameT0" + category).c_str(), "TanPhi firmware - TanPhi emulator; #Delta TanPhi (adim); Entries", nbinTan, -tanmin, tanmin);
  m_plots["h_tanPhi_whenGoodX"+category] = new TH1F (("h_tanPhi_whenGoodX" + category).c_str(), "TanPhi firmware - TanPhi emulator when goodX; #Delta TanPhi (adim); Entries", nbinTan, -tanmin, tanmin);
  m_plots["h_pos"+category] = new TH1F (("h_pos" + category).c_str(), "Position firmware - Position emulator; #Delta Position (cm); Entries / 0.01 cm", nbinx, -xmin , xmin );
  m_plots["h_pos_whenSameT0"+category] = new TH1F (("h_pos_whenSameT0" + category).c_str(), "Position firmware - Position emulator; #Delta Position (cm); Entries / 0.01 cm", nbinx, -xmin , xmin );
  m_plots["h_BXFW"+category] = new TH1F (("h_BXFW" + category).c_str(), "BX firmware - eventBX; #Delta BX (adim.); Entries / BX", nbinBX, -bxmin , bxmin );
  m_plots["h_BXEmul"+category] = new TH1F (("h_BXEmul" + category).c_str(), "BX emulator - eventBX; #Delta BX (adim.); Entries / BX", nbinBX, -bxmin , bxmin );
}
*/

/**************************************************************************************
 * 					BUCLE
 *************************************************************************************/ 					

//for (Int_t i = 0; i < 1000; i++){
for (Int_t i = 0; i < n1; i++){
  ntuple->GetEntry(i);

  for (std::size_t iTrigEm = 0; iTrigEm < nTrigsEm; iTrigEm++) {
    int bestTrigFw = -1; 
    int bestTimeFw = 99999;      
    int biggestNumOfHits = 0;    
    int totalNumOfHits = 0;    
    int biggestTrigFw = -1;    

    std::vector <std::string> outCategories = {};   
 
    if (qualityEm->at(iTrigEm) == 1) outCategories.push_back("Q1");
    else if (qualityEm->at(iTrigEm) == 2) outCategories.push_back("Q2");
    else if (qualityEm->at(iTrigEm) == 3) outCategories.push_back("Q3");
    else if (qualityEm->at(iTrigEm) == 4) outCategories.push_back("Q4");
    else if (qualityEm->at(iTrigEm) == 6) {outCategories.push_back("Q6"); outCategories.push_back("Correlated");}
    else if (qualityEm->at(iTrigEm) == 8) {outCategories.push_back("Q8"); outCategories.push_back("Correlated");}
    else if (qualityEm->at(iTrigEm) == 9) {outCategories.push_back("Q9"); outCategories.push_back("Correlated");}
    
    for (auto & outCat : outCategories) {
      m_plots["h_chi2"+outCat]->Fill(chi2Em->at(iTrigEm));  
    }
  } // loop over emul



} // loop over entries

TFile * file = new TFile ("outPlots_badchi2.root", "RECREATE");

/**************************************************************************************
 * 			           WRITE HISTOGRAMS
**************************************************************************************/

for (auto & category : qualityCategories) {

  m_plots["h_chi2"+category] -> Write();

} 

delete file; 



}//end macro



