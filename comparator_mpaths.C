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

const int numberOfHitsTh = 999999;

/***************************************************************************************
 *
 *			           SUBPROGRAMAS
 *
 * *************************************************************************************/

int numFromHits (std::vector <int> wires1){
  int numOfHits = (wires1[0]!=-1) + (wires1[1]!=-1) + (wires1[2]!=-1) + (wires1[3]!=-1);
  return numOfHits;
}

bool sameFit (std::vector <int> wires1, std::vector <int> wires2, std::vector <int> tdcs1, std::vector <int> tdcs2) {
  //cout << wires1[0] << " " << wires2[0] << "  " << tdcs1[0] << " " << tdcs2[0] << endl;
  //cout << wires1[1] << " " << wires2[1] << "  " << tdcs1[1] << " " << tdcs2[1] << endl;
  //cout << wires1[2] << " " << wires2[2] << "  " << tdcs1[2] << " " << tdcs2[2] << endl;
  //cout << wires1[3] << " " << wires2[3] << "  " << tdcs1[3] << " " << tdcs2[3] << endl << endl;

  if ((wires1[0]!=wires2[0] || abs(tdcs1[0]-tdcs2[0])>0)  && !(wires1[0]==-1 && wires2[0] == -1)) return false; 
  if ((wires1[1]!=wires2[1] || abs(tdcs1[1]-tdcs2[1])>0)  && !(wires1[1]==-1 && wires2[1] == -1)) return false; 
  if ((wires1[2]!=wires2[2] || abs(tdcs1[2]-tdcs2[2])>0)  && !(wires1[2]==-1 && wires2[2] == -1)) return false; 
  if ((wires1[3]!=wires2[3] || abs(tdcs1[3]-tdcs2[3])>0)  && !(wires1[3]==-1 && wires2[3] == -1)) return false; 
  return true; 
}

int sharedHits (std::vector <int> wires1, std::vector <int> wires2, std::vector <int> tdcs1, std::vector <int> tdcs2) {
  int numberOfSharedHits = 0;   
  for (int i=0; i<4; i++){
    numberOfSharedHits += (wires1[i]==wires2[i] && tdcs1[i]==tdcs2[i] && wires1[i]!=-1 && wires2[i]!=-1);
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

void comparator_mpaths(){



  gStyle->SetOptStat(1111111);

/**************************************************************************************
 * 			              NTUPLAS
**************************************************************************************/
// Cargo la ntupla de las primitivas
TFile *f = new TFile("./outputMPsEmul.root");
TTree *ntuple = (TTree*)f->Get("ntuple");

UInt_t nTrigsEm = 0;
std::vector<short> *wheelEm= 0;
std::vector<short> *sectorEm= 0;
std::vector<short> *stationEm= 0;
std::vector<short> *superlayerEm= 0;
std::vector<std::vector<int>> *wiresEm= 0;
std::vector<std::vector<int>> *tdcsEm= 0;

ntuple->SetBranchAddress("numberOfMPaths",  &nTrigsEm);
ntuple->SetBranchAddress("Wheel",&wheelEm);
ntuple->SetBranchAddress("Sector",&sectorEm);
ntuple->SetBranchAddress("Station",&stationEm);
ntuple->SetBranchAddress("Superlayer",&superlayerEm);
ntuple->SetBranchAddress("Cellnumber",&wiresEm);
ntuple->SetBranchAddress("tdcTime",&tdcsEm);

// Cargo la ntupla del firmware
TFile *f1 = new TFile("./outputMPsfw.root");
TTree *ntuple1 = (TTree*)f1->Get("ntuple");

UInt_t nTrigsFw = 0;
std::vector<short> *wheelFw= 0;
std::vector<short> *sectorFw= 0;
std::vector<short> *stationFw= 0;
std::vector<short> *superlayerFw= 0;
std::vector<std::vector<int>> *wiresFw= 0;
std::vector<std::vector<int>> *tdcsFw= 0;

ntuple1->SetBranchAddress("numberOfMPaths",&nTrigsFw);
ntuple1->SetBranchAddress("Wheel",&wheelFw);
ntuple1->SetBranchAddress("Sector",&sectorFw);
ntuple1->SetBranchAddress("Station",&stationFw);
ntuple1->SetBranchAddress("Superlayer",&superlayerFw);
ntuple1->SetBranchAddress("Cellnumber",&wiresFw);
ntuple1->SetBranchAddress("tdcTime",&tdcsFw);

 // Cargo la ntupla del firmware
TFile *f2 = new TFile("./hits.root");
TTree *ntuple2 = (TTree*)f2->Get("ntuple");

UInt_t nHits = 0;
std::vector<int> *hitTdcs=0;
std::vector<int> *hitLayers=0;
std::vector<int> *hitCells=0;
std::vector<short> *hitWheels=0;
std::vector<short> *hitSectors=0;
std::vector<short> *hitStations=0;
std::vector<short> *hitSuperlayers=0;

ntuple2->SetBranchAddress("numberOfHits",&nHits);
ntuple2->SetBranchAddress("hitTdc",&hitTdcs);
ntuple2->SetBranchAddress("hitLayer",&hitLayers);
ntuple2->SetBranchAddress("hitCell",&hitCells);
ntuple2->SetBranchAddress("hitWheel",&hitWheels);
ntuple2->SetBranchAddress("hitStation",&hitStations);
ntuple2->SetBranchAddress("hitSector",&hitSectors);
ntuple2->SetBranchAddress("hitSuperlayer",&hitSuperlayers);


int n1 = ntuple->GetEntries(); 
int n2 = ntuple1->GetEntries(); 
int n3 = ntuple2->GetEntries(); 

int nEntries = -1; 
if (n1 < n2) nEntries = n1; 
else nEntries = n2;  

cout << "Emulator segments: " << ntuple->GetEntries() << " AB7 outputs " << ntuple1->GetEntries() << " Hit events " << n3 << endl; 

/**************************************************************************************
 * 				Declaramos histogramas
**************************************************************************************/

std::map<std::string, TH1*> m_plots;
std::map<std::string, TH2*> m_plots2;
std::map<std::string, TEfficiency*> m_plotsEff;

std::vector <std::string> qualityCategories = {"3h", "4h", "All", "All3h"};  
//LABELS
std::vector <std::string> fitLabels = {"Same fit", "Different fit"}; 
std::vector <std::string> eventLabels = {"Events with same fits", "Events with at least one different fit"}; 


m_plots["h_goodEvents"] = new TH1F ("h_goodEvents","; ; Entries",2 , -0.5, 1.5);
for (unsigned int i = 0; i < eventLabels.size(); i++){
  m_plots["h_goodEvents"]->GetXaxis()->SetBinLabel(i+1, eventLabels.at(i).c_str());
}

for (auto & category : qualityCategories) {


  m_plots["h_sameFit"+category] = new TH1F (("h_sameFit" + category).c_str(), "; ; Entries",2 , -0.5, 1.5);
  for (unsigned int i = 0; i < fitLabels.size(); i++){
    m_plots["h_sameFit"+category]->GetXaxis()->SetBinLabel(i+1, fitLabels.at(i).c_str());
  }

  m_plotsEff["hMixerEfficiencyVsHits"+category] = new TEfficiency (("h_MixerEfficiencyVsHits" + category).c_str(), ("MixerEfficiency vs Hits " + category +  "; Number of hits; ").c_str(),38, 2.5, 40.5);
 
}

  m_plotsEff["hMixerEfficiency3hIn4h"] = new TEfficiency ("h_MixerEfficiency3hIn4h", "MixerEfficiency 3hIn4h vs Hits; Number of hits; " ,38, 2.5, 40.5);

/**************************************************************************************
 * 					BUCLE
 *************************************************************************************/ 					

//for (Int_t i = 1; i < 2; i++){
//for (Int_t i = 9994; i < 9995; i++){
for (Int_t i = 0; i < nEntries; i++){
//cout << "************************************************************************" << endl;  
  //if (i!=96) continue;
  //if (i!=277) continue;
  ntuple->GetEntry(i);
  ntuple1->GetEntry(i);
  ntuple2->GetEntry(i);

  short oldStation=-1, oldWheel=-1, oldSector=-1, oldSuperlayer=-1;
  bool gotSameFit = true;   
  int numberOfHits = 0; 

  
  for (std::size_t iTrigEm = 0; iTrigEm < nTrigsEm; iTrigEm++) {
    int bestTrigFw = -1; 
    int bestTimeFw = 99999;      
    int biggestNumOfHits = 0;    
    int totalNumOfHits = 0;    
    int biggestTrigFw = -1;    
    int mostHitsTrigFw = -1;    

    if (superlayerEm->at(iTrigEm) == 3) continue;

    if (oldStation == -1) { 
      oldStation =  stationEm->at(iTrigEm);
      oldWheel =  wheelEm->at(iTrigEm);
      oldSector =  sectorEm->at(iTrigEm);
      oldSuperlayer =  superlayerEm->at(iTrigEm);

      for (std::size_t iTrigHits = 0; iTrigHits < nHits; iTrigHits++) { 
        if (oldStation == hitStations->at(iTrigHits) && oldWheel == hitWheels->at(iTrigHits) && oldSector == hitSectors->at(iTrigHits) && oldSuperlayer == hitSuperlayers->at(iTrigHits)   ) 
          numberOfHits++;
      } 
    //  cout << numberOfHits << endl; 

    } else if  ( wheelEm->at(iTrigEm) != oldWheel || stationEm->at(iTrigEm) != oldStation || sectorEm->at(iTrigEm) != oldSector || superlayerEm->at(iTrigEm) != oldSuperlayer ) { 
      if (gotSameFit)  m_plots["h_goodEvents"]->Fill(0);
      else m_plots["h_goodEvents"]->Fill(1);
      gotSameFit = true; 
      oldStation =  stationEm->at(iTrigEm);
      oldWheel =  wheelEm->at(iTrigEm);
      oldSector =  sectorEm->at(iTrigEm);
      oldSuperlayer =  superlayerEm->at(iTrigEm);
      numberOfHits=0;
      for (std::size_t iTrigHits = 0; iTrigHits < nHits; iTrigHits++) { 
        if (oldStation == hitStations->at(iTrigHits) && oldWheel == hitWheels->at(iTrigHits) && oldSector == hitSectors->at(iTrigHits) && oldSuperlayer == hitSuperlayers->at(iTrigHits)   ) 
          numberOfHits++;
      }
    }

    if (numberOfHits > numberOfHitsTh ) continue; 

    for (std::size_t iTrigFw = 0; iTrigFw < nTrigsFw; iTrigFw++) {
      if ( (wheelEm->at(iTrigEm) == wheelFw->at(iTrigFw)) && (stationEm->at(iTrigEm) == stationFw->at(iTrigFw)) && (sectorEm->at(iTrigEm) == sectorFw->at(iTrigFw)) ) {
        if (sameFit(wiresEm->at(iTrigEm), wiresFw->at(iTrigFw), tdcsEm->at(iTrigEm), tdcsFw->at(iTrigFw))) {
          bestTrigFw = iTrigFw; break;
        } else {
          if (sharedHits(wiresEm->at(iTrigEm), wiresFw->at(iTrigFw), tdcsEm->at(iTrigEm), tdcsFw->at(iTrigFw)) != numFromHits(wiresEm->at(iTrigEm))) continue;
          else {
            mostHitsTrigFw = iTrigFw; 
            if (debug) {
              std::vector <int> wires1 = wiresEm->at(iTrigEm);
              std::vector <int> wires2 = wiresFw->at(iTrigFw);
              std::vector <int> tdcs1 = tdcsEm->at(iTrigEm);
              std::vector <int> tdcs2 = tdcsFw->at(iTrigFw);
              cout << wires1[0] << " " << wires2[0] << "  " << tdcs1[0] << " " << tdcs2[0] << endl;
              cout << wires1[1] << " " << wires2[1] << "  " << tdcs1[1] << " " << tdcs2[1] << endl;
              cout << wires1[2] << " " << wires2[2] << "  " << tdcs1[2] << " " << tdcs2[2] << endl;
              cout << wires1[3] << " " << wires2[3] << "  " << tdcs1[3] << " " << tdcs2[3] << endl << endl;
            }
          }
        }
      }
    } // loop over fw

/**************************************************************************************
 * 			    Select Quality categories
**************************************************************************************/
    // Filling plots
    std::vector <std::string> qualityCategories = {"All", "sameQuality", "emulBetterQuality", "fwBetterQuality"};  
    std::vector <std::string> outCategories = {"All"};
    if (numFromHits(wiresEm->at(iTrigEm)) == 4 ) outCategories.push_back("4h");
    else outCategories.push_back("3h");  
    if (bestTrigFw != -1) {

/**************************************************************************************
* 			           Fill Histograms
**************************************************************************************/
      for (auto & outCat : outCategories) {
        m_plotsEff["hMixerEfficiencyVsHits"+outCat] -> Fill(1,numberOfHits);
        m_plots["h_sameFit"+outCat]->Fill(0);
        if (outCat == "3h") m_plotsEff["hMixerEfficiencyVsHitsAll3h"] -> Fill(1,numberOfHits);
      } 
    } // if bestTrig!= -1
    else { 
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //                                                    Haven't found the same fit                                                         //
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      if (mostHitsTrigFw != -1) {
        m_plotsEff["hMixerEfficiency3hIn4h"] -> Fill(1,numberOfHits);          
         m_plotsEff["hMixerEfficiencyVsHitsAll3h"] -> Fill(1,numberOfHits);
      } else {
        m_plotsEff["hMixerEfficiency3hIn4h"] -> Fill(0,numberOfHits);  
        m_plotsEff["hMixerEfficiencyVsHitsAll3h"] -> Fill(0,numberOfHits);        
      }
      for (auto & outCat : outCategories) {
        if (numberOfHits < 5 && verLosMalos) { 
          cout <<"************** " << i << " **************" << endl; 
              std::vector <int> wires1 = wiresEm->at(iTrigEm);
              std::vector <int> tdcs1 = tdcsEm->at(iTrigEm);
              cout << wires1[0] << " " << tdcs1[0] << endl;
              cout << wires1[1] << " " << tdcs1[1] << endl;
              cout << wires1[2] << " " << tdcs1[2] << endl;
              cout << wires1[3] << " " << tdcs1[3] << endl << endl;
        }
        m_plotsEff["hMixerEfficiencyVsHits"+outCat] -> Fill(0,numberOfHits);
        m_plots["h_sameFit"+outCat]->Fill(1);
      }
	}
//	break;  
     

  } // loop over emul



} // loop over entries

TFile * file = new TFile ("outPlots_mpaths.root", "RECREATE");

/**************************************************************************************
 * 			           WRITE HISTOGRAMS
**************************************************************************************/

m_plots["h_goodEvents"] -> Write();
m_plotsEff["hMixerEfficiency3hIn4h"] -> Write();

for (auto & category : qualityCategories) {

  m_plotsEff["hMixerEfficiencyVsHits"+category] -> Write();
  m_plots["h_sameFit"+category] -> Write();
}/* 
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
  m_plots["h_tanPhi_whenSameT0"+category] -> Write();
  m_plots["h_pos"+category] -> Write();
  m_plots["h_pos_whenSameT0"+category] -> Write();
} */

delete file; 



}//end macro



