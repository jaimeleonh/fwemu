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

const int hitsToStudy = 3; 

const bool printHits = false; 
//const bool printHits = true; 


const int numberOfHitsTh = 12;


/*************************************************************************************
 *                                GLOBAL VARIABLES
 * **********************************************************************************/

std::map<std::string, TH1*> m_plots;
std::map<std::string, TH2*> m_plots2;
std::map<std::string, TEfficiency*> m_plotsEff;

struct hit
{
  short wheel; 
  short sector; 
  short station; 
  double shift; 
  int sl; 
  int la; 
  int cell; 
  int tdc;
};

/**************************************************************************************
 * 			                               NTUPLAS
**************************************************************************************/
// Cargo la ntupla de las primitivas
TFile *f = new TFile("./primitives.root");
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

// Cargo la ntupla del firmware
TFile *f1 = new TFile("./outputfw.root");
TTree *ntuple1 = (TTree*)f1->Get("ntuple");

UInt_t nTrigsFw = 0;
std::vector<double> *positionFw= 0;
std::vector<double> *directionFw= 0;
std::vector<int> *bxTimeFw= 0;
std::vector<int> *lastHitBXFw= 0;
std::vector<int> *arrivalBXFw= 0;
std::vector<short> *qualityFw= 0;
std::vector<double> *chi2Fw= 0;
std::vector<short> *wheelFw= 0;
std::vector<short> *sectorFw= 0;
std::vector<short> *stationFw= 0;
std::vector<std::vector<short>> *lateralitiesFw= 0;
std::vector<std::vector<int>> *wiresFw= 0;
std::vector<std::vector<int>> *tdcsFw= 0;


 // Cargo la ntupla del firmware
TFile *f2 = new TFile("./hits.root");
TTree *ntuple2 = (TTree*)f2->Get("ntuple");

UInt_t nHits = 0;
std::vector<int> *hitTdcs=0;
std::vector<int> *hitLayers=0;
std::vector<int> *hitCells=0;
std::vector<double> *hitShift=0;
std::vector<short> *hitWheels=0;
std::vector<short> *hitSectors=0;
std::vector<short> *hitStations=0;
std::vector<short> *hitSuperlayers=0;



/***************************************************************************************
 *
 *                    			           SUBPROGRAMAS
 *
 * *************************************************************************************/


int getFwSL (int iTrigFw) {
  if (wiresFw->at(iTrigFw)[0] == -1 && wiresFw->at(iTrigFw)[0] == -1 && wiresFw->at(iTrigFw)[0] == -1 && wiresFw->at(iTrigFw)[3] == -1 ) return 3; 
  else if (wiresFw->at(iTrigFw)[4] == -1 && wiresFw->at(iTrigFw)[5] == -1 && wiresFw->at(iTrigFw)[6] == -1 && wiresFw->at(iTrigFw)[7] == -1 ) return 1; 
  else return 0; 
}

int getEmSL (int iTrigEm) {
  if (wiresEm->at(iTrigEm)[0] == -1 && wiresEm->at(iTrigEm)[0] == -1 && wiresEm->at(iTrigEm)[0] == -1 && wiresEm->at(iTrigEm)[3] == -1 ) return 3; 
  else if (wiresEm->at(iTrigEm)[4] == -1 && wiresEm->at(iTrigEm)[5] == -1 && wiresEm->at(iTrigEm)[6] == -1 && wiresEm->at(iTrigEm)[7] == -1 ) return 1; 
  else return 0;
}

/***************************************** FILL LATERALITY PLOTS **************************************************/

std::vector <std::vector <int> > cellLayouts = {{0,-1,-2,-3},{0,-1,-2,-1},{0,-1,0,-1},{0,-1,0,1},{0,1,0,-1},{0,1,0,1},{0,1,2,1},{0,1,2,3}};

int getNonValidLayer (std::vector<bool> myValids) {
  for (int i = 0; i < 4; i++){
    if (myValids[i] == false ) return i + 1; 
  }
  return 0; 
}

bool compareLays (std::vector<int> refLayout, std::vector<int> myLayout, std::vector<bool> myValids) {
  for (int i = 0; i < 4; i++){
    if (myValids[i] && refLayout[i] != myLayout[i] ) return false; 
  } 
  return true; 
}

int getLatComb ( int iTrigEm, int sl, std::vector<bool> myValids ) {
  int sum4 = sl/2 * 4;
  int pot = 1;

  int myLat = 0; 
  for (int i = 0; i<4; i++){ 
    if (myValids[3-i]) { 
      myLat = myLat + lateralitiesEm->at(iTrigEm)[sum4+3-i]*pot;  
    } 
    pot = 2*pot; 
  }
  //cout << "myLat " << myLat << endl; 
  if (myLat < 0) cout << "ERROR iTrigEm " << iTrigEm << " myLat " << myLat << endl;  
  return myLat; 
}

void fillLats (int iTrigEm, int sl, std::string plot) {
 
  int sum4 = sl/2 * 4; 
  std::vector <bool> validLayout; 
  int nonValid; 


  if ( wiresEm->at(iTrigEm)[0+sum4] != -1 ) { // FILL 0***
    
    validLayout.push_back ( true );

    std::vector <int> cellLayout = {0};
    cellLayout.push_back ( 2*(wiresEm->at(iTrigEm)[1+sum4] - wiresEm->at(iTrigEm)[0+sum4]) - 1 ); validLayout.push_back ( wiresEm->at(iTrigEm)[1+sum4] != -1 );
    cellLayout.push_back ( 2*(wiresEm->at(iTrigEm)[2+sum4] - wiresEm->at(iTrigEm)[0+sum4])     ); validLayout.push_back ( wiresEm->at(iTrigEm)[2+sum4] != -1 );
    cellLayout.push_back ( 2*(wiresEm->at(iTrigEm)[3+sum4] - wiresEm->at(iTrigEm)[0+sum4]) - 1 ); validLayout.push_back ( wiresEm->at(iTrigEm)[3+sum4] != -1 );

    //cout << "CellLayout " << cellLayout[0] << " " << cellLayout[1] << " " << cellLayout[2] << " " << cellLayout[3] << " " << endl; 

    nonValid = getNonValidLayer(validLayout); 

    for (auto & layout : cellLayouts) {
      if (compareLays (layout, cellLayout, validLayout ) ) {
        //if (plot=="badLats" && layout[0] == 0 && layout[1] == -1 && layout[2] == 0 && layout[3] == -1) cout << "Aqui hay uno feo con" << getLatComb(iTrigEm, sl, validLayout) << endl; 
        if (plot == "badLats" && layout[0] == 0 && layout[1] == -1 && layout[2] == 0 && layout[3] == -1 && getLatComb(iTrigEm, sl, validLayout) == 4 && nonValid == 4 ) cout << "Aqui hay uno feo " << endl; 
        //if (plot=="badLats" && layout[0] == 0 && layout[1] == -1 && layout[2] == -2 && layout[3] == -3 && getLatComb(iTrigEm, sl, validLayout), nonValid == 5 ) cout << "Aqui hay uno feo " << endl; 
        m_plots2[plot+to_string(layout[0])+to_string(layout[1])+to_string(layout[2])+to_string(layout[3])] -> Fill (getLatComb(iTrigEm, sl, validLayout), nonValid );
      }
    }


  } else { //FILL X***
 
    nonValid = 1; 

    validLayout.push_back ( false );
    validLayout.push_back ( wiresEm->at(iTrigEm)[1+sum4] != -1 );
    validLayout.push_back ( wiresEm->at(iTrigEm)[2+sum4] != -1 );
    validLayout.push_back ( wiresEm->at(iTrigEm)[3+sum4] != -1 );


    // FILL X-1**
    std::vector <int> cellLayout1 = {-1};
    cellLayout1.push_back (    wiresEm->at(iTrigEm)[1+sum4] - wiresEm->at(iTrigEm)[1+sum4]  + 1     );
    cellLayout1.push_back ( 2*(wiresEm->at(iTrigEm)[2+sum4] - wiresEm->at(iTrigEm)[1+sum4]) + 2     );
    cellLayout1.push_back ( 2*(wiresEm->at(iTrigEm)[3+sum4] - wiresEm->at(iTrigEm)[1+sum4]) + 1     );
    //cout << "CellLayout1 " << cellLayout1[0] << " " << cellLayout1[1] << " " << cellLayout1[2] << " " << cellLayout1[3] << " " << endl; 
   
    for (auto & layout : cellLayouts) {
      if (compareLays (layout, cellLayout1, validLayout ) ) {
        m_plots2[plot+to_string(layout[0])+to_string(layout[1])+to_string(layout[2])+to_string(layout[3])] -> Fill (getLatComb(iTrigEm, sl, validLayout), nonValid );
      }
    }

    // FILL X1**
    std::vector <int> cellLayout2 = {-1};
    cellLayout2.push_back (    wiresEm->at(iTrigEm)[1+sum4] - wiresEm->at(iTrigEm)[1+sum4] - 1      );
    cellLayout2.push_back ( 2*(wiresEm->at(iTrigEm)[2+sum4] - wiresEm->at(iTrigEm)[1+sum4])         );
    cellLayout2.push_back ( 2*(wiresEm->at(iTrigEm)[3+sum4] - wiresEm->at(iTrigEm)[1+sum4]) - 1     );
    //cout << "CellLayout2 " << cellLayout2[0] << " " << cellLayout2[1] << " " << cellLayout2[2] << " " << cellLayout2[3] << " " << endl; 
  
    for (auto & layout : cellLayouts) {
      if (compareLays (layout, cellLayout2, validLayout ) ) {
    //    if ( layout[0] == 0 && layout[1] == -1 && layout[2] == 0 && layout[3] == -1) cout << "Aqui hay uno feo " << endl; 
        //if (layout[0] == 0 && layout[1] == -1 && layout[2] == 0 && layout[3] == -1 && getLatComb(iTrigEm, sl, validLayout), nonValid == 5 ) cout << "Aqui hay uno feo " << endl; 
        //if (plot == "badLats" && layout[0] == 0 && layout[1] == -1 && layout[2] == 0 && layout[3] == -1 && getLatComb(iTrigEm, sl, validLayout) == 5 && nonValid == 5 ) cout << "Aqui hay uno feo " << endl; 
        m_plots2[plot+to_string(layout[0])+to_string(layout[1])+to_string(layout[2])+to_string(layout[3])] -> Fill (getLatComb(iTrigEm, sl, validLayout), nonValid );
      }
    }
 
  }
}



/******************************************************************************************************************/


bool bigger (pair<short,size_t> i, pair<short, size_t> j) { return (i.first>j.first) ; } 


std::map < int, std::vector < pair < short, size_t>  > > orderPrimitives () {

  std::map < int,  std::vector < pair < short, size_t>  > > primitiveIndexes;

  for (std::size_t iTrigEm = 0; iTrigEm < nTrigsEm; iTrigEm++) {
    primitiveIndexes[ 100*sectorEm->at(iTrigEm) + 10*( wheelEm->at(iTrigEm) + 2 ) + stationEm->at(iTrigEm)].push_back(make_pair(qualityEm->at(iTrigEm),iTrigEm)); 
  }


  for (auto & par : primitiveIndexes) {
    std::sort ( par.second.begin(), par.second.end() , bigger);  
  }

  return primitiveIndexes; 

}

void seeFWPrims ( int wheel, int sector, int station, int &maxQuality, int &numFits, int *sl, bool *h4sl ) {
  maxQuality = -1; 
  numFits = 0; 
  sl[0] = 0; sl[1] = 0;   
  h4sl[0] = false; h4sl[1] = false;   

  for (std::size_t iTrigFw = 0; iTrigFw < nTrigsFw; iTrigFw++) {
    if (wheelFw->at(iTrigFw)!=wheel || sectorFw->at(iTrigFw)!=sector || stationFw->at(iTrigFw)!=station) continue;
    if (qualityFw->at(iTrigFw) > maxQuality ) maxQuality = qualityFw->at(iTrigFw);
    if (getFwSL(iTrigFw) != 0 ) {
      sl[getFwSL(iTrigFw) / 2] = iTrigFw; 
      if (qualityFw->at(iTrigFw) > 2) h4sl[getFwSL(iTrigFw) / 2] = true;
    }
    numFits++;
  }
}


void printFWPrims ( int wheel, int sector, int station ) {
  for (std::size_t iTrigFw = 0; iTrigFw < nTrigsFw; iTrigFw++) {
    if (wheelFw->at(iTrigFw)!=wheel || sectorFw->at(iTrigFw)!=sector || stationFw->at(iTrigFw)!=station) continue;
	  cout << "--------- Primitive in iTrigFW=" << iTrigFw << "-----------" << endl;
    for (int deb = 0; deb < 8; deb++)
        cout << wiresFw->at(iTrigFw)[deb] << " " << tdcsFw->at(iTrigFw)[deb] << " "<< lateralitiesFw->at(iTrigFw)[deb]  << endl;  
    cout << "ArrivalBX: "<< arrivalBXFw->at(iTrigFw) << endl; 
  }
}

void printJMHits (std::vector <hit> hits) {
  for (auto & hit : hits) 
    cout << hit.sl << " "
         << hit.la << " "
         << hit.cell << " " 
         << hit.tdc << endl;
}

void printFullHits (std::vector <hit> hits) {
  for (auto & hit : hits) 
    cout << hit.wheel << " " 
         << hit.sector << " " 
         << hit.station << " " 
         << hit.shift << " " 
         << hit.sl << " " 
         << hit.la << " " 
         << hit.cell << " " 
         << hit.tdc << endl;
}

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

void fillLatencies () {
  for (size_t iTrigFw = 0; iTrigFw < nTrigsFw; iTrigFw++){
    m_plots["h_Latencies"]->Fill(arrivalBXFw->at(iTrigFw) - lastHitBXFw->at(iTrigFw));
  }
}



/***************************************************************************************
 *										       
 *	               *****************PROGRAMA PRINCIPAL*******************		       
 *										       
 * *************************************************************************************/

void comparator_fw(){


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

ntuple1->SetBranchAddress("numberOfTrigs",&nTrigsFw);
ntuple1->SetBranchAddress("Position",&positionFw);
ntuple1->SetBranchAddress("Direction",&directionFw);
ntuple1->SetBranchAddress("Time",&bxTimeFw);
ntuple1->SetBranchAddress("lastHitBX",&lastHitBXFw);
ntuple1->SetBranchAddress("ArrivalBX",&arrivalBXFw);
ntuple1->SetBranchAddress("Quality",&qualityFw);
ntuple1->SetBranchAddress("chi2",&chi2Fw);
ntuple1->SetBranchAddress("Wheel",&wheelFw);
ntuple1->SetBranchAddress("Sector",&sectorFw);
ntuple1->SetBranchAddress("Station",&stationFw);
ntuple1->SetBranchAddress("Lateralities",&lateralitiesFw);
ntuple1->SetBranchAddress("Cellnumber",&wiresFw);
ntuple1->SetBranchAddress("tdcTime",&tdcsFw);

ntuple2->SetBranchAddress("numberOfHits",&nHits);
ntuple2->SetBranchAddress("hitTdc",&hitTdcs);
ntuple2->SetBranchAddress("hitLayer",&hitLayers);
ntuple2->SetBranchAddress("hitCell",&hitCells);
ntuple2->SetBranchAddress("hitShift",&hitShift);
ntuple2->SetBranchAddress("hitWheel",&hitWheels);
ntuple2->SetBranchAddress("hitStation",&hitStations);
ntuple2->SetBranchAddress("hitSector",&hitSectors);
ntuple2->SetBranchAddress("hitSuperlayer",&hitSuperlayers);

gStyle->SetOptStat(1111111);


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


std::vector <std::string> qualityCategories = {"All", "sameQuality", "emulBetterQuality", "fwBetterQuality","Q1", "Q2", "Q3", "Q4", "Q6", "Q8", "Q9", "3h", "4h"};  
//std::vector <std::string> qualityNumbers = {"Q1", "Q2", "Q3", "Q4", "Q6", "Q8", "Q9"};  

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

//LABELS
std::vector <std::string> fitLabels = {"Same fit", "Different fit"}; 
std::vector <std::string> eventLabels = {"Events with same fits", "Events with at least one different fit"}; 


m_plots["h_hitDistr"] = new TH1F ("h_hitDistr","; ; Entries", 38 , -2.5, 40.5);
m_plots["h_goodEvents"] = new TH1F ("h_goodEvents","; ; Entries",2 , -0.5, 1.5);
for (unsigned int i = 0; i < eventLabels.size(); i++){
  m_plots["h_goodEvents"]->GetXaxis()->SetBinLabel(i+1, eventLabels.at(i).c_str());
}


m_plots["h_Latencies"] = new TH1F ("h_Latencies","Distribution of latencies; BXs; Entries", 300 , 0, 300);



std::vector <std::string> q9Labels = {"Other correlated", "4h in SL1 and 4h in SL3", "Other" }; 
m_plots["h_badQ9"] = new TH1F ("h_badQ9","; ; Entries", 3, -0.5, 2.5);
for (unsigned int i = 0; i < q9Labels.size(); i++){
  m_plots["h_badQ9"]->GetXaxis()->SetBinLabel(i+1, q9Labels.at(i).c_str());
}

m_plots["h_badQ9_ArrBXDif"] = new TH1F ("h_badQ9_ArrBXDif","; #Delta BX ; Entries", 101, -0.5, 100.5 );


for (auto & category : qualityCategories) {


  m_plots["h_sameFit"+category] = new TH1F (("h_sameFit" + category).c_str(), "; ; Entries",2 , -0.5, 1.5);
  for (unsigned int i = 0; i < fitLabels.size(); i++){
    m_plots["h_sameFit"+category]->GetXaxis()->SetBinLabel(i+1, fitLabels.at(i).c_str());
  }

  m_plotsEff["hEfficiencyVsHits"+category] = new TEfficiency (("h_EfficiencyVsHits" + category).c_str(), ("EfficiencyVsHits " + category + "; Number of hits; ").c_str(),38, 2.5, 40.5);
  

  m_plots2["h_numOfHits"+category] = new TH2F (("h_numOfHits" + category).c_str(), "; Emulator number of Hits; Firmware number of Hits", 8 , 1.5, 9.5, 11 , -1.5, 9.5);

  m_plots["h_chi2"+category] = new TH1F (("h_chi2" + category).c_str(), "chi2 firmware - chi2 emulator; #Delta chi2 (cm2); Entries", nbinChi2, -chi2min, chi2min);
  m_plots["h_tanPhiDistEmul"+category] = new TH1F (("h_tanPhiDistEmul" + category).c_str(), "TanPhi emulator distribution; TanPhi (adim); Entries", nbinDistPsi, -psiDistMin, psiDistMin);
  m_plots["h_tanPhiDistFW"+category] = new TH1F (("h_tanPhiDistFW" + category).c_str(), "TanPhi firmware distribution; TanPhi (adim); Entries", nbinDistPsi, -psiDistMin, psiDistMin);
  m_plots["h_xDistEmul"+category] = new TH1F (("h_xDistEmul" + category).c_str(), "Position emulator distribution; Position (cm); Entries", nbinDistx, -xDistMin, xDistMin);
  m_plots["h_xDistFW"+category] = new TH1F (("h_xDistFW" + category).c_str(), "Position firmware distribution; Position (cm); Entries", nbinDistx, -xDistMin, xDistMin);
  m_plots["h_timeDistEmul"+category] = new TH1F (("h_timeDistEmul" + category).c_str(), "Time emulator distribution; Time (ns); Entries", nbinDistTime, -timeDistMin, timeDistMin);
  m_plots["h_timeDistFW"+category] = new TH1F (("h_timeDistFW" + category).c_str(), "Time firmware distribution; Time (ns); Entries", nbinDistTime, -timeDistMin, timeDistMin);
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
  m_plots["h_BXFW"+category] = new TH1F (("h_BXFW" + category).c_str(), "BX firmware - eventBX;  BX firmware - eventBX (adim.); Entries / BX", nbinBX, -bxmin , bxmin );
  m_plots["h_BXEmul"+category] = new TH1F (("h_BXEmul" + category).c_str(), "BX emulator - eventBX; BX emulator - eventBX (adim.); Entries / BX", nbinBX, -bxmin , bxmin );

}

std::vector <std::string> latCombs = {"0000","0001","0010","0011","0100","0101","0110","0111","1000","1001","1010","1011","1100","1101","1110","1111"};
std::vector <std::string> missing = {"4-hit", "Miss Lay 1", "Miss Lay 2", "Miss Lay 3", "Miss Lay 4"};

for (auto & layout : cellLayouts) {
  m_plots2["lat"+to_string(layout[0])+to_string(layout[1])+to_string(layout[2])+to_string(layout[3])] = new TH2F (  ("h_lat_" + to_string(layout[0])+to_string(layout[1])+to_string(layout[2])+to_string(layout[3])).c_str(), "", 16, -0.5, 15.5, 5, -0.5, 4.5 ) ;
  m_plots2["badLats"+to_string(layout[0])+to_string(layout[1])+to_string(layout[2])+to_string(layout[3])] = new TH2F (  ("h_bad_lat_" + to_string(layout[0])+to_string(layout[1])+to_string(layout[2])+to_string(layout[3])).c_str(), "", 16, -0.5, 15.5, 5, -0.5, 4.5 ) ;

  
  for (unsigned int i = 0; i < latCombs.size(); i++){
    m_plots2["lat"+to_string(layout[0])+to_string(layout[1])+to_string(layout[2])+to_string(layout[3])]->GetXaxis()->SetBinLabel(i+1, latCombs.at(i).c_str());
    m_plots2["badLats"+to_string(layout[0])+to_string(layout[1])+to_string(layout[2])+to_string(layout[3])]->GetXaxis()->SetBinLabel(i+1, latCombs.at(i).c_str());
  }
  for (unsigned int j = 0; j < missing.size(); j++){
    m_plots2["lat"+to_string(layout[0])+to_string(layout[1])+to_string(layout[2])+to_string(layout[3])]->GetYaxis()->SetBinLabel(j+1, missing.at(j).c_str());
    m_plots2["badLats"+to_string(layout[0])+to_string(layout[1])+to_string(layout[2])+to_string(layout[3])]->GetYaxis()->SetBinLabel(j+1, missing.at(j).c_str());
  }

}




/**************************************************************************************
 * 					BUCLE
 *************************************************************************************/ 					

//for (Int_t i = 0; i < 1; i++){
for (Int_t i = 0; i < nEntries; i++){
//  if (i!=9513) continue;
  cout << "**************************************"<< i <<"**********************************" << endl;  
  //if (i!=1805) continue;
  ntuple->GetEntry(i);
  ntuple1->GetEntry(i);
  ntuple2->GetEntry(i);


  fillLatencies();



  short oldStation=-1, oldWheel=-1, oldSector=-1;
  bool gotSameFit = true;   
  int numberOfHits = 0; 
  
  std::vector <hit> hitsInChamber;
  hitsInChamber.clear();


  bool newPrintHits = false;


  std::map < int,  std::vector < pair < short, size_t>  > > primitiveIndexes =  orderPrimitives();


  for (auto & chan : primitiveIndexes) {
    int vistaLaPrimera = false; 
    for (auto & pair : chan.second ) {
      size_t iTrigEm = pair.second; 

    //if ( round ( bxTimeEm -> at (iTrigEm) / 25. )  != eventBX   ) continue; 
    if ( vistaLaPrimera ) continue; 
    //vistaLaPrimera = true; 

  // for (std::size_t iTrigEm = 0; iTrigEm < nTrigsEm; iTrigEm++) {
    
    int bestTrigFw = -1; 
    int bestTimeFw = 99999;      
    int biggestNumOfHits = 0;    
    int totalNumOfHits = 0;    
    int biggestTrigFw = -1;    

    if ( getEmSL(iTrigEm) == 1 ) fillLats (iTrigEm, 1, "lat");
    else if ( getEmSL(iTrigEm) == 3 ) fillLats (iTrigEm, 3, "lat");
    else {
      fillLats (iTrigEm, 1, "lat");
      fillLats (iTrigEm, 3, "lat");
    }



    if (oldStation == -1) { 
      oldStation =  stationEm->at(iTrigEm);
      oldWheel =  wheelEm->at(iTrigEm);
      oldSector =  sectorEm->at(iTrigEm);
      hitsInChamber.clear();
      newPrintHits = printHits;

      for (std::size_t iTrigHits = 0; iTrigHits < nHits; iTrigHits++) { 
        if (oldStation == hitStations->at(iTrigHits) && oldWheel == hitWheels->at(iTrigHits) && oldSector == hitSectors->at(iTrigHits)   ) {
          numberOfHits++;
          hitsInChamber.push_back(hit{
              hitWheels->at(iTrigHits),  
              hitSectors->at(iTrigHits),  
              hitStations->at(iTrigHits), 
              hitShift->at(iTrigHits),  
              hitSuperlayers->at(iTrigHits) - 1,  
              hitLayers->at(iTrigHits), 
              hitCells->at(iTrigHits), 
              hitTdcs->at(iTrigHits) });  
        }
      } 
    //  cout << numberOfHits << endl; 

    } else if  ( wheelEm->at(iTrigEm) != oldWheel || stationEm->at(iTrigEm) != oldStation || (sectorEm->at(iTrigEm) != oldSector ) ) { 
      if (gotSameFit)  m_plots["h_goodEvents"]->Fill(0);
      else m_plots["h_goodEvents"]->Fill(1);
      gotSameFit = true; 
      oldStation =  stationEm->at(iTrigEm);
      oldWheel =  wheelEm->at(iTrigEm);
      oldSector =  sectorEm->at(iTrigEm);
      numberOfHits=0;
      hitsInChamber.clear();
      newPrintHits = printHits;
      for (std::size_t iTrigHits = 0; iTrigHits < nHits; iTrigHits++) { 
        if (oldStation == hitStations->at(iTrigHits) && oldWheel == hitWheels->at(iTrigHits) && oldSector == hitSectors->at(iTrigHits)   ) { 
          numberOfHits++;
          hitsInChamber.push_back(hit{
              hitWheels->at(iTrigHits),  
              hitSectors->at(iTrigHits),  
              hitStations->at(iTrigHits), 
              hitShift->at(iTrigHits),  
              hitSuperlayers->at(iTrigHits) - 1,  
              hitLayers->at(iTrigHits), 
              hitCells->at(iTrigHits), 
              hitTdcs->at(iTrigHits) });  
         // cout << "Insert hit " << endl;
        }
      }
    }

    m_plots["h_hitDistr"] -> Fill (numberOfHits);
    if (numberOfHits > numberOfHitsTh ) continue; 


    //cout << newPrintHits <<  " " << qualityEm->at(iTrigEm) << " " << hitsInChamber.size() << endl;
    if (newPrintHits && qualityEm->at(iTrigEm) >= 6) {
      printFullHits (hitsInChamber);
      newPrintHits = false;
    }




    for (std::size_t iTrigFw = 0; iTrigFw < nTrigsFw; iTrigFw++) {
      if ( (wheelEm->at(iTrigEm) == wheelFw->at(iTrigFw)) && (stationEm->at(iTrigEm) == stationFw->at(iTrigFw)) && (sectorEm->at(iTrigEm) == sectorFw->at(iTrigFw)) ) {
	int numberOfHits = sharedHits(wiresEm->at(iTrigEm), wiresFw->at(iTrigFw), tdcsEm->at(iTrigEm), tdcsFw->at(iTrigFw), lateralitiesEm->at(iTrigEm), lateralitiesFw->at(iTrigFw));
	if (biggestNumOfHits < numberOfHits) {
	  biggestNumOfHits= numberOfHits; //cout << biggestNumOfHits << endl;  
	  biggestTrigFw = iTrigFw; 
	} 
        if (!sameFit(wiresEm->at(iTrigEm), wiresFw->at(iTrigFw), tdcsEm->at(iTrigEm), tdcsFw->at(iTrigFw), lateralitiesEm->at(iTrigEm), lateralitiesFw->at(iTrigFw))) continue;
        int deltaT0 = abs (bxTimeEm->at(iTrigEm) - bxTimeFw->at(iTrigFw));
        if (deltaT0 <= bestTimeFw) {
          bestTrigFw = iTrigFw; 
          bestTimeFw = deltaT0; 
        }
      } // If same chamber
    } // loop over fw

/**************************************************************************************
 * 			    Select Quality categories
**************************************************************************************/
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
          if (qualityGroup(qualityEm->at(iTrigEm)) == 1 ) outCategories.push_back("3h");
          if (qualityGroup(qualityEm->at(iTrigEm)) == 2 ) outCategories.push_back("4h");
        } 	
      else if (qualityEm->at(iTrigEm) > qualityFw->at(bestTrigFw)) outCategories.push_back("emulBetterQuality"); 	
      else if (qualityEm->at(iTrigEm) < qualityFw->at(bestTrigFw)) outCategories.push_back("fwBetterQuality"); 

        if (debug) {
	  double myChi2;
	  if (qualityEm->at(iTrigEm)==6||qualityEm->at(iTrigEm)==8||qualityEm->at(iTrigEm)==9) myChi2 = chi2Fw->at(bestTrigFw)/(1024.*100);
	  else myChi2 = chi2Fw->at(bestTrigFw)/(1024.*100);
	  cout << "----------"<<i<<"----------" << endl;
	  cout << "Wh:" << wheelEm->at(iTrigEm) << " Se:" << sectorEm->at(iTrigEm) << " St:" << stationEm->at(iTrigEm) << endl; 
	  for (int deb = 0; deb < 8; deb++){
            cout << wiresEm->at(iTrigEm)[deb] << " " << tdcsEm->at(iTrigEm)[deb] << " "<< lateralitiesEm->at(iTrigEm)[deb]  << endl;  
	  } 
	    cout << "PosEmul "    << positionEm->at(iTrigEm) - shift->at(iTrigEm)  << " PosFW "  << positionFw->at(bestTrigFw) << endl;
	    cout << "TimeEmul "   << bxTimeEm->at(iTrigEm)    << " TimeFW " << bxTimeFw->at(bestTrigFw)                        << endl;
	    cout << "TanPsiEmul " << directionEm->at(iTrigEm) << " TanPsiFW "  << -directionFw->at(bestTrigFw)                 << endl;
	    cout << "chi2Emul " << chi2Em->at(iTrigEm) << " chi2FW "  << myChi2                                                << endl;
	} 
/**************************************************************************************
 * 			           Fill Histograms
**************************************************************************************/
      for (auto & outCat : outCategories) {
        m_plotsEff["hEfficiencyVsHits"+outCat] -> Fill(1,numberOfHits);

	double myChi2;
	if (qualityFw->at(bestTrigFw)==6||qualityFw->at(bestTrigFw)==8||qualityFw->at(bestTrigFw)==9) myChi2 = chi2Fw->at(bestTrigFw)/(1024.*100);
	else myChi2 = chi2Fw->at(bestTrigFw)/(1024.*100);
	//myChi2 = 4*chi2Fw->at(bestTrigFw)/(1024.*100);
        m_plots["h_sameFit"+outCat]->Fill(0);
        m_plots["h_tanPhiDistEmul"+outCat] -> Fill (directionEm->at(iTrigEm));
        m_plots["h_tanPhiDistFW"+outCat] -> Fill (-directionFw->at(bestTrigFw));
        m_plots["h_xDistEmul"+outCat] -> Fill (positionEm->at(iTrigEm));
        m_plots["h_xDistFW"+outCat] -> Fill (positionFw->at(bestTrigFw)+shift->at(iTrigEm));
        m_plots["h_timeDistFW"+outCat] -> Fill (bxTimeFw->at(bestTrigFw) - (int) (eventBX*25));
        m_plots["h_timeDistEmul"+outCat] -> Fill (bxTimeEm->at(iTrigEm) - (int) (eventBX*25));
	//cout << "bxTimeEm " << bxTimeEm->at(iTrigEm) << " EventBX*25 " <<eventBX*25 << " bxTimeEm->at(iTrigEm) - eventBX*25) " << bxTimeEm->at(iTrigEm) - (int)( eventBX*25 ) << endl; 

	m_plots["h_chi2"+outCat] -> Fill ( chi2Em->at(iTrigEm) - myChi2 );
      //  if ( fabs ( chi2Em->at(iTrigEm) - myChi2) > 1E-2) cout << chi2Em->at(iTrigEm) << " " << myChi2 << " " << i << endl;

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
        } else { 
//	  cout << i << endl; 
	}
        if (bxTimeFw->at(bestTrigFw) - bxTimeEm->at(iTrigEm) == 0 ){
	  m_plots["h_tanPhi_whenSameT0"+outCat]->Fill(-directionFw->at(bestTrigFw) - directionEm->at(iTrigEm));  
          m_plots["h_pos_whenSameT0"+outCat]->Fill(positionFw->at(bestTrigFw) + shift->at(iTrigEm) - positionEm->at(iTrigEm));
          if ( fabs(positionFw->at(bestTrigFw) + shift->at(iTrigEm) - positionEm->at(iTrigEm) ) >= xmin / ((double) nbinx )) {
//            cout << i << endl;  
          } else {
          }
        } else {
  //        cout << i << endl;  
	}
//        if ( fabs(directionFw->at(bestTrigFw) + directionEm->at(iTrigEm) ) >= 0.0001 && outCat == "Q9") cout << i << endl; 
        //if ( fabs(directionFw->at(bestTrigFw) + directionEm->at(iTrigEm) ) >= tanmin / ((double) nbinTan )) cout << i << endl; 
      } 
    } // if bestTrig!= -1
    else { 
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //                                                    Haven't found the same fit                                                         //
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      if ( getEmSL(iTrigEm) == 1 ) fillLats (iTrigEm, 1, "badLats");
      else if ( getEmSL(iTrigEm) == 3 ) fillLats (iTrigEm, 3, "badLats");
      else {
        fillLats (iTrigEm, 1, "badLats");
        fillLats (iTrigEm, 3, "badLats");
      }

      if (qualityEm->at(iTrigEm) == 9) {
        int numPrims, maxQual, sl[2];
        bool h4sl[2]; 
        seeFWPrims ( wheelEm->at(iTrigEm), sectorEm->at(iTrigEm), stationEm->at(iTrigEm), maxQual, numPrims, sl, h4sl );
        if ( maxQual > 4 ) m_plots["h_badQ9"]->Fill( 0 );
        else {
          if (numPrims == 2 && h4sl[0] && h4sl[1])  {
            m_plots["h_badQ9"]->Fill( 1 ); 
            m_plots["h_badQ9_ArrBXDif"]->Fill( abs(arrivalBXFw->at(sl[0]) - arrivalBXFw->at(sl[1])) );
            if ( abs(arrivalBXFw->at(sl[0]) - arrivalBXFw->at(sl[1])) < 20 )  cout << i << endl; 

          }
          else  m_plots["h_badQ9"]->Fill( 2 );
        }
      }

      gotSameFit = false; 

      //cout << biggestNumOfHits << endl; 
      std::vector <std::string> outCategories = {"All"};
      if (qualityEm->at(iTrigEm) == 1) outCategories.push_back("Q1");
      else if (qualityEm->at(iTrigEm) == 2) outCategories.push_back("Q2");
      else if (qualityEm->at(iTrigEm) == 3) outCategories.push_back("Q3");
      else if (qualityEm->at(iTrigEm) == 4) outCategories.push_back("Q4");
      else if (qualityEm->at(iTrigEm) == 6) outCategories.push_back("Q6");
      else if (qualityEm->at(iTrigEm) == 8) outCategories.push_back("Q8");
      else if (qualityEm->at(iTrigEm) == 9) outCategories.push_back("Q9");
      if (qualityGroup(qualityEm->at(iTrigEm)) == 1 ) outCategories.push_back("3h");
      if (qualityGroup(qualityEm->at(iTrigEm)) == 2 ) outCategories.push_back("4h");
      for (auto & outCat : outCategories) {
        m_plotsEff["hEfficiencyVsHits"+outCat] -> Fill(0,numberOfHits);
        m_plots["h_sameFit"+outCat]->Fill(1);
	m_plots2["h_numOfHits"+outCat]->Fill(numFromQuality(qualityEm->at(iTrigEm)),biggestNumOfHits);
      }
      if (numFromQuality(qualityEm->at(iTrigEm))==hitsToStudy && biggestNumOfHits!=0 && verLosMalos) {
        //if (numFromQuality(qualityFw->at(biggestTrigFw)) == 4) continue; 
        //if (numFromQuality(qualityEm->at(iTrigEm))==4 && biggestNumOfHits==4 &&  && biggestNumOfHits!=0 && verLosMalos) 
        //if (numFromQuality(qualityEm->at(iTrigEm))==4 && biggestNumOfHits == 3) 
        //if (false == true && numFromQuality(qualityEm->at(iTrigEm))==4 && biggestNumOfHits == 3) 
        //if (numFromQuality(qualityEm->at(iTrigEm)) == biggestNumOfHits && numFromQuality(qualityEm->at(iTrigEm))==4 && biggestNumOfHits == 3) 
	      cout << "----------- " << i <<" ---------" << endl;
	      cout << qualityEm->at(iTrigEm)<< " "<<numFromQuality(qualityEm->at(iTrigEm)) << " " << biggestNumOfHits << endl; 
	      cout << "Wh:" << wheelEm->at(iTrigEm) << " Se:" << sectorEm->at(iTrigEm) << " St:" << stationEm->at(iTrigEm) << endl; 
	      for (int deb = 0; deb < 8; deb++){
          cout << wiresEm->at(iTrigEm)[deb] << " " << tdcsEm->at(iTrigEm)[deb] << " "<< lateralitiesEm->at(iTrigEm)[deb]  << " " << wiresFw->at(biggestTrigFw)[deb] << " " << tdcsFw->at(biggestTrigFw)[deb] << " "<< lateralitiesFw->at(biggestTrigFw)[deb]  << endl;  
	      } 
	      cout << "PosEmul "    << positionEm->at(iTrigEm) - shift->at(iTrigEm)  << " PosFW "  << positionFw->at(biggestTrigFw) << endl;
	      cout << "TimeEmul "   << bxTimeEm->at(iTrigEm)    << " TimeFW " << bxTimeFw->at(biggestTrigFw)                        << endl;
	      cout << "TanPsiEmul " << directionEm->at(iTrigEm) << " TanPsiFW "  << -directionFw->at(biggestTrigFw)                 << endl;
	      cout << "chi2Emul " << chi2Em->at(iTrigEm) << " chi2FW "  << chi2Fw->at(biggestTrigFw)/(1024.*100)                          << endl;
        cout << "Hits in this chamber:" << endl;
        printFWPrims ( wheelEm->at(iTrigEm), sectorEm->at(iTrigEm), stationEm->at(iTrigEm) );
        if (true) printJMHits(hitsInChamber);
        if (true) printFullHits(hitsInChamber);
	    } else if (numFromQuality(qualityEm->at(iTrigEm))==hitsToStudy && biggestNumOfHits==0 && verLosMalos) {
	    //} else if (numFromQuality(qualityEm->at(iTrigEm))==8 && biggestNumOfHits==0 && verLosMalos && true==false) 
        //if (numFromQuality(qualityEm->at(iTrigEm))==4 && biggestNumOfHits == 3) 
        //if (false == true && numFromQuality(qualityEm->at(iTrigEm))==4 && biggestNumOfHits == 3) 
        //if (numFromQuality(qualityEm->at(iTrigEm)) == biggestNumOfHits && numFromQuality(qualityEm->at(iTrigEm))==4 && biggestNumOfHits == 3) 
	    cout << "----------- " << i <<" ---------" << endl;
	    //cout << qualityEm->at(iTrigEm)<< " "<<numFromQuality(qualityEm->at(iTrigEm)) << " " << biggestNumOfHits << endl; 
	    cout << "Wh:" << wheelEm->at(iTrigEm) << " Se:" << sectorEm->at(iTrigEm) << " St:" << stationEm->at(iTrigEm) << endl; 
	    for (int deb = 0; deb < 8; deb++){
        cout << wiresEm->at(iTrigEm)[deb] << " " << tdcsEm->at(iTrigEm)[deb] << " "<< lateralitiesEm->at(iTrigEm)[deb] << endl;  
	    } 
	    cout << "PosEmul "    << positionEm->at(iTrigEm) - shift->at(iTrigEm)  << " PosFW "  << -1 << endl;
	    cout << "TimeEmul "   << bxTimeEm->at(iTrigEm)    << " TimeFW " << -1                        << endl;
	    cout << "TanPsiEmul " << directionEm->at(iTrigEm) << " TanPsiFW "  << -1                 << endl;
	    cout << "chi2Emul " << chi2Em->at(iTrigEm) << " chi2FW "  << -1                         << endl;
      cout << "Hits in this chamber:" << endl;
      if (numberOfHits == 10) printJMHits(hitsInChamber);
      if (numberOfHits == 10) printFullHits(hitsInChamber);
  	}
//	break;  
  } 

  } // loop over emul


} 



  //for (int i = 0; i < 8; i++) cout << "-1 ";
  //cout << endl;

} // loop over entries

TFile * file = new TFile ("outPlots.root", "RECREATE");

/**************************************************************************************
 * 			           WRITE HISTOGRAMS
**************************************************************************************/
/*
file->cd();
file->Write();
file->Close();
*/


m_plots["h_goodEvents"] -> Write();
m_plots["h_Latencies"] -> Write();
m_plots["h_hitDistr"] -> Write();
m_plots["h_badQ9"] -> Write();
m_plots["h_badQ9_ArrBXDif"] -> Write();

for (auto & layout : cellLayouts) {
  m_plots2["lat"+to_string(layout[0])+to_string(layout[1])+to_string(layout[2])+to_string(layout[3])] -> Write();
  m_plots2["badLats"+to_string(layout[0])+to_string(layout[1])+to_string(layout[2])+to_string(layout[3])] -> Write();
}
  
  for (auto & category : qualityCategories) {

  m_plotsEff["hEfficiencyVsHits"+category] -> Write();
  m_plots2["h_numOfHits"+category] -> Write();
  m_plots["h_sameFit"+category] -> Write();
  m_plots["h_timeDistFW"+category] -> Write();
  m_plots["h_timeDistEmul"+category] -> Write();
  m_plots["h_tanPhiDistFW"+category] -> Write();
  m_plots["h_tanPhiDistEmul"+category] -> Write();
  m_plots["h_xDistFW"+category] -> Write();
  m_plots["h_xDistEmul"+category] -> Write();
  m_plots2["h_tanPhi2D"+category] -> Write();
  m_plots2["h_pos2D"+category] -> Write();
  m_plots2["h_time2D"+category] -> Write();
  m_plots["h_chi2"+category] -> Write();
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
} 

delete file; 



}//end macro



