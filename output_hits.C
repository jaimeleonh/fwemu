#include <stdio.h>
#include <stdlib.h>

void output_hits() {
/*       TString dir = gSystem->UnixPathName("./outputfw.C");
       dir.ReplaceAll("outputfw.C","");
       dir.ReplaceAll("/./","/");
*/
       UInt_t nHits = 0;
       std::vector<int> m_tdc;
       std::vector<int> m_layer;
       std::vector<int> m_cellNumber;
       std::vector<short> m_wheel;
       std::vector<short> m_sector;
       std::vector<short> m_station;
       std::vector<short> m_superlayer;
       
       m_tdc.clear();
       m_layer.clear();
       m_cellNumber.clear();
       m_wheel.clear();
       m_sector.clear();
       m_station.clear();
       m_superlayer.clear();

       TFile *hfile = new TFile("hits.root","RECREATE");
       TTree *m_tree = new TTree("ntuple", "Root tree");
       
       m_tree->Branch("numberOfHits",  &nHits);
       m_tree->Branch("hitTdc",  &m_tdc);
       m_tree->Branch("hitLayer",  &m_layer);
       m_tree->Branch("hitCell",  &m_cellNumber);
       m_tree->Branch("hitWheel",  &m_wheel);
       m_tree->Branch("hitSector",  &m_sector);
       m_tree->Branch("hitStation",  &m_station);
       m_tree->Branch("hitSuperlayer",  &m_superlayer);

       int tdc=0;
       int layer=0;
       int cellNumber=0;
       short wheel=0;
       short sector=0;
       short station=0;
       short superlayer=0;
       double dum=0.;

 
       //ifstream output("all6prims.txt");
       ifstream output("allHits_newPrecision.txt");
       //ifstream output("sl1Prims.txt");
       //ifstream output("allPrims_new.txt");
       if (output.is_open()){
         while(output>>wheel>>sector>>station>>dum>>superlayer>>layer>>cellNumber>>tdc) {
         //while(output>>quality>> position>>  direction>> time>> chi2>> shift>> wheel>> sector>>station>>wi1>>wi2>>wi3>>wi4>>wi5>>wi6>>wi7>>wi8>>tdc1>>tdc2>>tdc3>>tdc4>>tdc5>>tdc6>>tdc7>>tdc8>>lat1>>lat2>>lat3>>lat4>>lat5>>lat6>>lat7>>lat8>>BX) {
         if (sector != -1) {
	   nHits++;
           //cout << sector << endl; 
           m_tdc.push_back(tdc);
           m_layer.push_back(layer);
           m_cellNumber.push_back(cellNumber);
           m_wheel.push_back(wheel);
           m_sector.push_back(sector);
           m_station.push_back(station);
           m_superlayer.push_back(superlayer+1);

	 } else {
	   m_tree->Fill();
	   nHits=0;
           m_tdc.clear();
           m_layer.clear();
           m_cellNumber.clear();
           m_wheel.clear();
           m_sector.clear();
           m_station.clear();
           m_superlayer.clear();
	 }

         };//while
       };//if cablig
       //m_tree->Fill();
       m_tdc.clear();
       m_layer.clear();
       m_cellNumber.clear();
       m_wheel.clear();
       m_sector.clear();
       m_station.clear();
       m_superlayer.clear();
       
       hfile->Write(); 

/*
	TFile *f = new TFile ("outputfw.root","RECREATE");
	TTree *T = new TTree("ntuple", "data from ascii file");
	Long64_t nlines = T->ReadFile(Form("%soutputfw2.txt",dir.Data()),"event:position:tanphi:bxTime:quality:wheel:sector:station");
	T->Write();
*/





} // end void
