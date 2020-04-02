#include <stdio.h>
#include <stdlib.h>

void output_mpaths_firmware() {
/*       TString dir = gSystem->UnixPathName("./outputfw.C");
       dir.ReplaceAll("outputfw.C","");
       dir.ReplaceAll("/./","/");
*/



       ifstream output ("output_mixer_2303.txt");
       //ifstream output ("output_mixer_181219.txt");
       TFile *hfile = new TFile("outputMPsfw.root","RECREATE");








       UInt_t nMPaths = 0;
       std::vector<short> m_superlayer;
       std::vector<short> m_wheel;
       std::vector<short> m_sector;
       std::vector<short> m_station;
       std::vector<int> m_arrivalBX;
       std::vector<int> m_lastHitBX;
       std::vector<std::vector<int>> m_wires;
       std::vector<std::vector<int>> m_tdcs;
       
       m_superlayer.clear();
       m_wheel.clear();
       m_sector.clear();
       m_station.clear();
       m_arrivalBX.clear();
       m_lastHitBX.clear();
       m_wires.clear();
       m_tdcs.clear();

       TTree *m_tree = new TTree("ntuple", "Root tree");
       
       m_tree->Branch("numberOfMPaths",  &nMPaths);
       m_tree->Branch("Cellnumber",  &m_wires);
       m_tree->Branch("tdcTime",  &m_tdcs);
       m_tree->Branch("Superlayer",  &m_superlayer);
       m_tree->Branch("Wheel",  &m_wheel);
       m_tree->Branch("Sector",  &m_sector);
       m_tree->Branch("Station",  &m_station);
       m_tree->Branch("arrivalBX",  &m_arrivalBX);
       m_tree->Branch("lastHitBX",  &m_lastHitBX);

       int wi1, wi2, wi3, wi4, wi5, wi6, wi7, wi8;
       int tdc1, tdc2, tdc3, tdc4, tdc5, tdc6, tdc7, tdc8;
       std::vector<int> cells;
       std::vector<int> tdcs;
       short superlayer;
       short wheel;
       short sector;
       short station;
       int arrivalBX;
       int lastHitBX;
       

       if (output.is_open()){
         while(output>>wheel>>sector>>station>>superlayer>>arrivalBX>>lastHitBX>>wi1>>tdc1>>wi2>>tdc2>>wi3>>tdc3>>wi4>>tdc4) {
         //while(!output.eof()){
         //output>>index>>position>>direction>>time>>quality>>wheel>>sector>>station;
         if (wi1 != -1 || wi2 != -1) {
	       nMPaths++;
	   
           cells.clear();
	       cells.push_back(wi1);
	       cells.push_back(wi2);
           cells.push_back(wi3);
	       cells.push_back(wi4);

	       tdcs.clear();
           tdcs.push_back(tdc1);
           tdcs.push_back(tdc2);
           tdcs.push_back(tdc3);
           tdcs.push_back(tdc4);
           
           m_wheel.push_back(wheel);
           m_sector.push_back(sector);
           m_station.push_back(station);
           m_superlayer.push_back(superlayer);
           m_arrivalBX.push_back(arrivalBX);
           m_lastHitBX.push_back(lastHitBX);
           m_wires.push_back(cells);
           m_tdcs.push_back(tdcs);
	      } else {
	       m_tree->Fill();
	       nMPaths=0;
           m_wheel.clear();
           m_sector.clear();
           m_station.clear();
           m_superlayer.clear();
           m_arrivalBX.clear();
           m_lastHitBX.clear();
           m_wires.clear();
           m_tdcs.clear();
	 }

         };//while
       };//if cablig
       //m_tree->Fill();
       m_wires.clear();
       m_tdcs.clear();
       m_wheel.clear();
       m_sector.clear();
       m_station.clear();
       m_superlayer.clear();
       m_arrivalBX.clear();
       m_lastHitBX.clear();
       
       hfile->Write(); 

/*
	TFile *f = new TFile ("outputfw.root","RECREATE");
	TTree *T = new TTree("ntuple", "data from ascii file");
	Long64_t nlines = T->ReadFile(Form("%soutputfw2.txt",dir.Data()),"event:position:tanphi:bxTime:quality:wheel:sector:station");
	T->Write();
*/





}
