#include <stdio.h>
#include <stdlib.h>

void output_primitives() {
/*       TString dir = gSystem->UnixPathName("./outputfw.C");
       dir.ReplaceAll("outputfw.C","");
       dir.ReplaceAll("/./","/");
*/
       UInt_t nTrigs = 0;
       UInt_t eventBX = 0;
       std::vector<double> m_shift;
       std::vector<double> m_position;
       std::vector<double> m_direction;
       std::vector<double> m_time;
       std::vector<double> m_chi2;
       std::vector<short> m_quality;
       std::vector<short> m_wheel;
       std::vector<short> m_sector;
       std::vector<short> m_station;
       std::vector<std::vector<short>> m_lateralities;
       std::vector<std::vector<int>> m_wires;
       std::vector<std::vector<int>> m_tdcs;
       
       m_position.clear();
       m_direction.clear();
       m_time.clear();
       m_chi2.clear();
       m_quality.clear();
       m_wheel.clear();
       m_sector.clear();
       m_station.clear();
       m_lateralities.clear();
       m_wires.clear();
       m_tdcs.clear();

       TFile *hfile = new TFile("primitives.root","RECREATE");
       TTree *m_tree = new TTree("ntuple", "Root tree");
       
       m_tree->Branch("eventBX",  &eventBX);
       m_tree->Branch("numberOfTrigs",  &nTrigs);
       m_tree->Branch("Position",  &m_position);
       m_tree->Branch("Direction",  &m_direction);
       m_tree->Branch("chi2",  &m_chi2);
       m_tree->Branch("Time",  &m_time);
       m_tree->Branch("Quality",  &m_quality);
       m_tree->Branch("Shift",  &m_shift);
       m_tree->Branch("Wheel",  &m_wheel);
       m_tree->Branch("Sector",  &m_sector);
       m_tree->Branch("Station",  &m_station);
       m_tree->Branch("Lateralities",  &m_lateralities);
       m_tree->Branch("Cellnumber",  &m_wires);
       m_tree->Branch("tdcTime",  &m_tdcs);

       int BX; 
       double position;
       double shift;
       double direction;
       double time; 
       double chi2; 
       short quality;
       short wheel;
       short sector;
       short station;
       short lat1, lat2, lat3, lat4, lat5, lat6, lat7, lat8;
       int wi1, wi2, wi3, wi4, wi5, wi6, wi7, wi8;
       int tdc1, tdc2, tdc3, tdc4, tdc5, tdc6, tdc7, tdc8;
       std::vector<short> lateralities;
       std::vector<int> cells;
       std::vector<int> tdcs;
       
       //ifstream output("all6prims.txt");
       //ifstream output("all4Prims_primos.txt");
       //ifstream output("slprims.txt");
       //ifstream output("corPrims_newPrimos.txt");
       ifstream output("sl1Prims_newPrimos.txt");
       //ifstream output("allPrims_new.txt");
       if (output.is_open()){
         while(output>>quality>> position>>  direction>> time>> chi2>> shift>> wheel>> sector>>station>>wi1>>wi2>>wi3>>wi4>>wi5>>wi6>>wi7>>wi8>>tdc1>>tdc2>>tdc3>>tdc4>>tdc5>>tdc6>>tdc7>>tdc8>>lat1>>lat2>>lat3>>lat4>>lat5>>lat6>>lat7>>lat8>>BX) {
         if (quality != -1) {
	   nTrigs++;
/*
           if (lat1 == -1) lat1 = 0; 
           if (lat2 == -1) lat2 = 0; 
           if (lat3 == -1) lat3 = 0; 
           if (lat4 == -1) lat4 = 0; 
           if (lat5 == -1) lat5 = 0; 
           if (lat6 == -1) lat6 = 0; 
           if (lat7 == -1) lat7 = 0; 
           if (lat8 == -1) lat8 = 0; 

           if (wi1 == -1) wi1 = 0; 
           if (wi2 == -1) wi2 = 0; 
           if (wi3 == -1) wi3 = 0; 
           if (wi4 == -1) wi4 = 0; 
           if (wi5 == -1) wi5 = 0; 
           if (wi6 == -1) wi6 = 0; 
           if (wi7 == -1) wi7 = 0; 
           if (wi8 == -1) wi8 = 0; 
 	   
           if (tdc1 == -1) tdc1 = 0; 
           if (tdc2 == -1) tdc2 = 0; 
           if (tdc3 == -1) tdc3 = 0; 
           if (tdc4 == -1) tdc4 = 0; 
           if (tdc5 == -1) tdc5 = 0; 
           if (tdc6 == -1) tdc6 = 0; 
           if (tdc7 == -1) tdc7 = 0; 
           if (tdc8 == -1) tdc8 = 0; 
*/           
           lateralities.clear();
           lateralities.push_back(lat1);
           lateralities.push_back(lat2);
           lateralities.push_back(lat3);
           lateralities.push_back(lat4);
           lateralities.push_back(lat5);
           lateralities.push_back(lat6);
           lateralities.push_back(lat7);
           lateralities.push_back(lat8);
	  
           cells.clear(); 
	   cells.push_back(wi1);
	   cells.push_back(wi2);
	   cells.push_back(wi3);
	   cells.push_back(wi4);
	   cells.push_back(wi5);
	   cells.push_back(wi6);
	   cells.push_back(wi7);
	   cells.push_back(wi8);

	   tdcs.clear();
           tdcs.push_back(tdc1);
           tdcs.push_back(tdc2);
           tdcs.push_back(tdc3);
           tdcs.push_back(tdc4);
           tdcs.push_back(tdc5);
           tdcs.push_back(tdc6);
           tdcs.push_back(tdc7);
           tdcs.push_back(tdc8);
           
           m_position.push_back(position);
	   m_shift.push_back(shift);
           m_direction.push_back(direction);
           m_chi2.push_back(chi2);
           m_time.push_back(time);
           m_quality.push_back(quality);
           m_wheel.push_back(wheel);
           m_sector.push_back(sector);
           m_station.push_back(station);
           m_lateralities.push_back(lateralities);
           m_wires.push_back(cells);
           m_tdcs.push_back(tdcs);
	   eventBX = BX; 

	 } else {
	   m_tree->Fill();
	   nTrigs=0;
	   eventBX=0;
	   m_position.clear();
           m_direction.clear();
           m_time.clear();
           m_chi2.clear();
           m_shift.clear();
           m_quality.clear();
           m_wheel.clear();
           m_sector.clear();
           m_station.clear();
           m_lateralities.clear();
           m_wires.clear();
           m_tdcs.clear();
	 }

         };//while
       };//if cablig
       //m_tree->Fill();
       m_position.clear();
       m_direction.clear();
       m_time.clear();
       m_chi2.clear();
       m_shift.clear();
       m_quality.clear();
       m_wheel.clear();
       m_sector.clear();
       m_station.clear();
       m_lateralities.clear();
       m_wires.clear();
       m_tdcs.clear();
       
       hfile->Write(); 

/*
	TFile *f = new TFile ("outputfw.root","RECREATE");
	TTree *T = new TTree("ntuple", "data from ascii file");
	Long64_t nlines = T->ReadFile(Form("%soutputfw2.txt",dir.Data()),"event:position:tanphi:bxTime:quality:wheel:sector:station");
	T->Write();
*/





} // end void
