#include <stdio.h>
#include <stdlib.h>

void output_firmware() {
/*       TString dir = gSystem->UnixPathName("./outputfw.C");
       dir.ReplaceAll("outputfw.C","");
       dir.ReplaceAll("/./","/");
*/
       UInt_t nTrigs = 0;


       std::vector<double> m_position;
       std::vector<double> m_direction;
       std::vector<int> m_time;
       std::vector<int> m_arrivalBX;
       std::vector<double> m_chi2;
       std::vector<short> m_quality;
       std::vector<short> m_superlayer;
       std::vector<short> m_wheel;
       std::vector<short> m_sector;
       std::vector<short> m_station;
       std::vector<std::vector<short>> m_lateralities;
       std::vector<std::vector<int>> m_wires;
       std::vector<std::vector<int>> m_tdcs;
       
       m_position.clear();
       m_direction.clear();
       m_time.clear();
       m_arrivalBX.clear();
       m_chi2.clear();
       m_quality.clear();
       m_superlayer.clear();
       m_wheel.clear();
       m_sector.clear();
       m_station.clear();
       m_lateralities.clear();
       m_wires.clear();
       m_tdcs.clear();

       TFile *hfile = new TFile("outputfw.root","RECREATE");
       TTree *m_tree = new TTree("ntuple", "Root tree");
       
       m_tree->Branch("numberOfTrigs",  &nTrigs);
       m_tree->Branch("Position",  &m_position);
       m_tree->Branch("Direction",  &m_direction);
       m_tree->Branch("Time",  &m_time);
       m_tree->Branch("ArrivalBX",  &m_arrivalBX);
       m_tree->Branch("chi2",  &m_chi2);
       m_tree->Branch("Superlayer",  &m_superlayer);
       m_tree->Branch("Quality",  &m_quality);
       m_tree->Branch("Wheel",  &m_wheel);
       m_tree->Branch("Sector",  &m_sector);
       m_tree->Branch("Station",  &m_station);
       m_tree->Branch("Lateralities",  &m_lateralities);
       m_tree->Branch("Cellnumber",  &m_wires);
       m_tree->Branch("tdcTime",  &m_tdcs);

       int index; 
       double position;
       double direction;
       int time; 
       int arrivalBX; 
       double chi2;
       short quality;
       short superlayer;
       short wheel;
       short sector;
       short station;
       short lat1, lat2, lat3, lat4, lat5, lat6, lat7, lat8;
       int wi1, wi2, wi3, wi4, wi5, wi6, wi7, wi8;
       int tdc1, tdc2, tdc3, tdc4, tdc5, tdc6, tdc7, tdc8;
       std::vector<short> lateralities;
       std::vector<int> cells;
       std::vector<int> tdcs;
       

       //ifstream output ("output_prim_181219.txt");
       ifstream output ("output_primitives_2802_sl1.txt");
       //ifstream output ("output_100120.txt");
       //ifstream output ("output_primitives_cor3001.txt");
       //ifstream output ("output_primitives_sl13001.txt");
       //ifstream output ("output_last_new.txt");
       //ifstream output ("output_2perBX.txt");
       //ifstream output ("output_per4.txt");
       if (output.is_open()){
         //  while(output>>index>>position>>direction>>time>>quality>>wheel>>sector>>station>>wi1>>wi2>>wi3>>wi4>>wi5>>wi6>>wi7>>wi8>>tdc1>>tdc2>>tdc3>>tdc4>>tdc5>>tdc6>>tdc7>>tdc8>>lat1>>lat2>>lat3>>lat4>>lat5>>lat6>>lat7>>lat8) {
         //while(output>>index>>position>>direction>>time>>quality>>chi2>>wheel>>sector>>station>>wi1>>wi2>>wi3>>wi4>>wi5>>wi6>>wi7>>wi8>>tdc1>>tdc2>>tdc3>>tdc4>>tdc5>>tdc6>>tdc7>>tdc8>>lat1>>lat2>>lat3>>lat4>>lat5>>lat6>>lat7>>lat8) { // Before 28/02
         while(output>>index>>position>>direction>>time>>arrivalBX>>quality>>chi2>>wheel>>sector>>station>>wi1>>wi2>>wi3>>wi4>>wi5>>wi6>>wi7>>wi8>>tdc1>>tdc2>>tdc3>>tdc4>>tdc5>>tdc6>>tdc7>>tdc8>>lat1>>lat2>>lat3>>lat4>>lat5>>lat6>>lat7>>lat8) { // After 28/02
         //while(!output.eof()){
         //output>>index>>position>>direction>>time>>quality>>wheel>>sector>>station;
         if (index != -1) {
	  // cout << index << " " << position << endl; 
	   nTrigs++;
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
           m_direction.push_back(direction);
           m_time.push_back(time);
           m_arrivalBX.push_back(arrivalBX);
           m_quality.push_back(quality);
           m_chi2.push_back(chi2);
           m_wheel.push_back(wheel);
           m_sector.push_back(sector);
           m_station.push_back(station);
           m_lateralities.push_back(lateralities);
           m_wires.push_back(cells);
           m_tdcs.push_back(tdcs);
	 } else {
	   m_tree->Fill();
	   nTrigs=0;
	   m_position.clear();
           m_direction.clear();
           m_chi2.clear();
           m_time.clear();
           m_arrivalBX.clear();
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
       m_quality.clear();
       m_chi2.clear();
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





}
