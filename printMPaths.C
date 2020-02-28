//#include "tdrstyle.C"
#include <TLatex.h>


void printMPaths () {

  const bool fileOK = false; 

  //setTDRStyle(); 
  gStyle->SetOptStat(111111);
  std::vector<std::string> stuffTags = { "Time", "Pos", "Psi"};

  TString file = "outPlots_mpaths.root";
  gSystem->Exec("mkdir newPlots");

  TFile data1(file);
  TFile outPlots("finalPlots_mpaths.root","RECREATE");

  char name [128];

  for (auto&& keyAsObj : *data1.GetListOfKeys()){
    auto key = (TKey*) keyAsObj;
    cout << key->GetName() << " " << key->GetClassName() << endl;
   
    if (std::string(key->GetClassName()) == "TH1F") {
      TH1F * plott = (TH1F*) data1.Get(std::string(key->GetName()).c_str());
      plott->Draw();
      std::string name = key->GetName();
      gPad->SaveAs(("newPlots/" + name + ".png").c_str());
      gPad->SaveAs(("newPlots/" + name + ".pdf").c_str());
    } else if (std::string(key->GetClassName()) == "TH2F") {
      TH2F * plott = (TH2F*) data1.Get(std::string(key->GetName()).c_str());
      plott->Draw("colz");
      std::string name = key->GetName();
      gPad->SaveAs(("newPlots/" + name + ".png").c_str());
      gPad->SaveAs(("newPlots/" + name + ".pdf").c_str());
    } else if  (std::string(key->GetClassName()) == "TEfficiency") {
      TEfficiency * plott = (TEfficiency*) data1.Get(std::string(key->GetName()).c_str());
      plott->Draw();
      std::string name = key->GetName();
      gPad->SaveAs(("newPlots/" + name + ".png").c_str());
      gPad->SaveAs(("newPlots/" + name + ".pdf").c_str());
    }
    
  }


  /*std::string nameH = "h_goodEvents";
  sprintf(name,"%s",nameH.c_str());
  m_plots[name] = (TH1F*) data1.Get(name);
  m_plots[name]->Scale (1./ (double) m_plots[name]->GetEntries());
  m_plots[name]->Draw("histo");
  outPlots.cd(); 
  m_plots[name]->Write();
  sprintf(name,"newPlots/%s.png", nameH.c_str());
  gPad->SaveAs(name);
  sprintf(name,"newPlots/%s.pdf", nameH.c_str());
  gPad->SaveAs(name);
  */
} // end macro
