//#include "tdrstyle.C"
#include <TLatex.h>


void printPlots () {

  const bool fileOK = false; 

  //setTDRStyle(); 
  gStyle->SetOptStat(111111);
  std::vector<std::string> stuffTags = { "Time", "Pos", "Psi"};

  TString file = "outPlots.root";
  gSystem->Exec("mkdir newPlots");
  gSystem->Exec("rm newPlots/*");

  TFile data1(file);
  TFile outPlots("finalPlots.root","RECREATE");

  int a = 0;


  TLatex latex;
  latex.SetTextSize(0.03);

 

  std::vector <std::string> qualityCategories = {"All", "sameQuality", "emulBetterQuality", "fwBetterQuality"};  
  std::vector <std::string> qualityNumbers = {"Q1", "Q2", "Q3", "Q4", "Q6", "Q8", "Q9", "3h", "4h"};  
  std::vector <std::string> plots = {"h_chi2","h_sameFit","h_time", "h_time_whenGoodX","h_BX","h_BXFW","h_BXEmul","h_tanPhi","h_tanPhi_whenSameT0","h_tanPhi_whenGoodX","h_pos","h_pos_whenSameT0",};
  std::vector <std::string> plots2D = {"h_numOfHits","h_tanPhi2D","h_pos2D","h_time2D"};
  std::vector <std::string> plotsEff = {"h_EfficiencyVsHits"};
  std::vector <std::string> specialPlots = {"h_BX","h_xDist","h_tanPhiDist","h_timeDist"};
  std::map<std::string, TH1*> m_plots;
  std::map<std::string, TH2*> m_plots2;
  std::map<std::string, TEfficiency*> m_plotsEff;

  char name [128];

  if (a == 0) {
  std::string nameH = "h_goodEvents";
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
  }

  if (a==0){
  std::string nameH = "h_hitDistr";
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
  }

  for (const auto & category : qualityNumbers) {
    for (auto & plot : plotsEff) {
        std::string nameHisto = plot + category;
        sprintf(name,"%s",nameHisto.c_str());
        m_plotsEff[name] = (TEfficiency*) data1.Get(name);
        m_plotsEff[name]->Draw();
	outPlots.cd(); 
        m_plotsEff[name]->Write();
        sprintf(name,"newPlots/%s.png", nameHisto.c_str());
        gPad->SaveAs(name);
        sprintf(name,"newPlots/%s.pdf", nameHisto.c_str());
        gPad->SaveAs(name);
    }
  } 

  for (const auto & category : qualityCategories) {
      
      //gStyle->SetOptStat(111111);
      for (auto & plot : plotsEff) {
        std::string nameHisto = plot + category;
        sprintf(name,"%s",nameHisto.c_str());
        m_plotsEff[name] = (TEfficiency*) data1.Get(name);
        m_plotsEff[name]->Draw();
	outPlots.cd(); 
        m_plotsEff[name]->Write();
        sprintf(name,"newPlots/%s.png", nameHisto.c_str());
        gPad->SaveAs(name);
        sprintf(name,"newPlots/%s.pdf", nameHisto.c_str());
        gPad->SaveAs(name);
      }
      for (auto & plot : plots) {
        std::string nameHisto = plot + category;
        sprintf(name,"%s",nameHisto.c_str());
        m_plots[name] = (TH1F*) data1.Get(name);
        m_plots[name]->Scale (1./ (double) m_plots[name]->GetEntries());
        m_plots[name]->Draw("histo");
	outPlots.cd(); 
        m_plots[name]->Write();
        sprintf(name,"newPlots/%s.png", nameHisto.c_str());
	if (plot == "h_sameFit") gPad->SetLogy(0);
	else gPad->SetLogy();
        gPad->SaveAs(name);
        sprintf(name,"newPlots/%s.pdf", nameHisto.c_str());
	if (plot == "h_sameFit") gPad->SetLogy(0);
	else gPad->SetLogy();
        gPad->SaveAs(name);
	//outPlots.cd(); 
        //m_plots[name]->Write();
        //gPad->SetName(name);
        //gPad->Write();
        //m_plots[name]->Write();
      }
      for (auto & plot : plots2D) {
        std::string nameHisto = plot + category;
        sprintf(name,"%s",nameHisto.c_str());
        m_plots2[name] = (TH2F*) data1.Get(name);
        if (plot == "h_numOfHits") m_plots2[name]->Draw("colztext");
        else m_plots2[name]->Draw("colz");
	outPlots.cd(); 
        m_plots2[name]->Write();
	gPad->SetLogy(0);
        sprintf(name,"newPlots/%s.png", nameHisto.c_str());
        gPad->SaveAs(name);
        sprintf(name,"newPlots/%s.pdf", nameHisto.c_str());
        gPad->SaveAs(name);
	//outPlots.cd(); 
        //m_plots2[name]->Write();
        //gPad->SetName(name);
        //gPad->Write();
      }
      
  } // chamber
  for (const auto & category : qualityCategories) {
      for (const auto & plot : specialPlots) {
      gStyle->SetOptStat(0);

      std::string nameHist1 =  plot + "FW" + category;
      //std::string nameHist1 = "h_BXFW" + category;
      sprintf(name,"%s",nameHist1.c_str());
      TH1F * E1 = (TH1F*) data1.Get(name);
      E1->Scale (1./ (double) E1->GetEntries());
      

      nameHist1 = plot + "Emul" + category;
      //nameHist1 = "h_BXEmul" + category;
      sprintf(name,"%s",nameHist1.c_str());
      TH1F * E2 = (TH1F*) data1.Get(name);
      E2->Scale (1./ (double) E2->GetEntries());
      
      E1->SetMarkerStyle(20);
      E2->SetMarkerStyle(21);

      E1->SetMarkerColor(2);
      E2->SetMarkerColor(3);
      
      TLegend *leg = new TLegend(0.7,0.7,0.9,0.9);
      leg->AddEntry(E1,"Firmware");
      leg->AddEntry(E2,"Emulator");
      
      if (plot == "h_BX") {
        E1->SetTitle("BX FW/Emulator - EventBX");
        E1->GetXaxis()->SetTitle("BX FW/Emulator - EventBX");
      } else if (plot == "h_tanPsiDist") {
        E1->SetTitle("TanPsi FW/Emulator");
        E1->GetXaxis()->SetTitle("TanPsi FW/Emulator");
      } else if (plot == "h_xDist") {
        E1->SetTitle("Position FW/Emulator");
        E1->GetXaxis()->SetTitle("Position FW/Emulator (cm)");
      } else if (plot == "h_timeDist") {
        E1->SetTitle("Time FW/Emulator - L1A");
        E1->GetXaxis()->SetTitle("Time FW/Emulator - L1A");
      }

      E1->Draw("");
      E2->Draw("same");
      leg->Draw();

      gPad->SetLogy();
      sprintf(name,"newPlots/%sEvent%s.png", plot.c_str(), category.c_str());
      gPad->SaveAs(name);
      sprintf(name,"newPlots/%sEvent%s.pdf", plot.c_str(), category.c_str());
      gPad->SaveAs(name);
	outPlots.cd(); 
        gPad->SetName(name);
        gPad->Write();
    }
  }
/*
  for (auto & plot : plots) {
     
        sprintf(name,"%sQ1",plot.c_str());
        TH1F* E1 = (TH1F*) data1.Get(name);
        sprintf(name,"%sQ2",plot.c_str());
        TH1F* E2 = (TH1F*) data1.Get(name);
        sprintf(name,"%sQ3",plot.c_str());
        TH1F* E3 = (TH1F*) data1.Get(name);
        sprintf(name,"%sQ4",plot.c_str());
        TH1F* E4 = (TH1F*) data1.Get(name);
        sprintf(name,"%sQ6",plot.c_str());
        TH1F* E5 = (TH1F*) data1.Get(name);
        sprintf(name,"%sQ8",plot.c_str());
        TH1F* E6 = (TH1F*) data1.Get(name);
        sprintf(name,"%sQ9",plot.c_str());
        TH1F* E7 = (TH1F*) data1.Get(name);

        E1->Scale(1./(double) E1->GetEntries());
        E2->Scale(1./(double) E2->GetEntries());
        E3->Scale(1./(double) E3->GetEntries());
        E4->Scale(1./(double) E4->GetEntries());
        E5->Scale(1./(double) E5->GetEntries());
        E6->Scale(1./(double) E6->GetEntries());
        E7->Scale(1./(double) E7->GetEntries());



        E1->SetLineColor(kBlue);
        E2->SetLineColor(kRed);
        E3->SetLineColor(kGreen);
        E4->SetLineColor(kOrange);
        E5->SetLineColor(kBlack);
        E6->SetLineColor(kGray);
        E7->SetLineColor(kMagenta);
	
        TLegend *leg = new TLegend(0.6,0.6,0.80,0.8);
        leg->AddEntry(E1,"Q1","l");
        leg->AddEntry(E2,"Q2","l");
        leg->AddEntry(E3,"Q3","l");
        leg->AddEntry(E4,"Q4","l");
        leg->AddEntry(E5,"Q6","l");
        leg->AddEntry(E6,"Q8","l");
        leg->AddEntry(E7,"Q9","l");

        //E1->Draw("histo");
        //E2->Draw("histosame");
        //E3->Draw("histosame");
        //E4->Draw("histosame");
        //E5->Draw("histosame");
        //E6->Draw("histosame");
        //E7->Draw("histosame");
        //E1->Draw();
        E2->Draw("same");
        E3->Draw("same");
        E4->Draw("same");
        E5->Draw("same");
        E6->Draw("same");
        E7->Draw("same");
        leg->Draw();
     

     sprintf(name,"newPlots/%s.png", (plot + "_perQuality").c_str());
     gPad->SetLogy();
     gPad->SaveAs(name);
  }
*/      

  for (auto & plot : plots) {
      gStyle->SetOptStat(0);
     
        sprintf(name,"%sQ1",plot.c_str());
        TH1F* E1 = (TH1F*) data1.Get(name);
        sprintf(name,"%sQ2",plot.c_str());
        TH1F* E2 = (TH1F*) data1.Get(name);
        sprintf(name,"%sQ3",plot.c_str());
        TH1F* E3 = (TH1F*) data1.Get(name);
        sprintf(name,"%sQ4",plot.c_str());
        TH1F* E4 = (TH1F*) data1.Get(name);
        sprintf(name,"%sQ6",plot.c_str());
        TH1F* E5 = (TH1F*) data1.Get(name);
        sprintf(name,"%sQ8",plot.c_str());
        TH1F* E6 = (TH1F*) data1.Get(name);
        sprintf(name,"%sQ9",plot.c_str());
        TH1F* E7 = (TH1F*) data1.Get(name);

        E1->Add(E2);
        E3->Add(E4);
  
        std::string nameFile;  

        int entriesQ9 = E7->GetEntries();
        int entriesQ8 = E6->GetEntries();
        int entriesQ6 = E5->GetEntries();
        int entries4h = E3->GetEntries();
        int entries3h = E1->GetEntries();

        if (E1->GetEntries() != 0) E1->Scale(1./(double) E1->GetEntries());
        if (E3->GetEntries() != 0) E3->Scale(1./(double) E3->GetEntries());
        if (E5->GetEntries() != 0) E5->Scale(1./(double) E5->GetEntries());
        if (E6->GetEntries() != 0) E6->Scale(1./(double) E6->GetEntries());
        if (E7->GetEntries() != 0) E7->Scale(1./(double) E7->GetEntries());


        E1->SetFillColorAlpha(2,0.4);
        E3->SetFillColorAlpha(3,0.4);
        E5->SetFillColorAlpha(4,0.4);
        E6->SetFillColorAlpha(95,0.4);
        E7->SetFillColorAlpha(6,0.2);
   
        
        E1->SetMarkerStyle(20);
        E3->SetMarkerStyle(21);
        E5->SetMarkerStyle(22);
        E6->SetMarkerStyle(23);
        E7->SetMarkerStyle(47);


        E1->SetMarkerColor(2);
        E3->SetMarkerColor(3);
        E5->SetMarkerColor(4);
        E6->SetMarkerColor(95);
        E7->SetMarkerColor(6);

        E1->SetLineColor(2);
        E3->SetLineColor(3);
        E5->SetLineColor(4);
        E6->SetLineColor(95);
        E7->SetLineColor(6);
        /*E1->SetLineColor(kBlue);
        E3->SetLineColor(kGreen);
        E5->SetLineColor(kBlack);
        E6->SetLineColor(kOrange);
        E7->SetLineColor(kRed);*/
  
        //E1->SetFillStyle(3305);
        //E3->SetFillStyle(3144);
        //E5->SetFillStyle(3490);
        //E6->SetFillStyle(3244);
        //E7->SetFillStyle(3353);
        //E1->SetFillStyle(3444);
        //E3->SetFillStyle(3444);
        //E5->SetFillStyle(3444);
        //E6->SetFillStyle(3444);
        //E7->SetFillStyle(3444);
	

	
        TLegend *leg = new TLegend(0.7,0.7,0.9,0.9);
        leg->AddEntry(E1,"3h","fp");
        leg->AddEntry(E3,"4h","fp");
        leg->AddEntry(E5,"Q6","fp");
        leg->AddEntry(E6,"Q8","fp");
        leg->AddEntry(E7,"Q9","fp");

        if (plot == "h_sameFit") E7->GetYaxis()->SetRangeUser(0,1);
	else E7->GetYaxis()->SetRangeUser(0.0001,30);


        E7->Draw("bar");
        E6->Draw("barsame");
        E3->Draw("barsame");
        E5->Draw("barsame");
        E1->Draw("barsame");
        leg->Draw();
     
        double binx = E1->GetBinCenter (1) + E1->GetBinCenter (E1->GetNbinsX()) / 10; 
        nameFile = "Q9: "+std::to_string(entriesQ9);
        sprintf(name,"%s",nameFile.c_str());
        latex.DrawLatex(0,16,name);
        nameFile = "Q8: "+std::to_string(entriesQ8);
        sprintf(name,"%s",nameFile.c_str());
        latex.DrawLatex(0,10,name);
        nameFile = "Q6: "+std::to_string(entriesQ6);
        sprintf(name,"%s",nameFile.c_str());
        latex.DrawLatex(0,6,name);
        nameFile = "4h: "+std::to_string(entries4h);
        sprintf(name,"%s",nameFile.c_str());
        latex.DrawLatex(0,4,name);
        nameFile = "3h: "+std::to_string(entries3h);
        sprintf(name,"%s",nameFile.c_str());
        latex.DrawLatex(0,2.5,name);

     if (plot == "h_sameFit") gPad->SetLogy(0);
     else gPad->SetLogy();
     sprintf(name,"newPlots/%s.png", (plot + "_perQuality").c_str());
     gPad->SaveAs(name);
     sprintf(name,"newPlots/%s.pdf", (plot + "_perQuality").c_str());
     gPad->SaveAs(name);
     outPlots.cd(); 
     gPad->SetName(name);
     gPad->Write();
  }
} // end macro
