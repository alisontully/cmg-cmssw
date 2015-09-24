#include <iostream>
#include <fstream>

#include "TMultiGraph.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TFitResultPtr.h"
#include "TCanvas.h"
#include "TString.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TGraph.h"
#include "TStyle.h"

int main(int argc, char *argv[])
{
    std::string lepName = "electron";
    std::string detector_location;
    
    int nCuts = 50;
    int nCuts1 = 20;
    
    double scale_fac = 1;
    
    /*************************** Input fit results ****************************/
    
    TMultiGraph *mg = new TMultiGraph();
    TMultiGraph *mg1 = new TMultiGraph();
    TMultiGraph *mg2 = new TMultiGraph();
    TMultiGraph *mg3 = new TMultiGraph();
    TMultiGraph *mg_alpha = new TMultiGraph();
    TMultiGraph *mg_mean = new TMultiGraph();
    
    // Create vectors
    Double_t sigma_CB[nCuts], sigma_CB_err[nCuts], sigma_evt_by_evt[nCuts], sigma_evt_by_evt_err_l[nCuts], sigma_evt_by_evt_err_h[nCuts], sigma_eff[nCuts], sigma_eff_err[nCuts];
    Double_t sigma_eff_CB[nCuts], sigma_eff_CB_err[nCuts], sigma_eff_check[nCuts];
    Double_t alpha_CB[nCuts], alpha_CB_err[nCuts];
    Double_t mean_CB[nCuts], mean_CB_err[nCuts];
    
    // Create calibration lines
    Double_t line[nCuts], line_u[nCuts], line_l[nCuts];
    
    // Input data results
    std::ifstream input(TString("fitresult_Data_"+lepName+".txt"));
    for (int i=0; i<nCuts1; i++) {
        
        input >> sigma_CB[i];
        input >> sigma_CB_err[i];
        input >> sigma_evt_by_evt[i];
        input >> sigma_evt_by_evt_err_l[i];
        input >> sigma_evt_by_evt_err_h[i];
        input >> sigma_eff[i];
        input >> sigma_eff_err[i];
        input >> sigma_eff_check[i];
        input >> alpha_CB[i];
        input >> alpha_CB_err[i];
        input >> mean_CB[i];
        input >> mean_CB_err[i];
        
        sigma_evt_by_evt[i] = sigma_evt_by_evt[i]*scale_fac;
        sigma_evt_by_evt_err_l[i] = sigma_evt_by_evt_err_l[i]*scale_fac;
        sigma_evt_by_evt_err_h[i] = sigma_evt_by_evt_err_h[i]*scale_fac;
        
        sigma_eff_CB[i] = 1.0*sigma_eff_check[i]/sigma_CB[i];
        sigma_eff_CB_err[i] = sqrt(pow(sigma_eff_err[i]/sigma_CB[i], 2) + pow(sigma_eff_check[i]*sigma_CB_err[i]/(sigma_CB[i]*sigma_CB[i]), 2));
        
        line[i] = sigma_evt_by_evt[i];
        line_u[i] = 1.2*sigma_evt_by_evt[i];
        line_l[i] = 0.8*sigma_evt_by_evt[i];
        
//        std::cout << "sigma_evt_by_evt\t" << i << "\t" << sigma_evt_by_evt[i] << std::endl;
//        std::cout << "sigma_CB\t" << i << "\t" << sigma_CB[i] << "\t±\t" << sigma_CB_err[i] << std::endl;
//        std::cout << "sigma_eff\t" << i << "\t" << sigma_eff[i] << "\t±\t" << sigma_eff_err[i] << std::endl;
//        std::cout << "sigma_eff check\t" << i << "\t" << sigma_eff_check[i] << std::endl;
//        std::cout << "sigma_eff_CB\t" << i << "\t" << sigma_eff_CB[i] << "\t±\t" << sigma_eff_CB_err[i] << std::endl;
        
    }
    
    line[0] = 0;
    line_u[0] = 0;
    line_l[0] = 0;
    
    line[nCuts1-1] = 10.;
    line_u[nCuts1-1] = 1.2*10.;
    line_l[nCuts1-1] = 0.8*10.;
    
    input.close();
    
    
    
    
    // Create vectors for MC inputs
    Double_t sigma_CB1[nCuts], sigma_CB_err1[nCuts], sigma_evt_by_evt1[nCuts], sigma_evt_by_evt_err_l1[nCuts], sigma_evt_by_evt_err_h1[nCuts];
    Double_t sigma_eff_CB1[nCuts], sigma_eff_CB_err1[nCuts], sigma_eff1[nCuts], sigma_eff_err1[nCuts], sigma_eff_check1[nCuts];
    Double_t alpha_CB1[nCuts], alpha_CB_err1[nCuts];
    Double_t mean_CB1[nCuts], mean_CB_err1[nCuts];
    
    Double_t sigma_CB2[nCuts], sigma_CB_err2[nCuts], sigma_eff2[nCuts], sigma_eff_err2[nCuts];
    Double_t sigma_eff_CB2[nCuts], sigma_eff_CB_err2[nCuts], sigma_eff_check2[nCuts];
    Double_t alpha_CB2[nCuts], alpha_CB_err2[nCuts];
    Double_t mean_CB2[nCuts], mean_CB_err2[nCuts];
    
    // Input MC results
    std::ifstream input1(TString("fitresult_MC_"+lepName+".txt"));
    std::ifstream input2(TString("fitresult_zmass_diff_"+lepName+".txt"));
    for (int i=0; i<nCuts; i++) {
        input1 >> sigma_CB1[i];
        input1 >> sigma_CB_err1[i];
        input1 >> sigma_evt_by_evt1[i];
        input1 >> sigma_evt_by_evt_err_l1[i];
        input1 >> sigma_evt_by_evt_err_h1[i];
        input1 >> sigma_eff1[i];
        input1 >> sigma_eff_err1[i];
        input1 >> sigma_eff_check1[i];
        input1 >> alpha_CB1[i];
        input1 >> alpha_CB_err1[i];
        input1 >> mean_CB1[i];
        input1 >> mean_CB_err1[i];
        
        sigma_evt_by_evt1[i] = sigma_evt_by_evt1[i]*scale_fac;
        sigma_evt_by_evt_err_l1[i] = sigma_evt_by_evt_err_l1[i]*scale_fac;
        sigma_evt_by_evt_err_h1[i] = sigma_evt_by_evt_err_h1[i]*scale_fac;
        
        sigma_eff_CB1[i] = 1.0*sigma_eff_check1[i]/sigma_CB1[i];
        sigma_eff_CB_err1[i] = sqrt(pow(sigma_eff_err1[i]/sigma_CB1[i], 2) + pow(sigma_eff_check1[i]*sigma_CB_err1[i]/(sigma_CB1[i]*sigma_CB1[i]), 2));
        
//        std::cout << "sigma_CB 1\t" << i << "\t" << sigma_CB1[i] << "\t±\t" << sigma_CB_err1[i] << std::endl;
//        std::cout << "sigma_eff 1\t" << i << "\t" << sigma_eff1[i] << "\t±\t" << sigma_eff_err1[i] << std::endl;
//        std::cout << "sigma_eff check1\t" << i << "\t" << sigma_eff_check1[i] << std::endl;
//        std::cout << "sigma_eff_CB1\t" << i << "\t" << sigma_eff_CB1[i] << "\t±\t" << sigma_eff_CB_err1[i] << std::endl;
        
        input2 >> sigma_CB2[i];
        input2 >> sigma_CB_err2[i];
        input2 >> sigma_eff2[i];
        input2 >> sigma_eff_err2[i];
        input2 >> sigma_eff_check2[i];
        input2 >> alpha_CB2[i];
        input2 >> alpha_CB_err2[i];
        input2 >> mean_CB2[i];
        input2 >> mean_CB_err2[i];
        
        sigma_eff_CB2[i] = 1.0*sigma_eff2[i]/sigma_CB2[i];
        sigma_eff_CB_err2[i] = sqrt(pow(sigma_eff_err2[i]/sigma_CB2[i], 2) + pow(sigma_eff2[i]*sigma_CB_err2[i]/(sigma_CB2[i]*sigma_CB2[i]), 2));
        
//        std::cout << "sigma_CB\t" << i << "\t" << sigma_CB2[i] << "\t±\t" << sigma_CB_err2[i] << std::endl;
//        std::cout << "sigma_eff\t" << i << "\t" << sigma_eff2[i] << "\t±\t" << sigma_eff_err2[i] << std::endl;
//        std::cout << "sigma_eff check\t" << i << "\t" << sigma_eff_check2[i] << std::endl;
//        std::cout << "sigma_eff_CB\t" << i << "\t" << sigma_eff_CB2[i] << "\t±\t" << sigma_eff_CB_err2[i] << std::endl;
        
        
    }
    
    input1.close();
    input2.close();
    
    
    /******************************** Create Graphs ***************************/
    
    // Create data graphs
    TGraphAsymmErrors sigma_CB_evt_by_evt(nCuts1, sigma_evt_by_evt, sigma_CB, sigma_evt_by_evt_err_l, sigma_evt_by_evt_err_h, sigma_CB_err, sigma_CB_err);
    TGraphAsymmErrors alpha_CB_evt_by_evt(nCuts1, sigma_evt_by_evt, alpha_CB, sigma_evt_by_evt_err_l, sigma_evt_by_evt_err_h, alpha_CB_err, alpha_CB_err);
    TGraphAsymmErrors mean_CB_evt_by_evt(nCuts1, sigma_evt_by_evt, mean_CB, sigma_evt_by_evt_err_l, sigma_evt_by_evt_err_h, mean_CB_err, mean_CB_err);
    TGraphAsymmErrors sigma_eff_CB_evt_by_evt(nCuts1, sigma_evt_by_evt, sigma_eff_CB, sigma_evt_by_evt_err_l, sigma_evt_by_evt_err_h, sigma_eff_CB_err, sigma_eff_CB_err);
    TGraph line_middle(nCuts1, line, line);
    TGraph line_upper(nCuts1, line, line_u);
    TGraph line_lower(nCuts1, line, line_l);
    
    // Set graph style
    line_upper.SetLineStyle(2);
    line_lower.SetLineStyle(2);
    
    line_middle.SetLineColor(4);
    line_upper.SetLineColor(4);
    line_lower.SetLineColor(4);
    
    line_middle.SetLineWidth(2);
    line_upper.SetLineWidth(2);
    line_lower.SetLineWidth(2);
    
    sigma_CB_evt_by_evt.SetLineWidth(2);
    sigma_eff_CB_evt_by_evt.SetLineWidth(2);
    alpha_CB_evt_by_evt.SetLineWidth(2);
    mean_CB_evt_by_evt.SetLineWidth(2);
    
    gStyle->SetEndErrorSize(2);
    
    
    // Create MC graphs
    TGraphAsymmErrors sigma_CB_evt_by_evt1(nCuts, sigma_evt_by_evt1, sigma_CB1, sigma_evt_by_evt_err_l1, sigma_evt_by_evt_err_h1, sigma_CB_err1, sigma_CB_err1);
    TGraphAsymmErrors sigma_eff_CB_evt_by_evt1(nCuts, sigma_evt_by_evt1, sigma_eff_CB1, sigma_evt_by_evt_err_l1, sigma_evt_by_evt_err_h1, sigma_eff_CB_err1, sigma_eff_CB_err1);
    TGraphAsymmErrors alpha_CB_evt_by_evt1(nCuts, sigma_evt_by_evt1, alpha_CB1, sigma_evt_by_evt_err_l1, sigma_evt_by_evt_err_h1, alpha_CB_err1, alpha_CB_err1);
    TGraphAsymmErrors mean_CB_evt_by_evt1(nCuts, sigma_evt_by_evt1, mean_CB1, sigma_evt_by_evt_err_l1, sigma_evt_by_evt_err_h1, mean_CB_err1, mean_CB_err1);
    
    TGraphAsymmErrors sigma_CB_evt_by_evt2(nCuts, sigma_evt_by_evt1, sigma_CB2, sigma_evt_by_evt_err_l1, sigma_evt_by_evt_err_h1, sigma_CB_err2, sigma_CB_err2);
    TGraphAsymmErrors alpha_CB_evt_by_evt2(nCuts, sigma_evt_by_evt1, alpha_CB2, sigma_evt_by_evt_err_l1, sigma_evt_by_evt_err_h1, alpha_CB_err2, alpha_CB_err2);
    TGraphAsymmErrors sigma_eff_CB_evt_by_evt2(nCuts, sigma_evt_by_evt1, sigma_eff_CB2, sigma_evt_by_evt_err_l1, sigma_evt_by_evt_err_h1, sigma_eff_CB_err2, sigma_eff_CB_err2);
    
    
    // Set graph style
    sigma_CB_evt_by_evt1.SetLineWidth(2);
    sigma_eff_CB_evt_by_evt1.SetLineWidth(2);
    alpha_CB_evt_by_evt1.SetLineWidth(2);
    mean_CB_evt_by_evt1.SetLineWidth(2);
    
    sigma_CB_evt_by_evt2.SetLineWidth(2);
    sigma_eff_CB_evt_by_evt2.SetLineWidth(2);
    alpha_CB_evt_by_evt2.SetLineWidth(2);
    
    gStyle->SetEndErrorSize(2);
    
    sigma_CB_evt_by_evt1.SetLineColor(kRed);
    sigma_eff_CB_evt_by_evt1.SetLineColor(kRed);
    alpha_CB_evt_by_evt1.SetLineColor(kRed);
    mean_CB_evt_by_evt1.SetLineColor(kRed);
    
    sigma_CB_evt_by_evt1.SetMarkerColor(kRed);
    sigma_eff_CB_evt_by_evt1.SetMarkerColor(kRed);
    alpha_CB_evt_by_evt1.SetMarkerColor(kRed);
    mean_CB_evt_by_evt1.SetMarkerColor(kRed);
    
    
    sigma_CB_evt_by_evt2.SetLineColor(kGreen);
    sigma_eff_CB_evt_by_evt2.SetLineColor(kGreen);
    alpha_CB_evt_by_evt2.SetLineColor(kGreen);
    
    sigma_CB_evt_by_evt2.SetMarkerColor(kGreen);
    sigma_eff_CB_evt_by_evt2.SetMarkerColor(kGreen);
    alpha_CB_evt_by_evt2.SetMarkerColor(kGreen);
    
    
    // Add to multigraphs
    mg->Add(&sigma_CB_evt_by_evt, "*");
    mg->Add(&sigma_CB_evt_by_evt1, "*");
    mg->Add(&line_middle, "l");
    mg->Add(&line_upper, "l");
    mg->Add(&line_lower, "l");
    
    mg1->Add(&sigma_eff_CB_evt_by_evt, "*");
    mg1->Add(&sigma_eff_CB_evt_by_evt1, "*");
    
    mg2->Add(&sigma_CB_evt_by_evt1, "*");
    mg2->Add(&sigma_CB_evt_by_evt2, "*");
    mg2->Add(&line_middle, "l");
    mg2->Add(&line_upper, "l");
    mg2->Add(&line_lower, "l");
    
    mg3->Add(&sigma_eff_CB_evt_by_evt1, "*");
    mg3->Add(&sigma_eff_CB_evt_by_evt2, "*");
    
    mg_alpha->Add(&alpha_CB_evt_by_evt, "*");
    mg_alpha->Add(&alpha_CB_evt_by_evt1, "*");
    
    mg_mean->Add(&mean_CB_evt_by_evt1, "*");
    mg_mean->Add(&mean_CB_evt_by_evt, "*");
    
    
    
    /*********************************** Plot *********************************/
    
    Double_t sigma_evt_by_evt_err[nCuts];
    // Create graph to fit for linear fit
    TGraphErrors *tg_cal = new TGraphErrors(nCuts, sigma_evt_by_evt1, sigma_CB1, sigma_evt_by_evt_err, sigma_CB_err1);
    //TF1 func1("func1", "[0]*x", 0.7*scale_fac, 2.5*scale_fac);
    TF1 func1("func1", "[0]*x", 0.8*scale_fac, 2.5*scale_fac);
    TFitResultPtr fitresults1 = tg_cal->Fit("func1", "SR");
    tg_cal->SetTitle("");
    
    
    mg->SetTitle(TString("; #sigma_{event-by-event} [GeV/c^{2}];"));
    mg1->SetTitle(TString("; #sigma_{event-by-event} [GeV/c^{2}];"));
    mg2->SetTitle(TString("; #sigma_{event-by-event} [GeV/c^{2}];"));
    mg3->SetTitle(TString("; #sigma_{event-by-event} [GeV/c^{2}];"));
    mg_alpha->SetTitle(TString("; #sigma_{event-by-event} [GeV/c^{2}];"));
    mg_mean->SetTitle("; #sigma_{event-by-event} [GeV/c^{2}];");
    
    tg_cal->GetYaxis()->SetTitle("#sigma_{CB} [GeV/c^{2}]");
    tg_cal->GetXaxis()->SetTitle("#sigma_{event-by-event} [GeV/c^{2}]");
    
    
    TCanvas c("c", "c", 1000, 800);
    c.SetLeftMargin(0.15);
    c.SetBottomMargin(0.18);
    mg->Draw("AP");
    
    mg->GetYaxis()->SetTitle("#sigma_{CB} [GeV/c^{2}]");
    mg->GetYaxis()->SetTitleOffset(1.2);
    mg->GetYaxis()->SetTitleSize(0.06);
    mg->GetYaxis()->SetLabelSize(0.06);
    mg->GetYaxis()->SetTitleOffset(1.2);
    mg->GetXaxis()->SetTitleOffset(1.2);
    mg->GetXaxis()->SetLabelSize(0.06);
    mg->GetXaxis()->SetLimits(0.,5.5*scale_fac);
    mg->SetMinimum(0.);
    mg->SetMaximum(5.5*scale_fac);
    mg->GetXaxis()->SetTitleSize(0.06);
    
    
    TLegend *legend = new TLegend(0.2,0.7,0.4,0.9);
    legend->AddEntry(&sigma_CB_evt_by_evt,"Data","l");
    legend->AddEntry(&sigma_CB_evt_by_evt1,"MC","l");
    legend->AddEntry(&line_middle,"Perfect Calibration","l");
    legend->AddEntry(&line_upper,"+/- 20% envelope","l");
    legend->Draw();
    if (scale_fac != 1) {
        c.SaveAs(TString("Plots/sigma/calibrated/sigma_plot1"+detector_location+"_"+lepName+".png"));
    }
    else {
        c.SaveAs(TString("Plots/sigma/sigma_plot1"+detector_location+"_"+lepName+".png"));
    }
    
    
    
    
    TCanvas c1("c1", "c1", 1000, 800);
    c1.SetLeftMargin(0.15);
    c1.SetBottomMargin(0.18);
    
    mg1->Draw("AP");
    mg1->GetYaxis()->SetTitle("#sigma_{effective}/#sigma_{CB} [GeV/c^{2}]"); //"#alpha");
    mg1->GetYaxis()->SetTitleSize(0.06);
    mg1->GetYaxis()->SetLabelSize(0.06);
    mg1->GetYaxis()->SetTitleOffset(1.2);
    mg1->GetXaxis()->SetTitleOffset(1.2);
    mg1->GetXaxis()->SetLabelSize(0.06);
    mg1->GetXaxis()->SetTitleSize(0.06);
    mg1->SetMinimum(0.2);
    mg1->SetMaximum(3.5);
    
    TLegend *legend1 = new TLegend(0.7,0.7,0.9,0.9);
    legend1->AddEntry(&sigma_eff_CB_evt_by_evt,"Data","l");
    legend1->AddEntry(&sigma_eff_CB_evt_by_evt1,"MC","l");
    legend1->Draw();
    
    if (scale_fac != 1) {
        c1.SaveAs(TString("Plots/sigma/calibrated/sigma_plot2"+detector_location+"_"+lepName+".png"));
    }
    else {
        c1.SaveAs(TString("Plots/sigma/sigma_plot2"+detector_location+"_"+lepName+".png"));
    }
    
    TCanvas c2("c2", "c2", 1000, 800);
    c2.SetLeftMargin(0.15);
    c2.SetBottomMargin(0.18);
    
    mg2->Draw("AP");
    mg2->GetYaxis()->SetTitle("#sigma_{CB} [GeV/c^{2}]"); //"#alpha");
    mg2->GetYaxis()->SetTitleSize(0.06);
    mg2->GetYaxis()->SetLabelSize(0.06);
    mg2->GetYaxis()->SetTitleOffset(1.2);
    mg2->GetXaxis()->SetTitleOffset(1.2);
    mg2->GetXaxis()->SetLabelSize(0.06);
    mg2->GetXaxis()->SetTitleSize(0.06);
    mg2->GetXaxis()->SetLimits(0.,5.5*scale_fac);
    mg2->SetMinimum(0.);
    mg2->SetMaximum(5.5*scale_fac);
    
    TLegend *legend2 = new TLegend(0.2,0.7,0.4,0.9);
    legend2->AddEntry(&sigma_CB_evt_by_evt1,"MC","l");
    legend2->AddEntry(&sigma_CB_evt_by_evt2,"MC Z mass difference","l");
    legend2->AddEntry(&line_middle,"Perfect Calibration","l");
    legend2->AddEntry(&line_upper,"+/- 20% envelope","l");
    legend2->Draw();
    
    if (scale_fac != 1) {
        c2.SaveAs(TString("Plots/sigma/calibrated/sigma_MC_plot1"+detector_location+"_"+lepName+".png"));
    }
    else {
        c2.SaveAs(TString("Plots/sigma/sigma_MC_plot1"+detector_location+"_"+lepName+".png"));
    }
    
    TCanvas c3("c3", "c3", 1000, 800);
    c3.SetLeftMargin(0.15);
    c3.SetBottomMargin(0.18);
    
    mg3->Draw("AP");
    mg3->GetYaxis()->SetTitle("#sigma_{effective}/#sigma_{CB} [GeV/c^{2}]"); //"#alpha");
    mg3->GetYaxis()->SetTitleSize(0.06);
    mg3->GetYaxis()->SetLabelSize(0.06);
    mg3->GetYaxis()->SetTitleOffset(1.2);
    mg3->GetXaxis()->SetTitleOffset(1.2);
    mg3->GetXaxis()->SetLabelSize(0.06);
    mg3->GetXaxis()->SetTitleSize(0.06);
    mg3->SetMinimum(0.2);
    mg3->SetMaximum(3.5);
    
    TLegend *legend3 = new TLegend(0.5,0.7,0.9,0.9);
    legend3->AddEntry(&sigma_eff_CB_evt_by_evt1,"MC","l");
    legend3->AddEntry(&sigma_eff_CB_evt_by_evt2,"MC Z mass difference","l");
    legend3->Draw();
    
    if (scale_fac != 1) {
        c3.SaveAs(TString("Plots/sigma/calibrated/sigma_MC_plot2"+detector_location+"_"+lepName+".png"));
    }
    else {
        c3.SaveAs(TString("Plots/sigma/sigma_MC_plot2"+detector_location+"_"+lepName+".png"));
    }
    
    TCanvas c_alpha("c_alpha", "c_alpha", 1000, 800);
    c_alpha.SetLeftMargin(0.15);
    c_alpha.SetBottomMargin(0.18);
    
    mg_alpha->Draw("AP");
    mg_alpha->GetYaxis()->SetTitle("#alpha");
    mg_alpha->GetYaxis()->SetTitleSize(0.06);
    mg_alpha->GetYaxis()->SetLabelSize(0.06);
    mg_alpha->GetYaxis()->SetTitleOffset(1.2);
    mg_alpha->GetXaxis()->SetTitleOffset(1.2);
    mg_alpha->GetXaxis()->SetLabelSize(0.06);
    mg_alpha->GetXaxis()->SetTitleSize(0.06);
    
    
    TLegend *legend_a = new TLegend(0.7,0.7,0.9,0.9);
    legend_a->AddEntry(&alpha_CB_evt_by_evt,"Data","l");
    legend_a->AddEntry(&alpha_CB_evt_by_evt1,"MC","l");
    legend_a->Draw();
    
    if (scale_fac != 1) {
        c_alpha.SaveAs(TString("Plots/sigma/calibrated/alpha"+detector_location+"_"+lepName+".png"));
    }
    else {
        c_alpha.SaveAs(TString("Plots/sigma/alpha"+detector_location+"_"+lepName+".png"));
    }
    
    
    TCanvas c_mean("c_mean", "c_mean", 1000, 800);
    c_mean.SetLeftMargin(0.15);
    c_mean.SetBottomMargin(0.18);
    mg_mean->Draw("AP");
    
    mg_mean->GetYaxis()->SetTitle("mean [GeV/c^{2}]");
    mg_mean->GetYaxis()->SetTitleOffset(1.2);
    mg_mean->GetYaxis()->SetTitleSize(0.06);
    mg_mean->GetYaxis()->SetLabelSize(0.06);
    mg_mean->GetYaxis()->SetTitleOffset(1.2);
    mg_mean->GetXaxis()->SetTitleOffset(1.2);
    mg_mean->GetXaxis()->SetLabelSize(0.06);
    mg_mean->GetXaxis()->SetTitleSize(0.06);
    
    TLegend *legend_m = new TLegend(0.7,0.3,0.9,0.5);
    legend_m->AddEntry(&mean_CB_evt_by_evt,"Data","l");
    legend_m->AddEntry(&mean_CB_evt_by_evt1,"MC","l");
    legend_m->Draw();
    
    if (scale_fac != 1) {
        c_mean.SaveAs(TString("Plots/sigma/calibrated/mean"+detector_location+"_"+lepName+".png"));
    }
    else {
        c_mean.SaveAs(TString("Plots/sigma/mean"+detector_location+"_"+lepName+".png"));
    }
    
    
    TCanvas c_fit("c_fit", "c_fit", 1000, 800);
    c_fit.SetLeftMargin(0.15);
    c_fit.SetBottomMargin(0.18);
    
    tg_cal->Draw("AP*");
    
    tg_cal->GetYaxis()->SetTitleOffset(1.2);
    tg_cal->GetYaxis()->SetTitleSize(0.06);
    tg_cal->GetYaxis()->SetLabelSize(0.06);
    tg_cal->GetYaxis()->SetTitleOffset(1.2);
    tg_cal->GetXaxis()->SetTitleOffset(1.2);
    tg_cal->GetXaxis()->SetLabelSize(0.06);
    tg_cal->GetXaxis()->SetTitleSize(0.06);
    
    
    if (scale_fac != 1) {
        c_fit.SaveAs(TString("Plots/sigma/calibrated/calibration_fit"+detector_location+"_"+lepName+".png"));
    }
    else {
        c_fit.SaveAs(TString("Plots/sigma/calibration_fit"+detector_location+"_"+lepName+".png"));
    }
    
}
