#include <iostream>
#include <algorithm>
#include <fstream>

#include "RooRealVar.h"
#include "TFile.h"
#include "TTree.h"
#include "RooDataSet.h"
#include "RooArgSet.h"
#include "RooCBShape.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TH1.h"

RooFitResult* PlotMass(RooDataSet*, RooRealVar, RooCBShape, TString);
double sigma_eff_calc(RooRealVar, RooRealVar, RooCBShape, double, double);
Double_t effSigma(TH1*);

using namespace RooFit;

int main(int argc, char *argv[])
{
    std::string lepName = "electron";
    
    // Initialise variables
    RooRealVar zz1_mass("zz1_mass", "m_{4l} [GeV/c^{2}]",95, 145.);
    zz1_mass.setBins(50/0.01);
    RooRealVar zz1_delta_m("zz1_delta_m", "zz1_delta_m", 0.01, 5);
    RooRealVar Lep1_pdgId("Lep1_pdgId", "Lep1_pdgId", -100, 100);
    RooRealVar Lep2_pdgId("Lep2_pdgId", "Lep2_pdgId", -100, 100);
    RooRealVar Lep3_pdgId("Lep3_pdgId", "Lep3_pdgId", -100, 100);
    RooRealVar Lep4_pdgId("Lep4_pdgId", "Lep4_pdgId", -100, 100);
 
    TFile *file = new TFile("/afs/cern.ch/user/a/atully/CMSSW_7_4_7/src/CMGTools/HToZZ4L/cfg/Output_4l_newEleScale/GGHZZ4L_126_NewEleScale/fourLeptonTreeProducer/tree.root");
    TTree *tree = (TTree*)file->Get("tree");
    
    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus("zz1_mass", 1);
    tree->SetBranchStatus("zz1_delta_m", 1);
    tree->SetBranchStatus("Lep1_pdgId", 1);
    tree->SetBranchStatus("Lep2_pdgId", 1);
    tree->SetBranchStatus("Lep3_pdgId", 1);
    tree->SetBranchStatus("Lep4_pdgId", 1);

    
    // import data into RooDataSet
    RooDataSet *data = new RooDataSet("data", "data", RooArgSet(zz1_mass, zz1_delta_m, Lep1_pdgId, Lep2_pdgId, Lep3_pdgId, Lep4_pdgId), Import(*tree));
    data->Print("V");
    
    // Reduce data for lepton type
    if (lepName == "electron") {
	data = (RooDataSet*)data->reduce("abs(Lep1_pdgId)==11&&abs(Lep2_pdgId)==11&&abs(Lep3_pdgId)==11&&abs(Lep4_pdgId)==11");
    }

    if (lepName == "muon") {
        data = (RooDataSet*)data->reduce("abs(Lep1_pdgId)==13&&abs(Lep2_pdgId)==13&&abs(Lep3_pdgId)==13&&abs(Lep4_pdgId)==13");
    }
    
    data->Print("V");

    
    /********** Calculate bin boundaries for evenly filled bins ***************/
   
    int N = data->numEntries();
    double delta_m_array[N];
    
    // Fill array with delta_m values
    for (int i=0; i<N; i++) {
        delta_m_array[i] = data->get(i)->getRealValue("zz1_delta_m");
    }
    
    // Sort array in ascending order
    std::sort(delta_m_array, delta_m_array + N);
    
    int nCuts = 50;
    std::vector<double> cutVal(nCuts);
    
    for (int i=0; i<nCuts; i++) {
        int eventNum = round((i+1)*1.0*N/nCuts);
        cutVal[i] = delta_m_array[eventNum];
    }
    cutVal[nCuts-1] = delta_m_array[N-1];
    
    // Create vector of cut strings
    std::vector<std::string> cutString(nCuts);
    cutString[0] = "zz1_delta_m<"+std::to_string(cutVal[0]);
    for (int i=0; i<nCuts-2; i++) {
        cutString[i+1] = "zz1_delta_m<"+std::to_string(cutVal[i+1])+"&&zz1_delta_m>"+std::to_string(cutVal[i]);
    }
    cutString[nCuts-1] = "zz1_delta_m>"+std::to_string(cutVal[nCuts-2]);
    
    
    // Save cut string and values to file
    std::ofstream cutString_file;
    cutString_file.open(TString("Cuts/cutString_H_"+lepName+".txt"));
    for (int i=0; i<nCuts; i++) {
        cutString_file << cutString[i] << std::endl;
    }
    cutString_file.close();
    
    std::ofstream cutVal_file;
    cutVal_file.open(TString("Cuts/cutVal_H_"+lepName+".txt"));
    for (int i=0; i<nCuts; i++) {
        cutVal_file << cutVal[i] << std::endl;
    }
    cutVal_file << delta_m_array[0] << std::endl;
    cutVal_file.close();
    
    
    /******************************* Mass model *******************************/
    
    // Crystal-Ball
    RooRealVar mean("mean", "mean", 126., 124., 128.);
    RooRealVar sigma("sigma", "sigma", 2.6, 0.1, 10.);
    RooRealVar alpha("alpha", "alpha", 1.2, 0., 5.);
    RooRealVar n("n", "n", 0.81, 0., 100.);
    RooCBShape cb("cb", "cb", zz1_mass, mean, sigma, alpha, n);
    
    RooFitResult *fitresult_mass = cb.fitTo(*data, Save(true), Strategy(2));
    
    TCanvas c("c", "c", 600, 800);
    c.SetLeftMargin(0.15);
    RooPlot *frame = zz1_mass.frame();
    frame->GetYaxis()->SetTitle("Events / 0.5 GeV/c^{2}");
    frame->SetTitleOffset(1.4, "Y");
    frame->SetTitleOffset(0.9, "X");
    frame->SetTitle("");
    data->plotOn(frame);
    cb.plotOn(frame);
    frame->SetMaximum(1100);
    frame->SetTitleSize(0.05, "XY");
    frame->Draw();
    c.SaveAs(TString("Plots/Higgs/zzmass_"+lepName+".png"));
    
    
    /***************************** Fit H mass *********************************/
    
    Double_t sigma_CB[nCuts], sigma_CB_err[nCuts];
    Double_t sigma_evt_by_evt[nCuts], sigma_evt_by_evt_err_l[nCuts], sigma_evt_by_evt_err_h[nCuts];
    Double_t sigma_eff[nCuts], sigma_eff_err[nCuts], sigma_eff_check[nCuts];
    Double_t alpha_CB [nCuts], alpha_CB_err[nCuts];
    Double_t mean_CB [nCuts], mean_CB_err[nCuts];
    
    // Fit H  mass for each delta_m bin
    for (int i=0; i<nCuts; i++) {
        
        RooDataSet *data_bin = (RooDataSet*)data->reduce(TString(cutString[i]));
        data_bin->Print("V");
        
        // Set initial values
        if (i==0) {
            cb.getParameters(data_bin)->setRealValue("sigma", cutVal[i]*0.5);
        }
        else {
            cb.getParameters(data_bin)->setRealValue("sigma", (cutVal[i]+cutVal[i-1])*0.5);
        }
        
        // Fit and plot the mass
        RooFitResult *fitresult_bin = PlotMass(data_bin, zz1_mass, cb, TString("Higgs/zz1mass/zz1mass_"+lepName+"_bin_"+std::to_string(i)+".png"));
        
        // Save fit results to file
        cb.getParameters(*data_bin)->writeToFile(TString("FitResults_Higgs/zz1mass_"+lepName+"_bin_"+std::to_string(i)+".txt"));
        
        // Save sigma_CB, alpha and mean from the fit result
        sigma_CB[i] = cb.getParameters(*data_bin)->getRealValue("sigma");
        RooRealVar *sigma_CB_var = (RooRealVar*) fitresult_bin->floatParsFinal().find("sigma");
        sigma_CB_err[i] = sigma_CB_var->getError();
        
        alpha_CB[i] = cb.getParameters(*data_bin)->getRealValue("alpha");
        RooRealVar *alpha_CB_var = (RooRealVar*) fitresult_bin->floatParsFinal().find("alpha");
        alpha_CB_err[i] = alpha_CB_var->getError();
        
        mean_CB[i] = cb.getParameters(*data_bin)->getRealValue("mean");
        RooRealVar *mean_CB_var = (RooRealVar*) fitresult_bin->floatParsFinal().find("mean");
        mean_CB_err[i] = mean_CB_var->getError();
        
        
        /********************** Calculate effective sigma *********************/
        
        double stepsize = 0.01;
        sigma_eff[i] = sigma_eff_calc(zz1_mass, mean, cb, sigma_CB[i], stepsize);
        sigma_eff_err[i] = stepsize;
        
        // Check calculation of sigma effective
        TH1 *cb_hist = (TH1*) cb.createHistogram("cb_hist", zz1_mass, Binning(50/0.01, 95, 145));
        
        sigma_eff_check[i] = effSigma(cb_hist);
        
        
        // Calculate average delta_m in each bin
        double delta_m_sum = 0.;
        
        for (int j=0; j<data_bin->numEntries(); j++) {
            delta_m_sum = delta_m_sum + data_bin->get(j)->getRealValue("zz1_delta_m");
        }
        
        sigma_evt_by_evt[i] = delta_m_sum/data_bin->numEntries();
        sigma_evt_by_evt_err_h[i] = cutVal[i]-sigma_evt_by_evt[i];
        sigma_evt_by_evt_err_l[i] = sigma_evt_by_evt[i]-cutVal[i-1];
        std::cout << "sigma_evt_by_evt\t" << sigma_evt_by_evt[i] << std::endl;
 
    }
    
    sigma_evt_by_evt_err_l[0] = sigma_evt_by_evt[0] - delta_m_array[0];
    
    
    /************************** Save results to file **************************/
    
    std::ofstream fitresult;
    fitresult.open(TString("fitresult_Higgs_"+lepName+".txt"));
    
    // Print results
    for (int i=0; i<nCuts; i++) {
        
//        std::cout << "Bin number\t" << i << std::endl;
//        std::cout << "Crystal ball sigma\t" << sigma_CB[i] << "\t±\t" << sigma_CB_err[i] << std::endl;
//        std::cout << "Effective sigma\t" << sigma_eff[i] << "\t±\t" << sigma_eff_err[i] << std::endl;
//        std::cout << "Effective sigma check\t" << sigma_eff_check[i] << std::endl;
//        std::cout << "Average delta_m in bin\t" << sigma_evt_by_evt[i] << std::endl;
        
        fitresult << sigma_CB[i] << std::endl;
        fitresult << sigma_CB_err[i] << std::endl;
        fitresult << sigma_evt_by_evt[i] << std::endl;
        fitresult << sigma_evt_by_evt_err_l[i] << std::endl;
        fitresult << sigma_evt_by_evt_err_h[i] << std::endl;
        fitresult << sigma_eff[i] << std::endl;
        fitresult << sigma_eff_err[i] << std::endl;
        fitresult << sigma_eff_check[i] << std::endl;
        fitresult << alpha_CB[i] << std::endl;
        fitresult << alpha_CB_err[i] << std::endl;
        fitresult << mean_CB[i] << std::endl;
        fitresult << mean_CB_err[i] << std::endl;
    }
    
    fitresult.close();
    
    data->Delete();
    tree->Delete();
    file->Close();
}






RooFitResult* PlotMass(RooDataSet* data, RooRealVar var, RooCBShape pdf, TString name) {
    
    // Fit mass pdf to data
    RooFitResult *mass_fitresult = pdf.fitTo(*data, Extended(false), Save(true), Strategy(2)); // Strategy(0),
    mass_fitresult->Print();
    
    // Plot
    TCanvas c("c", "c", 800, 600);
    RooPlot *frame = var.frame(Bins(30));
    frame->SetTitle("");
    frame->GetXaxis()->SetTitle("m_{4l} [GeV/c^{2}]");
    //frame->GetXaxis()->SetTitleSize(1.2);
    frame->GetYaxis()->SetTitleOffset(1.2);
    data->plotOn(frame);
    pdf.plotOn(frame);
    frame->Draw();
    c.SaveAs(TString("Plots/"+name));
    
    return mass_fitresult;
    
    
}


double sigma_eff_calc(RooRealVar z1_mass_copy, RooRealVar mean, RooCBShape cb_copy, double sigma_CB, double stepsize)
{
    z1_mass_copy.setRange("integralrange", mean.getVal()-sigma_CB, mean.getVal()+sigma_CB);
    RooAbsReal* integral = cb_copy.createIntegral(z1_mass_copy,NormSet(z1_mass_copy), Range("integralrange"));
    double normalizedIntegralValue = integral->getVal();
    
    std::cout << "Integral value intial\t" << normalizedIntegralValue << std::endl;
    
    double check = 2, check_new = 2;
    if (normalizedIntegralValue < 0.683) {
        check = 1;
        check_new = 1;
    }
    
    double range_l = mean.getVal()-sigma_CB;
    double range_h = mean.getVal()+sigma_CB;
    
    while (check == check_new) {
        
        z1_mass_copy.setRange("integralrange", range_l, range_h);
        RooAbsReal* integral = cb_copy.createIntegral(z1_mass_copy,NormSet(z1_mass_copy), Range("integralrange"));
        normalizedIntegralValue = integral->getVal();
        
        //std::cout << "integral value\t" << normalizedIntegralValue << std::endl;
        
        if (normalizedIntegralValue<0.683) {
            range_l = range_l - stepsize;
            range_h = range_h + stepsize;
            check = 1;
        }
        else if (normalizedIntegralValue>0.683) {
            range_h = range_h - stepsize;
            range_l = range_l + stepsize;
            check = 2;
        }
    }
    
    std::cout << "Sigma effective\t" << range_h-mean.getVal() << std::endl;
    std::cout << "Integral value\t" << normalizedIntegralValue << std::endl;
    
    return range_h-mean.getVal();
}



//effsigma function from Chris
Double_t effSigma(TH1 * hist)
{
    
    TAxis *xaxis = hist->GetXaxis();
    Int_t nb = xaxis->GetNbins();
    if(nb < 10) {
        std::cout << "effsigma: Not a valid histo. nbins = " << nb << std::endl;
        return 0.;
    }
    
    Double_t bwid = xaxis->GetBinWidth(1);
    if(bwid == 0) {
        std::cout << "effsigma: Not a valid histo. bwid = " << bwid << std::endl;
        return 0.;
    }
    Double_t xmax = xaxis->GetXmax();
    Double_t xmin = xaxis->GetXmin();
    Double_t ave = hist->GetMean();
    Double_t rms = hist->GetRMS();
    
    Double_t total=0.;
    for(Int_t i=0; i<nb+2; i++) {
        total+=hist->GetBinContent(i);
    }
    //   if(total < 100.) {
    //     cout << "effsigma: Too few entries " << total << endl;
    //     return 0.;
    //   }
    Int_t ierr=0;
    Int_t ismin=999;
    
    Double_t rlim=0.683*total;
    Int_t nrms=rms/(bwid);    // Set scan size to +/- rms
    if(nrms > nb/10) nrms=nb/10; // Could be tuned...
    
    Double_t widmin=9999999.;
    for(Int_t iscan=-nrms;iscan<nrms+1;iscan++) { // Scan window centre
        Int_t ibm=(ave-xmin)/bwid+1+iscan;
        Double_t x=(ibm-0.5)*bwid+xmin;
        Double_t xj=x;
        Double_t xk=x;
        Int_t jbm=ibm;
        Int_t kbm=ibm;
        Double_t bin=hist->GetBinContent(ibm);
        total=bin;
        for(Int_t j=1;j<nb;j++){
            if(jbm < nb) {
                jbm++;
                xj+=bwid;
                bin=hist->GetBinContent(jbm);
                total+=bin;
                if(total > rlim) break;
            }
            else ierr=1;
            if(kbm > 0) {
                kbm--;
                xk-=bwid;
                bin=hist->GetBinContent(kbm);
                total+=bin;
                if(total > rlim) break;
            }
            else ierr=1;
        }
        Double_t dxf=(total-rlim)*bwid/bin;
        Double_t wid=(xj-xk+bwid-dxf)*0.5;
        if(wid < widmin) {
            widmin=wid;
            ismin=iscan;
        }
    }
    if(ismin == nrms || ismin == -nrms) ierr=3;
    if(ierr != 0) std::cout << "effsigma: Error of type " << ierr << std::endl;
    
    return widmin;
    
}


