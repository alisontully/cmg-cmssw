#include <iostream>
#include <algorithm>
#include <fstream>

#include "RooRealVar.h"
#include "TFile.h"
#include "TTree.h"
#include "RooDataSet.h"
#include "RooArgSet.h"
#include "RooCBShape.h"
#include "RooBreitWigner.h"
#include "RooFFTConvPdf.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TGraphAsymmErrors.h"
#include "RooDataHist.h"
#include "TH1.h"

using namespace RooFit;

RooFitResult* PlotMass(RooDataSet*, RooRealVar, RooAddPdf, RooFFTConvPdf, RooExponential, TString);
double sigma_eff_calc(RooRealVar, RooRealVar, RooCBShape, double, double);
Double_t effSigma(TH1*);


int main(int argc, char *argv[])
{
    std::string lepName = "muon";
    std::string dataType = "Data";
    std::string detector_location;// = "barrel";
    
    // Initialise variables
    RooRealVar z1_mass("z1_mass", "M_{Z} [GeV/c^{2}]",70, 110.);
    RooRealVar Lep1_eta("Lep1_eta", "Lep1_eta", -2.5, 2.5);
    RooRealVar Lep2_eta("Lep2_eta", "Lep2_eta", -2.5, 2.5);
    RooRealVar Lep1_pdgId("Lep1_pdgId", "Lepton 1 PDG ID", -100, 100);
    RooRealVar Lep2_pdgId("Lep2_pdgId", "Lepton 2 PDG ID", -100, 100);
    RooRealVar z1_delta_m("z1_delta_m", "z1_delta_m", 0.01, 5);
    
    RooRealVar z1_mass_copy("z1_mass_copy", "m_{2l} [GeV/c^{2}]", -50., 10.);
    z1_mass_copy.setBins(60/0.01);

    
    // Load file and tree
    TFile *file = new TFile();

    if (dataType == "Data") {
        if (lepName == "electron") {
	    file = new TFile("/afs/cern.ch/work/a/atully/private/Output_2l_new/DoubleEG_Run2015B_NewEleScale/twoLeptonTreeProducer/tree.root");
	}
        else if (lepName == "muon") {
            file = new TFile("/afs/cern.ch/work/a/atully/private/Output_data_new/DoubleMuon_merged.root");
        }
    }
    
    if (dataType == "MC") {
        file = new TFile("/afs/cern.ch/work/a/atully/private/Output_MC_2l/DYJetsToLL_LO_M50_50ns/twoLeptonTreeProducer/tree.root");
    }

    TTree *tree = (TTree*)file->Get("tree");
    
    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus("z1_mass", 1);
    tree->SetBranchStatus("Lep1_eta", 1);
    tree->SetBranchStatus("Lep2_eta", 1);
    tree->SetBranchStatus("Lep1_pdgId", 1);
    tree->SetBranchStatus("Lep2_pdgId", 1);
    tree->SetBranchStatus("z1_delta_m", 1);

    
    // import data into RooDataSet
    RooDataSet *data = new RooDataSet("data", "data", RooArgSet(z1_mass, Lep1_pdgId, Lep2_pdgId, Lep1_eta, Lep2_eta, z1_delta_m), Import(*tree));
    data->Print("V");
    
    
    // Cut for lepton type and detector location
    if (lepName == "muon") {
        if (detector_location == "barrel") {
            data = (RooDataSet*)data->reduce("Lep1_eta<1.1&&Lep2_eta<1.1&&abs(Lep1_pdgId)==13 && abs(Lep2_pdgId)==13");
        }
        else if (detector_location == "endcap") {
            data = (RooDataSet*)data->reduce("abs(Lep1_pdgId)==13 && abs(Lep2_pdgId)==13 && ((Lep1_eta>1.1&&Lep1_eta<2.4)||(Lep2_eta>1.1&&Lep2_eta<2.4))");
        }
        else {
            data = (RooDataSet*)data->reduce("abs(Lep1_pdgId)==13 && abs(Lep2_pdgId)==13");
            
        }
    }
    
    if (lepName == "electron") {
        if (detector_location == "barrel") {
            data = (RooDataSet*)data->reduce( "Lep1_eta<1.479&&Lep2_eta<1.479&&abs(Lep1_pdgId)==11 && abs(Lep2_pdgId)==11");
        }
        else if (detector_location == "endcap") {
            data = (RooDataSet*)data->reduce("abs(Lep1_pdgId)==11 && abs(Lep2_pdgId)==11 && ((Lep1_eta>1.479&&Lep1_eta<2.5)||(Lep2_eta>1.479&&Lep2_eta<2.5))");
        }
        else {
            data = (RooDataSet*)data->reduce("abs(Lep1_pdgId)==11 && abs(Lep2_pdgId)==11");
        }
    }

    

    /************ Calculate bin boundaries for evenly filled bins *************/
    
    int N = data->numEntries();
    double delta_m_array[N];
    
    // Fill array with delta_m values
    for (int i=0; i<N; i++) {
        delta_m_array[i] = data->get(i)->getRealValue("z1_delta_m");
    }
    
    // Sort array in ascending order
    std::sort(delta_m_array, delta_m_array + N);

    int nCuts;
    if (dataType == "Data") {
        nCuts = 20;
    }
    else if (dataType == "MC") {
        nCuts = 50;
    }
    
    std::vector<double> cutVal(nCuts);
    for (int i=0; i<nCuts; i++) {
        int eventNum = round((i+1)*1.0*N/nCuts);
        cutVal[i] = delta_m_array[eventNum];
    }
    cutVal[nCuts-1] = delta_m_array[N-1];
    
    // Create vector of cut strings
    std::vector<std::string> cutString(nCuts);
    cutString[0] = "z1_delta_m<"+std::to_string(cutVal[0]);
    for (int i=0; i<nCuts-2; i++) {
        cutString[i+1] = "z1_delta_m<"+std::to_string(cutVal[i+1])+"&&z1_delta_m>"+std::to_string(cutVal[i]);
    }
    cutString[nCuts-1] = "z1_delta_m>"+std::to_string(cutVal[nCuts-2]);
    
    
    // Save cut values and strings to file
    std::ofstream cutString_file;
    cutString_file.open(TString("Cuts/cutString"+detector_location+"_"+lepName+".txt"));
    for (int i=0; i<nCuts; i++) {
        cutString_file << cutString[i] << std::endl;
    }
    cutString_file.close();
    
    std::ofstream cutVal_file;
    cutVal_file.open(TString("Cuts/cutVal"+detector_location+"_"+lepName+".txt"));
    for (int i=0; i<nCuts; i++) {
        cutVal_file << cutVal[i] << std::endl;
    }
    
    cutVal_file << delta_m_array[0] << std::endl;

    cutVal_file.close();
    
    
    /************************** Create mass model *****************************/
    
    // Breit-Wigner
    // mean and width fixed to PDG values (http://pdg.lbl.gov/2014/listings/rpp2014-list-z-boson.pdf)
    RooRealVar m0("m0", "m0", 91.1876);
    RooRealVar width("width", "width", 2.4952);
    RooBreitWigner bw("bw", "bw", z1_mass, m0, width);
    
    // Crystal-Ball
    RooRealVar mean("mean", "mean", 0., -2., 1.);
    RooRealVar sigma("sigma", "sigma", 2.6, 0.1, 10.);
    RooRealVar alpha("alpha", "alpha", 1.2, 0., 5.);
    RooRealVar n("n", "n", 0.81, 0., 100.);
    RooCBShape cb("cb", "cb", z1_mass, mean, sigma, alpha, n);
    
    // Crystal-Ball copy for integration
    RooRealVar mean_copy("mean_copy", "mean_copy", 0., -2., 1.);
    RooRealVar sigma_copy("sigma_copy", "sigma_copy", 2.6, 0., 10.);
    RooRealVar alpha_copy("alpha_copy", "alpha_copy", 1.2, 0., 5.);
    RooRealVar n_copy("n_copy", "n_copy", 0.81, 0., 10000.);
    RooCBShape cb_copy("cb_copy", "cb_copy", z1_mass_copy, mean_copy, sigma_copy, alpha_copy, n_copy);
    
    // Convolution
    RooFFTConvPdf mass_sig("mass_sig", "mass_sig", z1_mass, bw, cb);
    
    // Background model
    RooRealVar bkg_exponent("bkg_exponent", "Exponent of Background", -0.006, -1, 1);
    RooExponential mass_bkg("mass_bkg", "mass_bkg", z1_mass, bkg_exponent);

    // Add signal and background mass models
    RooRealVar sf("sf", "sf", .9, 0.5, 1.);
    RooAddPdf mass_pdf("mass_pdf", "mass_pdf", RooArgList(mass_sig, mass_bkg), sf);
    
    
    
    /**************************** Fit z1 mass *********************************/
    
    Double_t sigma_CB[nCuts], sigma_CB_err[nCuts];
    Double_t sigma_evt_by_evt[nCuts], sigma_evt_by_evt_err_l[nCuts], sigma_evt_by_evt_err_h[nCuts];
    Double_t sigma_eff[nCuts], sigma_eff_err[nCuts], sigma_eff_check[nCuts];
    Double_t alpha_CB [nCuts], alpha_CB_err[nCuts];
    Double_t mean_CB [nCuts], mean_CB_err[nCuts];
    
    // Fit z1 mass for each delta_m bin
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
        RooFitResult *fitresult_bin = PlotMass(data_bin, z1_mass, mass_pdf, mass_sig, mass_bkg, TString("z1mass/"+dataType+"/"+lepName+"_z1mass"+detector_location+"_bin_"+std::to_string(i)+".png"));
        
        // Save fit results to file
        //mass_pdf.getParameters(*data_bin)->writeToFile(TString("FitResults/"+dataType+"/z1mass_"+lepName+""+detector_location+"_bin_"+std::to_string(i)+".txt"));
        
        cb_copy.getParameters(*data_bin)->setRealValue("mean_copy", cb.getParameters(*data_bin)->getRealValue("mean"));
        cb_copy.getParameters(*data_bin)->setRealValue("alpha_copy", cb.getParameters(*data_bin)->getRealValue("alpha"));
        cb_copy.getParameters(*data_bin)->setRealValue("sigma_copy", cb.getParameters(*data_bin)->getRealValue("sigma"));
        cb_copy.getParameters(*data_bin)->setRealValue("n_copy", cb.getParameters(*data_bin)->getRealValue("n"));
        
        // Save sigma_CB, alpha and mean from the fit result
        sigma_CB[i] = mass_pdf.getParameters(*data_bin)->getRealValue("sigma");
        RooRealVar *sigma_CB_var = (RooRealVar*) fitresult_bin->floatParsFinal().find("sigma");
        sigma_CB_err[i] = sigma_CB_var->getError();
        
        alpha_CB[i] = mass_pdf.getParameters(*data_bin)->getRealValue("alpha");
        RooRealVar *alpha_CB_var = (RooRealVar*) fitresult_bin->floatParsFinal().find("alpha");
        alpha_CB_err[i] = alpha_CB_var->getError();
        
        mean_CB[i] = mass_pdf.getParameters(*data_bin)->getRealValue("mean");
        RooRealVar *mean_CB_var = (RooRealVar*) fitresult_bin->floatParsFinal().find("mean");
        mean_CB_err[i] = mean_CB_var->getError();

 
        /********************** Calculate effective sigma *********************/
        
        double stepsize = 0.01;
        sigma_eff[i] = sigma_eff_calc(z1_mass_copy, mean, cb_copy, sigma_CB[i], stepsize);
        sigma_eff_err[i] = stepsize;
        
        
        // Check calculation of sigma effective
        TH1 *cb_hist = (TH1*) cb_copy.createHistogram("cb_hist", z1_mass_copy, Binning(60/0.01, -50, 10));
        sigma_eff_check[i] = effSigma(cb_hist);
        
        
        // Calculate average delta_m in each bin
        double delta_m_sum = 0.;

        for (int j=0; j<data_bin->numEntries(); j++) {
            delta_m_sum = delta_m_sum + data_bin->get(j)->getRealValue("z1_delta_m");
        }
        sigma_evt_by_evt[i] = delta_m_sum/data_bin->numEntries();
        sigma_evt_by_evt_err_h[i] = cutVal[i]-sigma_evt_by_evt[i];
        sigma_evt_by_evt_err_l[i] = sigma_evt_by_evt[i]-cutVal[i-1];

    }
    sigma_evt_by_evt_err_l[0] = sigma_evt_by_evt[0] - delta_m_array[0];
    
    
    /************************** Save results to file **************************/
    
    std::ofstream fitresult;
    fitresult.open(TString("fitresult"+detector_location+"_"+dataType+"_"+lepName+".txt"));
    
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












RooFitResult* PlotMass(RooDataSet* data, RooRealVar var, RooAddPdf pdf, RooFFTConvPdf pdf_sig, RooExponential pdf_bkg, TString name) {
    
    // Fit mass pdf to data
    RooFitResult *mass_fitresult = pdf.fitTo(*data, Extended(false), Save(true), Range(70.,110.), Strategy(2)); // Strategy(0),
    mass_fitresult->Print();

    // Plot
    TCanvas c("c", "c", 800, 600);
    RooPlot *frame = var.frame(Range(70., 110.));
    frame->SetTitle("");
    frame->GetXaxis()->SetTitle("m_{2l} [GeV/c^{2}]");
    frame->GetYaxis()->SetTitle("Events / 0.4 GeV/c^{2}");
    frame->SetTitleSize(0.05, "XY");
    frame->GetYaxis()->SetTitleOffset(1);
    frame->GetXaxis()->SetTitleOffset(0.83);
    data->plotOn(frame);
    pdf.plotOn(frame, LineColor(kBlack));
    pdf.plotOn(frame, Components(pdf_sig), LineStyle(2), LineColor(kBlue));
    pdf.plotOn(frame, Components(pdf_bkg), LineStyle(2), LineColor(kRed));
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


