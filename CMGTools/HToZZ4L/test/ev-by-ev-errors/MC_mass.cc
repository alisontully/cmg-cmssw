#include <iostream>
#include <fstream>

#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TPad.h"
#include "RooRealVar.h"
#include "TFile.h"
#include "TTree.h"
#include "RooArgSet.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooCBShape.h"
#include "RooCategory.h"
#include "RooPlot.h"
#include "RooFitResult.h"

using namespace RooFit;
double sigma_eff_calc(RooRealVar, RooRealVar, RooCBShape, double, double);
Double_t effSigma(TH1*);


int main(int argc, char *argv[])
{
    std::string lepName = "muon";
    std::string detector_location;// = "endcap";
    
    // Initialise variables
    RooRealVar z1_mass("z1_mass", "M_{Z} [GeV/c^{2}]",50, 130.);
    RooRealVar z1_delta_m("z1_delta_m", "z1_delta_m", 0.01, 5);
    RooRealVar Lep1_eta("Lep1_eta", "Lep1_eta", -2.5, 2.5);
    RooRealVar Lep2_eta("Lep2_eta", "Lep2_eta", -2.5, 2.5);
    RooRealVar Lep1_phi("Lep1_phi", "Lep1_phi", -3.14159, 3.14159);
    RooRealVar Lep2_phi("Lep2_phi", "Lep2_phi", -3.14159, 3.14159);
    RooRealVar Lep1_pdgId("Lep1_pdgId", "Lepton 1 PDG ID", -100, 100);
    RooRealVar Lep2_pdgId("Lep2_pdgId", "Lepton 2 PDG ID", -100, 100);
    RooRealVar Lep1_pt("Lep1_pt", "Lep1_pt", -100, 100);
    RooRealVar Lep2_pt("Lep2_pt", "Lep2_pt", -100, 100);
    
    // For MC data only
    RooRealVar Lep1_mcEta("Lep1_mcEta", "Lep1_mcEta", -2.5, 2.5);
    RooRealVar Lep2_mcEta("Lep2_mcEta", "Lep2_mcEta", -2.5, 2.5);
    RooRealVar Lep1_mcPhi("Lep1_mcPhi", "Lep1_mcPhi", -3.14159, 3.14159);
    RooRealVar Lep2_mcPhi("Lep2_mcPhi", "Lep2_mcPhi", -3.14159, 3.14159);
    RooRealVar Lep1_mcPt("Lep1_mcPt", "Lep1_mcPt", 0.01, 100);
    RooRealVar Lep2_mcPt("Lep2_mcPt", "Lep2_mcPt", 0.01, 100);
    
    RooRealVar z1_mass_diff("z1_mass_diff", "z1_mass_diff", -50, 10);
    z1_mass_diff.setBins(60/0.01);

    // Load file and tree
    TFile *file = new TFile("/afs/cern.ch/work/a/atully/private/Output_MC_2l/DYJetsToLL_LO_M50_50ns/twoLeptonTreeProducer/tree.root");
    TTree *tree = (TTree*)file->Get("tree");
    
    tree->SetBranchStatus("*", 0);
    // Enable interesting branches, e.g.:
    tree->SetBranchStatus("z1_mass", 1);
    tree->SetBranchStatus("Lep1_eta", 1);
    tree->SetBranchStatus("Lep2_eta", 1);
    tree->SetBranchStatus("Lep1_phi", 1);
    tree->SetBranchStatus("Lep2_phi", 1);
    tree->SetBranchStatus("Lep1_pdgId", 1);
    tree->SetBranchStatus("Lep2_pdgId", 1);
    tree->SetBranchStatus("Lep1_pt", 1);
    tree->SetBranchStatus("Lep2_pt", 1);
    tree->SetBranchStatus("Lep1_mcPhi", 1);
    tree->SetBranchStatus("Lep2_mcPhi", 1);
    tree->SetBranchStatus("Lep1_mcEta", 1);
    tree->SetBranchStatus("Lep2_mcEta", 1);
    tree->SetBranchStatus("Lep1_mcPt", 1);
    tree->SetBranchStatus("Lep2_mcPt", 1);
    tree->SetBranchStatus("z1_delta_m", 1);
    
    RooArgSet vars(z1_mass, Lep1_pdgId, Lep2_pdgId, Lep1_eta, Lep2_eta, Lep1_phi, Lep2_phi, Lep1_pt, Lep2_pt);
    vars.add(RooArgSet(Lep1_mcEta, Lep2_mcEta, Lep1_mcPhi, Lep2_mcPhi, Lep1_mcPt, Lep2_mcPt, z1_mass_diff, z1_delta_m));
    
    std::string cut_lep;
    if (lepName == "muon") {
        if (detector_location == "barrel") {
            cut_lep = "Lep1_eta<1.1&&Lep2_eta<1.1&&abs(Lep1_pdgId)==13 && abs(Lep2_pdgId)==13";
        }
        else if (detector_location == "endcap") {
            cut_lep = "abs(Lep1_pdgId)==13 && abs(Lep2_pdgId)==13 && ((Lep1_eta>1.1&&Lep1_eta<2.4)||(Lep2_eta>1.1&&Lep2_eta<2.4))";
        }
        else {
            cut_lep = "abs(Lep1_pdgId)==13 && abs(Lep2_pdgId)==13";
            
        }
    }
    
    if (lepName == "electron") {
        if (detector_location == "barrel") {
            cut_lep = "Lep1_eta<1.479&&Lep2_eta<1.479&&abs(Lep1_pdgId)==11 && abs(Lep2_pdgId)==11";
        }
        else if (detector_location == "endcap") {
            cut_lep = "abs(Lep1_pdgId)==11 && abs(Lep2_pdgId)==11 && ((Lep1_eta>1.479&&Lep1_eta<2.5)||(Lep2_eta>1.479&&Lep2_eta<2.5))";
        }
        else {
            cut_lep = "abs(Lep1_pdgId)==11 && abs(Lep2_pdgId)==11";
            
        }
    }

    int nCuts = 50;
    
    // Read vector of cut strings
    std::vector<std::string> cutString(nCuts);
    std::ifstream input_cuts(TString("Cuts/cutString"+detector_location+"_"+lepName+".txt"));
    for (int i=0; i<nCuts; i++) {
        std::getline(input_cuts, cutString[i]);
    }
    
    // Read vector of cut values
    std::vector<double> cutVal(nCuts);
    std::ifstream input_cutVal(TString("Cuts/cutVal"+detector_location+"_"+lepName+".txt"));
    for (int i=0; i<nCuts; i++) {
        input_cutVal >> cutVal[i];
    }
    
    double delta_m_min;
    input_cutVal >> delta_m_min;
    
    
    // Mass resolution model
    // Crystal-Ball
    RooRealVar mean("mean", "mean", 0., -2., 2.);
    RooRealVar sigma("sigma", "sigma", 2.6, 0., 5.);
    RooRealVar alpha("alpha", "alpha", 1.2, 0., 5.);
    RooRealVar n("n", "n", 0.81, 0., 100.);
    RooCBShape cb("cb", "cb", z1_mass_diff, mean, sigma, alpha, n);
    
    
    
    /************************* Fit to mass resolution *************************/
    
    // Create lepton vectors
    TLorentzVector *lep1 = new TLorentzVector();
    TLorentzVector *lep2 = new TLorentzVector();
    
    std::vector<double> m(2);
    
    std::vector<std::string> lep_pt(4);
    lep_pt[0] = "Lep1_pt";
    lep_pt[1] = "Lep2_pt";
    lep_pt[2] = "Lep1_mcPt";
    lep_pt[3] = "Lep2_mcPt";
    
    double lep_mass;
    if (lepName == "electron") {
        lep_mass = 0.000511;
    }
    else if (lepName == "muon") {
        lep_mass = 0.106;
    }
    
    Double_t sigma_CB[nCuts], sigma_CB_err[nCuts], sigma_eff[nCuts], sigma_eff_err[nCuts], sigma_eff_check[nCuts];
    Double_t alpha_CB[nCuts], alpha_CB_err[nCuts];
    Double_t mean_CB[nCuts], mean_CB_err[nCuts];

    for (int x=0; x<nCuts; x++) {
        // import data into RooDataSet
        RooDataSet *data = new RooDataSet("data", "data", vars, RooFit::Import(*tree), Cut(TString(cutString[x]+"&&"+cut_lep)));
        data = (RooDataSet*)data->reduce(TString(cutString[x]));
        data = (RooDataSet*)data->reduce(TString(cut_lep));
        data->Print("V");
        
        std::cout << "Cut String\t" << TString(cutString[x]+"&&"+cut_lep) << std::endl;
        
        RooDataSet data_new("data_new", "data_new", RooArgSet(z1_mass_diff));
        
        for (int i=0; i<data->numEntries(); i++) {
            for (int j=0; j<2; j++) {
                
                lep1->SetPtEtaPhiM(data->get(i)->getRealValue(TString(lep_pt[2*j])),
                                   data->get(i)->getRealValue("Lep1_eta"),
                                   data->get(i)->getRealValue("Lep1_phi"),
                                   lep_mass);
                
                lep2->SetPtEtaPhiM(data->get(i)->getRealValue(TString(lep_pt[j*2+1])),
                                   data->get(i)->getRealValue("Lep2_eta"),
                                   data->get(i)->getRealValue("Lep2_phi"),
                                   lep_mass);
                
                TLorentzVector z = lep1->operator+(*lep2);
                m[j] = z.M();
            }
            
            z1_mass_diff = m[0]-m[1];
            data_new.add(RooArgSet(z1_mass_diff));
        }
        
        data_new.Print("V");
        
        // Set initial values
        if (x==0) {
            cb.getParameters(data_new)->setRealValue("sigma", cutVal[x]*0.5);
        }
        else {
            cb.getParameters(data_new)->setRealValue("sigma", (cutVal[x]+cutVal[x-1])*0.5);
        }
        
        // Fit to mass resolution
        RooFitResult *fitresult_bin = cb.fitTo(data_new, Save(true), Strategy(2)); //, Range(-15,10),
        
        // Save sigma_CB, alpha and mean from the fit result
        sigma_CB[x] = cb.getParameters(data_new)->getRealValue("sigma");
        RooRealVar *sigma_CB_var = (RooRealVar*) fitresult_bin->floatParsFinal().find("sigma");
        sigma_CB_err[x] = sigma_CB_var->getError();
        
        alpha_CB[x] = cb.getParameters(data_new)->getRealValue("alpha");
        RooRealVar *alpha_CB_var = (RooRealVar*) fitresult_bin->floatParsFinal().find("alpha");
        alpha_CB_err[x] = alpha_CB_var->getError();
        
        mean_CB[x] = cb.getParameters(data_new)->getRealValue("mean");
        RooRealVar *mean_CB_var = (RooRealVar*) fitresult_bin->floatParsFinal().find("mean");
        mean_CB_err[x] = mean_CB_var->getError();

        
        /********************** Calculate effective sigma *********************/
        
        double stepsize = 0.01;
        sigma_eff[i] = sigma_eff_calc(z1_mass_diff, mean, cb, sigma_CB[i], stepsize);
        sigma_eff_err[i] = stepsize;

        
        // Check calculation of sigma effective
        TH1 *cb_hist = (TH1*) cb.createHistogram("cb_hist", z1_mass_diff, Binning(60/0.01, -15, 10));
        
        sigma_eff_check[x] = effSigma(cb_hist);
        
//        std::cout << "Sigma effective\t" << sigma_eff[x] << std::endl;
//        std::cout << "Integral value\t" << normalizedIntegralValue << std::endl;
//        std::cout << "Sigma effective check\t" << sigma_eff_check[x] << std::endl;

        
        /******************************** Plot ********************************/
        
        TCanvas c_cb("c_cb", "c_cb", 800, 600);
        RooPlot *frame = z1_mass_diff.frame(Range(-15.,10.));
        frame->SetTitle("");
        frame->GetXaxis()->SetTitle("M_{Z}(p^{rec}_{T}, #eta^{rec}, #phi^{rec}) - M_{Z}(p^{gen}_{T}, #eta^{rec}, #phi^{rec}) [GeV/c^{2}]"); //)/#deltam
        //frame->GetXaxis()->SetTitleSize(1.2);
        frame->GetYaxis()->SetTitleOffset(1.2);
        data_new.plotOn(frame);
        cb.plotOn(frame);
        frame->Draw();
        c_cb.SaveAs(TString("Plots/zmass_resolution/zmass_diff_"+lepName+"_"+detector_location+"_bin_"+std::to_string(x)+".png"));
                                          
        data_new.Delete();

        
    }
    
    /************************** Save results to file **************************/
    
    std::ofstream fitresult;
    fitresult.open(TString("fitresult"+detector_location+"zmass_diff_"+lepName+".txt"));
    
    // Print results
    for (int i=0; i<nCuts; i++) {
        
        fitresult << sigma_CB[i] << std::endl;
        fitresult << sigma_CB_err[i] << std::endl;
        fitresult << sigma_eff[i] << std::endl;
        fitresult << sigma_eff_err[i] << std::endl;
        fitresult << sigma_eff_check[i] << std::endl;
        fitresult << alpha_CB[i] << std::endl;
        fitresult << alpha_CB_err[i] << std::endl;
        fitresult << mean_CB[i] << std::endl;
        fitresult << mean_CB_err[i] << std::endl;
        
    }
    
    fitresult.close();
    
    file->Close();
    
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

