/*
* Anealing_Constants.cxx
*
*
*      Author: Gkougkousis Evangelos - Leonidas
*              egkougko@cern.ch
*               CERN-GENEVA
*/

#include "../RadDamage/Annealing_Constants.h"
#include "../RadDamage/Sensor.h"
#include "TFitResult.h"
#include "TVirtualFitter.h"
#include "TMath.h"

//#if !defined(__CINT__)
//ClassImp(Anealing_Constants);
//#endif
//
//#if !defined(__CLING__)
//ClassImp(Anealing_Constants);
//#endif
//
//#ifdef __CINT__
//#pragma link C++ class Anealing_Constants;
//#endif
//
//#ifdef __ROOTCLING__
//#pragma link C++ class Anealing_Constants;
//#endif

// Constructor
Annealing_constants::Annealing_constants(double a, double b, double c, double d, double e, double f, double g, double h, double i, Sensor* sensor)
{
    m_gA = a;
    m_gY = b;
    m_gC = c;
    m_ka = d;
    m_ky1 = e;
    m_Ea = f;
    m_Ey = g;
    m_cc = h;
    m_fAcRmv = i;
    m_ccErr = 0.0;
    m_gCErr = 0.0;
    m_fAcRmvErr = 0.0;
    c_sensor = sensor;
    c_base = sensor->Sensor::GetBase();
}
// --------------------------------------------------------------------------------------------------------------
Annealing_constants::Annealing_constants(Sensor* sensor, double a, double b, double c, double d, double e, double f, double g, double h, double i)
{
    m_gA = a;
    m_gY = b;
    m_gC = c;
    m_ka = d;
    m_ky1 = e;
    m_Ea = f;
    m_Ey = g;
    m_cc = h;
    m_fAcRmv = i;
    m_ccErr = 0.0;
    m_gCErr = 0.0;
    m_fAcRmvErr = 0.0;
    c_sensor = sensor;
    c_base = sensor->Sensor::GetBase();
}
// --------------------------------------------------------------------------------------------------------------
Annealing_constants::~Annealing_constants()
{
    
}
// --------------------------------------------------------------------------------------------------------------
Double_t Annealing_constants::DepVolt(double* x, double* par)
{
    double Q_e = 1.602176634*1e-19; // in qulomb
    double epsil_0 = 8.8541878128*1e-14;  // F/cm 
    double epsil_Si = 11.68; // relative 
    double thick = c_sensor->Sensor::GetThickness()*1e-4;
    double neff_0 = c_sensor->Sensor::GetNeff();
    double Volt = (Q_e*TMath::Power(thick, 2)/(2*epsil_0*epsil_Si))*(neff_0+par[0]*x[0]-par[1]*neff_0*(1-TMath::Exp(-par[2]*x[0])));
    return Volt;
}
// --------------------------------------------------------------------------------------------------------------
std::pair<double, double> Annealing_constants::EvalParam(TGraphErrors* &FitHist, int mode)
{
    std::vector<DataElement> data;
    if (!(c_base->RadBase::GetProfile(data, mode+1, c_base->RadBase::GetInputLim(mode), c_base->RadBase::GetInDataExt(mode),
          c_base->RadBase::GetInDataDir(mode), c_base->RadBase::GetInDataFile(mode))))
       {
        if (mode == 1) std::cout << __FUNCTION__ << " ERROR: Inporting acceptor removal data failed!" << std::endl; 
        else std::cout << __FUNCTION__ << " ERROR: Inporting depletion voltage evolution data failed!" << std::endl;
        return std::make_pair(-100,-100);
       }

    std::vector<double> vec;
    TGraphErrors *GraphData = new TGraphErrors(data.size());
    for (unsigned int k = 0; k < data.size(); k++)
        {
         vec.push_back((double)(data.at(k)).duration);
         GraphData->SetPoint(k, (double)(data.at(k)).duration, (double)(data.at(k)).dose);
         GraphData->SetPointError(k, (double)(data.at(k)).temperature, (double)(data.at(k)).fluence);
        }

    double wmax = *std::max_element(vec.begin(), vec.end());
    if (wmax != 0)
       {
        double fact1 = pow(10, fabs(floor(log10(fabs(wmax)))));
        wmax = ceil(wmax / fact1) * fact1;
       }
    else wmax = 0;
    double wmin = *std::min_element(vec.begin(), vec.end());
    if (wmin != 0)
       {
        double fact2 = pow(10, fabs(floor(log10(fabs(wmin)))));
        wmin = floor(wmin/fact2)*fact2;
       }
    else wmin = 0;
    if (c_base->RadBase::GetDebug()) std::cout << __FUNCTION__ << " INFO: Limnit estiamtion, min: " << wmin << ",\t max: " << wmax << std::endl;

    if (GraphData->Integral() == 0)
       {
        std::cout << __FUNCTION__ << " WARNING: Integral of histogram is 0 -> will not continue with linear fitting!" << std::endl;
        delete GraphData;
        return std::make_pair(-200, -200);
       }

    TF1* func;
    if (mode == 1) 
       {
        func = new TF1("AccRemoval", "[0]*TMath::Power(x,-[1])", wmin, wmax);
        GraphData->SetTitle("Acceptor Removal Coefficient");
        func->SetParameters(11.752e-7, 0.417);
        func->SetParLimits(0, 0, 1E-5);
        func->SetParLimits(1, 0, 1);
        func->SetParNames("Factor", "Exponent");
       }
    else {
          func = new TF1("Vdep", this, &Annealing_constants::DepVolt, wmin, wmax, 3, "Annealing_constants", "DepVolt");
          GraphData->SetTitle("Re-Introduction Rates");
          func->SetParameters(1.5E-2, 0.8, 3*m_cc);
          func->SetParLimits(0, -0.1, 0.1);
          func->SetParLimits(1, 0.1, 0.99);
          func->SetParLimits(2, 0.8*m_cc, 10*m_cc);
          func->SetParNames("g_{A}", "f", "c_{A}");
         }
    GraphData->SetMarkerColor(kBlack);  GraphData->SetMarkerStyle(20);  GraphData->SetMarkerSize(1.5);
    func->SetLineColor(kRed);           func->SetLineStyle(1);          func->SetLineWidth(2);
    TFitResultPtr myfit = GraphData->Fit(func, "S E M R");

    double chi2 = 0.0;
    std::pair <double, double> par0;
    std::pair <double, double> par1;
    if (!(myfit->IsValid()) || myfit->IsEmpty())
       {
        std::cout << __FUNCTION__ << ": WARNING : No fit result!" << std::endl;
        return std::make_pair(-300, -300);
       }
    else if (myfit->Ndf() == 0)
            {
             std::cout << __FUNCTION__ << ": WARNING : Zero degrees of freedom, fit unreliable!" << std::endl;
             return std::make_pair(-400, -400);
            }
    else {
          chi2 = fabs(1 - (myfit->MinFcnValue() / myfit->Ndf()));
          par0.first = myfit->Parameter(0);
          par0.second = myfit->ParError(0);
          par1.first = myfit->Parameter(1);
          par1.second = myfit->ParError(1);
         }
    
    if (mode == 1) FitHist = (TGraphErrors*)GraphData->Clone("Accep_Removal_Fit");
    else FitHist = (TGraphErrors*)GraphData->Clone("V_Dep_Fit");
    myfit = 0;
    delete func;
    delete GraphData;

    std::pair<double, double> res(0.0, 0.0);
    if (mode == 1) 
       {
        res.first = par0.first*pow(c_sensor->Sensor::GetNeff(), -1*par1.first);
        res.second = sqrt(pow(pow(c_sensor->Sensor::GetNeff(), -1*par1.first)*par0.second, 2) + 
                     pow(par0.first*log10(par1.first)*pow(c_sensor->Sensor::GetNeff(), -1*par1.first)*par1.second, 2));
        return res;
       }
    else {
          res.first = par0.first;
          res.second = par0.second;
          m_fAcRmv = par1.first;
          m_fAcRmvErr = par1.second;
          return res;
         }
}