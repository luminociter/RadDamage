/*
* Anealing_Constants.cxx
*
*
*      Author: Gkougkousis Evangelos - Leonidas
*              egkougko@cern.ch
*               CERN-GENEVA
*/

#include "../RadDamage/RadDamage.h"

//#if !defined(__CINT__)
//ClassImp(RadDamage);
//#endif
//
//#if !defined(__CLING__)
//ClassImp(RadDamage);
//#endif
//
//#ifdef __CINT__
//#pragma link C++ class RadDamage;
//#endif
//
//#ifdef __ROOTCLING__
//#pragma link C++ class RadDamage;
//#endif

//Constructor definition
RadDamage::RadDamage()
{
    m_base = new RadBase();
    m_sensor = new Sensor(m_base);
    m_AnCo = new Annealing_constants(1.81e-2,          // gA -> +/- 0.14e-2
                                     5.16e-2,          // gY -> +/- 0.09e-2
                                     1.0e-2,           // gC -> from fit
                                     2.4e13,           // ka -> +1.2 / -0.8 e+13
                                     1.5e15,           // ky -> +3.4/-1.1e+15
                                     1.09,             // Ea -> +/- 0.03
                                     1.33,             // Ey -> +/- 0.03
                                     6.4118e-14,       // cc -> from fit
                                     1.0,              // f -> from fit
                                     m_sensor);

    m_LeCo = { 1.23e-17,   // alpha_1
               7.07e-17,   // alpha_0_star
               3.3e-18,    // beta
               1.2e13,     // k01
               m_base->RadBase::GetBandGap(), // E1, 1.11 eV
               1.3 };       // E1_star

    m_CanvList = new TList();
    h_Canvas.clear();
    h_Canvas.resize(2, nullptr);
    m_RefThick = -1.0;
    m_RefAccept = -1.0;
    m_RefNDonor = -1.0;
    Initialize();
}
// --------------------------------------------------------------------------------------------------------------
//Detructor definition
RadDamage::~RadDamage()
{
    delete m_AnCo;
    delete m_sensor;
    delete m_base;
    delete m_CanvList;
}
// --------------------------------------------------------------------------------------------------------------
void RadDamage::Initialize()
{
    m_AlphaVecTRef.clear();
    m_LeakCurrVecTRef.clear();
    m_LeakCurrVecTUsr.clear();
    m_GiVecTRef.clear();
    m_GiVecTUsr.clear();
    m_PowerVecTRef.clear();
    m_PowerVecTUsr.clear();
    m_Vdep_vec.clear();
    m_Neff_vec.clear();
    m_Ndonor_vec.clear();
    m_Nacc_vec.clear();
    m_Temp_vec.clear();
    m_time_vec.clear();
    m_fluence_vec.clear();
    m_Accept = -1.0;
    m_Ndonor = -1.0;
    m_Accepr_Revers = -1.0;
    m_Neut_Revers = -1.0;
}
// --------------------------------------------------------------------------------------------------------------
void RadDamage::SetRefDate(int date) { m_base->RadBase::SetRefDate(date); }
// --------------------------------------------------------------------------------------------------------------
void RadDamage::SetInDataNames(int mode, TString DataDir, TString DataName, TString ext) 
{ 
    m_base->RadBase::SetInDataNames(mode, DataDir, DataName, ext); 
}
// --------------------------------------------------------------------------------------------------------------
void RadDamage::SetOutFileName(TString DataDir, TString DataName) { m_base->RadBase::SetOutFileName(DataDir, DataName); }
// --------------------------------------------------------------------------------------------------------------
void RadDamage::SetBandGap(double BandGap) 
{ 
    m_base->RadBase::SetBandGap(BandGap); 
    m_LeCo.E1 = BandGap;
}
// --------------------------------------------------------------------------------------------------------------
void RadDamage::SetEffBandGap(double EffBandGap) { m_base->RadBase::SetEffBandGap(EffBandGap); }
// --------------------------------------------------------------------------------------------------------------
void RadDamage::SetDebug(bool debug) { m_base->RadBase::SetDebug(debug); }
// --------------------------------------------------------------------------------------------------------------
void RadDamage::SetDoseRateScaling(float DoseRateScaling) { m_base->RadBase::SetDoseRateScaling(DoseRateScaling); }
// --------------------------------------------------------------------------------------------------------------
void RadDamage::SetTRef(double TRef) { m_base->RadBase::SetTRef(TRef); }
// --------------------------------------------------------------------------------------------------------------
void RadDamage::SetTUser(double TUser) { m_base->RadBase::SetTUser(TUser); }
// --------------------------------------------------------------------------------------------------------------
void RadDamage::SetInputLim(int mode, int inputlim) { m_base->RadBase::SetInputLim(mode, inputlim); }
// --------------------------------------------------------------------------------------------------------------
void RadDamage::SetRefStartCoord(double x, double y) { m_base->RadBase::SetRefStartCoord(x, y); }
// --------------------------------------------------------------------------------------------------------------
void RadDamage::SetRefEndCoord(double x, double y) { m_base->RadBase::SetRefEndCoord(x, y); }
// --------------------------------------------------------------------------------------------------------------
void RadDamage::SetInpTime(timeDr InpMin) { m_base->RadBase::SetInpTime(InpMin); }
// --------------------------------------------------------------------------------------------------------------
void RadDamage::SetInpCelc(bool InpCelc) { m_base->RadBase::SetInpCelc(InpCelc); }
// --------------------------------------------------------------------------------------------------------------
void RadDamage::SetDoPlotPDF(bool SetPlotPDF) { m_base->RadBase::SetDoPlotPDF(SetPlotPDF); }
// --------------------------------------------------------------------------------------------------------------
void RadDamage::SetConstant(bool doconst) { m_sensor->Sensor::SetConstant(doconst); }
// --------------------------------------------------------------------------------------------------------------
void RadDamage::SetThickness(double thickness) { m_sensor->Sensor::SetThickness(thickness); }
// --------------------------------------------------------------------------------------------------------------
void RadDamage::SetArea(double area) { m_sensor->Sensor::SetArea(area); }
// --------------------------------------------------------------------------------------------------------------
void RadDamage::SetNacceptor(double new_value) { m_Accept = new_value; }
// --------------------------------------------------------------------------------------------------------------
void RadDamage::SetNdonor(double new_value) { m_Ndonor = new_value; }
// --------------------------------------------------------------------------------------------------------------
void RadDamage::SetDonoRemovalFrac(double donoremovalFrac) { m_sensor->Sensor::SetDonoRemovalFrac(donoremovalFrac); }
// --------------------------------------------------------------------------------------------------------------
void RadDamage::SetNAcceptRevers(double NAcceptRevers) { m_Accepr_Revers = NAcceptRevers; }
// --------------------------------------------------------------------------------------------------------------
void RadDamage::SetAnCo(double gA, double gY, double gC, double ka, double ky, double Ea, double Ey, double cc, double f)
{ 
    m_AnCo->Annealing_constants::Set_gA(gA);
    m_AnCo->Annealing_constants::Set_gY(gY);
    m_AnCo->Annealing_constants::Set_gC(gC);
    m_AnCo->Annealing_constants::Set_ka(ka);
    m_AnCo->Annealing_constants::Set_ky1(ky);
    m_AnCo->Annealing_constants::Set_Ea(Ea);
    m_AnCo->Annealing_constants::Set_Ey(Ey);
    m_AnCo->Annealing_constants::Set_cc(cc);
    m_AnCo->Annealing_constants::Set_fAcRmv(f);
}
// --------------------------------------------------------------------------------------------------------------
void RadDamage::SetLeCo(double alpha1, double alpha0, double beta, double k1, double e1, double e1star)
{
    m_LeCo.alpha_1 = alpha1;
    m_LeCo.alpha_0_star = alpha0;
    m_LeCo.beta = beta;
    m_LeCo.k01 = k1;
    m_LeCo.E1 = e1;
    m_base->RadBase::SetBandGap(e1);
    m_LeCo.E1_star = e1star;
}
// --------------------------------------------------------------------------------------------------------------
bool RadDamage::Calculate()
{
    // Check that there is an input file
    if ((m_base->RadBase::GetInDataName(0)).empty() || m_base->RadBase::GetInDataName(0) == "")
       {
        std::cout << __FUNCTION__ << " ERROR: No input file set!" << std::endl;
        return false;
       }

    // Check that an output file has been set and create directory if not existing
    if ((m_base->RadBase::GetOutFileName()).empty() || m_base->RadBase::GetOutFileName() == "") m_base->SetOutFileName("", "");
    else m_base->RadBase::RecursMkDir((m_base->RadBase::GetOutDirName()).c_str());

    // Estimate the acceptor removal coefficient from fit if reference data provided
    if (m_base->RadBase::GetInDataFile(1) != "" || !(m_base->RadBase::GetInDataFile(1).empty()))
       {
        AccRemov = NULL;
        m_sensor->Sensor::SetNacceptorZero(m_Accept);
        m_sensor->Sensor::SetNdonorZero(m_Ndonor);
        std::pair<double, double> AccRem = m_AnCo->Annealing_constants::EvalParam(AccRemov, 1);
        m_AnCo->Annealing_constants::Set_cc(AccRem.first);
        m_AnCo->Annealing_constants::Set_ccErr(AccRem.second);
        if (RadBase::PrintFitInfo(AccRemov, &h_Canvas.at(0))) m_CanvList->Add(h_Canvas.at(0));
       }

    // Estimate the re-introduction function from fit if data is available
    if (m_base->RadBase::GetInDataFile(2) != "" || !(m_base->RadBase::GetInDataFile(2).empty()))
       {
        Sensor* r_sensor = new Sensor(m_base);
        if (m_RefThick > 0) r_sensor->Sensor::SetThickness(m_RefThick);
        else r_sensor->Sensor::SetThickness(m_sensor->Sensor::GetThickness());
        if (m_RefAccept > 0) r_sensor->Sensor::SetNacceptorZero(m_RefAccept);
        else r_sensor->Sensor::SetNacceptorZero(m_sensor->Sensor::GetNacceptor());
        if (m_RefNDonor > 0) r_sensor->Sensor::SetNdonorZero(m_RefNDonor);
        else r_sensor->Sensor::SetNdonorZero(m_sensor->Sensor::GetNdonor());
        Annealing_constants* r_AnCo = new Annealing_constants(r_sensor);
        // calculate the removal coeeficient for this sensor and set it at the annealing constants
        TGraphErrors* temp = NULL;
        std::pair<double, double> AccRemRef = r_AnCo->Annealing_constants::EvalParam(temp, 1);
        r_AnCo->Annealing_constants::Set_cc(AccRemRef.first);
        r_AnCo->Annealing_constants::Set_ccErr(AccRemRef.second);
        DepVltData = NULL;
        AccRemRef = r_AnCo->Annealing_constants::EvalParam(DepVltData, 2);
        m_AnCo->Annealing_constants::Set_gC(AccRemRef.first);
        m_AnCo->Annealing_constants::Set_gCErr(AccRemRef.second);
        m_AnCo->Annealing_constants::Set_fAcRmv(r_AnCo->Annealing_constants::Get_fAcRmv());
        m_AnCo->Annealing_constants::Set_fAcRmvErr(r_AnCo->Annealing_constants::Get_fAcRmvErr());
        if (RadBase::PrintFitInfo(DepVltData, &h_Canvas.at(1))) m_CanvList->Add(h_Canvas.at(1));
        delete r_AnCo;
        delete r_sensor;
        delete temp;
       }

    // Find out how many files there are to read
    unsigned int files = 0.0;
    std::vector<std::string> filenames;
    if (m_base->RadBase::GetInDataDir(0).empty() || m_base->RadBase::GetInDataDir(0) == "")
       {
        std::cout << __FUNCTION__ << " ERROR: Input data directory not set! Abording..." << std::endl; 
        return false;
       }
    else if (m_base->RadBase::GetInDataFile(0).empty() || m_base->RadBase::GetInDataFile(0) == "")
            {
             files = RadBase::CountFiles(m_base->RadBase::GetInDataDir(0).c_str(), m_base->RadBase::GetInDataExt(0).c_str());
             if (files <= 0 )
                {
                 std::cout << __FUNCTION__ << " ERROR: No imput scenario files detected in directory ";
                 std::cout << m_base->RadBase::GetInDataDir(0) << std::endl;
                 return false;
                }
             filenames = RadBase::ListFileNames(m_base->RadBase::GetInDataDir(0).c_str(), m_base->RadBase::GetInDataExt(0).c_str());
            }
    else {
          files = 1;
          filenames.push_back(m_base->RadBase::GetInDataFile(0).c_str());
         }

    // Initalize the variables and create the output file
    long double time = 0.0;
    unsigned int steps = 0;
    double Tref = m_base->RadBase::GetTRef();
    double Tusr = m_base->RadBase::GetTUser();
    double TrefGr = Tref - 273.15;
    double TusrGr = Tusr - 273.15;
    std::string output_file = (m_base->RadBase::GetOutFileName()).insert((m_base->RadBase::GetOutFileName()).size(), ".root");
    TFile* file = new TFile(output_file.c_str(), "UPDATE");
 
    for (unsigned int q = 0; q < files; q++)
        {
         // Read in the input data profile, a vector containing the individual data elements (duration, temperature, dose rate) is returned
         temp_rad_profile.clear();
         if (!(m_base->RadBase::GetProfile(temp_rad_profile, 1, m_base->RadBase::GetInputLim(0), 
             m_base->RadBase::GetInDataExt(0), m_base->RadBase::GetInDataDir(0), filenames.at(q))))
            {
             std::cout << __FUNCTION__ << " ERROR: Inporting radiation profile failed!" << std::endl; 
             continue;
            }
         steps = temp_rad_profile.size();
         std::cout << __FUNCTION__ << " INFO: Profile succesfully read from file: \"" << filenames.at(q) << "\". Profile legnth: " << steps << std::endl;

         // Reset the sensort object variables
         time = 0.0;
         m_totalDose.push_back(0);
         m_sensor->Sensor::Initialize();
         m_sensor->Sensor::SetNacceptorZero(m_Accept);
         m_sensor->Sensor::SetNdonorZero(m_Ndonor);
         m_sensor->Sensor::SetNAcceptRevers(m_Accepr_Revers);

         // Create the objects in the vectors to store properties
         m_Temp_vec.push_back(std::vector<double>());
         m_Neff_vec.push_back(std::vector<double>());
         m_Vdep_vec.push_back(std::vector<double>());
         m_Ndonor_vec.push_back(std::vector<double>());
         m_Nacc_vec.push_back(std::vector<double>());
         m_totalDose.push_back(0.0);
         m_time_vec.push_back(std::vector<double>());
         m_fluence_vec.push_back(std::vector<double>());

         // Iterate through the profile and irradiate the sensor at each step
         for (unsigned int t = 0; t < steps; t++)
             {
              m_sensor->Sensor::SetTemp((temp_rad_profile.at(t)).temperature);
              if (m_base->RadBase::GetDebug()) std::cout << (temp_rad_profile.at(t)).dose << "\t" << temp_rad_profile.at(t).duration
                                                         << "\t" << m_sensor->Sensor::GetTemp() << "\t" << m_totalDose.at(q) << std::endl;
              m_sensor->Sensor::Irradiate(m_LeCo, m_AnCo, (temp_rad_profile.at(t)).dose, (temp_rad_profile.at(t)).duration, m_totalDose.at(q));
              (m_Temp_vec.at(q)).push_back((temp_rad_profile.at(t)).temperature);
              (m_Neff_vec.at(q)).push_back(fabs(m_sensor->Sensor::GetNeff()));
              (m_Vdep_vec.at(q)).push_back(m_sensor->Sensor::DepletV(fabs(m_sensor->Sensor::GetNeff())));
              (m_Ndonor_vec.at(q)).push_back(m_sensor->Sensor::GetNdonor());
              (m_Nacc_vec.at(q)).push_back(m_sensor->Sensor::GetNacceptor());
              time += (temp_rad_profile.at(t)).duration; 
              m_totalDose.at(q) += (temp_rad_profile.at(t)).dose; // dose (neq/cm2)
              (m_time_vec.at(q)).push_back(time/((double)(24*3600)));  // time vector in days, assuming input time in seconds
              (m_fluence_vec.at(q)).push_back(m_totalDose.at(q));
              m_base->RadBase::ProgressBar(t, steps);
             }
         // Store the objects in vectors
         m_AlphaVecTRef.push_back(m_sensor->Sensor::GetAlphaVec(Tref));
         m_LeakCurrVecTRef.push_back(m_sensor->Sensor::GetLeakageCurrent(Tref));
         m_LeakCurrVecTUsr.push_back(m_sensor->Sensor::GetLeakageCurrent(Tusr));
         m_GiVecTRef.push_back(m_sensor->Sensor::GetGi(Tref));
         m_GiVecTUsr.push_back(m_sensor->Sensor::GetGi(Tusr));
         m_PowerVecTRef.push_back(m_sensor->Sensor::GetPowerConsumption(Tref));
         m_PowerVecTUsr.push_back(m_sensor->Sensor::GetPowerConsumption(Tusr));
        }

    // Data output, plotting and visualisation
    std::cout << __FUNCTION__ << " INFO: Processing finished, writing data..." << std::endl;
    // plots as function of time
    int norm;
    if (n_Normalise) norm = 11;
    else norm = 10;
    std::vector<std::vector<double>> minutes;
    for (unsigned int a = 0; a < m_time_vec.size(); a++)
        {
         minutes.push_back(std::vector<double>());
         for (unsigned int q = 0; q < (m_time_vec.at(a).size()); q++) (minutes.at(a)).push_back((m_time_vec.at(a)).at(q)*(24*60)-1);
        }
    m_base->RadBase::PlotVectors(m_time_vec, m_Neff_vec, "Date", "N_{eff} [1/cm^{3}]", "Effective Doping [N_{eff}]", &filenames, file, norm);
    m_base->RadBase::PlotVectors(m_time_vec, m_Ndonor_vec, "Date", "N_{don.} [1/cm^{3}]", "Donors [N_{don.}]", &filenames, file, norm);
    m_base->RadBase::PlotVectors(m_time_vec, m_Nacc_vec, "Date", "N_{acc.} [1/cm^{3}]", "Acceptors [N_{acc.}]", &filenames, file, norm);
    m_base->RadBase::PlotVectors(m_time_vec, m_Vdep_vec, "Date", "U_{dep} [V]", "Depletion Voltage", &filenames, file, norm);
    m_base->RadBase::PlotVectors(minutes, m_AlphaVecTRef, "Minutes at 60^{o}C", Form("#alpha (%.1f^{o}C) [A/cm]", TrefGr), Form("#alpha Factors @ %.1f^{o}C", TrefGr), &filenames, file, 9);
    m_base->RadBase::PlotVectors(minutes, m_LeakCurrVecTRef, "Minutes at 60^{o}C", Form("I_{leak} (%.1f^{o}C) [A/module]", TrefGr), Form("Leakage Current per module @ %.1f^{o}C", TrefGr), &filenames, file, 9);
    m_base->RadBase::PlotVectors(minutes, m_LeakCurrVecTUsr, "Minutes at 60^{o}C", Form("I_{leak} (%.1f^{o}C) [A/module]", TusrGr), Form("Leakage Current per module @ %.1f^{o}C", TusrGr), &filenames, file, 9);
    m_base->RadBase::PlotVectors(minutes, m_GiVecTRef, "Minutes at 60^{o}C", Form("G_{i} (%.1f^{o}C) [A/cm^{3}]", TrefGr), Form("Leakage Current @ %.1f^{o}C", TrefGr), &filenames, file, 9);
    m_base->RadBase::PlotVectors(minutes, m_GiVecTUsr, "Minutes at 60^{o}C", Form("G_{i} (%.1f^{o}C) [A/cm^{3}]", TusrGr), Form("Leakage Current @ %.1f^{o}C", TusrGr), &filenames, file, 9);
    m_base->RadBase::PlotVectors(m_time_vec, m_PowerVecTRef, "Date",  Form("P (%.1f^{o}C) [W/cm^{2}]", TrefGr), Form("Power Consumption @ %.1f^{o}C", TrefGr), &filenames, file, norm);
    m_base->RadBase::PlotVectors(m_time_vec, m_PowerVecTUsr, "Date", Form("P (%.1f^{o}C) [W/cm^{2}]", TusrGr), Form("Power Consumption @ %.1f^{o}C", TusrGr), &filenames, file, norm);
    m_base->RadBase::PlotVectors(m_time_vec, m_Temp_vec, "Date", "Temperature [K]", "Temperature", &filenames, file, 12);
    m_base->RadBase::PlotVectors(m_time_vec, m_fluence_vec, "Date", "Fluence [n_{eq}/cm^{2}]", "Fluence", &filenames, file, 10);

    file->cd();
    m_CanvList->Write();
    file->Close();
    return true;
}