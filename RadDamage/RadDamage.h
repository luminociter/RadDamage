/*
* RadDamage.h
*
*
*      Author: Gkougkousis Evangelos - Leonidas
*              egkougko@cern.ch
*               CERN-GENEVA
*/

#ifndef __RadDamage__
#define __RadDamage__

#include "Sensor.h"

class RadDamage : public RadBase
{
private:
    std::vector<DataElement> temp_rad_profile;
    std::vector<std::vector<double>> m_Vdep_vec;
    std::vector<std::vector<double>> m_Neff_vec;
    std::vector<std::vector<double>> m_Ndonor_vec;
    std::vector<std::vector<double>> m_Nacc_vec;
    std::vector<std::vector<double>> m_Temp_vec;
    std::vector<std::vector<double>> m_time_vec;
    std::vector<std::vector<double>> m_fluence_vec;
    std::vector<std::vector<double>> m_AlphaVecTRef;
    std::vector<std::vector<double>> m_LeakCurrVecTRef;
    std::vector<std::vector<double>> m_LeakCurrVecTUsr;
    std::vector<std::vector<double>> m_GiVecTRef;
    std::vector<std::vector<double>> m_GiVecTUsr;
    std::vector<std::vector<double>> m_PowerVecTRef;
    std::vector<std::vector<double>> m_PowerVecTUsr;
    std::vector<long double> m_totalDose;
    double m_RefThick;
    double m_RefAccept;
    double m_RefNDonor;
    double m_Accept;
    double m_Ndonor;
    double m_Accepr_Revers;
    double m_Neut_Revers;
    bool n_Normalise;
    Sensor* m_sensor;
    RadBase* m_base;
    Annealing_constants* m_AnCo;
    leakage_current_consts m_LeCo;
    TGraphErrors* AccRemov;
    TGraphErrors* DepVltData;
    TList* m_CanvList;
    std::vector<TCanvas* > h_Canvas;

public:
    RadDamage();
    virtual ~RadDamage(); // default Destructor
    bool Calculate();
    void Initialize();

    // Set functions RadBase
    void SetRefDate(int date);
    void SetInDataNames(int mode = 0, TString DataDir = "", TString DataName = "", TString ext = "");
    void SetOutFileName(TString DataDir, TString DataName = "");
    void SetBandGap(double BandGap);
    void SetEffBandGap(double EffBandGap);
    void SetDebug(bool debug = false);
    void SetDoseRateScaling(float DoseRateScaling);
    void SetTRef(double TRef);
    void SetTUser(double TUser);
    void SetInputLim(int mode = 0, int inputlim = 0);
    void SetRefStartCoord(double x = -999.0, double y = -999.0);
    void SetRefEndCoord(double x = -999.0, double y = -999.0);
    void SetInpTime(timeDr InpMin = minutes);
    void SetInpCelc(bool InpCelc = true);
    void SetDoPlotPDF(bool SetPlotPDF = true);
    void SetDoNormalize(bool norm = true) { n_Normalise = norm; };
    // Set functions Sensor
    void SetConstant(bool doconst = false);
    void SetThickness(double thickness = 200);
    void SetArea(double area = 5.95);
    void SetNacceptor(double new_value);
    void SetNdonor(double new_value);
    void SetDonoRemovalFrac(double donoremovalFrac);
    void SetNAcceptRevers(double NAcceptRevers);
    void SetRefThick(double RefThick = -1) { m_RefThick = RefThick; };
    void SetRefAccept(double RefAccept = -1) { m_RefAccept = RefAccept; };
    void SetRefNDonor(double RefNDonor = -1) { m_RefNDonor = RefNDonor; };
    // Set Annealing and Leckage Current constants
    void SetAnCo(double gA, double gY, double gC, double ka, double ky, double Ea, double Ey, double cc, double f);
    void SetLeCo(double alpha1, double alpha0, double beta, double k1, double e1, double e1star);
};

#endif
