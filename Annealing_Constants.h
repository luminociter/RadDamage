/*
* Annealing_constants.h
*
*
*      Author: Gkougkousis Evangelos - Leonidas
*              egkougko@cern.ch
*               CERN-GENEVA
*/

#ifndef __Annealing_constants__
#define __Annealing_constants__

#include "RadBase.h"

class Sensor;
class Annealing_constants : public RadBase
{
private:
    double m_gA;
    double m_gY;
    double m_gC;
    double m_gCErr;
    double m_ka;
    double m_ky1;
    double m_Ea;
    double m_Ey;
    double m_cc;
    double m_ccErr;
    double m_fAcRmv;
    double m_fAcRmvErr;
    Sensor* c_sensor;
    RadBase* c_base;

public:
    double DepVolt(double* x, double* par);
    Annealing_constants(Sensor* sensor, double a = 0, double b = 0, double c = 0, double d = 0, 
                        double e = 0, double f = 0, double g = 0, double h = 0, double i = 0);
    Annealing_constants(double a, double b, double c, double d, double e, 
                        double f, double g, double h, double i, Sensor* sensor);
    virtual ~Annealing_constants(); // default Destructor
    std::pair<double, double> EvalParam(TGraphErrors*& FitHist, int mode = 1);

    double Get_gA() { return m_gA; };
    double Get_gY() { return m_gY; };
    double Get_gC() { return m_gC; };
    double Get_gCErr() { return m_gCErr; };
    double Get_ka(double T = 0.0) { return m_ka*exp(-m_Ea/(8.6173303e-5*T)); };
    double Get_ky1(double T = 0.0) { return m_ky1*exp(-m_Ey/(8.6173303e-5*T)); };                     
    double Get_Ea() { return m_Ea; };
    double Get_Ey() { return m_Ey; };
    double Get_cc() { return m_cc; };
    double Get_ccErr() { return m_ccErr; };
    double Get_fAcRmv() { return m_fAcRmv; };
    double Get_fAcRmvErr() { return m_fAcRmvErr; };
    void Set_gA(double gA) { m_gA = gA; };
    void Set_gY(double gY) { m_gY = gY; };
    void Set_gC(double gC) { m_gC = gC; };
    void Set_gCErr(double gCErr) { m_gCErr = gCErr; };
    void Set_ka(double ka) { m_ka = ka; };
    void Set_ky1(double ky1) { m_ky1 = ky1; };
    void Set_Ea(double Ea) { m_Ea = Ea; };
    void Set_Ey(double Ey) { m_Ey = Ey; };
    void Set_cc(double cc) { m_cc = cc; };
    void Set_ccErr(double ccErr) { m_ccErr = ccErr; };
    void Set_fAcRmv(double f) { m_fAcRmv = f; };
    void Set_fAcRmvErr(double fErr) { m_fAcRmvErr = fErr; };
};

#endif
