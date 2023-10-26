/*
* Anealing_Constants.cxx
*
*
*      Author: Gkougkousis Evangelos - Leonidas
*              egkougko@cern.ch
*               CERN-GENEVA
*/

#include "../RadDamage/Sensor.h"
#include "../RadDamage/Annealing_Constants.h"

//#if !defined(__CINT__)
//ClassImp(Sensor);
//#endif
//
//#if !defined(__CLING__)
//ClassImp(Sensor);
//#endif
//
//#ifdef __CINT__
//#pragma link C++ class Sensor;
//#endif
//
//#ifdef __ROOTCLING__
//#pragma link C++ class Sensor;
//#endif

//Constructor definition
Sensor::Sensor(double a, double b, int c, double d, double e, RadBase* base)
{
  m_Nacceptor_zero = a;
  m_Ndonor_zero = b;
  Temperature = c;
  m_NAccept_BenefAnnealing = d;
  m_donoRemovalFrc = e;
  s_base = base;
  m_area = 5.95;
  m_thickness = 200;
  m_doConstant = false;
  Initialize();
}
// --------------------------------------------------------------------------------------------------------------
Sensor::Sensor(RadBase* base)
{
  m_Nacceptor_zero = 0;
  m_Ndonor_zero = 1.7e12; // IBL, balk concentration in cm-3
  Temperature = 0;
  m_NAccept_BenefAnnealing = 0;
  m_donoRemovalFrc = 0.99;
  s_base = base;
  m_area = 5.95;
  m_thickness = 200;
  m_doConstant = false;
  Initialize();
}
// --------------------------------------------------------------------------------------------------------------
void Sensor::Initialize()
{
  m_leakage_current_alpha.clear();
  time_history.clear();
  G_i.clear();
  m_leakage_current.clear();
  alpha_vec.clear();
  m_powerconsumption.clear();
  inverse_weighted_temp.clear();
  m_donoRemovalFrc = 0.99;
  m_volume = m_area * m_thickness * 1e-4; // this is in centimeters, LHCb sensor volume
  m_N_Eff0 = m_Nacceptor_zero - m_Ndonor_zero;
  m_Nacceptor = m_Nacceptor_zero;
  m_Nacceptor_stable_constdamage = 0;
  m_Nacceptor_stable_acceptorremoval = 0;
  m_NAccept_BenefAnnealing = 0;
  m_max_long_term = 0;
  m_NAcceptor_reverseannealing = 0;
  m_Ndonor = m_Ndonor_zero;
  m_Ndonor_stable_donorremoval = 0;
  m_Ndonor_stable_constdamage = 0;
}
// --------------------------------------------------------------------------------------------------------------
// default Destructor
Sensor::~Sensor()
{

}
// --------------------------------------------------------------------------------------------------------------
double Sensor::GetTemp()
{
  return Temperature;
}
// --------------------------------------------------------------------------------------------------------------
void Sensor::SetTemp(double value)
{
  Temperature = value;
}
// --------------------------------------------------------------------------------------------------------------
double Sensor::GetNacceptor()
{
  return m_Nacceptor;
}
// --------------------------------------------------------------------------------------------------------------
double Sensor::GetNdonor()
{
  return m_Ndonor;
}
// --------------------------------------------------------------------------------------------------------------
void Sensor::SetNacceptorZero(double new_value)
{
  m_Nacceptor_zero = new_value; 
  m_N_Eff0 = m_Nacceptor_zero - m_Ndonor_zero;
  m_Nacceptor = m_Nacceptor_zero;
}
// --------------------------------------------------------------------------------------------------------------
void Sensor::SetNdonorZero(double new_value)
{
  m_Ndonor_zero = new_value;
  m_N_Eff0 = m_Nacceptor_zero - m_Ndonor_zero;
  m_Ndonor = m_Ndonor_zero;
}
// --------------------------------------------------------------------------------------------------------------
void Sensor::SetDonoRemovalFrac(double donoremovalFrac)
{
  if (donoremovalFrac > 0 && donoremovalFrac < 1) donoremovalFrac = m_donoRemovalFrc;
  else {
        std::cout << __FUNCTION__ << " ERROR: Invalid donor removal fraction set, setting to 0.99!" << std::endl;
        m_donoRemovalFrc = 0.99;
       }
}
// --------------------------------------------------------------------------------------------------------------
void Sensor::SetNAcceptRevers(double NAcceptRevers)
{
  if (NAcceptRevers >= 0) m_NAccept_BenefAnnealing = NAcceptRevers;
  else {
        std::cout << __FUNCTION__ << " ERROR: Invalid initial reverseble acceptor number, setting to 0!" << std::endl;
        m_NAccept_BenefAnnealing = 0;
       }
}
// --------------------------------------------------------------------------------------------------------------
std::vector<double> Sensor::GetGi(double T)
{
  if (T != s_base->RadBase::GetTRef())
     {
      std::vector <double> G_imod(G_i.size());
      for (unsigned int k = 0; k < G_i.size(); k++)
          {
           G_imod.at(k) = G_i.at(k) * (pow(T/s_base->RadBase::GetTRef(), 2) *
                          exp(-(s_base->RadBase::GetEffBandGap()/(2*8.6173303e-5))*(1/T - 1/s_base->RadBase::GetTRef())));
          }
      return G_imod;
     }
  else return G_i;
}
// --------------------------------------------------------------------------------------------------------------
std::vector<double> Sensor::GetLeakageCurrent(double T)
{
  if (T != s_base->RadBase::GetTRef()) 
     { 
      std::vector <double> leakage(m_leakage_current.size());
      for (unsigned int k = 0; k < m_leakage_current.size(); k++)
          {
           leakage.at(k) = m_leakage_current.at(k) * (pow(T/s_base->RadBase::GetTRef(), 2) *
                            exp(-(s_base->RadBase::GetEffBandGap()/(2*8.6173303e-5))*(1/T - 1/s_base->RadBase::GetTRef())));
          }
      return leakage;
     }
  else return m_leakage_current;
}
// --------------------------------------------------------------------------------------------------------------
std::vector<double> Sensor::GetAlphaVec(double T)
{
  if (T != s_base->RadBase::GetTRef())
     { 
      std::vector <double> alpha(alpha_vec.size());
      for (unsigned int k = 0; k < alpha_vec.size(); k++)
          {
           alpha.at(k) = alpha_vec.at(k) * (pow(T/s_base->RadBase::GetTRef(), 2) *
                        exp(-(s_base->RadBase::GetEffBandGap()/(2*8.6173303e-5))*(1/T - 1/s_base->RadBase::GetTRef())));
          }
      return alpha;
     }
  else return alpha_vec;
}
// --------------------------------------------------------------------------------------------------------------
std::vector<double> Sensor::GetPowerConsumption(double T)
{
  if (T != s_base->RadBase::GetTRef())
     { 
      std::vector <double> consumpt(m_powerconsumption.size());
      for (unsigned int k = 0; k < m_powerconsumption.size(); k++)
          {
           consumpt.at(k) = m_powerconsumption.at(k) * (pow(T/s_base->RadBase::GetTRef(), 2) *
                            exp(-(s_base->RadBase::GetEffBandGap()/(2*8.6173303e-5))*(1/T - 1/s_base->RadBase::GetTRef())));
          }
      return consumpt;
     }
  else return m_powerconsumption;
}
// --------------------------------------------------------------------------------------------------------------.
// Calculate depletion voltage from given effective doping concentration
double Sensor::DepletV(double doping_conc)
{
    // q/2*epsilon*epsilon_0*Neff*D^2, units are correct
    return (1.6021766208e-13*abs(doping_conc)*pow(m_thickness, 2))/(2*11.68*8.854187817);
}
// --------------------------------------------------------------------------------------------------------------
void Sensor::Irradiate(leakage_current_consts leconsts, Annealing_constants *constants, long double fluence, float duration, long double totalDose)
{
  double limit = 1e-30;

  // calculating the effective doping concentration
  if (s_base->RadBase::GetDebug()) cout << duration << " " << constants->Annealing_constants::Get_gA() << " " << fluence << " " << constants->Annealing_constants::Get_ka(GetTemp())
                                        << " " << m_Nacceptor << " " << constants->Annealing_constants::Get_gC() << " " << constants->Annealing_constants::Get_gY()
                                        << " " << constants->Annealing_constants::Get_ky1(GetTemp()) << " " << m_Ndonor << endl;

  //================================= Effective doping estimation for depletion voltage calculations =================================
  // Constant damage terms
  m_Nacceptor_stable_acceptorremoval = (constants->Annealing_constants::Get_fAcRmv())*fabs(m_N_Eff0)
                                       *(1-TMath::Exp(-constants->Annealing_constants::Get_cc()*(fluence+totalDose)));
  m_Nacceptor_stable_constdamage += (constants->Annealing_constants::Get_gC())*(fluence);
  //--------------------------------------------------------------------------------------------
  // Beneficial annealing - Introduction of acceptors that slwoly decreases
  if (!m_doConstant)
     {
      m_NAccept_BenefAnnealing = (constants->Annealing_constants::Get_gA()*(fluence/duration)/(constants->Annealing_constants::Get_ka(GetTemp())))*
                                 (1-TMath::Exp(-constants->Annealing_constants::Get_ka(GetTemp())*duration)) +
                                 m_NAccept_BenefAnnealing*TMath::Exp(-constants->Annealing_constants::Get_ka(GetTemp())*duration);
     }
  else m_NAccept_BenefAnnealing = (constants->Annealing_constants::Get_gA()*fluence+ m_NAccept_BenefAnnealing)*(TMath::Exp(-constants->Annealing_constants::Get_ka(GetTemp())*duration));
 //--------------------------------------------------------------------------------------------
 // Reverse anneling - acceptor intorduction
  if (!m_doConstant)
     {
      // Estimate availabele number of defects to be annealed
      m_max_long_term = (constants->Annealing_constants::Get_gY()*(fluence/duration)/(constants->Annealing_constants::Get_ky1(GetTemp())))*
                        (1-TMath::Exp(-constants->Annealing_constants::Get_ky1(GetTemp())*duration))+
                        m_max_long_term*TMath::Exp(-constants->Annealing_constants::Get_ky1(GetTemp())*duration);
      // Estimate the number of acceptors to introduce from the annealable defercts that were calculated above
      m_NAcceptor_reverseannealing += (constants->Annealing_constants::Get_gY()*(fluence/duration)/(constants->Annealing_constants::Get_ky1(GetTemp())))*
                                      (constants->Annealing_constants::Get_ky1(GetTemp())*duration+TMath::Exp(-constants->Annealing_constants::Get_ky1(GetTemp())*duration)-1)+
                                      m_max_long_term *(1-TMath::Exp(-constants->Annealing_constants::Get_ky1(GetTemp()) *duration));
     }
  else {
        double temp = (constants->Annealing_constants::Get_gY() * fluence + m_max_long_term) * (1 - 1 / (1 + constants->Annealing_constants::Get_ky1(GetTemp()) * duration));
        m_NAcceptor_reverseannealing += temp;
        m_max_long_term += constants->Annealing_constants::Get_gY() * fluence - temp;
       }
  // Total Dopant concentration
  m_Nacceptor = m_Nacceptor_zero - m_Nacceptor_stable_acceptorremoval + m_Nacceptor_stable_constdamage + m_NAcceptor_reverseannealing + m_NAccept_BenefAnnealing;
  m_Ndonor = m_Ndonor_zero - m_Ndonor_stable_donorremoval + m_Ndonor_stable_constdamage;
  //================================= Leackage current calculations =================================
  vector<double> tmp(3);
  vector<double> tmp2(2);
  double G_i_tmp = 0;   

  if (use_CMS_paper_alpha0)
     {
      double CMS_alpha0_c1 = -8.9e-17; 
      double CMS_alpha0_c2 = 4.6e-14;
      double CMS_beta = 2.9e-18;
      double inverse_weighted_temp_tmp = 0;
      inverse_weighted_temp_tmp = duration / GetTemp();
      inverse_weighted_temp.push_back(inverse_weighted_temp_tmp);
      // needs to be updated at every itaration in the while loop since the averaged temperature needs to be recalculated
      tmp.at(0) = fluence * duration;
      tmp.at(1) = fluence * duration * leconsts.alpha_1;
      tmp.at(2) = fluence * duration * CMS_beta;
      m_leakage_current_alpha.push_back(tmp);
      // time comes in seconds and k01 is seconds-1 so units are fine here
      tmp2.at(0) = duration * leconsts.k01* exp(-leconsts.E1/(8.6173303e-5*GetTemp()));
      // t/60 to convert time from seconds to minutes since this is the proposed unit here
      tmp2.at(1) = duration / 60.0; 
      if (time_history.size() > 0)
         {
          tmp2.at(0) += (time_history.back()).at(0);
          tmp2.at(1) += (time_history.back()).at(1);
         }
      time_history.push_back(tmp2);
      unsigned int i = 0;
      double temperature_average_i = 0;
      while (i < m_leakage_current_alpha.size())
            {
             temperature_average_i = 0;
             // calculate the average inverse temperature weighted with the individual time from the moment of irradiation (i) to now (m_leakage_current_alpha.size())
             for (unsigned int j = i; j < m_leakage_current_alpha.size(); j++) temperature_average_i += inverse_weighted_temp.at(j);
             // for i=0 the time is not the difference but just the total time
             if (i >= 1) 
                {
                 temperature_average_i /= (double)(time_history.back().at(1)-time_history.at(i-1).at(1));
                 G_i_tmp += m_leakage_current_alpha.at(i).at(0) * (CMS_alpha0_c1 + CMS_alpha0_c2 * temperature_average_i) +
                            m_leakage_current_alpha.at(i).at(1) * exp(-(time_history.back().at(0)-time_history.at(i-1).at(0))) - 
                            m_leakage_current_alpha.at(i).at(2) * log(time_history.back().at(1)-time_history.at(i-1).at(1));
                }
             else {
                   temperature_average_i /= (double)(time_history.back().at(1));
                   G_i_tmp += m_leakage_current_alpha.at(i).at(0) * (CMS_alpha0_c1 + CMS_alpha0_c2 * temperature_average_i) +
                              m_leakage_current_alpha.at(i).at(1) * exp(-time_history.back().at(0)) - 
                              m_leakage_current_alpha.at(i).at(2) * log(time_history.back().at(1));
                  }
             i++;
            }
     }
  else {
        // std::cout << "=================================================================" << std::endl;
        // std::cout << "Time History elements: " << time_history.size() << " Leckage current constants " << m_leakage_current_alpha.size() << std::endl;
        double fact = 1.0; // CMS Paper additional factor
        // temperature dependence of alpha0 is in the theta function of the time calculation - high temp = longer times and vice versa
        tmp.at(0) = fluence * leconsts.alpha_0_star;
        // std::cout << fluence << " " << leconsts.alpha_0_star << std::endl;
        tmp.at(1) = fluence * leconsts.alpha_1;
        tmp.at(2) = fluence * leconsts.beta;
        // std::cout << fluence << " " << duration << " " << GetTemp() << std::endl;
        m_leakage_current_alpha.push_back(tmp);
        // k01 is in seconds-1 so time units are fine here
        tmp2.at(0) = duration * leconsts.k01 * exp(-leconsts.E1/((8.6173303e-5)*GetTemp()));
        // T-ref is hardcoded and fixed, it cannot be changed in a trivial way since the other
        // parameters, especially alpha 0 star were evaluated at this reference temperature!!!
        // t/60 to convert time from seconds to minutes since this is the units for the logarithmic term
        tmp2.at(1) = ((double)duration/(double)(60.))*exp((-leconsts.E1_star/(8.6173303e-5))*(((double)(1.0)/(double)GetTemp()) - ((double)1.0/s_base->RadBase::GetTRef())));
        if (time_history.size() > 0) 
           {
            tmp2.at(0) += (time_history.back()).at(0);
            tmp2.at(1) += (time_history.back()).at(1);
           }
        time_history.push_back(tmp2);
        for (unsigned int i =0; i < m_leakage_current_alpha.size(); i++)
            {
             // std::cout << "------------" << std::endl; 
             // std::cout << m_leakage_current_alpha.size() << " " << i << " " << G_i_tmp << " " << G_i_tmp / (totalDose + fluence) << std::endl;
             if (i > 0) 
                {
                 // std::cout << "Here" << std::endl;
                 G_i_tmp += (m_leakage_current_alpha.at(i)).at(0) * fact +
                            (m_leakage_current_alpha.at(i)).at(1) * exp(-((time_history.back()).at(0) - (time_history.at(i-1)).at(0))) -                                  
                            (m_leakage_current_alpha.at(i)).at(2) * log((time_history.back()).at(1) - (time_history.at(i-1)).at(1));
                }
             else {
                   // std::cout << "There" << std::endl;
                   G_i_tmp = (m_leakage_current_alpha.at(i)).at(0) * fact +
                             (m_leakage_current_alpha.at(i)).at(1) * exp(-(time_history.back()).at(0)) -
                             (m_leakage_current_alpha.at(i)).at(2) * log((time_history.back()).at(1));
                  }
             // std::cout << m_leakage_current_alpha.size() << " " << i << " " << G_i_tmp << " " << G_i_tmp/(totalDose+fluence) << std::endl;
            }
       }

  // Populate a, lecakge current and power consumption vectors (this is at T = Tref, might be usefull to do user scaling)
  G_i.push_back(G_i_tmp);
  if ( G_i_tmp/(totalDose+fluence) > 5e-16) alpha_vec.push_back(1e-17);
  else alpha_vec.push_back(G_i_tmp/(totalDose+fluence));
  // std::cout << G_i_tmp / (totalDose + fluence) << " " << G_i_tmp << std::endl;
  m_leakage_current.push_back(G_i_tmp*m_volume); 
  m_powerconsumption.push_back(G_i_tmp*DepletV(fabs(GetNeff())));

  return;
}
// --------------------------------------------------------------------------------------------------------------