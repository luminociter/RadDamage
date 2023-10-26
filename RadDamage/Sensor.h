/*
* Sensor.h
*
*
*      Author: Gkougkousis Evangelos - Leonidas
*              egkougko@cern.ch
*               CERN-GENEVA
* 
* All atributes of a sensor (like irradiation, doping concentration, leakage currents, etc) will be stored in this class
* Class to handle Sensors
* Content of a sensor
*/

#ifndef __Sensor__
#define __Sensor__

#include "RadBase.h"
#include "Annealing_Constants.h"

class Annealing_constants;
class Sensor : public RadBase
{
private:							

  double m_Nacceptor;
  double m_Nacceptor_zero;
  double m_Nacceptor_stable_constdamage;
  double m_Nacceptor_stable_acceptorremoval;
  double m_NAccept_BenefAnnealing;
  double m_max_long_term;
  double m_NAcceptor_reverseannealing;
  double m_Ndonor;
  double m_Ndonor_zero;
  double m_Ndonor_stable_donorremoval;
  double m_Ndonor_stable_constdamage;

  double m_N_Eff0;
  double Temperature;
  double m_thickness = 200;                 // sensor thickness in um      
  double m_volume;                          // sensor volume  
  double m_area;                            // sensor surface area in cm^2
  const bool use_CMS_paper_alpha0 = false;  // if false: temperature averaging for leakage current calculation a la Moll (theta function for time scale change), 
                                            // if true: use CMS approach of temperature average (attention: constants can only be changed directly in the code in this case)
  double m_donoRemovalFrc;                  // fraction of donors which can be removed by donor removal
  bool m_doConstant;
  // initial donor concentration (right now the code only works for originally n-type bulk material!)
  vector<vector<double>> m_leakage_current_alpha;
  vector<vector<double>> time_history;
  vector<double> G_i;
  vector<double> m_leakage_current;
  vector<double> alpha_vec;
  vector<double> m_powerconsumption;
  vector<double> inverse_weighted_temp;
  RadBase* s_base;

public:
  // Constructors
  Sensor(double a, double b, int c, double d, double e, RadBase* base); //Constructor
  Sensor(RadBase* base);
  virtual ~Sensor(); // default Destructor

  // Set ang Get functions
  RadBase* GetBase() { return s_base; };
  void SetThickness(double thickness = 200) { m_thickness = thickness; };
  double GetThickness() { return m_thickness; };
  void SetArea(double area = 5.95) { m_area = area; };
  void SetConstant(bool doconst = false) { m_doConstant = doconst; };
  double GetArea() { return m_area; };
  double GetNeff() { return m_Nacceptor - m_Ndonor; };
  void SetNacceptorZero(double new_value);
  double GetNacceptor();
  void SetNdonorZero(double new_value);
  double GetNdonor();
  void SetDonoRemovalFrac(double donoremovalFrac);
  void SetNAcceptRevers(double NAcceptRevers);
  void SetTemp(double value);
  double GetTemp();
  vector<double> GetGi(double T);
  vector<double> GetLeakageCurrent(double T);
  vector<double> GetAlphaVec(double T);
  vector<double> GetPowerConsumption(double T);

  void Initialize();
  double DepletV(double doping_conc);
  //Calcualte irradiated step parameters
  void Irradiate(leakage_current_consts leconsts, Annealing_constants* constants, 
                 long double fluence, float duration, long double totalDose);
};

#endif