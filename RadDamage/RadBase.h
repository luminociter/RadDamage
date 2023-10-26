/*
* RadBase.h
*
*
*      Author: Gkougkousis Evangelos - Leonidas
*              egkougko@cern.ch
*               CERN-GENEVA
*/

#ifndef __RadBase__
#define __RadBase__

#include <iostream>  // Basic IO
#include <fstream>   // Read and write files
#include <sstream>   // to get string into stream
#include <cstdlib>   // convert arguments of function call into int, float,...
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>     // for fabs()
#include <math.h>
#include <time.h>
#include <ctime>
#include "TFile.h"  
#include "TCanvas.h"
#include "TROOT.h"
#include "TMultiGraph.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TVector.h"
#include "TAxis.h"
#include "TFunction.h"
#include "TPaveText.h"
#include "TPad.h"
#include "TText.h"
#include "TGaxis.h"
#include "TDirectory.h"
#include "TSystemDirectory.h"
#include "TF1.h"
#include "TMath.h"
#include "TLine.h"
#include "TH1.h"
#include "TDatime.h"
#include "TLegend.h"

#ifdef _WIN32
#include <Windows4Root.h>
#include <windows.h>
#elif defined __GNUC__
#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h> 
#include <regex>
#endif

// leakage current constants will be stored in this struct
struct leakage_current_consts { double alpha_1; double alpha_0_star; double beta; double k01; double E1; double E1_star; };
struct DataElement { long double duration; long double temperature; long double dose; long double fluence;};
enum timeDr {days, hours, minutes, seconds};

class RadBase
{
private:
    bool m_debug;                          // additional debugging console outputs
    bool m_plot_pdf;                       // plots are by default saved as root file, but can also be directly exported as pdf
    timeDr m_InpMin;                       // Input file in days, hours, min, sec
    bool m_InpCelc;                        // Input file in C or K
    double m_TRef;                         // set a reference temperature for the volume corrected leakage current plot (it will only effect this one plot!) Now: implemented!
    double m_TUser;                        // Usee temperature for plot rescaling
    double m_bandGap;                      // Hamburg leackage current mid band-gap energy
    double m_EffBandGap;                   // Effective mid band gap value in eV used for scaling with temperature
    float m_DoseRateScaling;               // This parameter is multiplied to all fluence rates caclulated from the input profile
    int m_limit_file_size_input[3];        // How many lines should be read from the profile, 0 for everything
    TDatime m_RefDate;                     // set a date for the plots to begin (to be correct it has to be equal to the beginning of your (!) temp/irr profile)   
    std::pair<float, float> m_start_coord; // Coordiates of the bottom right angle of the ladder (the closest to the beam) at fully closed position
    std::pair<float, float> m_end_coord;   // Coordiates of the top left corner of the module at fully closed position
    TString m_InDataName[3];
    TString m_InDataDir[3];
    TString m_InExt[3];
    TString m_ofname;
    TString m_ofdir;

    std::pair<float, float> CalcCoordinates(std::pair<float, float> initial, float dist);
    int CreateDir(const char* path);

public:

    RadBase(); // default constructor
    virtual ~RadBase(); // default Destructor
    bool ProgressBar(Long64_t evt, Long64_t total);
    bool GetProfile(vector<DataElement>& Temperature_Profile_vector, int mode, int limit_file_size_input = 0, 
                                        TString InExt = "", TString InDataDir = "", TString InDataName = "");
    void PlotVectors(std::vector<std::vector<double>> vecx, std::vector<std::vector<double>> vecy,
                     const std::string& xName, const std::string& yName, const std::string& plotname, 
                     std::vector<std::string>* series, TFile* file, int mode);
    double EstFluence(double lumi, std::pair<float, float> start, std::pair<float, float> stop, float z = 1.0);
    int DirExists(const char* path);
    int RecursMkDir(const char* dirname);
    bool PrintFitInfo(TGraphErrors* histo, TCanvas** ca, std::string funcName = "none");
    unsigned int CountFiles(const char* dir, const char* ext = "");
    std::vector<std::string> ListFileNames(const char* path, const char* ext = "");

    // Get and Set functions
    void SetRefDate(int date = 0);
    void SetInDataNames(int mode = 0, TString DataDir = "", TString DataName = "", TString ext = "");
    void SetOutFileName(TString DataDir = "", TString DataName = "");
    std::string GetInDataName(int mode = 0) { return (std::string)(m_InDataName[mode] + m_InDataDir[mode] + m_InExt[mode]); };
    std::string GetOutFileName() { return (std::string)(m_ofdir + m_ofname); };
    std::string GetInDataExt(int mode = 0) { return (std::string)m_InExt[mode]; };
    std::string GetInDataDir(int mode = 0) { return (std::string)m_InDataDir[mode]; };
    std::string GetInDataFile(int mode = 0) { return (std::string)m_InDataName[mode]; };
    std::string GetOutDirName() { return (std::string)m_ofdir; };
    double GetBandGap() { return m_bandGap; };
    void SetBandGap(double BandGap) { m_bandGap = BandGap; };
    double GetEffBandGap() { return m_EffBandGap; };
    void SetEffBandGap(double EffBandGap) { m_EffBandGap = EffBandGap; };
    bool GetDebug() { return m_debug; };
    void SetDebug(bool debug = false) { m_debug = debug; };
    void SetDoPlotPDF(bool doPlotPDF = true) { m_plot_pdf = doPlotPDF; };
    bool GetDoPlotPDF() { return m_plot_pdf; };
    void SetDoseRateScaling(float DoseRateScaling) { m_DoseRateScaling = DoseRateScaling; };
    float GetDoseRateScaling() { return m_DoseRateScaling; };
    void SetTRef(double TRef);
    double GetTRef() { return m_TRef; };
    void SetTUser(double TUser);
    double GetTUser() { return m_TUser; };
    void SetInputLim(int mode, int inputlim = 0) { m_limit_file_size_input[mode] = inputlim; };
    int GetInputLim(int mode) { return m_limit_file_size_input[mode]; };
    bool SetRefStartCoord(double x = -999.0, double y = -999.0);
    bool SetRefEndCoord(double x = -999.0, double y = -999.0);
    std::pair<double, double> GetRefStartCoord() { return m_start_coord; };
    std::pair<double, double> GetRefEndCoord() { return m_end_coord; };
    void SetInpTime(timeDr timeD = hours);
    timeDr GetInpTime() { return m_InpMin; };
    void SetInpCelc(bool InpCelc = true) { m_InpCelc = InpCelc; };
    bool GetInpCelc() { return m_InpCelc; };

    // Remember to put them back to private members
    unsigned int GetGlobalDayFromDate(std::string date);
    Int_t GetDateFromGlobalDay(Int_t date);
};

#endif