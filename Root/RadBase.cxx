/*
* RadBase.cxx
*
*
*      Author: Gkougkousis Evangelos - Leonidas
*              egkougko@cern.ch
*               CERN-GENEVA
*/

#include "../RadDamage/RadBase.h"

//#if !defined(__CINT__)
//ClassImp(RadBase);
//#endif
//
//#if !defined(__CLING__)
//ClassImp(RadBase);
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
RadBase::RadBase()
{
   m_debug = false;
   m_plot_pdf = true;
   m_InpCelc = true;
   m_TRef = 294.15;
   m_TUser = 253.15;
   m_bandGap = 1.11;
   m_EffBandGap = 1.214;
   m_DoseRateScaling = 1.;
   for (unsigned int i = 0; i < 3; i++)
       {
        m_limit_file_size_input[i] = 0;
        m_InDataName[i] = "";
        m_InDataDir[i] = "";
        m_InExt[i] = "";
       }
   m_ofname = "";
   m_ofdir = "";
   m_start_coord = std::make_pair(-999.0, -999.0);
   m_end_coord = std::make_pair(-999.0, -999.0);
   m_InpMin = minutes;
}
// --------------------------------------------------------------------------------------------------------------
RadBase::~RadBase()
{
}
// --------------------------------------------------------------------------------------------------------------
void RadBase::SetInpTime(timeDr timeD)
{
    if (timeD == days || timeD == hours || timeD == minutes || timeD == seconds) m_InpMin = timeD;
    else std::cout << __FUNCTION__ << " ERROR: Invalid optino for inout file time units!!" << std::endl;
}
// --------------------------------------------------------------------------------------------------------------
bool RadBase::SetRefStartCoord(double x, double y)
{
    if (x != -999 && y != -999) m_start_coord = std::make_pair(x, y);
    else {
          std::cout << __FUNCTION__ << " ERROR: Ladder start coordinates not set!!" << std::endl;
          return false;
         }
    return true;
}
// --------------------------------------------------------------------------------------------------------------
bool RadBase::SetRefEndCoord(double x, double y)
{
    if (x != -999 && y != -999) m_end_coord = std::make_pair(x, y);
    else {
          std::cout << __FUNCTION__ << " ERROR: Ladder end coordinates not set!!" << std::endl;
          return false;
         }
    return true;
}
// --------------------------------------------------------------------------------------------------------------
void RadBase::SetInDataNames(int mode, TString DataDir, TString DataName, TString ext)
{
    bool valid = false;
    if(DataDir.Length() > 1 && DataDir != "")
      {
       valid = true;
       TString sLast = DataDir[DataDir.Length() - 1];
       if (!sLast.Contains("/")) DataDir = DataDir + "/";
      }
    if (mode == 0) m_InDataDir[0] = DataDir;
    else { 
          if (!valid) m_InDataDir[mode] = m_InDataDir[0]; 
          else m_InDataDir[mode] = DataDir;
         }

    if (DataName.Length() > 4 && ext == "")
       {
        if (DataName.Contains(".txt")) { m_InDataName[mode] = DataName(0, DataName.Length() - 4); m_InExt[mode] = ".txt"; }
        else if (DataName.Contains(".csv")) { m_InDataName[mode] = DataName(0, DataName.Length() - 4); m_InExt[mode] = ".csv"; }
        else if (ext != "") { m_InDataName[mode] = DataName; m_InExt[mode] = ext;}
        else { m_InDataName[mode] = DataName; m_InExt[mode] = ".txt";}
       }
    else if (DataName != "" && ext != "") { m_InDataName[mode] = DataName; m_InExt[mode] = ext; }
    else { m_InDataName[mode] = ""; m_InExt[mode] = ""; }
}
// --------------------------------------------------------------------------------------------------------------
void RadBase::SetOutFileName(TString DataDir, TString DataName)
{
  if (DataDir == "" && m_InDataDir[0] != "") m_ofdir = m_InDataDir[0];
  else if (DataDir != "") 
          {
           TString sLast = DataDir[DataDir.Length() - 1];
           if (!sLast.Contains("/")) DataDir = DataDir + "/";
           m_ofdir = DataDir;
          }
  if (DataName == "" && m_InDataName[0] != "") m_ofname = m_InDataName[0];
  else if (DataName.Contains(".root")) m_ofname = DataName(0, DataName.Length() - 5);
  else m_ofname = DataName;
}
// --------------------------------------------------------------------------------------------------------------
void RadBase::SetRefDate(int date)
{
    if (date != 0) 
       {
        if (((float)date/100000000) > 1 || date <= 0)
           {
            std::cout << __FUNCTION__ << " ERROR: Incorrect number of digits for reference day, using system date!" << std::endl;
            m_RefDate.Set(date, 0);
           }
        else { 
              Int_t dy = date / 10000;
              Int_t dm = (date - dy * 10000) / 100;
              Int_t dd = (date - dy * 10000 - dm * 100);
              if (dm > 12)
                 {
                  std::cout << __FUNCTION__ << " ERROR: Invalid month for reference date, using system date!" << std::endl;
                  m_RefDate.Set(date, 0);
                 }
              else if (dd > 31)
                      {
                       std::cout << __FUNCTION__ << " ERROR: Invalid day for reference date, using system date!" << std::endl;
                       m_RefDate.Set(date, 0);
                      }
              else m_RefDate.Set(date, 0);
             }
       }
    else m_RefDate.Set();
}
// --------------------------------------------------------------------------------------------------------------
void RadBase::SetTRef(double TRef)
{ 
    if (!m_InpCelc) 
       { 
         if (TRef >= 213.15 && TRef <= 773.15) m_TRef = TRef;
         else {
               m_TRef = 294.15;
               std::cout << __FUNCTION__ << " ERROR: Reference temperature out of supported range, setting to 294.12 K!" << std::endl;
              }
       }
    else {        
          if (TRef >= -60 && TRef <= 500) m_TRef = TRef + 273.15;
          else {
                m_TRef = 294.15;
                std::cout << __FUNCTION__ << " ERROR: Reference temperature out of supported range, setting to 21C!" << std::endl;
               }
         }
}
// --------------------------------------------------------------------------------------------------------------
void RadBase::SetTUser(double TUser)
{
    if (!m_InpCelc)
       { 
        if (TUser >= 213.15 && TUser <= 773.15) m_TUser = TUser;
        else {
              m_TUser = 294.15;
              std::cout << __FUNCTION__ << " ERROR: Reference temperature out of supported range, setting to 294.12 K!" << std::endl;
             }
       }
    else {        
          if (TUser >= -60 && TUser <= 500) m_TUser = TUser + 273.15;
          else {
                m_TUser = 294.15;
                std::cout << __FUNCTION__ << " ERROR: Reference temperature out of supported range, setting to 21C!" << std::endl;
               }
         }
}
// --------------------------------------------------------------------------------------------------------------
// function to read in the temp/rad profile
// reads in a file of format (timestep/temperature/radiation_dose) and gives back a vector 
// of DataElements containing each the timestep/temperature/dose information of one line of the input file
bool RadBase::GetProfile(vector<DataElement>& Temperature_Profile_vector, int mode, int limit_file_size_input,
                         TString InExt, TString InDataDir, TString InDataName)
{
    char ifname[2048];
    static FILE* fd_file = NULL;  // internal
    size_t size = 0;
    size_t lCurPos = 0;
    std::pair<float, float> StrtCrd = std::make_pair(-999., -999.);
    std::pair<float, float> EndCrd = std::make_pair(-999., -999.);

    // Check and create filename
    if (!(InDataDir.IsNull())) strncpy(ifname, InDataDir.Data(), sizeof(ifname));
    if (InDataName.IsNull()) { std::cout << __FUNCTION__ << " ERROR: Input data file name not set! Abording..." << std::endl; return false; }
    strcat(ifname, InDataName.Data());
    strcat(ifname, InExt.Data());

    // Check if file exists
    fd_file = fopen(ifname, "r");
    if (fd_file == NULL)
       {
        std::cout << __FUNCTION__ << " ERROR: Failed to open txt file: " << ifname << std::endl; 
        return false; 
       }

    // Check file size
#ifdef _WIN32
    lCurPos = _ftelli64(fd_file);
    _fseeki64(fd_file, 0, 2); // carefull the fseek third argument is as follows: 0 -> start of file, 1-> current position, 2-> end of file
    size = _ftelli64(fd_file);
    _fseeki64(fd_file, lCurPos, 0); // carefull the fseek third argument is as follows: 0 -> start of file, 1-> current position, 2-> end of file
#else
    lCurPos = ftell(fd_file);
    fseek(fd_file, 0, 2);
    size = ftell(fd_file); // carefull the fseek third argument is as follows: 0 -> start of file, 1-> current position, 2-> end of file
    fseek(fd_file, lCurPos, 0); // carefull the fseek third argument is as follows: 0 -> start of file, 1-> current position, 2-> end of file
#endif
    if (size <= 0) 
       { 
        std::cout << __FUNCTION__ << " ERROR: txt file is zero size, skipping..." << std::endl; 
        return false; 
       }

    unsigned int nel, ntime;
    double lumi, temp;
    float dist;
    double par1, par2, par3, par4;
    DataElement tmp_data;
    char* arr_ptr = &ifname[0];
    while (fgets(ifname, sizeof(ifname), fd_file))
          {
           if (mode == 1) // Read the anealing scenario in
              {
               nel = 0;
               ntime = 0;
               lumi = 0;
               temp = 0;
               dist = 0;
               nel = sscanf(ifname, "%i\t%lg\t%f\t%lg\n", &ntime, &temp, &dist, &lumi);
               if (nel != 4 && strlen(arr_ptr) == 1) continue;
               else if (!(nel == 4 && (ntime > 0 && lumi >= 0 && dist >= 0)))
                       {
                        std::cout << __FUNCTION__ << " ERROR: Invalid parameters from input file, Time: " << ntime 
                                  << " Temperature: " << temp << " Luminocity: " << lumi << " Distace: " << dist << std::endl;
                        return false;
                       }
               tmp_data.fluence = (long double)lumi;
               if (m_InpMin == days) tmp_data.duration = (long double)(ntime*3600*24);
               else if (m_InpMin == hours) tmp_data.duration = (long double)(ntime*3600);
               else if (m_InpMin == minutes) tmp_data.duration = (long double)(ntime*60);
               else tmp_data.duration = (long double)(ntime);
               if (m_InpCelc) tmp_data.temperature = (long double)(temp + 273.15);
               else tmp_data.temperature = (long double)temp;
               StrtCrd = CalcCoordinates(m_start_coord, dist);
               if (StrtCrd.first == -2000 || StrtCrd.first == -1000 || StrtCrd.second == -2000 || StrtCrd.second == -1000)
                  {
                   std::cout << __FUNCTION__ << " WARNING:Invalid start coordiantes, skipping point.." << std::endl;
                   continue;
                  }
                   if (EndCrd.first == -2000 || EndCrd.first == -1000 || EndCrd.second == -2000 || EndCrd.second == -1000)
               if (m_end_coord.first != -999 && m_end_coord.second != 999)
                  {
                   EndCrd = CalcCoordinates(m_end_coord, dist);
                   if (EndCrd.first == -2000 || EndCrd.first == -1000 || EndCrd.second == -2000 || EndCrd.second == -1000)
                      {
                       std::cout << __FUNCTION__ << " WARNING:Invalid end coordiantes, skipping point.." << std::endl;
                       continue;
                      }
                  }
               tmp_data.dose = (long double)(EstFluence(lumi, StrtCrd, EndCrd));
               if (m_debug)
                  {
                   std::cout << ntime << "\t\t" << temp << "\t" << lumi << ":" << dist << std::endl;
                   std::cout << tmp_data.duration << "\t\t" << tmp_data.temperature << "\t" << tmp_data.dose << std::endl;
                  }
              }
            else if (mode == 2 || mode == 3) // Read the acceptor removal dataset in
                    { 
                     par1 = 0; 
                     par2 = 0;
                     par3 = 0;
                     par4 = 0;
                     nel = sscanf(ifname, "%lg\t%lg\t%lg\t%lg\n", &par1, &par2, &par3, &par4);
                     if (nel != 4 && strlen(arr_ptr) == 1) continue;
                     else if (!(nel == 4 && (par1 >= 0 && par2 >= 0 && par3 > 0 && par4 >= 0)))
                             {
                              std::cout << __FUNCTION__ << " ERROR: Invalid parameters from input file,"; 
                              if (mode == 3) std::cout << " Fluence: " << par1;
                              else std::cout << " N_eff(0): " << par1;
                              if (mode == 3) std::cout << ", Fluence uncertainty: " << par2;
                              else std::cout << ", N_eff(0) uncertainty: " << par2;
                              if (mode == 3) std::cout << ", Depletion voltage: " << par3;
                              else std::cout << ", Removal coefficient: " << par3;
                              if (mode == 3) std::cout << ", Depletion voltage uncertainty: " << par4 << std::endl;
                              else std::cout << ", Removal coefficient uncertainty: " << par4 << std::endl;
                              return false;
                             }
                     tmp_data.duration = par1;
                     tmp_data.temperature = par2;
                     tmp_data.dose = par3;
                     tmp_data.fluence = par4;
                    }
            Temperature_Profile_vector.push_back(tmp_data);
            limit_file_size_input--;
            if (limit_file_size_input == 0) break;
          }

  fclose(fd_file);
  return true;
}
// --------------------------------------------------------------------------------------------------------------
// Calculate coordiates of the ladder edge in relative reference assuming 45 degree lateral movement
std::pair<float, float> RadBase::CalcCoordinates(std::pair<float, float> initial, float dist)
{
    if (initial.first != -999 && initial.second != -999)
       { 
        if (dist > 0)
           {
            float x = initial.first - (dist / sqrt(2));
            float y = initial.first + initial.second - x;
            return std::make_pair(x, y);
           }
        else if (dist == 0) return initial;
        else {
              std::cout << __FUNCTION__ << " ERROR: Negative distnace, cannot calculate coordinates!!" << std::endl;
              return std::make_pair(-1000, -1000);
             }        
       }
    else {
          std::cout << __FUNCTION__ << " ERROR: Invalid initial coordiantes!!" << std::endl;
          return std::make_pair(-2000, -2000);
         }
}
// --------------------------------------------------------------------------------------------------------------
// Clculate the avaerage fluence along the ladder using the coordiantes
double RadBase::EstFluence(double lumi, std::pair<float, float> start, std::pair<float, float> stop, float z)
{
    double alpa[7] = { 3.63, 2.72E-04, -3.30E-06, -3.43E-09, 6.36E-12, 1.94E-15, 1.63E-18 };
    double k[7] = { 2.30, -1.22E-04, -5.11E-06, 6.22E-09, 2.79E-11, -6.87E-14, 4.18E-17 };
    double alp = alpa[0] + alpa[1] * z + alpa[2] * pow(z, 2) + alpa[3] * pow(z, 3) + alpa[4] * pow(z, 4) + alpa[5] * pow(z, 5) + alpa[6] * pow(z, 6);
    double kappa = k[0] + k[1]*z + k[2]*pow(z, 2) + k[3]*pow(z, 3) + k[4]*pow(z, 4) + k[5]*pow(z, 5) + k[6]*pow(z, 6);
    // Calculate fluence fror a single pixel
    if (stop.first == -999.0 || stop.second == -999.0)
       {
        if (start.first != -999 || start.second != -999) return m_DoseRateScaling*alp*lumi*(pow(sqrt(pow(start.first*0.1,2)+pow(start.second*0.1,2)),-1*kappa))*1e13;
        else {
              std::cout << __FUNCTION__ << " ERROR: Invalid pixel coordiantes, cannot calculate fluence!!" << std::endl;
              return -1.0;
             }
       }
    // Caase where mean fluence across ladder is to be calculated
    else if (start.first != -999.0 || start.second != -999.0)
            {
             double nPix = 256*256*3;
             double area = (start.first-stop.first)*(stop.second-start.second);
             double dS = area/nPix;
             double ddist = sqrt(dS);
             // Numberical integration calculating fluence for each elementary volume
             double x = start.first;
             double y = start.second;
             double fluence = 0.0;
             while(x > stop.first)
                  {
                   while (y < stop.second)
                         {
                          fluence += alp*lumi*(pow(sqrt(pow(x*0.1,2)+pow(y*0.1,2)),-1*kappa))*1e13;
                          y = y + ddist;
                         }
                   y = start.second;
                   x = x - ddist;
                  }
             return m_DoseRateScaling*(fluence/nPix);
            }
    else {
          std::cout << __FUNCTION__ << " ERROR: Invalid start coordiantes, cannot calculate ladder fluence!!" << std::endl;
          return -1.0;
         }
}
// --------------------------------------------------------------------------------------------------------------
unsigned int RadBase::GetGlobalDayFromDate(std::string date)
{
    // date is in form yyyymmdd
    unsigned int ui_date = atoi(date.c_str());
    unsigned int dy = ui_date / 10000;                     // Extract year
    unsigned int dm = (ui_date - dy * 10000) / 100;        // Extract month
    unsigned int dd = (ui_date - dy * 10000 - dm * 100);   // Extract day of the month
    dm = (dm + 9) % 12;                                    // Shift month to start with 0 for March and add up to 
                                                           // 11 for February to correctly acount for leap years
    unsigned int df = dy;                                  // Calculate the corrections year
    if (dm / 10 > 0) df = df - dm/10;
    // else dm--;
    unsigned int cumulday[12] = {59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 0, 31};
    return (dy-1)*365 + df/4 - df/100 + df/400 + cumulday[dm] + dd;
}
// --------------------------------------------------------------------------------------------------------------
Int_t RadBase::GetDateFromGlobalDay(Int_t day)
{
    int cumulday[12] = {31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365};
    unsigned int y = (unsigned int)(day/365) + 1;
    int d1 = day - ((y-1)*365 + (y-1)/4 - (y-1)/100 + (y-1)/400);
    while (d1 < 0)
          {
           y--;
           d1 = day - ((y-1)*365 + (y-1)/4 - (y-1)/100 + (y-1)/400);
          }
    unsigned int m = 0;
    if (d1 != 0)
       {
        unsigned int y1 = y;
        for (unsigned int k = 0; k < 12; k++)
            {
             if (d1 <= cumulday[k])
                {
                 m = k + 1;
                 if (m < 3) y1 = y - 1;
                 if (k > 0) d1 = day - ((y-1)*365 + y1/4 - y1/100 + y1/400) - cumulday[k-1];
                 else d1 = day - ((y-1)*365 + y1/4 - y1/100 + y1/400);
                 if (d1 == 0) // we are at the last day of the month
                    {
                     unsigned int mday[12] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
                     d1 = mday[m - 1];
                     if (m == 2)  // correction for february
                        {
                         if (y1 % 4) d1++;
                         if (y1 % 100) d1--;
                         if (y1 % 400) d1++;
                        }
                    }
                 break;
                }
            }
       }
    else { m = 12; d1 = 31; y--; } // we are at the last day of the year
    if (GetDebug()) std::cout << __FUNCTION__ << " INFO: " << day << "\t" << d1 
                              << "\t" << y << "\t" << m << "\t" << d1 << std::endl;
    return y * 10000 + m * 100 + d1;
}
// --------------------------------------------------------------------------------------------------------------
// function to plot different informations in root and pdf files
void RadBase::PlotVectors(std::vector<std::vector<double>> vecx, std::vector<std::vector<double>> vecy, const std::string& xName, const std::string &yName,
                          const std::string &plotname, std::vector<std::string>* series, TFile* file, int mode)
{
  std::string output_file_name_pdf = GetOutFileName();
  output_file_name_pdf.insert(output_file_name_pdf.size(),"_");
  output_file_name_pdf.insert(output_file_name_pdf.size(),plotname);
  output_file_name_pdf.insert(output_file_name_pdf.size(),".pdf");

  // Fix temperature period values
  if (mode == 12)
     {
      std::vector<double> tempx;
      std::vector<double> tempy;
      for (unsigned int m = 0; m < vecx.size(); m++)
          {
           tempx.clear();
           tempy.clear();
           for (unsigned int k = 0; k < (vecx.at(m)).size(); k++)
               {
                if (k > 0)  
                   {
                    tempx.push_back((vecx.at(m)).at(k-1));
                    tempy.push_back((vecy.at(m)).at(k));
                   }
                else {
                      tempx.push_back(0);
                      tempy.push_back((vecy.at(m)).at(0));
                     }
                tempy.push_back((vecy.at(m)).at(k));
                tempx.push_back((vecx.at(m)).at(k));
               }
           (vecx.at(m)) = tempx;
           (vecy.at(m)) = tempy;
          }
     }

  // In case of plot as a function of time, convert the days from simulation into a more handy date
  if (mode >= 10) 
     {
      unsigned long int globalday = GetGlobalDayFromDate(std::to_string(m_RefDate.GetDate()));
      TDatime* da = new TDatime();
      double day = 0.0;
      double intday = 0;
      double floatday = 0;
      for (unsigned int m = 0; m < vecx.size(); m++)
          {
           for (unsigned int k = 0; k < (vecx.at(m)).size(); k++)
               {
                day = (vecx.at(m)).at(k) + globalday;
                floatday = modf(day, &intday);
                intday = GetDateFromGlobalDay(intday);
                da->Set(intday, (int)(ceil(floatday*24*3600)));
                (vecx.at(m)).at(k) = da->Convert();
               }
          }
      delete da; 
     }

  // Normalization part for the y axis
  if (mode == 11)
     {
      std::vector<double> maxima;
      for (unsigned int k = 0; k < vecy.size(); k++)
          {
           maxima.push_back(*std::max_element((vecy.at(k)).begin(), (vecy.at(k)).end()));
          }
      double max = *std::max_element(maxima.begin(), maxima.end());
      for (unsigned int k = 0; k < vecy.size(); k++)
          {
           for (unsigned int m = 0; m < (vecy.at(k)).size(); m++)
               {
                (vecy.at(k)).at(m) = (vecy.at(k)).at(m) / max;
               }
          }
      }

  std::vector<TGraphErrors*> grs;
  TMultiGraph* mg = new TMultiGraph();
  mg->SetTitle(Form("%s; %s; %s", plotname.c_str(), xName.c_str(), yName.c_str()));
  grs.resize(vecx.size(), NULL);
  for (unsigned int d = 0; d < vecx.size(); d++)
      {
       grs.at(d) = new TGraphErrors((vecx.at(d)).size(), &(vecx.at(d))[0], &(vecy.at(d))[0], NULL, NULL);
       grs.at(d)->SetTitle((series->at(d)).c_str());
       int k = d%4;
       int m = d/4;
       if (k == 0) grs.at(d)->SetLineColor(kBlue+2*m);
       else if (k == 1)grs.at(d)->SetLineColor(kRed+2*m);
       else if (k == 2)grs.at(d)->SetLineColor(kGreen+2*m);
       else grs.at(d)->SetLineColor(kOrange+2*m);
       grs.at(d)->SetLineWidth(2);
       mg->Add(grs.at(d),"");
      }

  gROOT->SetBatch(kTRUE);
  TGaxis::SetMaxDigits(3);
  gStyle->SetOptStat(kFALSE);
  TCanvas* c1 = new TCanvas(plotname.c_str(), plotname.c_str(), 0, 0, 1024, 668);
  c1->SetFillStyle(4000);    c1->SetFillColor(4000); c1->SetFrameFillColor(4000); c1->SetFrameFillStyle(4000);
  c1->SetFrameBorderMode(0); c1->SetTicks(1, 1);     c1->SetGrid(1, 1);
  c1->cd();
  gPad->SetLeftMargin(0.125); gPad->SetBottomMargin(0.125);
  mg->Draw("APL");
  mg->GetXaxis()->SetTitleOffset(1.25);
  mg->GetYaxis()->SetTitleOffset(1.25);
  // In case of plot as a function of time, convert the days from simulation into a more handy date
  if (mode >= 10) 
     {
      mg->GetXaxis()->SetTimeDisplay(1);
      mg->GetXaxis()->SetNdivisions(6,2,0);
      mg->GetXaxis()->SetTimeFormat("%d/%m/%Y");
      mg->GetXaxis()->SetTimeOffset(0,"gmt");
     }
  gPad->Update();

  TLegend* legend = new TLegend(0.1, 0.7, 0.48, 0.9);
  for (unsigned int d = 0; d < vecx.size(); d++) legend->AddEntry(grs.at(d), (series->at(d)).c_str(), "l");
  legend->Draw();
  c1->Update();
  gROOT->SetBatch(kFALSE);
  file->cd();
  c1->Write();

  if (mode != 66) if (m_plot_pdf) c1->Print(output_file_name_pdf.c_str());

  if (mode == 1) cout << "Final " << plotname << " is: " << (vecy.at((vecy.size()-1))).at((vecy.at((vecy.size() - 1))).size()-1) << " " << yName << endl;
  else cout << "Final " << plotname << " is: " << (vecy.at((vecy.size()-1))).at((vecy.at((vecy.size() - 1))).size()-1) << endl;

  for (int m = grs.size()-1; m > -1; m--) delete grs.at(m);
  delete mg;
  delete c1;
}
// --------------------------------------------------------------------------------------------------------------
bool RadBase::ProgressBar(Long64_t evt, Long64_t total)
{
    float prd = (float)total/(float)100;
    if ((remainder((float)evt, prd) <= 0  && ceil(remainder((float)evt, prd)) == 0) && (evt + 1) < total)
       {
        std::cout << "<" << std::setfill('=') << std::setw(floor((0.2*((float)evt/prd)))) << "";
        std::cout << std::setfill(' ') << std::setw(20 - floor((0.2*((float)evt / prd))) + 2) << std::right << "> :";
        std::cout << std::left << round((float)evt / prd) << "%" << ", Processed entries: " << evt + 1 << " / " << total << "\r" << std::setfill(' ');
        return false;
       }
    else if ((evt + 1) == total) 
            {
             std::cout << "<" << std::setfill('=') << std::setw(19) << "=" << std::setfill(' ') 
                       << "> :" << 100 << "%" << ", Processed entries: "  << evt+1 << " / " 
                       << total << "\r" << std::setfill(' ') << std::endl;
             return true;
            }
    return false;
}
// --------------------------------------------------------------------------------------------------------------
// Function to create directories in linux and windows
int RadBase::RecursMkDir(const char* dirname)
{
    int check = 0;
    const size_t len = strlen(dirname);
#ifdef _WIN32
    char _path[MAX_PATH];
#else
    char _path[PATH_MAX];
#endif
    char *p;

    /* Copy string so its mutable */
    if (len > sizeof(_path)-1) 
       {
        std::cout << __FUNCTION__ << " ERROR: Path too long!" << std::endl;
        return -1;
       }
    strcpy(_path, dirname);

    /* Iterate the string */
    unsigned int t = 0;
    for (p = _path + 1; *p; p++) 
        {
         if (*p == '/') 
            {
             t++;
             if (t == 1) continue;
             *p = '\0'; // Temporarily truncate 
             if (CreateDir(_path) != 1) return -1;
             *p = '/';
            }
        }
    if (CreateDir(_path) != 1) return -1;
    return 1;
}
// --------------------------------------------------------------------------------------------------------------
int RadBase::CreateDir(const char* path)
{
    int check = 0;

    if (DirExists(path) != 1)
       {
#ifdef _WIN32
       check = 0;
       check = CreateDirectory(path, NULL);
#else
       check = mkdir(path, 0777);
       if (check == -1) check = 0;
       else check = 1;
#endif
       if (check == 1 && DirExists(path) != 1)
          {
           std::cout << __FUNCTION__ << " ERROR: Folder creation failed: " << path << "!" << std::endl;
           return -1;
          }
       }
    else check = 1;

    return check;
}
// --------------------------------------------------------------------------------------------------------------
// Function to check if directory exists
int RadBase::DirExists(const char* path)
{
    struct stat info;
    if (stat(path, &info) != 0) return 0;
    else if (info.st_mode & S_IFDIR) return 1;
    else return -1;
}
// --------------------------------------------------------------------------------------------------------------
unsigned int RadBase::CountFiles(const char* dir, const char* ext)
{
    // Count the number of files in the idrectory
    unsigned int nfiles = 0;
    TSystemDirectory* directory = new TSystemDirectory("", dir);
    TList* files = directory->GetListOfFiles();
    unsigned int size = 0;
    while (ext[size] != '\0') size++;
    if (files)
       {
        TSystemFile* sfile;
        TIter next(files);
        TString fname;
        while ((sfile = (TSystemFile*)next()))
              {
               fname = sfile->GetName();
               if (!sfile->IsDirectory())
                  {
                   if (strcmp(ext, "") == 0 || size == 0) nfiles++;
                   else if (fname.EndsWith(ext)) nfiles++;               
                  }
              }
       }
    return nfiles;
}
// --------------------------------------------------------------------------------------------------------------
// Function to return vector of filenames in a directory                                                                                                                        
std::vector<std::string> RadBase::ListFileNames(const char* path, const char* ext)
{
    std::vector<std::string> filenames;
    std::string search_path = path;
    std::string search_ext;
    unsigned int size = 0;
    while (ext[size] != '\0') size++;
    if (size != 0 && strcmp(ext, "") != 0) search_ext = ext;
    else search_ext = "*";
#ifdef _WIN32
    search_path += "*.";
    search_path += search_ext;
    std::cout << "1: " << search_path << std::endl;
    WIN32_FIND_DATA fd;
    HANDLE hFind = ::FindFirstFile(search_path.c_str(), &fd);
    if (hFind != INVALID_HANDLE_VALUE)
       {
        do {
            if (!(fd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)) 
               {
                filenames.push_back(fd.cFileName);
                if (size != 0 && strcmp(ext, "") != 0) filenames.back() = (filenames.back()).substr(0, (filenames.back()).size()-(size+1));
               }
           } 
        while (::FindNextFile(hFind, &fd));
        ::FindClose(hFind);
       }
#else
    struct dirent *pdir;
    DIR *path_dir = opendir(path);
    while ((pdir = readdir(path_dir)))
          {
           search_path = pdir->d_name;
           if ((size != 0 && strcmp(ext, "") != 0) && search_path.find(search_ext) != std::string::npos) filenames.push_back(pdir->d_name);
           else if (strcmp(pdir->d_name, ".") !=0 && strcmp(pdir->d_name, "..") != 0) filenames.push_back(pdir->d_name);
          }
#endif

    return filenames;
}
// --------------------------------------------------------------------------------------------------------------
bool RadBase::PrintFitInfo(TGraphErrors* histo, TCanvas** ca, std::string funcName)
{
  TIter nextfunc(histo->GetListOfFunctions());
  TFunction* func = 0;
  std::vector <std::string> FitNames;
  while ((func = (TFunction*)nextfunc())) FitNames.push_back(func->GetName());
  bool found = false;
  if (m_debug)
     {
      std::cout << __FUNCTION__ << " INFO: Found functions in histogram " << histo->GetName() << ": ";
      for (unsigned int k = 0; k < FitNames.size(); k++) std::cout << FitNames.at(k) << " ";
      std::cout << std::endl;
     }
  if (funcName == "none" && FitNames.size() == 1) { found = true; funcName = FitNames.at(0); }
  else if (funcName != "none" && FitNames.size() > 1)
          {
           for (unsigned int k = 0; k < FitNames.size(); k++)
               {
                if (FitNames.at(k) == funcName) { found = true; break; }
               }
          }
  if (found)
     {
      std::string fitname;
      std::string par0;
      std::string par1;
      if (funcName == "AccRemoval") 
         {
          fitname = "Acceptor Removal Fit";
          par0 = Form("A: %e +/- %e", histo->GetFunction(funcName.c_str())->GetParameter(0), histo->GetFunction(funcName.c_str())->GetParError(0));
          par1 = Form("B: %e +/- %e", histo->GetFunction(funcName.c_str())->GetParameter(1), histo->GetFunction(funcName.c_str())->GetParError(1));
         }
      else if (funcName == "Vdep") 
              {
               fitname = "Depletion Voltage Fit";
               par0 = Form("g_{A}: %e +/- %e", histo->GetFunction(funcName.c_str())->GetParameter(0), histo->GetFunction(funcName.c_str())->GetParError(0));
               par1 = Form("f: % e + / -% e", histo->GetFunction(funcName.c_str())->GetParameter(1), histo->GetFunction(funcName.c_str())->GetParError(1));
              }
      else {
            std::cout << __FUNCTION__ << " WARNING: Unsupported fit type for histogram " << histo->GetName() << ": " << funcName << std::endl;
            return false;
           }    
      gROOT->SetBatch(kTRUE);
      TGaxis::SetMaxDigits(3);
      gStyle->SetOptStat(kFALSE);
      (*ca) = new TCanvas(Form("%s", histo->GetName()), histo->GetTitle(), 0, 0, 600, 600);
      (*ca)->SetFillStyle(4000);    (*ca)->SetFillColor(4000); (*ca)->SetFrameFillColor(4000); (*ca)->SetFrameFillStyle(4000);
      (*ca)->SetFrameBorderMode(0); (*ca)->SetTicks(1, 1);     (*ca)->SetGrid(1, 1);            
      if (funcName == "AccRemoval") { (*ca)->SetLogy(1); (*ca)->SetLogx(1); }
      (*ca)->cd();
      gPad->SetLeftMargin(0.125); gPad->SetBottomMargin(0.125);
      histo->Draw("AP");
      histo->GetXaxis()->SetTitleOffset(1.25);
      if (funcName == "AccRemoval") histo->GetYaxis()->SetTitleOffset(1.75);
      else histo->GetYaxis()->SetTitleOffset(1.25);
      if (funcName == "AccRemoval")
         {
          histo->GetXaxis()->SetTitle("N_{eff}(0) [cm^{-3}]"); 
          histo->GetYaxis()->SetTitle("c_{A} [cm^{-2}]");
         }
      else if (funcName == "Vdep")
              {
               histo->GetXaxis()->SetTitle("#Phi [n_{eq}/cm^{2}]");
               histo->GetYaxis()->SetTitle("Bias Voltage [V]");
              }
      double wmax = TMath::MaxElement(histo->GetN(), histo->GetX());
      if (wmax != 0)
        {
         double fact1 = pow(10, fabs(floor(log10(fabs(wmax)))));
         wmax = ceil(wmax / fact1) * fact1;
        }
      else wmax = 0;
      double wmin = TMath::MinElement(histo->GetN(), histo->GetX());
      if (wmin != 0)
         {
          double fact2 = pow(10, fabs(floor(log10(fabs(wmin)))));
          wmin = floor(wmin/fact2)*fact2;
         }
      else wmin = 0;
      histo->GetXaxis()->SetRangeUser(wmin, wmax);
      gPad->Update();
      TPaveText* paveA = new TPaveText(0.50, 0.60, 0.90, 0.75, "NDCNB");
      paveA->SetTextAlign(11);
      paveA->SetFillStyle(0);
      paveA->SetBorderSize(0);
      paveA->AddText(fitname.c_str());
      ((TText*)paveA->GetListOfLines()->Last())->SetTextColor(1);
      ((TText*)paveA->GetListOfLines()->Last())->SetTextFont(32);
      ((TText*)paveA->GetListOfLines()->Last())->SetTextSize(0.035);
      paveA->AddText(Form("Entires: %lu", (long unsigned int)(histo->GetN())));
      ((TText*)paveA->GetListOfLines()->Last())->SetTextColor(histo->GetLineColor());
      paveA->AddText(par0.c_str());
      ((TText*)paveA->GetListOfLines()->Last())->SetTextColor(histo->GetLineColor());
      paveA->AddText(par1.c_str());
      ((TText*)paveA->GetListOfLines()->Last())->SetTextColor(histo->GetLineColor());
      paveA->Draw();
      (*ca)->Update();
      gROOT->SetBatch(kFALSE);
      return true;
     }
  else {
        std::cout << __FUNCTION__ << " WARNING: Requested function not found in histogram " << histo->GetName() << std::endl;
        return false;
       }
}
// --------------------------------------------------------------------------------------------------------------