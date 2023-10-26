{
  gInterpreter->Reset();
  gInterpreter->LoadMacro("C:/Users/Vagelis/Desktop/VELOU1Radiation/RadDamage/Backup/Root/RadBase.cxx+");
  gInterpreter->LoadMacro("C:/Users/Vagelis/Desktop/VELOU1Radiation/RadDamage/Backup/Root/Annealing_Constants.cxx+");
  gInterpreter->LoadMacro("C:/Users/Vagelis/Desktop/VELOU1Radiation/RadDamage/Backup/Root/Sensor.cxx+");
  gInterpreter->LoadMacro("C:/Users/Vagelis/Desktop/VELOU1Radiation/RadDamage/Backup/Root/RadDamage.cxx+");

  gROOT->ProcessLine("RadDamage* Damage = new RadDamage()");

  // General Options
  gROOT->ProcessLine("Damage->SetDoPlotPDF(false)");                    // Generates PDF printouts of the plots, not yet implemented
  gROOT->ProcessLine("Damage->SetDebug(false)");						// Debug variable
  gROOT->ProcessLine("Damage->SetInputLim(0, 0)");		 				// Ignore lines starting from the end of the input file of each the annealing scenario
  gROOT->ProcessLine("Damage->SetInputLim(1, 0)");						// Ignore lines starting from the end of the input file for acceptor removal fit
  gROOT->ProcessLine("Damage->SetInputLim(2, 0)");						// Ignore lines starting from the end of the input file for re-introduction fit
  gROOT->ProcessLine("Damage->SetInpTime(minutes)");                	// Unit of the time variable at the input file (days, hours, minutes, seconds)
  gROOT->ProcessLine("Damage->SetInpCelc(true)");                   	// Temperature units in Celcius (true) or Kelvin (false)
  // gROOT->ProcessLine("Damage->SetRefStartCoord(5.2, 5.1)");	    	// Inner edge coordinates of VELO on fully closed position
  gROOT->ProcessLine("Damage->SetRefStartCoord(0, 5.1)");				// Inner edge coordinates of VELO on fully closed position
  // gROOT->ProcessLine("Damage->SetRefEndCoord(-37.3, 19.1)");			// Farthest edge coordinates of VELO on fully closed position
  gROOT->ProcessLine("Damage->SetRefDate(20220603)");					// Reference day for x-axis on plots
  gROOT->ProcessLine("Damage->SetDoseRateScaling(1.0)");				// Rescale factor for estimated dose from luminocity
  gROOT->ProcessLine("Damage->SetTRef(21.0)");                      	// Reference temperature for the Hamburg Leackage current model
  gROOT->ProcessLine("Damage->SetTUser(-20.0)");                    	// User temperture for all calcualtions to be returned
  gROOT->ProcessLine("Damage->SetEffBandGap(1.30)");         			// Set the effective bandgap energy for temperature scaling of leackage current (eV)
  gROOT->ProcessLine("Damage->SetBandGap(1.11)");             			// Band-Gap energy for Hamburg leackage current model calculations (eV)
  gROOT->ProcessLine("Damage->SetDoNormalize(false)");
  // Sensor Options
  gROOT->ProcessLine("Damage->SetConstant(false)");                 	// Use constant model for Neff calculations, not full differential equation approximation
  gROOT->ProcessLine("Damage->SetThickness(200)");						// Sensor thickness in um
  gROOT->ProcessLine("Damage->SetArea(5.95)");                     		// Sensor area in cm^2
  gROOT->ProcessLine("Damage->SetNacceptor(1.7e12)");               	// Initial acceptor concentration
  gROOT->ProcessLine("Damage->SetNdonor(0.0)");                    		// Initial donor concentration
  gROOT->ProcessLine("Damage->SetDonoRemovalFrac(0.99)");           	// Donor removel fraction (removable donors/initial donors)
  gROOT->ProcessLine("Damage->SetNAcceptRevers(0.0)");              	// Initial Reversible Acceptors
  // Reference Sensor Oprions
  gROOT->ProcessLine("Damage->SetRefThick(300)");                   	// Thickness of the reference sensor for the introduction rate fit
  gROOT->ProcessLine("Damage->SetRefAccept(1e12)");                 	// Initial acceptor concentration for the reference sensor
  gROOT->ProcessLine("Damage->SetRefNDonor(0)");                    	// Initail donor concentration for the reference sensor

  // Additional Available Options
  // gROOT->ProcessLine("Damage->SetBandGap(double BandGap)");			// Set the BandGap Eeff
  // Set leackage Current constnsts
  // gROOT->ProcessLine("Damage->SetLeCo(double alpha1, double alpha0, double beta, double k1, double e1, double e1star)");
  // Set Annealing Constants
  // gROOT->ProcessLine("Damage->SetAnCo(double gA, double gY, double gC, double ta, double ty, double Ea, double Ey, double cc, double f)");

  gROOT->ProcessLine("Damage->SetInDataNames(0, \"C:/Users/Vagelis/Desktop/VELOU1Radiation/RadDamage/Backup/Scenarios\")");  
  gROOT->ProcessLine("Damage->SetInDataNames(1, \"C:/Users/Vagelis/Desktop/VELOU1Radiation/RadDamage/Backup/Share\", \"AcceptorRemoval\")");
  gROOT->ProcessLine("Damage->SetInDataNames(2, \"C:/Users/Vagelis/Desktop/VELOU1Radiation/RadDamage/Backup/Share\", \"VDepVelo1\")");
  gROOT->ProcessLine("Damage->SetOutFileName(\"C:/Users/Vagelis/Desktop/VELOU1Radiation/RadDamage/Backup/Results\", \"Scenarios\")");
  gROOT->ProcessLine("Damage->Calculate()");
}