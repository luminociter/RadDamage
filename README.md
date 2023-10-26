# RadDamage
Calculation of the Leakage Current and Neff with respect to fluence and annealing for silicon sensors

The code is subdivided in several classes written in c++. It is controlled by a steering macro (Run.C), though which all preferences can be adjusted and only requires root to run via the following command:
.x <path_to_steering_macro>/Run.C
 
Several switches and options are available at the Run.C file with in-lined comments explaining each one of them. 
 
The code accepts as input 4 column txt files, containing the duration of each step, the temperature, the distance from the velo closed position (that being references as zero) and the luminosity in fb-1. once you set the folder where the tabulated text scenario files are located, the code will automatically loop over all txt files contained in the location, calculate, and plot the different scenarios, using the name of each txt file as the scenario's name for the legends. The program will start by performing the receptor removal fit using the data from bibliography, located at the share folder. Then it will then continue and perform the depletion voltage fit for the reference sensor, to determine the re-introduction coefficient. Data for this step are also located in the share file and can be updated accordingly. Finlay, once these constants have been evaluated, the program will perform the calculation for each of the available scenario files.
 
Once the code runs, it will perform the acceptor removal fit and process the scenarios. The resulting root file will be created at the results subfolder. By default, it will process all scenario files it finds in the Scenario folder and use the filename to create the legends in the plots.
Scenario files are txt tab separated files with four columns.
- First is the duration of the step (in seconds by default, but can be changed to hours or days, provided the appropriate option for the time units at the Run.C file is correctly adjusted).
- The second column corresponds to the temperature in Celsius, but again it can be changed to kelvin provided that at the Run.C file the unit for temperature is set to kelvin.
- The third column is the distance of the velo from its closed position laterally. This means that at open it is at 30 and at fully closed at 1 because of the presence of the shims. Once the shims are removed it goes to 0. The geometry of the modules is already programmed in as well at the fact that they are not orthogonal to the translation plane.
- The fourth and final colon is the amount of luminosity recorded in each step.

Concerning single pixels or entire modules, the code works as follows:
-	If only the gROOT->ProcessLine("Damage->SetRefStartCoord(0, 5.1)"); line in the run file is set, then the leakage current for a single pixel with these coordinates is calculated.
-	If you both the gROOT->ProcessLine("Damage->SetRefStartCoord(0, 5.1)");  and gROOT->ProcessLine("Damage->SetRefEndCoord(-37.3, 19.1)"); lines are set, then the leakage current for a rectangular module forming a 45* angle with the horizontal displacement axis having the start and stop coordinates specified is calculated. The coordinates are assumed to be with VELO is fully closed (no shims). Negative numbers can be used (to refer to the A or C side) to set the coordinates of the modules of interest, instead of the default which the closest one. In module mode, a trapezoidal integration with numerical assisted uncertainty estimation is implemented to correctly evaluate the mW/cm2 and power dissipation.
-	All calculations assume z = 0. In principle it’s not difficult to integrate the full geometry, meaning generating a 3D plot with the leakage current evolution per pixel for the full velo, and one with the per module, per side and per station evolution. The idea was to do a worst-case scenario and evaluate if extrapolations for different replacement cases actually modify the end working point. The way the framework is written though, I do not see why with minimal coding (<10 lines) this cannot be expanded to a full 3D simulation.

To run on lxplux:
•	Clone the entire contents of RadDamage folder to a local directory
•	source /cvmfs/sft.cern.ch/lcg/releases/ROOT/v6.20.06-3f7fd/x86_64-centos7-gcc9-opt/ROOT-env.sh
•	do: root -l Run.C+
