# Chi2_optimization_2022




How to run the script

make

./analyzeHGCOctTB <file_list> outFileName.root <dataset> <configuration> <chi2-method> 
<file_list.txt> : contains files to be analyzed (can be found in file_list folder)
<configuration> : alpha (always use alpha for config-1)
<chi2-method> : 0,1,2,3
	      case 0 : rechit energy sum as input to chi2",
	      case 1 : detector scale with No offset for H hadrons as input to chi2",
              case 2 : detector scale with offset (0.4 GeV) for H hadrons as input to chi2",
  	      case 3 : detector scale with offset (0.4 GeV) for H hadrons AND events around the core of beam energy as input to chi2",

alwyas use chi2 method -case-1


Example:
`./analyzeHGCOctTB file_list/skimmed_ntuples_0to1M.txt outFileName.root data alpha 1 
How to run the script for sim:

./analyzeHGCOctTB file_list/v16_v8/<sim_file.txt> outFileName.root sim alpha 1 



Descrition of scripts:

AnalyzeHGCOctTB.cc : Main analysis code
AnalyzeHGCOctTB.h : Initialize histos here
HGCNtupleVariables.h: Tree variable initialization
Others are helping classes
