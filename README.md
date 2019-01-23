# This is GreensboroCorrelations, a PHENIX analysis code for doing "simple" two-particle correlations

## These are the main analysis files

File name | description
--------- | -----------
GreensboroCorrelations.C | class implementation file, has only the constructor and destructor
GreensboroCorrelations.h | class header file, has all the global level variables, function declarations, etc
GreensboroCorrelationsFunctions.C | class file for various functions to compute the various multi-particle correlations
GreensboroCorrelationsInitializations.C | class file for the initialization functions, particularly the TTree and histograms
GreensboroCorrelationsLinkDef.h | standard ROOT link def file
GreensboroCorrelationsOffsets.C | class file for storing the Q-vector offsets
GreensboroCorrelationsOffsetsRBR.C | class file for storing the run-by-run offsets, not currently used
GreensboroCorrelationsProcessEvent.C | class file for process_event, this is the main part of the analysis code
GreensboroCorrelationsSpecialEventCuts.C | class file for the special event cut that matches the BBC charge and NFVTXtracks

## Additional basic PHENIX universe files

File name | description
--------- | -----------
Makefile.am | standard makefile
autogen.sh | standard autogen file
configure.in | standard configure file

## Macros to test the code

File name | description
--------- | -----------
RunMyMacro.C | from offline/AnalysisTrain/pat/macro
Run_GreensboroCorrelations.C | user defined macro to run the code
Run_GreensboroCorrelationsRun14.C | specific settings for Run14 AuAu/HeAu
Run_GreensboroCorrelationsRun15.C | specific settings for Run15 pAu
Run_GreensboroCorrelationsRun16.C | specific settings for Run14 dAu

## Additional files

File name | description
--------- | -----------
.gitattributes | get the linguistics right

