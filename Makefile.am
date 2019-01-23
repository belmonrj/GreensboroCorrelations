AUTOMAKE_OPTIONS = foreign

INCLUDES = \
	-I$(includedir) \
	-I$(OFFLINE_MAIN)/include \
	-I$(ROOTSYS)/include

lib_LTLIBRARIES = \
  libGreensboroCorrelations.la

noinst_HEADERS = \
  GreensboroCorrelations.h \
  GreensboroCorrelationsLinkDef.h

libGreensboroCorrelations_la_SOURCES = \
  GreensboroCorrelations.C \
  GreensboroCorrelationsInitializations.C \
  GreensboroCorrelationsProcessEvent.C \
  GreensboroCorrelationsFunctions.C \
  GreensboroCorrelationsSpecialEventCuts.C \
  GreensboroCorrelations_Dict.C

libGreensboroCorrelations_la_LIBADD = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -lSubsysReco \
  -lTOAD

# Rule for generating CINT dictionaries from class headers.
GreensboroCorrelations_Dict.C: \
  GreensboroCorrelations.h \
  GreensboroCorrelationsLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

clean-local:
	rm -f *Dict*
