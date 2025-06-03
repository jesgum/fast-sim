#!/bin/bash

NEVENTS=10000
CONFIGURATION="-b --configuration json://configuration.json "
TIMEFRAME=300
GENERATOR_CONFIG="--configFile ../external_generator/ini/pythia8_pp.ini"

o2-sim-dpl-eventgen --generator external --nEvents ${NEVENTS} ${CONFIGURATION} ${GENERATOR_CONFIG} --aggregate-timeframe ${TIMEFRAME} \
| o2-sim-mctracks-to-aod ${CONFIGURATION} \
| o2-analysis-onthefly-tracker ${CONFIGURATION} \
| o2-analysis-onthefly-tofpid ${CONFIGURATION} \
| o2-analysis-onthefly-richpid ${CONFIGURATION}