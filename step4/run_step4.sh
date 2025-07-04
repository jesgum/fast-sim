#!/bin/bash

NEVENTS=10000
CONFIGURATION="-b --configuration json://configuration.json "
TIMEFRAME=300
GENERATOR_CONFIG="--configFile $O2DPG_MC_CONFIG_ROOT/MC/config/ALICE3/ini/xicc_pp.ini"

o2-sim-dpl-eventgen --generator external --nEvents ${NEVENTS} ${CONFIGURATION} ${GENERATOR_CONFIG} --aggregate-timeframe ${TIMEFRAME} \
| o2-sim-mctracks-to-aod ${CONFIGURATION} \
| o2-analysis-onthefly-tracker ${CONFIGURATION} \
| o2-analysis-onthefly-tofpid ${CONFIGURATION} \
| o2-analysis-onthefly-richpid ${CONFIGURATION} \
| o2-analysis-alice3-decaypreselector ${CONFIGURATION} \
| o2-analysis-alice3-multicharm-table ${CONFIGURATION}