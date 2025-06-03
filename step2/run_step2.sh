#!/bin/bash

NEVENTS=10000
CONFIGURATION="-b --configuration json://configuration.json "
TIMEFRAME=300


o2-sim-dpl-eventgen --generator pythia8pp --nEvents ${NEVENTS} ${CONFIGURATION} --aggregate-timeframe ${TIMEFRAME} \
| o2-sim-mctracks-to-aod ${CONFIGURATION} \
| o2-analysis-onthefly-tracker ${CONFIGURATION} \
| o2-analysis-alice3-decaypreselector ${CONFIGURATION} \
| o2-analysis-alice3-decayfinder ${CONFIGURATION}