#!/bin/bash

#source /afs/cern.ch/cms/LCG/LCG-2/UI/cms_ui_env.sh
eval `scramv1 runtime -sh`
#source /afs/cern.ch/cms/ccs/wm/scripts/Crab/crab.sh
source /cvmfs/cms.cern.ch/crab3/crab.sh

#source /afs/cern.ch/cms/PPD/PdmV/tools/wmclient/current/etc/wmclient.sh
## source wmclient script without the shared libraries for slc5
source /afs/cern.ch/cms/PPD/PdmV/tools/wmclient/current/etc/wmclient_testful.sh
export PATH=/afs/cern.ch/cms/PPD/PdmV/tools/wmcontrol:${PATH}
export PYTHONPATH=/afs/cern.ch/cms/PPD/PdmV/tools/wmcontrol:${PYTHONPATH}
voms-proxy-init -voms cms
export X509_USER_PROXY=` voms-proxy-info -path`
