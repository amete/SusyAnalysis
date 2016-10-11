#!/bin/bash
generator="sherpa"
mkdir -p /data/uclhc/uci/user/amete/analysis_n0226_run/EWK2L/figures/${generator};
cd /data/uclhc/uci/user/amete/analysis_n0226_run/EWK2L/figures/${generator};
#plotFigures --regions SR2L-preMT250        --channels em    --variables mT2lep --binValues 30,0,300  --plotLog;
plotFigures --regions CR2L-VVDF,VR2L-VVDF,CR2L-Top        --channels em    --variables ptL1,mT2lep --binValues 20,0,200  --plotLog;
plotFigures --regions CR2L-VVSF,VR2L-VVSF                 --channels sf    --variables ptL1,mT2lep --binValues 20,0,200  --plotLog;
#plotFigures --regions SR2L-MT290,SR2L-MT2120,SR2L-MT2150  --channels sf,em --variables mT2lep --binValues 20,0,200  --plotLog;
##plotFigures --regions SR2L-preMT2  --channels sf,em --variables mT2lep,met,mll,ptL0,ptL1 --binValues 20,0,200  --plotLog;
##plotFigures --regions SR2L-preMT2  --channels sf,em --variables etaL0,etaL1 --binValues 30,-3.0,3.0  --plotLog;
##plotFigures --regions SR2L-preMT2  --channels sf,em --variables phiL0,phiL1 --binValues 35,-3.5,3.5  --plotLog;
##plotFigures --regions CR2L-VVDF,VR2L-VVDF,CR2L-Top                    --channels em    --variables mT2lep,met,mll,ptL0,ptL1 --binValues 20,0,200  --plotLog;
##plotFigures --regions CR2L-VVDF,VR2L-VVDF,CR2L-Top                    --channels em    --variables etaL0,etaL1 --binValues 30,-3.0,3.0 --plotLog;
##plotFigures --regions CR2L-VVDF,VR2L-VVDF,CR2L-Top                    --channels em    --variables phiL0,phiL1 --binValues 35,-3.5,3.5 --plotLog;
##plotFigures --regions CR2L-VVDF,VR2L-VVDF,CR2L-Top                    --channels em    --variables nCentralLJets,nCentralLJets30,nCentralBJets,nForwardJets --binValues 20,0,20 --plotLog;
##plotFigures --regions CR2L-VVSF,VR2L-VVSF                             --channels sf    --variables mT2lep,met,mll,ptL0,ptL1 --binValues 20,0,200  --plotLog;
##plotFigures --regions CR2L-VVSF,VR2L-VVSF                             --channels sf    --variables etaL0,etaL1 --binValues 30,-3.0,3.0 --plotLog;
##plotFigures --regions CR2L-VVSF,VR2L-VVSF                             --channels sf    --variables phiL0,phiL1 --binValues 35,-3.5,3.5 --plotLog;
##plotFigures --regions CR2L-VVSF,VR2L-VVSF                             --channels sf    --variables nCentralLJets,nCentralLJets30,nCentralBJets,nForwardJets --binValues 20,0,20 --plotLog;
cd -;
