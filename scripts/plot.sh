#!/bin/bash
generator="sherpa"
mkdir -p /data/uclhc/uci/user/amete/analysis_n0224_run/EWK2L/figures_6/${generator};
cd /data/uclhc/uci/user/amete/analysis_n0224_run/EWK2L/figures_6/${generator};
plotFigures --regions SR2L-preMT2  --channels sf,em --variables mT2lep,met --binValues 20,0,200  --plotLog;
#plotFigures --regions SR2L-MT290,SR2L-MT2120,SR2L-MT2150  --channels sf,em --variables mT2lep --binValues 20,0,200  --plotLog;
#plotFigures --regions SR2L-preMT2,SR2L-MT290,SR2L-MT2120,SR2L-MT2150  --channels sf,em --variables mT2lep,mll,ptL0,ptL1 --binValues 20,0,200  --plotLog;
#plotFigures --regions SR2L-preMT2,SR2L-MT290,SR2L-MT2120,SR2L-MT2150  --channels sf,em --variables nCentralLJets,nCentralLJets30,nCentralBJets,nForwardJets --binValues 20,0,20 --plotLog;
#plotFigures --regions CR2L-VVDF,VR2L-VVDF,CR2L-Top                    --channels em    --variables mT2lep --binValues 20,0,200  --plotLog;
#plotFigures --regions CR2L-VVDF,VR2L-VVDF,CR2L-Top                    --channels em    --variables mT2lep,mll,ptL0,ptL1 --binValues 20,0,200  --plotLog;
#plotFigures --regions CR2L-VVDF,VR2L-VVDF,CR2L-Top                    --channels em    --variables nCentralLJets,nCentralLJets30,nCentralBJets,nForwardJets --binValues 20,0,20 --plotLog;
#plotFigures --regions CR2L-VVSF,VR2L-VVSF                             --channels sf    --variables mT2lep --binValues 20,0,200  --plotLog;
#plotFigures --regions CR2L-VVSF,VR2L-VVSF                             --channels sf    --variables mT2lep,mll,ptL0,ptL1 --binValues 20,0,200  --plotLog;
#plotFigures --regions CR2L-VVSF,VR2L-VVSF                             --channels sf    --variables nCentralLJets,nCentralLJets30,nCentralBJets,nForwardJets --binValues 20,0,20 --plotLog;
cd -;
