#!/bin/bash

#inputdir=/afs/cern.ch/work/i/idas/public/SimOut/DeltaPt/Extended2026D83
inputdir=/eos/user/p/psuryade/SimOut/DeltaPt/Extended2026D86_r1
outdir=/eos/user/p/psuryade/SimOut/DeltaPt
pydir=$PWD/ReadSimResult/SimTrackAna/python
for i in `seq 0 9`
do
  echo Processing loop $i with file $inputdir/step1_${i}.root
  if [ -f $outdir/step1.root ] ; then
      rm $outdir/step1.root
  fi
  ln -s $inputdir/step1_${i}.root $outdir/step1.root 
  cmsRun $pydir/CellHitSum_cfg.py -n 6
  mv geantoutput.root $outdir/geantoutput_${i}.root
done

cd $outdir
mkdir geant_r1r
mv geantoutput_*.root geant_r1r/
cd geant_r1r
hadd geantoutput_r1r.root geantoutput_0.root geantoutput_1.root geantoutput_2.root geantoutput_3.root geantoutput_4.root geantoutput_5.root geantoutput_6.root geantoutput_7.root geantoutput_8.root geantoutput_9.root
