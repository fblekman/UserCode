source /afs/cern.ch/cms/sw/cmsset_default.sh

cd /afs/cern.ch/user/f/fblekman/scratch0/metsignif/CMSSW_3_4_2/src/Private/METSignifFirstData/test
cmsenv
nsls Jan10MetSignif/data | grep FREYA | awk '{print "rfio:/castor/cern.ch/user/f/fblekman/Jan10MetSignif/data/"$1}'  > filelist_data.txt
cp filelist_data.txt filelist.txt
root run.C 
mv outputhistos.root outputhistos_data.root
nsls Jan10MetSignif/MC/900 | grep FREYA | awk '{print "rfio:/castor/cern.ch/user/f/fblekman/Jan10MetSignif/MC/900/"$1}'  > filelist_mc.txt
cp filelist_mc.txt filelist.txt
root run.C
mv outputhistos.root outputhistos_mc.root
