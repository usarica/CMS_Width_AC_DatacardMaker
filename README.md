CreateWidthDatacards
====================

Datacards for Higgs width analysis

Just clone the repository version you need (see wiki for milestones/working versions)

./makeEverything.sh \<outputdir\> \<lowmassbound\> \<dimensions\> does everything.  
lowmassbound must be 220 (any other value gives wrong results as bkg ratios are fixed)  
dimension=2: 2D  
dimension=1: 1D(m4l)  
dimension=0: 1D(Dgg)  

Results are saved in outputdir/HCG/220 (with syst) and outputdir/HCG/220_noSyst (wo syst). The macro that produces the plot is utils/plotScan1D.C
