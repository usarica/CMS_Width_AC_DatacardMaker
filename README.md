# CMS Width/AC Data Card Maker

## Checkout instructions

```
git clone https://github.com/JHUGen/JHUGenMELA.git
(cd JHUGenMELA; git checkout -b from-v232 v2.3.2; ./setup.sh -j;)
git clone https://github.com/MELALabs/MelaAnalytics.git
(cd MelaAnalytics; git checkout -b from-v22 v2.2)
git clone https://github.com/usarica/CMSDataTools.git
git clone https://github.com/usarica/CMS_Width_AC_DatacardMaker.git
```

The dependencies on MELA and MelaAnalytics comes from CMSDataTools.
We use this repository for convenience in Condor submission; it does not affect data card compilation.
