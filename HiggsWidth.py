from HiggsAnalysis.CombinedLimit.PhysicsModel import *

### This is the base python class to study the Higgs width

class Higgswidth(PhysicsModel):
    def __init__(self):
        self.mHRange = []
        self.GGsmfixed = False
        self.GGsmval = 1.
        self.poiMap = []
        self.pois = {}
        self.verbose = False
    def setModelBuilder(self, modelBuilder):
        PhysicsModel.setModelBuilder(self,modelBuilder)
        self.modelBuilder.doModelBOnly = False

    def getYieldScale(self,bin,process):
         return 1
            
    def setPhysicsOptions(self,physOptions):
        for po in physOptions:
            if po.startswith("GGsmVal="):   
                self.GGsmval = po.replace("GGsmVal=","")
                print 'G/G_SM set to :', self.GGsmval
            if po == "GGsmfixed":
                print "Will fix CMS_zz4l_GGsm to 1 and float mu"
                self.GGsmfixed = True
            
    def doParametersOfInterest(self):
        """Create POI and other parameters, and define the POI set."""
        if self.GGsmfixed:
            self.modelBuilder.out.var("CMS_zz4l_GGsm").setVal(1.)
            self.modelBuilder.out.var("CMS_zz4l_GGsm").setConstant(True)
            self.modelBuilder.out.var("CMS_zz4l_mu").setVal(0.01,6.)
            self.modelBuilder.out.var("CMS_zz4l_mu").setConstant(False)
            print "Fixing CMS_zz4l_GGsm"
            poi = "CMS_zz4l_mu"
        else:
            self.modelBuilder.out.var("CMS_zz4l_GGsm").setRange(0.0001,50.)
            self.modelBuilder.out.var("CMS_zz4l_GGsm").setVal(float(self.GGsmval))
            self.modelBuilder.out.var("CMS_zz4l_GGsm").setConstant(False)
            self.modelBuilder.out.var("CMS_zz4l_mu").setVal(1.)
            self.modelBuilder.out.var("CMS_zz4l_mu").setConstant(True)
            poi = "CMS_zz4l_GGsm"
        
        self.modelBuilder.doSet("POI",poi)
        
higgswidth = Higgswidth()
