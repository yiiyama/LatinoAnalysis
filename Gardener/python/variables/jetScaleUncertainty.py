#
#
#      |     ___ \  
#      |        ) | 
#      |       __/  
#     _____| _____| 
#                                
#
#


from LatinoAnalysis.Gardener.gardening import TreeCloner
import numpy
import ROOT
import math
import sys
import optparse
import re
import warnings
import os.path
from collections import OrderedDict
from array import array;

class JESTreeMaker(TreeCloner):
    def __init__(self):
       pass

    def help(self):
        return '''Apply id/iso and filter lepton collection'''

    def addOptions(self,parser):
        pass

    def checkOptions(self,opts):
        pass

    def changeOrder(self, vectorname, vector, jetOrderList) :
        # vector is already linked to the otree branch
        # vector name is the "name" of that vector to be modified
        
        #for i in range( len(getattr(self.otree, vectorname)) ) :
          #pass
          #print " --> before ", vectorname, "[", i, "] = ", getattr(self.otree, vectorname)[i]

        # take vector and clone vector
        # equivalent of: temp_vector = itree."vector"
        temp_vector = getattr(self.itree, vectorname)
        # remix the order of vector picking from the clone
        for i in range( len(jetOrderList) ) :
          #print " --> [", i, " :: ", len(jetOrderList) ,"] :::>> ", len(temp_vector), " --> ", jetOrderList[i]      
          # otree."vectorname"[i] = temp_vector[jetOrderList[i]] <--- that is the "itree" in the correct position
          # setattr(self.otree, vector + "[" + str(i) + "]", temp_vector[ jetOrderList[i] ])
          vector.push_back ( temp_vector[ jetOrderList[i] ] )
          #vector.push_back ( 10000. )
        # set the default value for the remaining
        for i in range( len(temp_vector) - len(jetOrderList) ) :
          vector.push_back ( -9999. )
          
        #for i in range( len(getattr(self.otree, vectorname)) ) :
          #pass
          #print " --> after[ " , len(jetOrderList), "] ", vectorname, "[", i, "] = ", getattr(self.otree, vectorname)[i]

    def process(self,**kwargs):
        tree  = kwargs['tree']
        input = kwargs['input']
        output = kwargs['output']
        
        # Make two output directories
        outputSplit = os.path.split(output)
        outputUp = os.path.join( outputSplit[0] + 'Up' , outputSplit[1] )
        if outputUp and not os.path.exists(outputUp):
           os.system('mkdir -p '+outputUp)
        outputDown = os.path.join( outputSplit[0] + 'Down' , outputSplit[1] )
        if outputDown and not os.path.exists(outputDown):
           os.system('mkdir -p '+outputDown)
           
        # does that work so easily and give new variable itree and otree?
        self.connect(tree,input)

        nentries = self.itree.GetEntries()
        print 'Total number of entries: ',nentries 

        #
        # create branches for otree, the ones that will be modified!
        # see: https://root.cern.ch/phpBB3/viewtopic.php?t=12507
        # this is the list of variables to be modified
        #
        self.namesOldBranchesToBeModifiedVector = [
            'std_vector_jet_pt',
            'std_vector_jet_NumberSoftMu',
            'std_vector_jet_bjpb',
            'std_vector_jet_cmva',
            'std_vector_jet_csvv2ivf',
            'std_vector_jet_mass',
            'std_vector_jet_pfcsv',
            'std_vector_jet_softMuEta',
            'std_vector_jet_softMuIso',
            'std_vector_jet_softMuPhi',
            'std_vector_jet_softMuPt',
            #'std_vector_jet_softMuD0',  # from next latino production
            #'std_vector_jet_softMuDz',  # from next latino production
            'std_vector_jet_NumberSoftMu',
            'std_vector_jet_ssvhb',
            'std_vector_jet_ssvhe',
            'std_vector_jet_tche',
            'std_vector_jet_tchp',
            'std_vector_jet_eta',
            'std_vector_jet_QGRmax',
            'std_vector_jet_puid',
            'std_vector_jet_phi',
            'std_vector_jet_QGlikelihood',            
            'std_vector_jet_QGaxis2',
            #'std_vector_jet_HadronFlavour',            
            'std_vector_jet_QGaxis1',
            'std_vector_jet_QGRMScand',
            #'std_vector_jet_PartonFlavour',
           ]
        
        # and these variables NEED to be defined as functions in WWVar.C
        # e.g. mll, dphill, ...
        self.namesOldBranchesToBeModifiedSimpleVariable = [           
           'mjj',
           'njet'
           ]
        
        # jet variables with the structure "std_vector_jet_"NAME to be migrated to "jet"NAME"+number.
        # e.g. jetpt1, jeteta1, jetpt2, jeteta2, ...
        self.jetVariables = [
            'pt',
            'eta',
            'phi',
            'mass',
            'tche'
            ]
        
        self.jetVarList = []
        # maximum number of "single jet" variables to be saved
        maxnjets = 2 # 7 --> everything is available in form of std::vector -> these will be deprecated
        for jetVar in self.jetVariables:
          for i in xrange(maxnjets):
            self.jetVarList.append("jet"+jetVar+str(i+1))

        # clone the tree
        self.clone(outputUp,self.namesOldBranchesToBeModifiedVector + self.namesOldBranchesToBeModifiedSimpleVariable + self.jetVarList)
        self.upTree = self.otree
        self.upFile = self.ofile
        self.clone(outputDown,self.namesOldBranchesToBeModifiedVector + self.namesOldBranchesToBeModifiedSimpleVariable + self.jetVarList)
        self.downTree = self.otree
        self.downFile = self.ofile

        self.oldBranchesToBeModifiedVector = {}
        for bname in self.namesOldBranchesToBeModifiedVector:
          bvector =  ROOT.std.vector(float) ()
          self.oldBranchesToBeModifiedVector[bname] = bvector

        # now actually connect the branches
        for bname, bvector in self.oldBranchesToBeModifiedVector.iteritems():
            self.upTree.Branch(bname,bvector)
            self.downTree.Branch(bname,bvector)

        self.oldBranchesToBeModifiedSimpleVariable = {}
        for bname in self.namesOldBranchesToBeModifiedSimpleVariable:
          bvariable = numpy.ones(1, dtype=numpy.float32)
          self.oldBranchesToBeModifiedSimpleVariable[bname] = bvariable

        # now actually connect the branches
        for bname, bvariable in self.oldBranchesToBeModifiedSimpleVariable.iteritems():
            self.upTree.Branch(bname,bvariable,bname+'/F')
            self.downTree.Branch(bname,bvariable,bname+'/F')

        #self.jetVarDic = OrderedDict()
        self.jetVarDic = {}
        for bname in self.jetVarList:
          bvariable = numpy.ones(1, dtype=numpy.float32)
          self.jetVarDic[bname] = bvariable

        # now actually connect the branches
        for bname, bvariable in self.jetVarDic.iteritems():
            #print " bname   = ", bname
            #print " bvariable = ", bvariable
            self.upTree.Branch(bname,bvariable,bname+'/F')
            self.downTree.Branch(bname,bvariable,bname+'/F')
         
        # input tree  
        itree = self.itree

        # change this part into correct path structure... 
        cmssw_base = os.getenv('CMSSW_BASE')
        try:
            ROOT.gROOT.LoadMacro(cmssw_base+'/src/LatinoAnalysis/Gardener/python/variables/WWVar.C+g')

        except RuntimeError:
            ROOT.gROOT.LoadMacro(cmssw_base+'/src/LatinoAnalysis/Gardener/python/variables/WWVar.C++g')

        # Load jes uncertainty
        jecUnc = ROOT.JetCorrectionUncertainty(os.path.expandvars("${CMSSW_BASE}/src/LatinoAnalysis/Gardener/input/Summer15_25nsV6_MC_Uncertainty_AK4PFchs.txt"))
        #----------------------------------------------------------------------------------------------------
        print '- Starting eventloop'
        step = 5000

        # to be used later on in the code ...
        new_std_vector_jet_pt = ROOT.std.vector(float) ()

        #for i in xrange(10000):
        #for i in xrange(2000):
        for i in xrange(10):
            itree.GetEntry(i)

            if i > 0 and i%step == 0.:
              print i,'events processed :: ', nentries
                
            # scale jet pt
            # Scale Up
            print 'jet pt: ', itree.std_vector_jet_pt[0]
            jetPtUp = []
            for i in range(itree.std_vector_jet_pt.size()):
                if itree.std_vector_jet_pt[i] > 0:
                    jecUnc.setJetEta(itree.std_vector_jet_eta[i])
                    jecUnc.setJetPt(itree.std_vector_jet_pt[i])
                    jetPtUp.append(itree.std_vector_jet_pt[i]*(1+jecUnc.getUncertainty(True)))
                else:
                    break
                
            jetOrderUp = sorted(range(len(jetPtUp)), key=lambda k: jetPtUp[k], reverse=True)
            #print 'jet uncertainty: ', jecUnc.getUncertainty(True)
            print 'jet pt: ', jetPtUp
            print 'jet order: ', jetOrderUp
                           
            for bname, bvector in self.oldBranchesToBeModifiedVector.iteritems():
                bvector.clear()
                if 'jet_pt' in bname:
                    for i in range( len(jetOrderUp) ) :
                        bvector.push_back ( jetPtUp[jetOrderUp[i]] )
                    for i in range( len(getattr(self.itree, bname)) - len(jetOrderUp) ) :
                        bvector.push_back ( -9999. )
                else:
                    self.changeOrder( bname, bvector, jetOrderUp)
            
            jetpt1 = jetPtUp[jetOrderUp[0]]
            jetpt2 = jetPtUp[jetOrderUp[1]]
            jeteta1 = itree.std_vector_jet_eta[jetOrderUp[0]]
            jeteta2 = itree.std_vector_jet_eta[jetOrderUp[1]]
            jetphi1 = itree.std_vector_jet_phi[jetOrderUp[0]]
            jetphi2 = itree.std_vector_jet_phi[jetOrderUp[1]]
            jetmass1 = itree.std_vector_jet_mass[jetOrderUp[0]]
            jetmass2 = itree.std_vector_jet_mass[jetOrderUp[1]]
            WWUp = ROOT.WW(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, jetpt1, jetpt2, jeteta1, jeteta2, jetphi1, jetphi2, jetmass1, jetmass2)
            
            # set the list of jets into the object "WW"
            new_std_vector_jet_pt.clear()
            for iGoodJet in jetOrderUp :
                new_std_vector_jet_pt.push_back(jetPtUp[ iGoodJet ])
            WWUp.setJets(new_std_vector_jet_pt)
        
            # now fill the variables like "mjj", "njets" ...
            for bname, bvariable in self.oldBranchesToBeModifiedSimpleVariable.iteritems():
                bvariable[0] = getattr(WWUp, bname)()   
                
            # refill the single jet variables
            counter = 0
            varCounter = 0
            for bname, bvariable in self.jetVarDic.iteritems():
                if counter < len(jetOrderUp):
                    if 'jetpt' in bname:
                        bvariable[0] = jetPtUp[ jetOrderUp[counter] ]
                    else:
                        bvariable[0] = (getattr(self.itree, 'std_vector_jet_'+self.jetVariables[varCounter] ))[ jetOrderUp[counter] ]
                    counter += 1
                else:
                    bvariable[0] = -9999.
                if counter == maxnjets:
                    varCounter += 1
                    counter = 0

            self.upTree.Fill()
            
            # Scale Down
            jetPtDown = []
            for i in range(itree.std_vector_jet_pt.size()):
                if itree.std_vector_jet_pt[i] > 0:
                    jecUnc.setJetEta(itree.std_vector_jet_eta[i])
                    jecUnc.setJetPt(itree.std_vector_jet_pt[i])
                    jetPtDown.append(itree.std_vector_jet_pt[i]*(1-jecUnc.getUncertainty(False)))
                else:
                     break
            jetOrderDown = sorted(range(len(jetPtDown)), key=lambda k: jetPtDown[k], reverse=True)
            print 'jet pt down: ', jetPtDown          
                    
            for bname, bvector in self.oldBranchesToBeModifiedVector.iteritems():
                bvector.clear()
                if 'jet_pt' in bname:
                    for i in range( len(jetOrderDown) ) :
                        bvector.push_back ( jetPtDown[jetOrderDown[i]] )
                    for i in range( len(getattr(self.itree, bname)) - len(jetOrderDown) ) :
                        bvector.push_back ( -9999. )
                else:
                    self.changeOrder( bname, bvector, jetOrderDown)
            
            jetpt1 = jetPtDown[jetOrderDown[0]]
            jetpt2 = jetPtDown[jetOrderDown[1]]
            jeteta1 = itree.std_vector_jet_eta[jetOrderDown[0]]
            jeteta2 = itree.std_vector_jet_eta[jetOrderDown[1]]
            jetphi1 = itree.std_vector_jet_phi[jetOrderDown[0]]
            jetphi2 = itree.std_vector_jet_phi[jetOrderDown[1]]
            jetmass1 = itree.std_vector_jet_mass[jetOrderDown[0]]
            jetmass2 = itree.std_vector_jet_mass[jetOrderDown[1]]
            WWDown = ROOT.WW(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, jetpt1, jetpt2, jeteta1, jeteta2, jetphi1, jetphi2, jetmass1, jetmass2)
            
            # set the list of jets into the object "WW"
            new_std_vector_jet_pt.clear()
            for iGoodJet in jetOrderDown :
                new_std_vector_jet_pt.push_back(jetPtDown[ iGoodJet ])
            WWDown.setJets(new_std_vector_jet_pt)
        
            # now fill the variables like "mjj", "njets" ...
            for bname, bvariable in self.oldBranchesToBeModifiedSimpleVariable.iteritems():
                bvariable[0] = getattr(WWDown, bname)()  
                
            # refill the single jet variables
            counter = 0
            varCounter = 0
            for bname, bvariable in self.jetVarDic.iteritems():
                if counter < len(jetOrderDown):
                    if 'jetpt' in bname:
                        bvariable[0] = jetPtDown[ jetOrderDown[counter] ]
                    else:
                        bvariable[0] = (getattr(self.itree, 'std_vector_jet_'+self.jetVariables[varCounter] ))[ jetOrderDown[counter] ]
                    counter += 1
                else:
                    bvariable[0] = -9999.
                if counter == maxnjets:
                    varCounter += 1
                    counter = 0
                        
            self.downTree.Fill()

        self.otree = self.upTree
        self.ofile = self.upFile
        self.disconnect(True,False)
        self.otree = self.downTree
        self.ofile = self.downFile
        self.disconnect(True,True)
        print '- Eventloop completed'

