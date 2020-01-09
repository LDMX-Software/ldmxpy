
from __future__ import division

from numpy import linalg as la
from rootpy.tree import Tree

import AnalysisUtils as au

from EventModels import EcalEvent


class EcalAnalysis(object): 

    def __init__(self): 
        self.tree = None
    
    def initialize(self, params):

        self.tree = Tree('ecal_ntuple', model=EcalEvent)

    def process(self, event):
    
        #
        # Hit Level
        #
        
        # Get the collection of digitized Ecal hits
        ecal_dhits = event.get_collection('ecalDigis_recon')

        # counters
        layer_sum = 0
        max_layer_hit = 0
        max_cell_energy = 0
        ecal_hit_x = []
        ecal_hit_y = []

        for ecal_dhit in ecal_dhits: 
        
            # Get the digitized hit energy
            denergy = ecal_dhit.getEnergy()
            
            # Currently, if the energy of a digitized hit is below threshold,
            # its energy is set to zero.  These hits aren't considered.
            if denergy == 0: continue 
            self.tree.ecal_dhit_count += 1

            layer = ecal_dhit.getLayer()
            self.tree.ecal_dhit_layer.push_back(layer)
            layer_sum += layer
            max_layer_hit = max(max_layer_hit, layer)

            self.tree.total_ecal_denergy += denergy
            self.tree.ecal_dhit_energy.push_back(denergy)
            max_cell_energy = max(max_cell_energy, denergy)

            if layer == 1:
                self.tree.ecal_layer1_hit_count += 1
                self.tree.ecal_layer1_energy_sum += denergy

            if layer <= 20: 
                self.tree.trigger_energy_sum += denergy
        
        if self.tree.ecal_dhit_count > 0: 
            self.tree.average_ecal_layer_hit = layer_sum/self.tree.ecal_dhit_count
        
        self.tree.ecal_max_layer_hit = max_layer_hit
        self.tree.ecal_max_denergy_cell = max_cell_energy
        
        # Check if the ECal veto collection exist. If it does, use it to get the
        # BDT result. 
        if event.collection_exist('EcalVeto_recon'): 
       
            # Get the collection of BDT results from the event.
            ecal_veto_results = event.get_collection('EcalVeto_recon')

            self.tree.vecal_dhit_count        = ecal_veto_results[0].getNReadoutHits()
            self.tree.vtotal_ecal_denergy     = ecal_veto_results[0].getSummedDet()
            self.tree.vecal_summed_tight_iso  = ecal_veto_results[0].getSummedTightIso()
            self.tree.vecal_max_denergy_cell  = ecal_veto_results[0].getMaxCellDep()
            self.tree.vecal_shower_rms        = ecal_veto_results[0].getShowerRMS()
            self.tree.vecal_x_pos_std         = ecal_veto_results[0].getXStd()
            self.tree.vecal_y_pos_std         = ecal_veto_results[0].getYStd()
            self.tree.vaverage_ecal_layer_hit = ecal_veto_results[0].getAvgLayerHit()
            self.tree.vecal_max_layer_hit     = ecal_veto_results[0].getDeepestLayerHit()
            self.tree.vecal_layer_std         = ecal_veto_results[0].getStdLayerHit()
 
            # Get the BDT probability.
            self.tree.bdt_prob = ecal_veto_results[0].getDisc()
            
            # Check if the veto passed.
            if ecal_veto_results[0].passesVeto(): self.tree.passes_ecal_veto =  1

        # Check if the ECal veto collection exist. If it does, use it to get the
        # BDT result. 
        if event.collection_exist('MultiElectronVeto_recon'):
            
            print "Here"

            me_veto_results = event.get_collection('MultiElectronVeto_recon')
            
            if me_veto_results[0].cylinder_0_1_layer_0_0.size() != 0:
                for i in xrange(0, me_veto_results[0].cylinder_0_1_layer_0_0.size()):
                    print me_veto_results[0].cylinder_0_1_layer_0_0[i]
                    setattr(self.tree, "cyl%s_0_1_layer_0_0"  % i, me_veto_results[0].cylinder_0_1_layer_0_0[i])
                    setattr(self.tree, "cyl%s_0_1_layer_1_2"  % i, me_veto_results[0].cylinder_0_1_layer_1_2[i])
                    setattr(self.tree, "cyl%s_0_1_layer_3_6"  % i, me_veto_results[0].cylinder_0_1_layer_3_6[i])
                    setattr(self.tree, "cyl%s_0_1_layer_7_14" % i, me_veto_results[0].cylinder_0_1_layer_7_14[i])
                    setattr(self.tree, "cyl%s_0_1_layer_15"   % i, me_veto_results[0].cylinder_0_1_layer_15[i])

                    setattr(self.tree, "cyl%s_1_3_layer_0_0"  % i, me_veto_results[0].cylinder_1_3_layer_0_0[i])
                    setattr(self.tree, "cyl%s_1_3_layer_1_2"  % i, me_veto_results[0].cylinder_1_3_layer_1_2[i])
                    setattr(self.tree, "cyl%s_1_3_layer_3_6"  % i, me_veto_results[0].cylinder_1_3_layer_3_6[i])
                    setattr(self.tree, "cyl%s_1_3_layer_7_14" % i, me_veto_results[0].cylinder_1_3_layer_7_14[i])
                    setattr(self.tree, "cyl%s_1_3_layer_15"   % i, me_veto_results[0].cylinder_1_3_layer_15[i])

                    setattr(self.tree, "cyl%s_3_5_layer_0_0"  % i, me_veto_results[0].cylinder_3_5_layer_0_0[i])
                    setattr(self.tree, "cyl%s_3_5_layer_1_2"  % i, me_veto_results[0].cylinder_3_5_layer_1_2[i])
                    setattr(self.tree, "cyl%s_3_5_layer_3_6"  % i, me_veto_results[0].cylinder_3_5_layer_3_6[i])
                    setattr(self.tree, "cyl%s_3_5_layer_7_14" % i, me_veto_results[0].cylinder_3_5_layer_7_14[i])
                    setattr(self.tree, "cyl%s_3_5_layer_15"   % i, me_veto_results[0].cylinder_3_5_layer_15[i])
                
                    setattr(self.tree, "cyl%s_5_layer_0_0"  % i, me_veto_results[0].cylinder_5_layer_0_0[i])
                    setattr(self.tree, "cyl%s_5_layer_1_2"  % i, me_veto_results[0].cylinder_5_layer_1_2[i])
                    setattr(self.tree, "cyl%s_5_layer_3_6"  % i, me_veto_results[0].cylinder_5_layer_3_6[i])
                    setattr(self.tree, "cyl%s_5_layer_7_14" % i, me_veto_results[0].cylinder_5_layer_7_14[i])
                    setattr(self.tree, "cyl%s_5_layer_15"   % i, me_veto_results[0].cylinder_5_layer_15[i])
    
        
        self.tree.fill(reset=True)

    def finalize(self): 

        self.tree.write()

    def get_noise_rms(self):
        return (900 + 22*27.56)*(0.130/33000)
