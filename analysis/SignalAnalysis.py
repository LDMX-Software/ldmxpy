
from __future__ import division

import math

from numpy import linalg as la
from rootpy.tree import Tree

import AnalysisUtils as au
from EventModels import SignalEvent

class SignalAnalysis:

    def __init__(self): 
        self.tree = None

    def initialize(self):
        self.tree = Tree('signal_ntuple', model=SignalEvent)

    def process(self, event):
       
        # Get the collection of MC particles from the event
        particles = event.get_collection('SimParticles_sim')
       
        # Search the list of sim particles for the recoil electron. If it isn't 
        # found, throw an exception
        recoils = au.get_recoil_electrons(particles)
        self.tree.n_electrons = len(recoils)

        # Calculate the e- recoil truth momentum
        for i, recoil in enumerate(recoils): 
            recoil_e_tpvec  = recoil.getMomentum()
            recoil_e_tp     = la.norm(recoil_e_tpvec)
            recoil_e_tpt    = math.sqrt(recoil_e_tpvec[0]*recoil_e_tpvec[0] 
                                      + recoil_e_tpvec[1]*recoil_e_tpvec[1])
        
            self.tree.recoil_e_truth_p.push_back(recoil_e_tp)
            self.tree.recoil_e_truth_pt.push_back(recoil_e_tpt)
            self.tree.recoil_e_truth_px.push_back(recoil_e_tpvec[0])
            self.tree.recoil_e_truth_py.push_back(recoil_e_tpvec[1])
            self.tree.recoil_e_truth_pz.push_back(recoil_e_tpvec[2])

            self.tree.recoil_e_vertex_x.push_back(recoil.getVertex()[0])
            self.tree.recoil_e_vertex_y.push_back(recoil.getVertex()[1])
            self.tree.recoil_e_vertex_z.push_back(recoil.getVertex()[2])


        # Get the A' mass 
        aprime = au.get_ap(particles)
        
        if not aprime:
            raise RuntimeError("A' was not found.")

        self.tree.ap_mass =  aprime.getMass()

        '''
        #
        # Ecal
        #
        ecal_hits = event.get_collection('EcalSimHits_sim')

        for ecal_hit in ecal_hits:

            is_recoil_e_hit = 0
            for icontrib in xrange(0, ecal_hit.getNumberOfContribs()):
                contrib = ecal_hit.getContrib(icontrib)
                if contrib.particle == recoil_e: 
                    is_recoil_e_hit = 1
                    break

            self.tree.ecal_is_recoil_e_hit.push_back(is_recoil_e_hit)
            self.tree.ecal_hit_energy.push_back(ecal_hit.getEdep())
            self.tree.ecal_hit_x.push_back(ecal_hit.getPosition()[0])
            self.tree.ecal_hit_y.push_back(ecal_hit.getPosition()[1])

        # Get the collection of digitized Ecal hits
        ecal_dhits = event.get_collection('ecalDigis_recon')
        
        ecal_dhit_count = 0
        layer_sum = 0
        layers = []
        ecal_hit_x = []
        ecal_hit_y = []
        denergy_sum = 0
        max_denergy_cell = 0
        max_layer_hit = 0
        for ecal_dhit in ecal_dhits: 

            self.tree.ecal_dhit_energy.push_back(denergy)
            denergy_sum += denergy

            layer = ecal_dhit.getLayer()
            layers.append(layer)
            self.tree.ecal_dhit_layer.push_back(layer)
            layer_sum += layer

            # Find the hex cell with the largest energy deposition
            max_denergy_cell = max(denergy, max_denergy_cell)
            
            max_layer_hit = max(max_layer_hit, layer)

        self.tree.ecal_dhit_count = ecal_dhit_count
        self.tree.average_ecal_layer_hit = -9999
        if ecal_dhit_count != 0: 
            self.tree.average_ecal_layer_hit = layer_sum/ecal_dhit_count
        self.tree.total_ecal_denergy = denergy_sum
        self.tree.ecal_max_denergy_cell = max_denergy_cell
        self.tree.ecal_max_layer_hit = max_layer_hit
        self.ecal_layer_std = np.std(np.array(layers)) 

        # Get the collection of Ecal veto results from the event
        ecal_veto_results = event.get_collection('EcalVeto_recon')
        
        self.tree.ecal_summed_tight_iso = ecal_veto_results[0].getSummedTightIso()
        self.tree.ecal_shower_rms = ecal_veto_results[0].getShowerRMS()
        self.tree.ecal_x_pos_std = ecal_veto_results[0].getXStd()
        self.tree.ecal_y_pos_std = ecal_veto_results[0].getYStd()

        # Get the BDT probability.
        self.tree.bdt_prob = ecal_veto_results[0].getDisc()

        target_sp_hits = event.get_collection('TargetScoringPlaneHits_sim')
        self.tree.single_trk_p = -9999
        self.tree.single_trk_pdg = -9999
        if len(findable_dic) == 1:
            #print '[ TargetPhotoNuclearAnalysis ]: Only have a single findable track.'
            findable_track = findable_dic.itervalues().next()
            sim_particle = findable_track.getSimParticle()
            self.tree.single_trk_pdg = findable_track.getSimParticle().getPdgID()
            p_find = la.norm(sim_particle.getMomentum())
            #print '[ TargetPhotoNuclearAnalysis ]: p of findable track: %s' % p_find
            #if is_recoil: print '[ TargetPhotoNuclearAnalysis ]: Findable track is recoil e.'
            #print 'Sim particle:'
            #sim_particle.Print()
            min = 10000
            for target_sp_hit in target_sp_hits: 
                if target_sp_hit.getSimParticle() == sim_particle: 
                    #print '[ TargetPhotoNuclearAnalysis ]: Iterating over target sp sim particle.'
                    #target_sp_hit.getSimParticle().Print()
                    pvec = target_sp_hit.getMomentum()
                    if pvec[2] < 0: 
                        #print '[ TargetPhotoNuclearAnalysis ]: particle is going backwards, pz = %s' % pvec[2] 
                        continue
                    p = la.norm(pvec)
                    #print '[ TargetPhotoNuclearAnalysis ]: Momentum at scoring plane: %s' % p
                    diff = p_find - p
                    
                    if diff < min: 
                        #print '[ TargetPhotoNuclearAnalysis ]: Found best match'
                        self.tree.single_trk_p = p
                        min = diff
        '''

        self.tree.fill(reset=True)

    def finalize(self): 
      
        self.tree.write()


