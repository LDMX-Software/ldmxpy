
from __future__ import division

import math
import numpy as np
import Plotter
import ROOT as r

from numpy import linalg as la

class SignalAnalysis:

    def __init__(self): 
        self.initialize()

    def initialize(self):
        
        self.events = []

        self.ap_mass = []
        
        self.recoil_e_truth_p  = []
        self.recoil_e_truth_pt = []
        self.recoil_e_truth_px = []
        self.recoil_e_truth_py = []
        self.recoil_e_truth_pz = []

        #self.recoil_e_theta = 

        self.recoil_is_findable = []
        self.recoil_is_3p1_findable = []
        self.recoil_is_2p2_findable = []
        
        self.missing_layer = []
        self.missing_layer_apmass = []

        self.triggered = []

    def process(self, event):
        
        #print 'Event: %s' % event.get_event_number()

        # Get the collection of MC particles from the event
        particles = event.get_collection('SimParticles')
        
        if particles.GetEntriesFast() == 0: return

        # Recoil electron
        recoil_e = None

        # A'
        aprime = None
       
        for particle in particles: 

            #print 'PDG ID: %s' % particle.getPdgID()
            #print 'Status: %s' % particle.getGenStatus()
            if (particle.getPdgID() == 11) & (particle.getGenStatus() == 1):
                recoil_e = particle

            if particle.getPdgID() == 622: 
                aprime = particle

        self.ap_mass.append(aprime.getMass())

        # Calculate the truth momentum of the recoil
        recoil_e_pvec = recoil_e.getMomentum()
        self.recoil_e_truth_p.append(la.norm(recoil_e_pvec))
        self.recoil_e_truth_pt.append(
                math.sqrt(recoil_e_pvec[0]*recoil_e_pvec[0] + recoil_e_pvec[1]*recoil_e_pvec[1]))
        self.recoil_e_truth_px.append(recoil_e_pvec[0]) 
        self.recoil_e_truth_py.append(recoil_e_pvec[1]) 
        self.recoil_e_truth_pz.append(recoil_e_pvec[2]) 

        recoil_sim_hits = event.get_collection('RecoilSimHits')
        recoil_e_hit_count = np.zeros(10)
        for recoil_sim_hit in recoil_sim_hits:
            if recoil_sim_hit.getSimParticle() == recoil_e:
                recoil_e_hit_count[recoil_sim_hit.getLayerID() - 1] += 1 
        
        n_stereo_hits = 0
        is_3p1_findable = False
        is_2p2_findable = False
        for layer_n in xrange(0, 8):
            if recoil_e_hit_count[layer_n]*recoil_e_hit_count[layer_n + 1] != 0: 
                n_stereo_hits += 1

        if n_stereo_hits > 3:
            is_3p1_findable = True
            is_2p2_findable = True
        if (n_stereo_hits == 3) & ((recoil_e_hit_count[8] > 0) or (recoil_e_hit_count[9] > 0)):
            is_3p1_findable = True
            is_2p2_findable = True
        if (n_stereo_hits == 2) & ((recoil_e_hit_count[8] > 0) or (recoil_e_hit_count[9] > 0)):
            is_2p2_findable = True

        if is_3p1_findable or is_2p2_findable: 
            self.recoil_is_findable.append(1)
        else: self.recoil_is_findable.append(0)

        if is_3p1_findable: self.recoil_is_3p1_findable.append(1)
        else: self.recoil_is_3p1_findable.append(0)
        
        if is_2p2_findable: self.recoil_is_2p2_findable.append(1)
        else: self.recoil_is_2p2_findable.append(0)

        if not is_3p1_findable and not is_2p2_findable: 
            for layer_n in xrange(0, 10):
                if recoil_e_hit_count[layer_n] == 0:
                    self.missing_layer.append(layer_n + 1)
                    self.missing_layer_apmass.append(aprime.getMass())


        trigger_results = event.get_collection('Trigger')
        for trigger_result in trigger_results: 
            if trigger_result.passed(): self.triggered.append(1)
            else: self.triggered.append(0)

    def finalize(self): 
       
        masses = np.unique(self.ap_mass)
        self.recoil_is_findable = np.array(self.recoil_is_findable)
        self.recoil_is_3p1_findable = np.array(self.recoil_is_3p1_findable)
        self.recoil_is_2p2_findable = np.array(self.recoil_is_2p2_findable)
        self.ap_mass = np.array(self.ap_mass)
        self.recoil_e_truth_p = np.array(self.recoil_e_truth_p)
        self.recoil_e_truth_pt = np.array(self.recoil_e_truth_pt)
        self.triggered = np.array(self.triggered)
        self.missing_layer = np.array(self.missing_layer)
        self.missing_layer_apmass = np.array(self.missing_layer_apmass)

        labels = []
        
        accep = []
        accep_pcut = []

        trk_accep = []
        trk_accep_3p1 = []
        trk_accep_2p2 = []
        trk_accep_pcut = []
        trk_accep_energy_cut_err = []
        trk_accep_err = np.zeros(len(masses))

        recoil_e_truth_p_mass  = []
        recoil_e_truth_pt_mass = []
        
        missing_layer_mass = []

        trigger_accep = []
        trigger_accep_pcut = []

        for mass in masses:
            recoil_is_findable_mass = self.recoil_is_findable[self.ap_mass == mass]
            recoil_is_findable_mass_pcut = self.recoil_is_findable[(self.ap_mass == mass) & (self.recoil_e_truth_p < 1200)]
            print len(recoil_is_findable_mass)
            print len(recoil_is_findable_mass_pcut)
            trigger_mass = self.triggered[self.ap_mass == mass]
            trigger_mass_pcut = self.triggered[(self.ap_mass == mass) & (self.recoil_e_truth_p < 1200)]
           
            accep.append(
                    (len(recoil_is_findable_mass[(recoil_is_findable_mass == 1) 
                        & (trigger_mass == 1)])/len(recoil_is_findable_mass))*100)

            accep_pcut.append(
                    (len(recoil_is_findable_mass_pcut[(recoil_is_findable_mass_pcut == 1) 
                        & (trigger_mass_pcut == 1)])/len(recoil_is_findable_mass_pcut))*100)
            trk_accep.append(
                    (len(recoil_is_findable_mass[recoil_is_findable_mass == 1])/len(recoil_is_findable_mass))*100)
        
            recoil_is_findable_mass = self.recoil_is_3p1_findable[self.ap_mass == mass]
            trk_accep_3p1.append(
                    (len(recoil_is_findable_mass[recoil_is_findable_mass == 1])/len(recoil_is_findable_mass))*100)
            
            recoil_is_findable_mass = self.recoil_is_2p2_findable[self.ap_mass == mass]
            trk_accep_2p2.append(
                    (len(recoil_is_findable_mass[recoil_is_findable_mass == 1])/len(recoil_is_findable_mass))*100)

            trk_accep_pcut.append(
                    (len(recoil_is_findable_mass_pcut[recoil_is_findable_mass_pcut == 1])/len(recoil_is_findable_mass_pcut))*100)
            
            recoil_e_truth_p_mass.append(self.recoil_e_truth_p[self.ap_mass == mass])
            recoil_e_truth_pt_mass.append(self.recoil_e_truth_pt[self.ap_mass == mass])
            
            missing_layer_mass.append(self.missing_layer[self.missing_layer_apmass == mass])

            trigger_accep.append((len(trigger_mass[trigger_mass == 1])/len(trigger_mass))*100)
            trigger_accep_pcut.append((len(trigger_mass_pcut[trigger_mass_pcut == 1])/len(trigger_mass_pcut))*100)
            
            labels.append("$A'$ mass: %s MeV" % mass) 
        
        trk_accep_arr = [trk_accep, trk_accep_3p1, trk_accep_2p2, trk_accep_pcut]
        mass_arr = [masses, masses, masses, masses]

        plt = Plotter.Plotter('signal_analysis.pdf')

        plt.plot_hists(recoil_e_truth_p_mass, 
                       np.linspace(0, 4000, 120), 
                       labels=labels,
                       norm=True,
                       x_label='$p(e^{-})$ (GeV)',
                       ylog=True)
        
        plt.plot_hists(recoil_e_truth_pt_mass, 
                       np.linspace(0, 500, 90), 
                       labels=labels,
                       norm=True,
                       ylog=True, 
                       x_label='$p_{t}(e^{-})$ (GeV)')
       
        plt.plot_graphs([masses, masses], [trigger_accep, trigger_accep_pcut], 
                       np.zeros(len(masses)), np.zeros(len(masses)),
                       labels=['All', 'p < 1.2 GeV'],
                       xlog = True, 
                       x_label="$A'$ Mass (MeV)",
                       y_label="Trigger acceptance (%)", 
                       ylim=[0, 100])

        plt.plot_graphs([masses, masses], [accep, accep_pcut],
                       np.zeros(len(masses)), np.zeros(len(masses)), 
                       labels=['All', 'p < 1.2 GeV'],
                       xlog=True,
                       x_label="$A'$ Mass (MeV)",
                       y_label="Total Acceptance (%)", 
                       ylim=[0, 100])
        
        plt.plot_graphs(mass_arr, trk_accep_arr, 
                        np.zeros(len(masses)), np.zeros(len(masses)),
                        labels=['3p1+2p2', '3p1', '2p2', '3p1+2p2, p < 1.2 GeV'],
                        xlog=True,
                        x_label="$A'$ Mass (MeV)",
                        y_label="Tracker Acceptance (%)",
                        ylim=[0, 100])

        plt.plot_hists(missing_layer_mass, 
                       np.linspace(0, 11, 12), 
                       labels=labels,
                       norm=True,
                       x_label='Missing Recoil $e^-$ Hit')

        plt.close()
