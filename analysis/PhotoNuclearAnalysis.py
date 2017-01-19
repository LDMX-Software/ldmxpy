
from __future__ import division

import copy
import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import Plotter
import root_numpy as rnp
import ROOT as r
import rootpy.plotting.root2matplotlib as rplt


from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LogNorm
from scipy.stats import norm

class PhotoNuclearAnalysis(object) : 

    def __init__(self) : 

        self.initialize()

    def is_recoil(self, particle) : 
        return (particle.getPdgID() == 11) & (particle.getParentCount() == 0)

    def created_within_target(self, particle) : 
        if abs(particle.getVertex()[2]) < 0.550 : return True 
        return False
    
    def track_is_findable(self, sim_particle, sim_hits):
        hits = [0,0,0,0,0,0,0,0,0,0]
        print sim_hits.GetEntriesFast()
        for sim_hit_n in xrange(0,sim_hits.GetEntriesFast()):
            sim_hit = sim_hits.At(sim_hit_n)

            if sim_hit.getSimParticle().GetObject() == sim_particle:
                hits[sim_hit.getLayerID()] += 1

        print hitanalysis/s
        n_3d_hits = 0
        for layer_n in xrange(0,8,2):
            if hits[layer_n]*hits[layer_n+1] != 0: n_3d_hits += 1

        if n_3d_hits > 3: return True
        elif (n_3d_hits == 3) & (hits[8] > 0 or hits[9] > 0): return True

        return False

    def initialize(self) :

        self.events              = []
        self.pn_gamma_energy     = []
        self.pn_particle_mult    = []
        self.lead_neutron_ke     = []
        self.lead_neutron_pdgid  = []
        self.lead_neutron_theta  = []
        self.lead_neutron_phi    = []
        self.lead_em_ke     = []
        self.lead_em_pdgid  = []
        self.lead_em_theta  = []
        self.lead_em_phi    = []
        
        self.n_tracks = []
        self.n_charged = []

        self.total_energy_ecal = []
        self.total_energy_hcal = []
        self.energy_hcal = []
        self.hcal_z_pos = []
       
        self.target_total_energy = []
        
        self.event_count = -1

    def process(self, event) :
        
        print 'Event: %s' % event.get_event_number()
        
        # Get the collection of MC particles from the event
        particles = event.get_collection('SimParticles')
        
        recoil = None
        pn_gamma = None

        for particle in particles:

            #print 'PDG ID: %s' % particle.getPdgID()
            #print 'Parent count: %s' % particle.getParentCount()
            #print 'Daughtercount: %s' % particle.getDaughterCount()
            
            if self.is_recoil(particle):
                recoil = particle
                break

        # Use the recoil electron to find the PN gamma
        for daughter_count in xrange(0, recoil.getDaughterCount()):
            daughter = recoil.getDaughter(daughter_count)
            #print 'Daughter PDG ID: %s' % daughter.getPdgID()
            #print 'Daughter Parent count: %s' % daughter.getParentCount()
            #print 'Daughter Daughtercount: %s' % daughter.getDaughterCount()
            if (daughter.getPdgID() == 22) & self.created_within_target(daughter):
                pn_gamma = daughter
         
        self.pn_gamma_energy.append(pn_gamma.getEnergy())
        self.pn_particle_mult.append(pn_gamma.getDaughterCount())

    def finalize(self) : 

        plt = Plotter.Plotter("photo_nuclear_analysis.pdf")

        plt.plot_hist(self.pn_gamma_energy, 
                      np.linspace(2500, 4000, 41),
                      ylog=True,
                      x_label='$E(\gamma)$ (MeV)')
                    
        plt.plot_hist(self.pn_particle_mult,
                      np.linspace(0, 100, 101),
                      ylog=True,
                      x_label='PN Multiplicity')

        plt.close()
