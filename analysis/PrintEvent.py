from __future__ import division

import math
import ROOT as r
import numpy as np
import Plotter

from numpy import linalg as la
from scipy.stats import norm

class PrintEvent(object) :

    def calculate_w(self, particle): 
        pvec = particle.getMomentum()
        p = la.norm(pvec)
        ke = particle.getEnergy() - particle.getMass()
        theta = math.acos(pvec[2]/p)*180.0/3.14159 
        return 0.5*(p + ke)*(1.12 - 0.5*(pvec[2]/p)), theta, p

    def is_recoil(self, particle) : 
        return (particle.getPdgID() == 11) & (particle.getParentCount() == 0)
   
    def created_within_target(self, particle) :
        if abs(particle.getVertex()[2]) <= 0.550 : return True 
        return False
    
    def process(self, event):

        if event.get_event_number() != 55918: 
            return

        # Get the collection of MC particles from the even
        particles = event.get_collection('SimParticles_sim')
        
        # Loop through all of the sim particles in the event and find the recoil
        # electron.  The recoil electron can then be used to obtain associated
        # brem gamma involved in a PN reaction.
        recoil_e = None
        for particle in particles:
            
            # We only care to find the recoil for now
            if particle.getPdgID() != 11: continue

            # If the particle is found to be consistent with the recoil 
            # electron, then tag it and skip searching
            if self.is_recoil(particle):
                recoil_e = particle
                break

        # Search for the PN gamma and use it to get the PN daughters
        pn_gamma = None
        for daughter_count in xrange(0, recoil_e.getDaughterCount()):
            daughter = recoil_e.getDaughter(daughter_count)
            if daughter.getDaughterCount() == 0: continue
            
            if ((daughter.getPdgID() == 22) 
                    & self.created_within_target(daughter)
                    & self.created_within_target(daughter.getDaughter(0))):
                pn_gamma = daughter
                break
        
        for pn_daughter_count in xrange(0, pn_gamma.getDaughterCount()): 
            pn_daughter = pn_gamma.getDaughter(pn_daughter_count)
            
            ke = pn_daughter.getEnergy() - pn_daughter.getMass()

            wp, theta, p = self.calculate_w(pn_daughter)

            print 'PDG ID: %s, KE: %s, theta: %s, p: %s w: %s' % (pn_daughter.getPdgID(), ke, theta, p, wp)
    
    def finalize(self) :
        print "Done:"
