
import math

import ROOT as r

from numpy import linalg as la
from rootpy.tree import Tree

import AnalysisUtils as au
from EventModels import ElectroNuclearEvent

class ElectroNuclearAnalysis(object): 

    def __init__(self): 
        self.tree = None

    def initialize(self): 
        self.tree = Tree('en_ntuple', model=ElectroNuclearEvent)

        self.lfit = r.TF1('lfit', 'exp([1]-[0]*x)', 950, 1150)
        self.lfit.SetParameters(0.01093, 2.766)

        func = 'exp([1]-[0]*x)*([2]+[3]*x + [4]*pow(x,2)'
        func += ' + [5]*pow(x,3) + [6]*pow(x,4) + [7]*pow(x,5) + [8]*pow(x,6))'
        self.hfit = r.TF1('hfit', func, 1150, 3700)
        self.hfit.SetParameters(0.004008, 11.23, 1.242e-6, -1.964e-9, 
                                8.243e-13, 9.333e-17, -7.584e-20, -1.991e-23,
                                9.757e-27);

    def process(self, event): 

        # Get the collection of MC particles from the event.
        particles = event.get_collection('SimParticles_sim')

        # Search the list of sim particles for the recoil electron. If it isn't 
        # found, throw an exception
        recoil_e = au.get_recoil_electron(particles)

        # Persist the electron vertex
        self.tree.recoil_e_vx = recoil_e.getVertex()[0]
        self.tree.recoil_e_vy = recoil_e.getVertex()[1]
        self.tree.recoil_e_vz = recoil_e.getVertex()[2]
      
        # Get the momentum of the recoil electron.
        self.tree.recoil_e_tp = la.norm(recoil_e.getMomentum())

        # Calculate the q2 (beam energy - recoil energy)
        self.tree.q2 = (4000.0 - la.norm(recoil_e.getMomentum()))
        
        

        '''
        # Get all of the tracker hits in the event
        recoil_hits = event.get_collection('SiStripHits_recon')

        # Loop through of all the hits and find the truth momentum of the recoil 
        for recoil_hit in recoil_hits: 
            recoil_sim_hit = recoil_hit.getSimTrackerHits().At(0)
            if recoil_sim_hit.getLayerID() == 1:
                if recoil_sim_hit.getSimParticle() == recoil_e:
                    self.tree.recoil_e_tp = la.norm(recoil_sim_hit.getMomentum())
                    break

        # Use the recoil electron to retrieve the electronuclear daughters
        en_daughters = []
        en_weight = 1
        highest_w_nucleon_w = -9999
        highest_w_theta = -9999
        for idaughter in xrange(0, recoil_e.getDaughterCount()): 
            daughter = recoil_e.getDaughter(idaughter)

            if daughter.getProcessType() == 4:
                en_daughters.append(daughter)
       
                if (daughter.getPdgID() != 2212) and (daughter.getPdgID() != 2112): continue 
                ke = au.get_kinetic_energy(daughter)
                p_en = la.norm(daughter.getMomentum())
                theta = au.get_theta_z(daughter) 
                w = au.calculate_w(daughter, 0.5)  
                if w > highest_w_nucleon_w: 
                    highest_w_nucleon_w = w
                    highest_w_theta = theta

        enumbers = [33183, 39590]
        if event.get_event_number() in enumbers:
            print '###### Event: %s ########' % event.get_event_number()
            for daughter in en_daughters:
                print 'PDG ID: %s, ke: %s, theta: %s, endpoint: [%s, %s, %s], w: %s' % (
                        daughter.getPdgID(), 
                        au.get_kinetic_energy(daughter),
                        au.get_theta_z(daughter),
                        daughter.getEndPoint()[0],
                        daughter.getEndPoint()[1],
                        daughter.getEndPoint()[2], 
                        au.calculate_w(daughter, 0.5))  



        if (highest_w_nucleon_w > 1500) & (highest_w_theta > 100):
            en_weight = self.calculate_weight(highest_w_nucleon_w) 

        self.tree.en_particle_mult = len(en_daughters)
        self.tree.en_reaction_z = en_daughters[0].getVertex()[2]
        self.tree.en_weight = en_weight

        # Get the FindableTracks collection from the event
        findable_tracks = event.get_collection('FindableTracks_recon')
        findable_dic, stub_dic, axial_dic = au.get_findable_tracks_map(findable_tracks)
      
        if len(findable_dic) == 1: 
            
            findable_track = findable_dic.itervalues().next()
            sim_particle = findable_track.getSimParticle()
            self.tree.single_trk_pdg = findable_track.getSimParticle().getPdgID()
            
            if sim_particle == recoil_e: 
                for recoil_hit in recoil_hits: 
                    recoil_sim_hit = recoil_hit.getSimTrackerHits().At(0)
                    if recoil_sim_hit.getLayerID() == 1: 
                        if recoil_sim_hit.getSimParticle() == recoil_e:
                            p_find = la.norm(recoil_sim_hit.getMomentum())
                            break
            else: p_find = la.norm(sim_particle.getMomentum())
            
            self.tree.single_trk_p = p_find
        
        # Get the lead hadron in the event
        lead_hadron, lead_proton, lead_neutron, lead_pion = self.get_lead_hadrons(en_daughters)
        self.tree.lead_hadron_ke = au.get_kinetic_energy(lead_hadron)
        self.tree.lead_hadron_pdg_id = lead_hadron.getPdgID()
        self.tree.lead_hadron_theta_z = au.get_theta_z(lead_hadron) 

        # Determine if an event can be classified as a single neutron event   
        fneutron_count, fproton_count = self.hard_forward_hadron_count(en_daughters)
        if fneutron_count == 1: 
            self.tree.is_single_neutron = 1
        elif fneutron_count == 2: 
            self.tree.is_dineutron = 1
        elif fproton_count == 2: 
            self.tree.is_diproton = 1
        elif (fneutron_count == 1) & (fproton_count == 1): 
            self.tree.is_pn = 1

        # Check if the event is unphysical
        if self.is_ghost_event(lead_hadron): self.tree.is_ghost = 1
        '''

        self.tree.fill(reset=True)
    
    def finalize(self): 

        self.tree.write()

    def calculate_weight(self, w): 
        return self.lfit.Eval(w)/self.hfit.Eval(w)
   
    def get_lead_hadrons(self, en_daughters): 

        lead_ke         = -9999
        lead_proton_ke  = -9999
        lead_neutron_ke = -9999
        lead_pion_ke    = -9999

        lead_hadron  = None
        lead_proton  = None
        lead_neutron = None
        lead_pion    = None
        
        for en_daughter in en_daughters: 
            
            ke = au.get_kinetic_energy(en_daughter) 
            if lead_ke < ke: 
                lead_ke = ke
                lead_hadron = en_daughter

            if en_daughter.getPdgID() == 2112: 
                if lead_neutron_ke < ke: 
                    lead_neutron_ke = ke
                    lead_neutron = en_daughter
            
            if en_daughter.getPdgID() == 2212: 
                if lead_proton_ke < ke: 
                    lead_proton_ke = ke
                    lead_proton = en_daughter
            
            if (abs(en_daughter.getPdgID()) == 211) or (en_daughter.getPdgID() == 111):
                if lead_pion_ke < ke: 
                    lead_pion_ke = ke
                    lead_pion = en_daughter

        return lead_hadron, lead_proton, lead_neutron, lead_pion

    def is_ghost_event(self, lead_hadron): 
        
        ke = au.get_kinetic_energy(lead_hadron)
        theta_z = au.get_theta_z(lead_hadron)

        theta_z_theo = math.exp(-0.000575249*ke + 5.3643)
        
        #print 'KE: %s, Theta: %s, Theta (exp): %s' % (ke, theta_z, theta_z_theo)
        
        if theta_z > theta_z_theo: return 1
        else: return 0

    def hard_forward_hadron_count(self, daughters): 

        fneutron_count = 0
        fproton_count = 0

        for en_daughter in daughters: 

            ke = au.get_kinetic_energy(en_daughter)
            theta_z = au.get_theta_z(en_daughter)
            if ke > 100 and theta_z < 90:
                if en_daughter.getPdgID() == 2112: 
                    fneutron_count += 1
                elif en_daughter.getPdgID() == 2212: 
                    fproton_count += 1
        
        return fneutron_count, fproton_count
