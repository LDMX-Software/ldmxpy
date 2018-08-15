
import math

from rootpy.tree import Tree

from numpy import linalg as la
import AnalysisUtils as au
from EventModels import PhotoNuclearEvent

class PhotoNuclearAnalysis(object): 

    def __init__(self): 
   
        self.tree = None

    def initialize(self): 
        self.tree = Tree('pn_ntuple', model=PhotoNuclearEvent)

        self.event_count = 0
        self.tfile = open('events.log', 'w')

    def process(self, event):

        self.event_count += 1
        #print "################################################\n"
        #print "#                   Event %s                   #\n" % self.event_count
        #print "################################################\n"

        # Get the collection of MC particles from the event.
        particles = event.get_collection('SimParticles_sim')

        # Search the list of sim particles for the recoil electron.  If it isn't
        # found, throw an exception.
        #recoil_e = au.get_recoil_electron(particles)

        # Search the list of sim particles for the recoil electron. If it isn't 
        # found, throw an exception
        recoils = au.get_recoil_electrons(particles)
        self.tree.n_electrons = len(recoils)

        
        # Get the e- recoil truth information
        for i, recoil in enumerate(recoils): 
            self.tree.recoil_e_vertex_x.push_back(recoil.getVertex()[0])
            self.tree.recoil_e_vertex_y.push_back(recoil.getVertex()[1])
            self.tree.recoil_e_vertex_z.push_back(recoil.getVertex()[2])

        # Use the recoil electron to retrieve the gamma that underwent a 
        # photonuclear reaction.
        pn_gamma = au.get_pn_gamma(recoils)

        self.tree.pn_particle_mult = pn_gamma.getDaughterCount()
        self.tree.pn_gamma_energy  = pn_gamma.getEnergy()
        self.tree.pn_gamma_int_z   = pn_gamma.getEndPoint()[2]
        self.tree.pn_gamma_vertex_z   = pn_gamma.getVertex()[2]

        # Get the lead hadron in the event
        lead_hadron, lead_proton, lead_neutron, lead_pion = self.get_lead_hadrons(pn_gamma)
        self.tree.lead_hadron_ke = au.get_kinetic_energy(lead_hadron)
        self.tree.lead_hadron_pdg_id = lead_hadron.getPdgID()
        self.tree.lead_hadron_theta_z = au.get_theta_z(lead_hadron) 

        if lead_neutron:
            self.tree.lead_neutron_ke = au.get_kinetic_energy(lead_neutron)
            self.tree.lead_neutron_theta_z = au.get_theta_z(lead_neutron)

        if lead_proton: 
            self.tree.lead_proton_ke = au.get_kinetic_energy(lead_proton)
            self.tree.lead_proton_theta_z = au.get_theta_z(lead_proton)

        if lead_pion: 
            self.tree.lead_pion_ke = au.get_kinetic_energy(lead_pion)
            self.tree.lead_pion_theta_z = au.get_theta_z(lead_pion)

        # If the collection of PN weights exist, retrieve it from the event. 
        if event.collection_exist('PNweight_recon'): 
            weights = event.get_collection('PNweight_recon')
            
            # Retrieve the PN weight from the event
            pn_weight = float(weights[0].getWeight())
            self.tree.pn_weight = pn_weight
    
        # Determine if an event can be classified as a single neutron event   
        fneutron_count, fproton_count = self.hard_forward_hadron_count(pn_gamma)
        if fneutron_count == 1: 
            self.tree.is_single_neutron = 1
        elif fneutron_count == 2: 
            self.tree.is_dineutron = 1
        elif fproton_count == 2: 
            self.tree.is_diproton = 1
        elif (fneutron_count == 1) & (fproton_count == 1): 
            self.tree.is_pn = 1

        if self.is_ghost_event(lead_hadron):
            #print 'Is ghost'
            self.tree.is_ghost = 1

        # If the collection of findable tracks exist, retrieve it from the 
        # event. 
        #if event.collection_exist('FindableTracks_recon'): 

        #    findable_tracks = event.get_collection('FindableTracks_recon')
        #    findable_dic, stub_dic, axial_dic = au.get_findable_tracks_map(findable_tracks)    

        self.tree.fill(reset=True)

    def finalize(self): 

        self.tree.write()
        self.tfile.close()

    def get_lead_hadrons(self, pn_gamma): 

        lead_ke         = -9999
        lead_proton_ke  = -9999
        lead_neutron_ke = -9999
        lead_pion_ke    = -9999

        lead_hadron  = None
        lead_proton  = None
        lead_neutron = None
        lead_pion    = None
        
        for pn_daughter_count in xrange(0, pn_gamma.getDaughterCount()): 
            pn_daughter = pn_gamma.getDaughter(pn_daughter_count)
            
            ke = au.get_kinetic_energy(pn_daughter) 
            if lead_ke < ke: 
                lead_ke = ke
                lead_hadron = pn_daughter

            if pn_daughter.getPdgID() == 2112: 
                if lead_neutron_ke < ke: 
                    lead_neutron_ke = ke
                    lead_neutron = pn_daughter
            
            if pn_daughter.getPdgID() == 2212: 
                if lead_proton_ke < ke: 
                    lead_proton_ke = ke
                    lead_proton = pn_daughter
            
            if (abs(pn_daughter.getPdgID()) == 211) or (pn_daughter.getPdgID() == 111):
                if lead_pion_ke < ke: 
                    lead_pion_ke = ke
                    lead_pion = pn_daughter

        return lead_hadron, lead_proton, lead_neutron, lead_pion

    def hard_forward_hadron_count(self, pn_gamma): 

        #lh, lp, ln, lp = self.get_lead_hadrons(pn_gamma)

        fneutron_count = 0
        fproton_count = 0

        for pn_daughter_count in xrange(0, pn_gamma.getDaughterCount()): 
            pn_daughter = pn_gamma.getDaughter(pn_daughter_count)
            #if pn_daughter == ln: continue

            ke = au.get_kinetic_energy(pn_daughter)
            theta_z = au.get_theta_z(pn_daughter)
            if ke > 100 and theta_z < 90:
                if pn_daughter.getPdgID() == 2112: 
                    fneutron_count += 1
                elif pn_daughter.getPdgID() == 2212: 
                    fproton_count += 1
        
        return fneutron_count, fproton_count

    def is_ghost_event(self, lead_hadron): 
        
        ke = au.get_kinetic_energy(lead_hadron)
        theta_z = au.get_theta_z(lead_hadron)

        #if theta_z > 80: return 1
        #return 0
        theta_z_theo = math.exp(-0.000575249*ke + 5.3643)
        
        #print 'KE: %s, Theta: %s, Theta (exp): %s' % (ke, theta_z, theta_z_theo)
        
        if theta_z > theta_z_theo: 
            #print 'Ghost event found!'
            return 1
        else: return 0
       
        
