
import math

from rootpy.tree import Tree

from numpy import linalg as la
import AnalysisUtils as au
from EventModels import PhotoNuclearEvent

class PhotoNuclearAnalysis(object): 

    def __init__(self): 
   
        self.tree = None
        # Open the ROOT file that will be used to write the ntuple to.
        #self.file = root_open('pn_ntuple.root', 'recreate')

    def initialize(self): 
        self.tree = Tree('pn_ntuple', model=PhotoNuclearEvent)

        self.event_count = 0
        self.tfile = open('events.log', 'w')

    def process(self, event):

        self.event_count += 1

        # Get the collection of MC particles from the event.
        particles = event.get_collection('SimParticles_sim')

        # Search the list of sim particles for the recoil electron.  If it isn't
        # found, throw an exception.
        recoil_e = au.get_recoil_electron(particles)

        self.tree.recoil_e_vertex_x = recoil_e.getVertex()[0]
        self.tree.recoil_e_vertex_y = recoil_e.getVertex()[1]
        self.tree.recoil_e_vertex_z = recoil_e.getVertex()[2]

        # Use the recoil electron to retrieve the gamma that underwent a 
        # photonuclear reaction.
        pn_gamma = au.get_pn_gamma(recoil_e)

        self.tree.pn_particle_mult = pn_gamma.getDaughterCount()
        self.tree.pn_gamma_energy  = pn_gamma.getEnergy()
        self.tree.pn_gamma_int_z   = pn_gamma.getEndPoint()[2]
        self.tree.pn_gamma_vertex_z   = pn_gamma.getVertex()[2]

        # Get the lead hadron in the event
        hadron = self.get_lead_hadron(pn_gamma)
        self.tree.lead_hadron_ke = au.get_kinetic_energy(hadron)
        self.tree.lead_hadron_pdg_id = hadron.getPdgID()
        self.tree.lead_hadron_theta_z = au.get_theta_z(hadron) 

        # Get the PN weight from the event
        weights = event.get_collection('PNweight_recon') 
        pn_weight = float(weights[0].getWeight())
        self.tree.pn_weight = pn_weight

        # Get the FindableTracks collection from the event
        findable_tracks = event.get_collection('FindableTracks_recon')
        findable_dic, stub_dic, axial_dic = au.get_findable_tracks_map(findable_tracks)

        recoil_hits = event.get_collection('SiStripHits_recon')
        if len(findable_dic) == 1: 
           
            #print 'Event has single track'

            findable_track = findable_dic.itervalues().next()
            sim_particle = findable_track.getSimParticle()
            self.tree.single_trk_pdg = findable_track.getSimParticle().getPdgID()
            
            if sim_particle == recoil_e:
                p_max = 0
                pt_max = 0
                for recoil_hit in recoil_hits: 
                    recoil_sim_hit = recoil_hit.getSimTrackerHits().At(0)
                    if recoil_sim_hit.getLayerID() == 1: 
                        if recoil_sim_hit.getSimParticle() == recoil_e:
                            p = la.norm(recoil_sim_hit.getMomentum())
                            if p_max < p: 
                                p_max = p
                                pt_max = au.get_pt(recoil_sim_hit.getMomentum())
                    
                p_find = p_max
                pt_find = pt_max
                print 'Single track momentum: %s' % p_find
                print 'Single track momentum: %s' % pt_find
            else: 
                p_find = la.norm(sim_particle.getMomentum())
                pt_find = au.get_pt(sim_particle.getMomentum())
            
            self.tree.single_trk_p = p_find
            self.tree.single_trk_pt = pt_find
        
        #events = [3, 6, 73, 92, 109, 113, 121, 124, 129, 134, 145, 146, 205,
        #          208, 215, 238, 253, 290, 307, 335, 344, 355, 359, 366, 396,
        #          397, 419, 428, 442, 463, 505, 534, 560, 567, 587, 591, 625, 637]

        #events = [36, 57, 83]
        if self.event_count in events: 
            self.tfile.write("################################################\n")
            self.tfile.write("#                   Event %s                   #\n" % self.event_count)
            self.tfile.write("################################################\n")
            for pn_daughter_count in xrange(0, pn_gamma.getDaughterCount()): 
                pn_daughter = pn_gamma.getDaughter(pn_daughter_count)
            
                ke = pn_daughter.getEnergy() - pn_daughter.getMass()
                pvec = pn_daughter.getMomentum()
                p = la.norm(pvec)
                theta = math.acos(pvec[2]/p)*180.0/3.14159
                ep = pn_daughter.getEndPoint()

                self.tfile.write('\tPDG ID: %s, KE: %s, theta: %s, p: %s, end point: [%s, %s, %s] \n' % (
                    pn_daughter.getPdgID(), ke, theta, p, ep[0], ep[1], ep[2]))
            

        self.tree.fill(reset=True)

    def finalize(self): 

        self.tree.write()
        self.tfile.close()

    def get_lead_hadron(self, pn_gamma): 

        lead_ke = -9999
        lead_proton = -9999
        lead_neutron = -9999
        lead_pion = -9999
        lead_hadron = None
        for pn_daughter_count in xrange(0, pn_gamma.getDaughterCount()): 
            pn_daughter = pn_gamma.getDaughter(pn_daughter_count)
            
            ke = au.get_kinetic_energy(pn_daughter) 
            if lead_ke < ke: 
                lead_ke = ke
                lead_hadron = pn_daughter

        return lead_hadron


