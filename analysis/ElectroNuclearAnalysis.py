
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

        self.tree.recoil_e_vertex_x = recoil_e.getVertex()[0]
        self.tree.recoil_e_vertex_y = recoil_e.getVertex()[1]
        self.tree.recoil_e_vertex_z = recoil_e.getVertex()[2]

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
      
        recoil_hits = event.get_collection('SiStripHits_recon')
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
        hadron = self.get_lead_hadron(en_daughters)
        self.tree.lead_hadron_ke = au.get_kinetic_energy(hadron)
        self.tree.lead_hadron_pdg_id = hadron.getPdgID()
        self.tree.lead_hadron_theta_z = au.get_theta_z(hadron) 

        self.tree.fill(reset=True)
    
    def finalize(self): 

        self.tree.write()

    def calculate_weight(self, w): 
        return self.lfit.Eval(w)/self.hfit.Eval(w)
    
    def get_lead_hadron(self, daughters): 

        lead_ke = -9999
        lead_proton = -9999
        lead_neutron = -9999
        lead_pion = -9999
        lead_hadron = None
        for pn_daughter in daughters: 
            
            ke = au.get_kinetic_energy(pn_daughter) 
            if lead_ke < ke: 
                lead_ke = ke
                lead_hadron = pn_daughter

        return lead_hadron
