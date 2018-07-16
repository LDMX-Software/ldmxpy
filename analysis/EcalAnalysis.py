

from numpy import linalg as la
from rootpy.tree import Tree

import AnalysisUtils as au
#import Line as l

from EventModels import EcalEvent


class EcalAnalysis(object): 

    def __init__(self): 
        self.tree = None
    
    def initialize(self):

        self.tree = Tree('ecal_ntuple', model=EcalEvent)

        #self.tfile = open('events.log', 'w')

        #self.event_count = 0
        #self.cfile_name = ''

    def process(self, event):
    
        #if self.cfile_name != event.get_file_name(): 
        #    self.event_count = 0
        #    self.cfile_name = event.get_file_name()
        
        #self.event_count += 1

        # Get the collection of MC particles from the event
        particles = event.get_collection('SimParticles_sim')

        # Search the list of sim particles for the recoil electron.  If it isn't
        # found, throw an exception.
        recoil_e = au.get_recoil_electron(particles)

        # Get the collection of EcalVeto objects from the event.  These objects
        # contain the results from the BDT.
        ecal_veto_results = event.get_collection('EcalVeto_recon')

        # Get the BDT probability.
        self.tree.bdt_prob = ecal_veto_results[0].getDisc()
        self.tree.passes_ecal_veto =  0
        if ecal_veto_results[0].passesVeto(): self.tree.passes_ecal_veto =  1

        # Get the scoring plane hits at the Ecal and search for any hits
        # associted with the recoil electron
        ecal_sp_hits = event.get_collection('EcalScoringPlaneHits_sim')
        self.tree.recoil_e_ecal_sp_x = -9999
        self.tree.recoil_e_ecal_sp_y = -9999
        self.tree.recoil_e_ecal_sp_z = -9999
        max_p = 0
        recoil_e_ecal_sp_hit = None
        for ecal_sp_hit in ecal_sp_hits: 
            if ecal_sp_hit.getSimParticle() == recoil_e:
                if ecal_sp_hit.getMomentum()[2] <= 0: continue
                #print 'Ecal scoring plane: x: %s y: %s z: %s ' % (ecal_sp_hit.getPosition()[0],
                #        ecal_sp_hit.getPosition()[1], ecal_sp_hit.getPosition()[2])
                #print 'Ecal scoring plane: px: %s py: %s pz: %s' % (ecal_sp_hit.getMomentum()[0], 
                #        ecal_sp_hit.getMomentum()[1], ecal_sp_hit.getMomentum()[2])
                pvec = ecal_sp_hit.getMomentum()
                p = la.norm(pvec)
                if max_p < p: 
                    #print 'Max p found'
                    max_p = p
                    recoil_e_ecal_sp_hit = ecal_sp_hit
       
        #line = l.Line(recoil_e_ecal_sp_hit.getPosition(), recoil_e_ecal_sp_hit.getMomentum())

        #if self.tree.recoil_e_ecap_sp_x != -9999: 
            

        #
        # Hit Level
        #
        
        # Get the collection of digitized Ecal hits
        ecal_dhits = event.get_collection('ecalDigis_recon')
      
        ecal_dhit_count = 0
        ecal_dhit_l1_count = 0
        ecal_dhit_l1_energy_sum = 0
        denergy_sum = 0

        for ecal_dhit in ecal_dhits: 

            # Get the digitized hit energy
            denergy = ecal_dhit.getEnergy()

            # Currently, if the energy of a digitized hit is below threshold,
            # its energy is set to zero.  These hits aren't considered.
            if denergy == 0: continue 
            ecal_dhit_count += 1
            
            if ecal_dhit.getLayer() == 1:
                if ecal_dhit.getAmplitude() > 5*self.get_noise_rms():
                    ecal_dhit_l1_count += 1
                    ecal_dhit_l1_energy_sum += denergy

            denergy_sum += denergy
       
            #print "Distance to: %s" % line.distance_to_point(ecal_dhit.getPosition())

        self.tree.total_ecal_denergy = denergy_sum
        self.tree.ecal_layer1_hit_count = ecal_dhit_l1_count
        self.tree.ecal_layer1_energy_sum = ecal_dhit_l1_energy_sum
        '''

        if recoil_e_ecal_sp_hit: 
            self.tree.recoil_e_ecal_sp_x = recoil_e_ecal_sp_hit.getPosition()[0]
            self.tree.recoil_e_ecal_sp_y = recoil_e_ecal_sp_hit.getPosition()[1]
            self.tree.recoil_e_ecal_sp_z = recoil_e_ecal_sp_hit.getPosition()[2]

            self.tree.recoil_e_ecal_sp_p  = max_p
            self.tree.recoil_e_ecal_sp_px = recoil_e_ecal_sp_hit.getMomentum()[0]
            self.tree.recoil_e_ecal_sp_py = recoil_e_ecal_sp_hit.getMomentum()[1]
            self.tree.recoil_e_ecal_sp_pz = recoil_e_ecal_sp_hit.getMomentum()[2]
        '''
        ''' 
        if self.tree.bdt_prob > 0: 
            self.tfile.write('%s : %s : %s : %s\n' % (event.get_file_name(),
                                                 event.get_event_number(),
                                                 self.event_count, 
                                                 pn_weight))
            
            for pn_daughter_count in xrange(0, pn_gamma.getDaughterCount()): 
                pn_daughter = pn_gamma.getDaughter(pn_daughter_count)
            
                ke = pn_daughter.getEnergy() - pn_daughter.getMass()
                pvec = pn_daughter.getMomentum()
                p = la.norm(pvec)
                theta = math.acos(pvec[2]/p)*180.0/3.14159
                ep = pn_daughter.getEndPoint()

                self.tfile.write('\tPDG ID: %s, KE: %s, theta: %s, p: %s, end point: [%s, %s, %s] \n' % (
                    pn_daughter.getPdgID(), ke, theta, p, ep[0], ep[1], ep[2]))

        '''
        
        self.tree.fill(reset=True)

    def finalize(self): 

        self.tree.write()
        #self.tfile.close()

    def get_noise_rms(self):
        return (900 + 22*27.56)*(0.130/33000)
