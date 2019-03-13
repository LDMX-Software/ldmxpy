
import math
import ROOT as r

from rootpy.tree import Tree

from numpy import linalg as la
import AnalysisUtils as au
from EventModels import PhotoNuclearEvent

class PhotoNuclearAnalysis(object): 

    def __init__(self): 
   
        self.tree = None

    def initialize(self, params): 
        self.tree = Tree('pn_ntuple', model=PhotoNuclearEvent)

    def process(self, event):

        #self.event_count += 1
        #print "################################################\n"

        #print "#                   Event %s                   #\n" % self.event_count
        #print "################################################\n"

        # Get the collection of MC particles from the event.
        particles = event.get_collection('SimParticles_sim')

        # Search the list of sim particles for the recoil electron. If it isn't 
        # found, throw an exception
        recoil_e = au.get_recoil_electrons(particles)[0]

        # Persist the electron vertex
        self.tree.recoil_e_vx = recoil_e.getVertex()[0]
        self.tree.recoil_e_vy = recoil_e.getVertex()[1]
        self.tree.recoil_e_vz = recoil_e.getVertex()[2]

        # Retrieve the collection ot target scoring plane hits.  
        sp_hits = event.get_collection('TargetScoringPlaneHits_sim')

        # Create a list of all scoring plane hits associated with the recoil
        # electron.  Only hits created in the downstream scoring plane and 
        # with a position pz momentum are considered.
        recoil_sp_hits = [sp_hit for sp_hit in sp_hits if (
                            (sp_hit.getSimParticle() == recoil_e) & 
                            (sp_hit.getLayerID() == 2) & 
                            (sp_hit.getMomentum()[2] > 0)
                          )]
       
        if len(recoil_sp_hits) == 0: 
            print 'Missing recoil electron. Skipping event.'
            return

        pvec = recoil_sp_hits[0].getMomentum()
        pvec = r.TVector3(pvec[0], pvec[1], pvec[2])
        recoil_e_p = pvec.Mag()
        
        self.tree.recoil_e_p = recoil_e_p
        self.tree.recoil_e_px = pvec.Px()
        self.tree.recoil_e_py = pvec.Py()
        self.tree.recoil_e_pz = pvec.Pz()
        self.tree.recoil_e_pt = pvec.Pt()


        # Calculate the energy of the recoil electron
        m_e = 0.5109989461
        recoil_e_energy = math.sqrt(pvec.Mag()*pvec.Mag() + m_e*m_e)

        self.tree.recoil_e_energy = recoil_e_energy

        # Recoil electron polar angle 
        self.tree.recoil_e_theta = pvec.Theta()*(180/math.pi)

        # Create a vector denoting the incident electron momentum
        in_e_pvec = r.TVector3(0, 0, math.sqrt(4000*4000 - m_e*m_e))

        # Calculate the momentum transfer
        self.tree.q = (in_e_pvec - pvec).Mag()

        # Calculate the energy transfer
        self.tree.omega = 4000 - recoil_e_energy

        # Use the recoil electron to retrieve the gamma that underwent a 
        # photonuclear reaction.
        pn_gamma = au.get_pn_gamma([recoil_e])

        pn_gamma_pvec = pn_gamma.getMomentum()

        self.tree.pn_gamma_energy     = pn_gamma.getEnergy()
        self.tree.pn_gamma_int_z      = pn_gamma.getEndPoint()[2]
        self.tree.pn_gamma_vertex_z   = pn_gamma.getVertex()[2]

        pn_particles = [particle for particle in particles if (
            (particle.getPdgID() != 22) &
            (particle.getPdgID() < 10000) &
            (particle.getProcessType() == 9))]

        self.tree.pn_particle_mult = len(pn_particles)

        var = { 'h_nucleon_ke':-9999, 'h_nucleon_theta':-9999, 
                'h_nucleon_pdg':-9999,'h_p_ke':-9999, 'h_p_theta':-9999, 
                'h_n_ke':-9999, 'h_n_theta':-9999, 'h_pi_ke':-9999, 
                'h_pi_theta':-9999, 'h_pi0_ke':-9999, 'h_pi0_theta':-9999
        }
        
        for particle in pn_particles:
           
            #particle.Print()

            pvec = particle.getMomentum()
            vec = r.TVector3(pvec[0], pvec[1], pvec[2])

            # Calculate the polar angle
            theta = vec.Theta()*180/math.pi
            self.tree.hadron_theta.push_back(theta)
            
            self.tree.hadron_omega.push_back(self.tree.omega)
            self.tree.hadron_recoil_pt.push_back(self.tree.recoil_e_pt)
            self.tree.hadron_q.push_back(self.tree.q)

            # Calculate the kinetic energy
            ke = au.get_kinetic_energy(particle)
            self.tree.hadron_ke.push_back(ke)

            pdg = abs(particle.getPdgID())
            self.tree.hadron_pdgid.push_back(particle.getPdgID())

            if var['h_nucleon_ke'] < ke: 
                var['h_nucleon_ke'] = ke
                var['h_nucleon_theta'] = theta
                var['h_nucleon_pdg'] = particle.getPdgID()

            # Calculate the weight
            weight = particle.getEnergy()

            if pdg == 2112: 
                weight = ke
                if var['h_n_ke'] < ke: 
                    var['h_n_ke'] = ke
                    var['h_n_theta'] = theta
            elif pdg == 2212: 
                weight = ke
                if var['h_p_ke'] < ke: 
                    var['h_p_ke'] = ke
                    var['h_p_theta'] = theta
            elif pdg == 211: 
                if var['h_pi_ke'] < ke: 
                    var['h_pi_ke'] = ke
                    var['h_pi_theta'] = theta
            elif pdg == 111: 
                if var['h_pi0_ke'] < ke: 
                    var['h_pi0_ke'] = ke
                    var['h_pi0_theta'] = theta

            
        self.tree.lead_hadron_ke     = var['h_nucleon_ke']
        self.tree.lead_hadron_theta  = var['h_nucleon_theta']
        self.tree.lead_hadron_pdg_id = var['h_nucleon_pdg']
        self.tree.lead_p_ke          = var['h_p_ke']
        self.tree.lead_p_theta       = var['h_p_theta']
        self.tree.lead_n_ke          = var['h_n_ke']
        self.tree.lead_n_theta       = var['h_n_theta']
        self.tree.lead_pi_ke         = var['h_pi_ke']
        self.tree.lead_pi_theta      = var['h_pi_theta']
        self.tree.lead_pi0_ke        = var['h_pi0_ke']
        self.tree.lead_pi0_theta     = var['h_pi0_theta']

        event_type = au.classify_event(pn_particles, 200)
        self.tree.event_type = event_type

        self.tree.fill(reset=True)

    def finalize(self): 

        self.tree.write()

