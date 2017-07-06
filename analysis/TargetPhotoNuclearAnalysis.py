from __future__ import division

import math
import ROOT as r
import numpy as np
import Plotter

from numpy import linalg as la
from scipy.stats import norm

class TargetPhotoNuclearAnalysis(object) : 

    def __init__(self) : 
        self.initialize()
        self.file_name = ''

    def calc_sig_eff(self, sig, min, max, step, invert) :
        
        cut = min
        sig_eff = []
        cuts    = []

        while cut <= max : 
            if invert : 
                sig_integrated = len(sig[sig > cut])
            else : 
                sig_integrated = len(sig[sig < cut])
            sig_eff.append((sig_integrated/len(sig))*100)
            cuts.append(cut)
            cut += step
        
        return cuts, sig_eff

    def is_recoil(self, particle) : 
        return (particle.getPdgID() == 11) & (particle.getParentCount() == 0)

    def created_within_target(self, particle) :
        if abs(particle.getVertex()[2]) <= 0.550 : return True 
        return False
  
    def calculate_w(self, particle): 
        pvec = particle.getMomentum()
        p = la.norm(pvec)
        ke = particle.getEnergy() - particle.getMass()
        theta = math.acos(pvec[2]/p)*180.0/3.14159 
        return 0.5*(p + ke)*(1.12 - 0.5*(pvec[2]/p)), theta

    def initialize(self) :

        self.ntuple = {}
        self.variables = [
                'events',
                'pn_gamma_energy', 'pn_particle_mult', 'pn_interaction_z', 
                'track_count', 'stub_count', 'axial_count', 
                'ecal_hit_energy', 'total_ecal_energy', 'passes_ecal_veto',
                'hcal_hit_energy', 'total_hcal_energy', 'total_hcal_pe',
                'passes_hcal_veto', 'up_tp_energy', 'down_tp_energy',
                'lead_hadron_ke', 'lead_hadron_theta', 
                'lead_hadron_p', 'lead_hadron_pdgid', 'max_w', 'max_w_theta', 
                'lead_pion_ke', 'lead_pion_theta', 'lead_pion_p', 
                'lead_neutron_ke', 'lead_neutron_theta', 'lead_neutron_p', 
                'lead_proton_ke', 'lead_proton_theta', 'lead_proton_p',
                'total_recoil_hits'
        ]
        for layer_n in xrange(0, 10):
            self.variables.append('total_recoil_hits_l%s' % (layer_n + 1))
            self.variables.append('total_charge_l%s' % (layer_n + 1))

        for variable in self.variables: 
            self.ntuple[variable] = []
    
        self.colors = [r.kAzure + 2, r.kRed + 2, r.kGreen + 2, r.kViolet + 6, r.kOrange + 7]

        self.target_total_energy = []
        
        self.event_count = -1
        
        self.file_prefix = None

    def process(self, event) :
        
        self.event_count += 1
        self.ntuple['events'].append(event.get_event_number())
        #print "Event: %s" % self.event_count

        if self.event_count == 1: 
            self.file_prefix = event.get_file_name()[
                    event.get_file_name().rfind('/') + 1:-5]
            print self.file_prefix

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
        
        self.ntuple['pn_gamma_energy'].append(pn_gamma.getEnergy())
        self.ntuple['pn_particle_mult'].append(pn_gamma.getDaughterCount())
        self.ntuple['pn_interaction_z'].append(pn_gamma.getEndPoint()[2])

        lead_ke = -9999
        lead_hadron = None
        lead_proton_ke = -9999
        lead_proton = None
        lead_neutron_ke = -9999
        lead_neutron = None
        lead_pion_ke = -9999
        lead_pion = None
        max_w = -9999
        max_w_theta = -9999
        for pn_daughter_count in xrange(0, pn_gamma.getDaughterCount()): 
            pn_daughter = pn_gamma.getDaughter(pn_daughter_count)
            
            ke = pn_daughter.getEnergy() - pn_daughter.getMass()
            if lead_ke < ke: 
                lead_ke = ke
                lead_hadron = pn_daughter
            
            if ((abs(pn_daughter.getPdgID()) == 2212) & (lead_proton_ke < ke)):
                lead_proton_ke = ke
                lead_proton = pn_daughter
            elif ((abs(pn_daughter.getPdgID()) == 2112) & (lead_neutron_ke < ke)):
                lead_neutron_ke = ke
                lead_neutron = pn_daughter
            elif ((abs(pn_daughter.getPdgID()) == 211) & (lead_neutron_ke < ke)):
                lead_pion_ke = ke
                lead_pion = pn_daughter

            wp, theta = self.calculate_w(pn_daughter)
            if wp > max_w: 
                max_w = wp
                max_w_theta = theta
        
        self.ntuple['lead_hadron_ke'].append(lead_ke)
        self.ntuple['lead_proton_ke'].append(lead_proton_ke)
        self.ntuple['lead_neutron_ke'].append(lead_neutron_ke)
        self.ntuple['lead_pion_ke'].append(lead_pion_ke)
        
        self.ntuple['max_w'].append(max_w)
        self.ntuple['max_w_theta'].append(max_w_theta)

        pdgid = -9999
        if lead_hadron == None: 
            theta = -9999
            p = -9999
        else:
            pvec = lead_hadron.getMomentum()
            p = la.norm(pvec)
            theta = math.acos(pvec[2]/p)*180.0/3.14159 
            pdgid = lead_hadron.getPdgID()

        self.ntuple['lead_hadron_theta'].append(theta)
        self.ntuple['lead_hadron_p'].append(p)
        self.ntuple['lead_hadron_pdgid'].append(pdgid)
        
        if lead_proton == None: 
            theta = -9999
            p = -9999
        else:
            pvec = lead_proton.getMomentum()
            p = la.norm(pvec)
            theta = math.acos(pvec[2]/p)*180.0/3.14159 

        self.ntuple['lead_proton_theta'].append(theta)
        self.ntuple['lead_proton_p'].append(p)

        if lead_neutron == None: 
            theta = -9999
            p = -9999
        else:
            pvec = lead_neutron.getMomentum()
            p = la.norm(pvec)
            theta = math.acos(pvec[2]/p)*180.0/3.14159 

        self.ntuple['lead_neutron_theta'].append(theta)
        self.ntuple['lead_neutron_p'].append(p)


        if lead_pion == None: 
            theta = -9999
            p = -9999
        else:
            pvec = lead_pion.getMomentum()
            p = la.norm(pvec)
            theta = math.acos(pvec[2]/p)*180.0/3.14159 

        self.ntuple['lead_pion_theta'].append(theta)
        self.ntuple['lead_pion_p'].append(p)
        
        #
        # Trigger Pads
        #
        trigger_pad_hits = event.get_collection('TriggerPadSimHits_sim')
        total_up_tp_energy = 0
        total_down_tp_energy = 0
        for trigger_pad_hit in trigger_pad_hits:
            if (trigger_pad_hit.getPosition()[2] > 0): 
                total_down_tp_energy += trigger_pad_hit.getEdep()
            if (trigger_pad_hit.getPosition()[2] < 0): 
                total_up_tp_energy += trigger_pad_hit.getEdep()
    
        self.ntuple['down_tp_energy'].append(total_down_tp_energy)
        self.ntuple['up_tp_energy'].append(total_up_tp_energy)

        #
        # 'Tracking'
        #

        # Get the FindableTracks collection from the event
        findable_tracks = event.get_collection('FindableTracks_recon')
        findable_dic = {}
        stub_dic = {}
        axial_dic = {}
        
        # Create a map between a sim particle and a findable track. This 
        # particular map will only contain tracks that can be found using the
        # 4s and 3s1a strategies.
        for findable_track in findable_tracks:
            if (findable_track.is4sFindable() or findable_track.is3s1aFindable()):
                findable_dic[findable_track.getSimParticle()] = findable_track
            elif findable_track.is2sFindable() & (findable_track.getSimParticle() != recoil_e):
                stub_dic[findable_track.getSimParticle()] = findable_track
            elif (findable_track.is2aFindable() 
                    & (not findable_track.is4sFindable())
                    & (not findable_track.is3s1aFindable()) 
                    & (not findable_track.is2sFindable())):
                axial_dic[findable_track.getSimParticle()] = findable_track

        self.ntuple['track_count'].append(len(findable_dic))
        self.ntuple['stub_count'].append(len(stub_dic))
        self.ntuple['axial_count'].append(len(axial_dic))

        #
        # Ecal
        #

        ecal_hits = event.get_collection('ecalDigis_recon')
        energy_counter = 0
        for ecal_hit in ecal_hits: 
            self.ntuple['ecal_hit_energy'].append(ecal_hit.getEnergy())
            energy_counter += ecal_hit.getEnergy()

        self.ntuple['total_ecal_energy'].append(energy_counter)

        #
        # Hcal
        #

        hcal_hits = event.get_collection('hcalDigis_recon')
        energy_counter = 0
        pe_counter = 0
        for hcal_hit in hcal_hits: 
            self.ntuple['hcal_hit_energy'].append(hcal_hit.getEnergy())
            energy_counter += hcal_hit.getEnergy()
            pe_counter += hcal_hit.getPE()

        self.ntuple['total_hcal_energy'].append(energy_counter)
        self.ntuple['total_hcal_pe'].append(pe_counter)

        # Get the collection of Ecal veto results from the event
        ecal_veto_results = event.get_collection('EcalVeto_recon')
        
        # Check if the event would have passed the Ecal veto
        for ecal_veto_result in ecal_veto_results: 
            if ecal_veto_result.passesVeto(): 
                self.ntuple['passes_ecal_veto'].append(1)
            else: self.ntuple['passes_ecal_veto'].append(0)
        
        # Get the collection of Ecal veto results from the event
        hcal_veto_results = event.get_collection('HcalVeto_recon')
        
        # Check if the event would have passed the Ecal veto
        for hcal_veto_result in hcal_veto_results: 
            if hcal_veto_result.passesVeto(): 
                self.ntuple['passes_hcal_veto'].append(1)
            else: self.ntuple['passes_hcal_veto'].append(0)
        
        #
        # Hit Level
        #
        
        recoil_sim_hits = event.get_collection('RecoilSimHits_sim')
        self.ntuple['total_recoil_hits'].append(recoil_sim_hits.GetEntriesFast())
        hit_counter = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        charge_counter = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

        for recoil_sim_hit in recoil_sim_hits:
            hit_counter[recoil_sim_hit.getLayerID() - 1] += 1
            charge_counter[recoil_sim_hit.getLayerID() - 1] += recoil_sim_hit.getEdep()
        
        for layer_n in xrange(0, 10): 
            self.ntuple['total_recoil_hits_l%s' % (layer_n + 1)].append(hit_counter[layer_n])
            self.ntuple['total_charge_l%s' % (layer_n + 1)].append(charge_counter[layer_n])

    def finalize(self) :

        for variable in self.variables: 
            self.ntuple[variable] = np.array(self.ntuple[variable])

        plt = Plotter.Plotter(self.file_prefix + '_plots')

        single_track = self.ntuple['track_count'] == 1
        ecal_veto = self.ntuple['passes_ecal_veto'] == 1
        hcal_veto = self.ntuple['passes_hcal_veto'] == 1
        basic_veto = single_track & ecal_veto & hcal_veto
        cuts = { 'Single track':single_track, 'Ecal veto':ecal_veto, 'Hcal veto':hcal_veto, 'Basic veto':basic_veto }

        plt.plot_hists([
                        self.ntuple['pn_gamma_energy'],
                        self.ntuple['pn_gamma_energy'][single_track],
                        self.ntuple['pn_gamma_energy'][hcal_veto],
                        self.ntuple['pn_gamma_energy'][ecal_veto],
                        self.ntuple['pn_gamma_energy'][basic_veto]
                       ],
                      np.linspace(0, 4000, 161),
                      labels=['All', 'Single track', 'Hcal veto', 'Ecal veto', 'Basic veto'], 
                      ylog=True,
                      label_loc=2, 
                      x_label='$E(\gamma)$ (MeV)')
       
        plt.create_root_hist('pn_gamma_energy', 
                             self.ntuple['pn_gamma_energy'], 
                             160, 0, 4000,
                             'E(#gamma) (MeV)', 
                             color=self.colors[0])
        index = 1
        for name, cut in cuts.iteritems():
            plt.create_root_hist('pn_gamma_energy - %s' % name, 
                                 self.ntuple['pn_gamma_energy'][cut], 
                                 160, 0, 4000,
                                 'E(#gamma) (MeV)', 
                                 color=self.colors[index])
            index += 1

        plt.plot_hists([
                        self.ntuple['pn_particle_mult'],
                        self.ntuple['pn_particle_mult'][single_track],
                        self.ntuple['pn_particle_mult'][hcal_veto],
                        self.ntuple['pn_particle_mult'][ecal_veto],
                        self.ntuple['pn_particle_mult'][basic_veto],
                       ],
                      np.linspace(0, 120, 121),
                      labels=['All', 'Single track', 'Hcal veto', 'Ecal veto', 'Basic veto'], 
                      ylog=True,
                      x_label='PN Multiplicity')
         
        plt.create_root_hist('pn_particle_mult', 
                             self.ntuple['pn_particle_mult'], 
                             120, 0, 120,
                             'PN Multiplicity', 
                             color=self.colors[0])
        index = 1
        for name, cut in cuts.iteritems():
            plt.create_root_hist('pn_particle_mult - %s' % name, 
                                 self.ntuple['pn_particle_mult'][cut], 
                                 120, 0, 120,
                                 'PN Multiplicity', 
                                 color=self.colors[index])
            index += 1

        plt.plot_hists([
                        self.ntuple['pn_interaction_z'],
                        self.ntuple['pn_interaction_z'][single_track],
                        self.ntuple['pn_interaction_z'][hcal_veto],
                        self.ntuple['pn_interaction_z'][ecal_veto],
                        self.ntuple['pn_interaction_z'][basic_veto]
                       ],
                      np.linspace(-100, 350, 451),
                      labels=['All', 'Single track', 'Hcal veto', 'Ecal veto', 'Basic veto'], 
                      ylog=True,
                      x_label='PN $\gamma$ Interaction Point z (mm)')

        plt.plot_hists([
                        self.ntuple['pn_interaction_z'],
                        self.ntuple['pn_interaction_z'][single_track],
                        self.ntuple['pn_interaction_z'][hcal_veto],
                        self.ntuple['pn_interaction_z'][ecal_veto],
                        self.ntuple['pn_interaction_z'][basic_veto]
                       ],
                      np.linspace(-1, 1, 451),
                      labels=['All', 'Single track', 'Hcal veto', 'Ecal veto', 'Basic veto'], 
                      ylog=True,
                      x_label='PN $\gamma$ Interaction Point z (mm)')

        plt.create_root_hist('pn_interaction_z', 
                             self.ntuple['pn_interaction_z'], 
                             450, -1, 1,
                             'PN $\gamma$ Interaction Point z (mm)', 
                             color=self.colors[0])
        index = 1
        for name, cut in cuts.iteritems():
            plt.create_root_hist('pn_interaction_z - %s' % name, 
                                 self.ntuple['pn_interaction_z'][cut], 
                                 450, -1, 1,
                                 'PN $\gamma$ Interaction Point z (mm)',
                                 color=self.colors[index])
            index += 1

        plt.plot_hists([
                        self.ntuple['lead_hadron_ke'],
                        self.ntuple['lead_hadron_ke'][single_track],
                        self.ntuple['lead_hadron_ke'][hcal_veto],
                        self.ntuple['lead_hadron_ke'][ecal_veto],
                        self.ntuple['lead_hadron_ke'][basic_veto]
                       ],
                      np.linspace(0, 4000, 161),
                      labels=['All', 'Single track', 'Hcal veto', 'Ecal veto', 'Basic veto'], 
                      ylog=True,
                      x_label='Leading Hadron Kinetic Energy (MeV)')

        plt.create_root_hist('lead_hadron_ke', 
                             self.ntuple['lead_hadron_ke'], 
                             160, 0, 4000,
                             'Leading Hadron Kinetic Energy (MeV)',
                             color=self.colors[0])
        index = 1
        for name, cut in cuts.iteritems():
            plt.create_root_hist('lead_hadron_ke - %s' % name, 
                                 self.ntuple['lead_hadron_ke'][cut], 
                                 160, 0, 4000,
                                 'Leading Hadron Kinetic Energy (MeV)', 
                                 color=self.colors[index])
            index += 1

        plt.plot_hists([
                        self.ntuple['lead_hadron_theta'],
                        self.ntuple['lead_hadron_theta'][single_track],
                        self.ntuple['lead_hadron_theta'][hcal_veto],
                        self.ntuple['lead_hadron_theta'][ecal_veto],
                        self.ntuple['lead_hadron_theta'][basic_veto]
                       ],
                      np.linspace(0, 180, 361),
                      labels=['All', 'Single track', 'Hcal veto', 'Ecal veto', 'Basic veto'], 
                      ylog=True,
                      x_label='Leading Hadron Theta (degrees)')

        plt.create_root_hist('lead_hadron_theta', 
                             self.ntuple['lead_hadron_theta'], 
                             360, 0, 180,
                             'Leading Hadron Theta (degrees)',
                             color=self.colors[0])
        index = 1
        for name, cut in cuts.iteritems():
            plt.create_root_hist('lead_hadron_theta - %s' % name, 
                                 self.ntuple['lead_hadron_theta'][cut], 
                                 360, 0, 180,
                                 'Leading Hadron Theta (degrees)',
                                 color=self.colors[index])
            index += 1

        plt.plot_hist2d(self.ntuple['lead_hadron_ke'], 
                        self.ntuple['lead_hadron_theta'], 
                        np.linspace(0, 4000, 161),
                        np.linspace(0, 180, 361),
                        x_label='Leading Hadron Kinetic Energy (MeV)',
                        y_label='Leading Hadron Theta (degrees)')

        plt.plot_hists([
                        self.ntuple['lead_hadron_p'],
                        self.ntuple['lead_hadron_p'][single_track],
                        self.ntuple['lead_hadron_p'][hcal_veto],
                        self.ntuple['lead_hadron_p'][ecal_veto],
                        self.ntuple['lead_hadron_p'][basic_veto]
                       ],
                      np.linspace(0, 4000, 161),
                      labels=['All', 'Single track', 'Hcal veto', 'Ecal veto', 'Basic veto'], 
                      ylog=True,
                      x_label='Leading Hadron Momentum (MeV)')

        plt.create_root_hist('lead_hadron_p', 
                             self.ntuple['lead_hadron_p'], 
                             160, 0, 4000,
                             'Leading Hadron Momentum (MeV)',
                             color=self.colors[0])
        index = 1
        for name, cut in cuts.iteritems():
            plt.create_root_hist('lead_hadron_p - %s' % name, 
                                 self.ntuple['lead_hadron_p'][cut], 
                                 160, 0, 4000,
                                 'Leading Hadron Momentum (MeV)',
                                 color=self.colors[index])
            index += 1

        plt.plot_hists([
                        self.ntuple['lead_hadron_pdgid'],
                        self.ntuple['lead_hadron_pdgid'][single_track],
                        self.ntuple['lead_hadron_pdgid'][hcal_veto],
                        self.ntuple['lead_hadron_pdgid'][ecal_veto],
                        self.ntuple['lead_hadron_pdgid'][basic_veto]
                       ],
                      np.linspace(0, 4000, 161),
                      labels=['All', 'Single track', 'Hcal veto', 'Ecal veto', 'Basic veto'], 
                      ylog=True,
                      x_label='Leading Hadron Momentum (MeV)')

        plt.create_root_hist('lead_hadron_pdgid', 
                             self.ntuple['lead_hadron_pdgid'], 
                             160, 0, 4000,
                             'Leading Hadron Momentum (MeV)',
                             color=self.colors[0])
        index = 1
        for name, cut in cuts.iteritems():
            plt.create_root_hist('lead_hadron_pdgid - %s' % name, 
                                 self.ntuple['lead_hadron_pdgid'][cut], 
                                 160, 0, 4000,
                                 'Leading Hadron Momentum (MeV)',
                                 color=self.colors[index])
            index += 1

        plt.plot_hists([
                        self.ntuple['lead_proton_ke'],
                        self.ntuple['lead_proton_ke'][single_track],
                        self.ntuple['lead_proton_ke'][hcal_veto],
                        self.ntuple['lead_proton_ke'][ecal_veto],
                        self.ntuple['lead_proton_ke'][basic_veto]
                       ],
                      np.linspace(0, 4000, 161),
                      labels=['All', 'Single track', 'Hcal veto', 'Ecal veto', 'Basic veto'], 
                      ylog=True,
                      x_label='Leading $p$ Kinetic Energy (MeV)')

        plt.create_root_hist('lead_proton_ke', 
                             self.ntuple['lead_proton_ke'], 
                             160, 0, 4000,
                             'Leading #p Kinetic Energy (MeV)',
                             color=self.colors[0])
        index = 1
        for name, cut in cuts.iteritems():
            plt.create_root_hist('lead_proton_ke - %s' % name, 
                                 self.ntuple['lead_proton_ke'][cut], 
                                 160, 0, 4000,
                                 'Leading #p Kinetic Energy (MeV)', 
                                 color=self.colors[index])
            index += 1

        plt.plot_hists([
                        self.ntuple['lead_proton_theta'],
                        self.ntuple['lead_proton_theta'][single_track],
                        self.ntuple['lead_proton_theta'][hcal_veto],
                        self.ntuple['lead_proton_theta'][ecal_veto],
                        self.ntuple['lead_proton_theta'][basic_veto]
                       ],
                      np.linspace(0, 180, 361),
                      labels=['All', 'Single track', 'Hcal veto', 'Ecal veto', 'Basic veto'], 
                      ylog=True,
                      x_label='Leading $p$ Theta (degrees)')

        plt.create_root_hist('lead_proton_theta', 
                             self.ntuple['lead_proton_theta'], 
                             360, 0, 180,
                             'Leading #p Theta (degrees)', 
                             color=self.colors[0])
        index = 1
        for name, cut in cuts.iteritems():
            plt.create_root_hist('lead_proton_theta - %s' % name, 
                                 self.ntuple['lead_proton_theta'][cut], 
                                 360, 0, 180,
                                 'Leading #p Theta (degrees)',
                                 color=self.colors[index])
            index += 1

        plt.plot_hist2d(self.ntuple['lead_proton_ke'], 
                        self.ntuple['lead_proton_theta'], 
                        np.linspace(0, 4000, 161),
                        np.linspace(0, 180, 361),
                        x_label='Leading $p$ Kinetic Energy (MeV)',
                        y_label='Leading $p$ Theta (degrees)')

        plt.plot_hists([
                        self.ntuple['lead_proton_p'],
                        self.ntuple['lead_proton_p'][single_track],
                        self.ntuple['lead_proton_p'][hcal_veto],
                        self.ntuple['lead_proton_p'][ecal_veto],
                        self.ntuple['lead_proton_p'][basic_veto]
                       ],
                      np.linspace(0, 4000, 161),
                      labels=['All', 'Single track', 'Hcal veto', 'Ecal veto', 'Basic veto'], 
                      ylog=True,
                      x_label='Leading $p$ Momentum (MeV)')

        plt.create_root_hist('lead_proton_p', 
                             self.ntuple['lead_proton_p'], 
                             160, 0, 4000,
                             'Leading #p Momentum (MeV)', 
                             color=self.colors[0])
        index = 1
        for name, cut in cuts.iteritems():
            plt.create_root_hist('lead_proton_p - %s' % name, 
                                 self.ntuple['lead_proton_p'][cut], 
                                 160, 0, 4000,
                                 'Leading #p Momentum (MeV)', 
                                 color=self.colors[index])
            index += 1

        plt.plot_hists([
                        self.ntuple['lead_neutron_ke'],
                        self.ntuple['lead_neutron_ke'][single_track],
                        self.ntuple['lead_neutron_ke'][hcal_veto],
                        self.ntuple['lead_neutron_ke'][ecal_veto],
                        self.ntuple['lead_neutron_ke'][basic_veto]
                       ],
                      np.linspace(0, 4000, 161),
                      labels=['All', 'Single track', 'Hcal veto', 'Ecal veto', 'Basic veto'], 
                      ylog=True,
                      x_label='Leading $n$ Kinetic Energy')

        plt.create_root_hist('lead_neutron_ke', 
                             self.ntuple['lead_neutron_ke'], 
                             160, 0, 4000,
                             'Leading #n Kinetic Energy (MeV)', 
                             color=self.colors[0])
        index = 1
        for name, cut in cuts.iteritems():
            plt.create_root_hist('lead_neutron_ke - %s' % name, 
                                 self.ntuple['lead_neutron_ke'][cut], 
                                 160, 0, 4000,
                                 'Leading #n Kinetic Energy (MeV)',
                                 color=self.colors[index])
            index += 1

        plt.plot_hists([
                        self.ntuple['lead_neutron_theta'],
                        self.ntuple['lead_neutron_theta'][single_track],
                        self.ntuple['lead_neutron_theta'][hcal_veto],
                        self.ntuple['lead_neutron_theta'][ecal_veto],
                        self.ntuple['lead_neutron_theta'][basic_veto]
                       ],
                      np.linspace(0, 180, 361),
                      labels=['All', 'Single track', 'Hcal veto', 'Ecal veto', 'Basic veto'], 
                      ylog=True,
                      x_label='Leading $n$ Theta (degrees)')
        
        plt.create_root_hist('lead_neutron_theta', 
                             self.ntuple['lead_neutron_theta'], 
                             360, 0, 180,
                             'Leading #n Theta (degrees)', 
                             color=self.colors[0])
        index = 1
        for name, cut in cuts.iteritems():
            plt.create_root_hist('lead_neutron_theta - %s' % name, 
                                 self.ntuple['lead_neutron_theta'][cut], 
                                 360, 0, 180,
                                 'Leading #n Theta (degrees)',
                                 color=self.colors[index])
            index += 1


        plt.plot_hist2d(self.ntuple['lead_neutron_ke'], 
                        self.ntuple['lead_neutron_theta'], 
                        np.linspace(0, 4000, 161),
                        np.linspace(0, 180, 361),
                        x_label='Leading $n$ Kinetic Energy (MeV)',
                        y_label='Leading $n$ Theta (degrees)')
        
        plt.plot_hists([
                        self.ntuple['lead_neutron_p'],
                        self.ntuple['lead_neutron_p'][single_track],
                        self.ntuple['lead_neutron_p'][hcal_veto],
                        self.ntuple['lead_neutron_p'][ecal_veto],
                        self.ntuple['lead_neutron_p'][basic_veto]
                       ],
                      np.linspace(0, 4000, 161),
                      labels=['All', 'Single track', 'Hcal veto', 'Ecal veto', 'Basic veto'], 
                      ylog=True,
                      x_label='Leading $n$ Momentum (MeV)')

        plt.create_root_hist('lead_neutron_p', 
                             self.ntuple['lead_neutron_p'], 
                             160, 0, 4000,
                             'Leading #n Momentum (MeV)',
                             color=self.colors[0])
        index = 1
        for name, cut in cuts.iteritems():
            plt.create_root_hist('lead_neutron_p - %s' % name, 
                                 self.ntuple['lead_neutron_p'][cut], 
                                 160, 0, 4000,
                                 'Leading #n Momentum (MeV)',
                                 color=self.colors[index])
            index += 1

        plt.plot_hists([
                        self.ntuple['lead_pion_ke'],
                        self.ntuple['lead_pion_ke'][single_track],
                        self.ntuple['lead_pion_ke'][hcal_veto],
                        self.ntuple['lead_pion_ke'][ecal_veto],
                        self.ntuple['lead_pion_ke'][basic_veto]
                       ],
                      np.linspace(0, 4000, 161),
                      labels=['All', 'Single track', 'Hcal veto', 'Ecal veto', 'Basic veto'], 
                      ylog=True,
                      x_label='Leading $\pi$ Kinetic Energy')
       
        plt.create_root_hist('lead_pion_ke', 
                             self.ntuple['lead_pion_ke'], 
                             160, 0, 4000,
                             'Leading #pi Kinetic Energy (MeV)',
                             color=self.colors[0])
        index = 1
        for name, cut in cuts.iteritems():
            plt.create_root_hist('lead_pion_ke - %s' % name, 
                                 self.ntuple['lead_pion_ke'][cut], 
                                 160, 0, 4000,
                                 'Leading #pi Kinetic Energy (MeV)',
                                 color=self.colors[index])
            index += 1

        plt.plot_hists([
                        self.ntuple['lead_pion_theta'],
                        self.ntuple['lead_pion_theta'][single_track],
                        self.ntuple['lead_pion_theta'][hcal_veto],
                        self.ntuple['lead_pion_theta'][ecal_veto],
                        self.ntuple['lead_pion_theta'][basic_veto]
                       ],
                      np.linspace(0, 180, 361),
                      labels=['All', 'Single track', 'Hcal veto', 'Ecal veto', 'Basic veto'], 
                      ylog=True,
                      x_label='Leading $\pi$ Theta (degrees)')

        
        plt.create_root_hist('lead_pion_theta', 
                             self.ntuple['lead_pion_theta'], 
                             360, 0, 180,
                             'Leading #pi Theta (degrees)', 
                             color=self.colors[0])
        index = 1
        for name, cut in cuts.iteritems():
            plt.create_root_hist('lead_pion_theta - %s' % name, 
                                 self.ntuple['lead_pion_theta'][cut], 
                                 360, 0, 180,
                                 'Leading #pi Theta (degrees)', 
                                 color=self.colors[index])
            index += 1

        plt.plot_hist2d(self.ntuple['lead_pion_ke'], 
                        self.ntuple['lead_pion_theta'], 
                        np.linspace(0, 4000, 161),
                        np.linspace(0, 180, 361),
                        x_label='Leading $\pi$ Kinetic Energy (MeV)',
                        y_label='Leading $\pi$ Theta (degrees)')

        plt.plot_hists([
                        self.ntuple['lead_pion_p'],
                        self.ntuple['lead_pion_p'][single_track],
                        self.ntuple['lead_pion_p'][hcal_veto],
                        self.ntuple['lead_pion_p'][ecal_veto],
                        self.ntuple['lead_pion_p'][basic_veto]
                       ],
                      np.linspace(0, 4000, 161),
                      labels=['All', 'Single track', 'Hcal veto', 'Ecal veto', 'Basic veto'], 
                      ylog=True,
                      x_label='Leading $\pi$ Momentum (MeV)')

        plt.create_root_hist('lead_pion_p', 
                             self.ntuple['lead_pion_p'], 
                             160, 0, 4000,
                             'Leading #pi Momentum (MeV)', 
                             color=self.colors[0])
        index = 1
        for name, cut in cuts.iteritems():
            plt.create_root_hist('lead_pion_p - %s' % name, 
                                 self.ntuple['lead_pion_p'][cut], 
                                 160, 0, 4000,
                                 'Leading #pi Momentum (MeV)',
                                 color=self.colors[index])
            index += 1

        plt.plot_hists([
                        self.ntuple['max_w'],
                        self.ntuple['max_w'][single_track],
                        self.ntuple['max_w'][hcal_veto],
                        self.ntuple['max_w'][ecal_veto],
                        self.ntuple['max_w'][basic_veto]
                       ],
                      np.linspace(0, 5000, 251),
                      labels=['All', 'Single track', 'Hcal veto', 'Ecal veto', 'Basic veto'], 
                      ylog=True,
                      x_label='Max W (MeV)')

        plt.create_root_hist('max_w', 
                             self.ntuple['max_w'], 
                             250, 0, 5000,
                             'Max W (MeV)', 
                             color=self.colors[0])
        
        plt.plot_hists([
                        self.ntuple['max_w_theta'],
                        self.ntuple['max_w_theta'][single_track],
                        self.ntuple['max_w_theta'][hcal_veto],
                        self.ntuple['max_w_theta'][ecal_veto],
                        self.ntuple['max_w_theta'][basic_veto]
                       ],
                      np.linspace(0, 180, 361),
                      labels=['All', 'Single track', 'Hcal veto', 'Ecal veto', 'Basic veto'], 
                      ylog=True,
                      x_label='Max W (MeV)')

        plt.create_root_hist('max_w_theta', 
                             self.ntuple['max_w_theta'], 
                             360, 0, 180,
                             'Max W (MeV)', 
                             color=self.colors[0])
        index = 1
        for name, cut in cuts.iteritems():
            plt.create_root_hist('max_w_theta - %s' % name, 
                                 self.ntuple['max_w_theta'][cut], 
                                 360, 0, 180,
                                 'Max W (MeV)',
                                 color=self.colors[index])
            index += 1


        plt.plot_hists([
                        self.ntuple['max_w'][self.ntuple['max_w_theta'] > 100]
                       ],
                      np.linspace(0, 5000, 251),
                      labels=['All'],
                      ylog=True,
                      x_label='Max W($\theta >$ 100) (MeV)')

        plt.create_root_hist('max_w_theta_cut', 
                             self.ntuple['max_w'][self.ntuple['max_w_theta'] > 100], 
                             250, 0, 5000,
                             'Max W(#theta > 100) (MeV)', 
                             color=self.colors[0])
        index = 1
        for name, cut in cuts.iteritems():
            plt.create_root_hist('max_w_theta_cut - %s' % name, 
                                 self.ntuple['max_w'][(self.ntuple['max_w_theta'] > 100) & cut], 
                                 360, 0, 180,
                                'Max W(#theta > 100) (MeV)', 
                                 color=self.colors[index])

        plt.plot_hists([
                        self.ntuple['track_count'],
                        self.ntuple['track_count'][single_track],
                        self.ntuple['track_count'][hcal_veto],
                        self.ntuple['track_count'][ecal_veto],
                        self.ntuple['track_count'][basic_veto]
                       ], 
                       np.linspace(0, 10, 11),
                       labels=['All', 'Single track', 'Hcal veto', 'Ecal veto', 'Basic veto'], 
                       ylog=True,
                       x_label='Track Multiplicity')

        plt.create_root_hist('track_count', 
                             self.ntuple['track_count'], 
                             10, 0, 10,
                             'Track Multiplicity', 
                             color=self.colors[0])
        index = 1
        for name, cut in cuts.iteritems():
            plt.create_root_hist('track_count - %s' % name, 
                                 self.ntuple['track_count'][cut], 
                                 10, 0, 10,
                                 'Track Multiplicity', 
                                 color=self.colors[index])
            index += 1

        plt.plot_hists([
                        self.ntuple['stub_count'],
                        self.ntuple['stub_count'][single_track],
                        self.ntuple['stub_count'][hcal_veto],
                        self.ntuple['stub_count'][ecal_veto],
                        self.ntuple['stub_count'][basic_veto]
                       ],  
                       np.linspace(0, 10, 11),
                       labels=['All', 'Single track', 'Hcal veto', 'Ecal veto', 'Basic veto'], 
                       ylog=True,
                       x_label='Stub Multiplicity')

        plt.create_root_hist('stub_count', 
                             self.ntuple['stub_count'], 
                             10, 0, 10,
                             'Stub Multiplicity', 
                             color=self.colors[0])
        index = 1
        for name, cut in cuts.iteritems():
            plt.create_root_hist('stub_count - %s' % name, 
                                 self.ntuple['stub_count'][cut], 
                                 10, 0, 10,
                                 'Stub Multiplicity', 
                                 color=self.colors[index])
            index += 1

        plt.plot_hists([
                        self.ntuple['axial_count'],
                        self.ntuple['axial_count'][single_track],
                        self.ntuple['axial_count'][hcal_veto],
                        self.ntuple['axial_count'][ecal_veto],
                        self.ntuple['axial_count'][basic_veto]
                       ],  
                       np.linspace(0, 10, 11),
                       labels=['All', 'Single track', 'Hcal veto', 'Ecal veto', 'Basic veto'], 
                       ylog=True,
                       x_label='Axial Multiplicity')

        plt.create_root_hist('axial_count', 
                             self.ntuple['axial_count'], 
                             10, 0, 10,
                             'Axial Multiplicity', 
                             color=self.colors[0])
        index = 1
        for name, cut in cuts.iteritems():
            plt.create_root_hist('axial_count - %s' % name, 
                                 self.ntuple['axial_count'][cut], 
                                 10, 0, 10,
                                 'Axial Multiplicity', 
                                 color=self.colors[index])
            index += 1

        plt.plot_hists([
                        self.ntuple['down_tp_energy'],
                        self.ntuple['down_tp_energy'][single_track],
                        self.ntuple['down_tp_energy'][hcal_veto],
                        self.ntuple['down_tp_energy'][ecal_veto],
                        self.ntuple['down_tp_energy'][basic_veto]
                       ], 
                       np.linspace(0, 100, 201),
                       labels=['All', 'Single track', 'Hcal veto', 'Ecal veto', 'Basic veto'], 
                       ylog=True, 
                       x_label='Energy Deposited in Downstream Trigger Pad (MeV)')
        
        plt.create_root_hist('down_tp_energy', 
                             self.ntuple['down_tp_energy'], 
                             200, 0, 100,
                             'Energy Deposited in Downstream Trigger Pad (MeV)', 
                             color=self.colors[0])
        index = 1
        for name, cut in cuts.iteritems():
            plt.create_root_hist('down_tp_energy - %s' % name, 
                                 self.ntuple['down_tp_energy'][cut], 
                                 200, 0, 100,
                                 'Energy Deposited in Downstream Trigger Pad (MeV)', 
                                 color=self.colors[index])
            index += 1
        
        plt.plot_hists([
                        self.ntuple['up_tp_energy'],
                        self.ntuple['up_tp_energy'][single_track],
                        self.ntuple['up_tp_energy'][hcal_veto],
                        self.ntuple['up_tp_energy'][ecal_veto],
                        self.ntuple['up_tp_energy'][basic_veto]
                       ],
                       np.linspace(0, 100, 201),
                       labels=['All', 'Single track', 'Hcal veto', 'Ecal veto', 'Basic veto'], 
                       ylog=True, 
                       x_label='Energy Deposited in Upstream Trigger Pad (MeV)')
        
        plt.create_root_hist('up_tp_energy', 
                             self.ntuple['up_tp_energy'], 
                             200, 0, 100,
                             'Energy Deposited in Upstream Trigger Pad (MeV)', 
                             color=self.colors[0])
        index = 1
        for name, cut in cuts.iteritems():
            plt.create_root_hist('up_tp_energy - %s' % name, 
                                 self.ntuple['up_tp_energy'][cut], 
                                 200, 0, 100,
                                 'Energy Deposited in Upstream Trigger Pad (MeV)', 
                                 color=self.colors[index])
            index += 1
       
        plt.plot_hist2d(self.ntuple['up_tp_energy'], 
                        self.ntuple['down_tp_energy'], 
                        np.linspace(0, 100, 201),
                        np.linspace(0, 100, 201),
                        x_label='Energy Deposited in Upstream Trigger Pad (MeV)', 
                        y_label='Energy Deposited in Downstream Trigger Pad (MeV)') 

        plt.plot_hists([
                        self.ntuple['total_ecal_energy'],
                        self.ntuple['total_ecal_energy'][single_track],
                        self.ntuple['total_ecal_energy'][hcal_veto],
                        self.ntuple['total_ecal_energy'][ecal_veto],
                        self.ntuple['total_ecal_energy'][basic_veto]
                       ], 
                       np.linspace(0, 140, 141),
                       labels=['All', 'Single track', 'Hcal veto', 'Ecal veto', 'Basic veto'], 
                       ylog=True,
                       x_label='Total Energy Deposited in Ecal Si (MeV)')
        
        plt.create_root_hist('total_ecal_energy', 
                             self.ntuple['total_ecal_energy'], 
                             140, 0, 140,
                             'Energy Deposited in Upstream Trigger Pad (MeV)', 
                             color=self.colors[0])
        index = 1
        for name, cut in cuts.iteritems():
            plt.create_root_hist('total_ecal_energy - %s' % name, 
                                 self.ntuple['total_ecal_energy'][cut], 
                                 140, 0, 140,
                                 'Energy Deposited in Upstream Trigger Pad (MeV)',
                                 color=self.colors[index])
            index += 1
        
        plt.plot_hists([
                        self.ntuple['ecal_hit_energy'],
                        self.ntuple['ecal_hit_energy'][single_track],
                        self.ntuple['ecal_hit_energy'][hcal_veto],
                        self.ntuple['ecal_hit_energy'][ecal_veto],
                        self.ntuple['ecal_hit_energy'][basic_veto]
                       ], 
                       np.linspace(0, 100, 201),
                       labels=['All', 'Single track', 'Hcal veto', 'Ecal veto', 'Basic veto'], 
                       ylog=True,
                       x_label='Readout Hit Energy Deposited in Ecal Si (MeV)')
        
        plt.create_root_hist('ecal_hit_energy', 
                             self.ntuple['ecal_hit_energy'], 
                             200, 0, 100,
                             'Readout Hit Energy Deposited in Ecal Si (MeV)',
                             color=self.colors[0])
        index = 1
        for name, cut in cuts.iteritems():
            plt.create_root_hist('ecal_hit_energy - %s' % name, 
                                 self.ntuple['ecal_hit_energy'][cut], 
                                 200, 0, 100,
                                 'Readout Hit Energy Deposited in Ecal Si (MeV)', 
                                 color=self.colors[index])
            index += 1

        plt.plot_hists([
                        self.ntuple['hcal_hit_energy'],
                        self.ntuple['hcal_hit_energy'][single_track],
                        self.ntuple['hcal_hit_energy'][hcal_veto],
                        self.ntuple['hcal_hit_energy'][ecal_veto],
                        self.ntuple['hcal_hit_energy'][basic_veto]
                       ], 
                       np.linspace(0, 150, 301),
                       labels=['All', 'Single track', 'Hcal veto', 'Ecal veto', 'Basic veto'], 
                       ylog=True,
                       x_label='Readout Hit Energy Deposited in Hcal Scint (MeV)')
        
        plt.create_root_hist('hcal_hit_energy', 
                             self.ntuple['hcal_hit_energy'], 
                             300, 0, 150,
                             'Readout Hit Energy Deposited in Ecal Si (MeV)', 
                             color=self.colors[0])
        index = 1
        for name, cut in cuts.iteritems():
            plt.create_root_hist('hcal_hit_energy - %s' % name, 
                                 self.ntuple['hcal_hit_energy'][cut], 
                                 300, 0, 150,
                                 'Readout Hit Energy Deposited in Ecal Si (MeV)',
                                 color=self.colors[index])
            index += 1

        plt.plot_hists([
                        self.ntuple['total_hcal_energy'],
                        self.ntuple['total_hcal_energy'][single_track],
                        self.ntuple['total_hcal_energy'][hcal_veto],
                        self.ntuple['total_hcal_energy'][ecal_veto],
                        self.ntuple['total_hcal_energy'][basic_veto]
                       ],
                       np.linspace(0, 400, 401),
                       labels=['All', 'Single track', 'Hcal veto', 'Ecal veto', 'Basic veto'], 
                       ylog=True,
                       x_label='Total Energy Deposited in Hcal Scint (MeV)')

        plt.create_root_hist('total_hcal_energy', 
                             self.ntuple['total_hcal_energy'], 
                             400, 0, 400,
                             'Total Energy Deposited in Hcal Scint (MeV)', 
                             color=self.colors[0])
        index = 1
        for name, cut in cuts.iteritems():
            plt.create_root_hist('total_hcal_energy - %s' % name, 
                                 self.ntuple['total_hcal_energy'][cut], 
                                 400, 0, 400,
                                 'Total Energy Deposited in Hcal Scint (MeV)',
                                 color=self.colors[index])
            index += 1

        
        plt.plot_hists([
                        self.ntuple['total_recoil_hits'],
                        self.ntuple['total_recoil_hits'][single_track],
                        self.ntuple['total_recoil_hits'][hcal_veto],
                        self.ntuple['total_recoil_hits'][ecal_veto],
                        self.ntuple['total_recoil_hits'][basic_veto]
                       ],
                       np.linspace(0, 150, 151),
                       labels=['All', 'Single track', 'Hcal veto', 'Ecal veto', 'Basic veto'], 
                       ylog=True,
                       x_label='Recoil Hit Multiplicity')

        plt.create_root_hist('total_recoil_hits', 
                             self.ntuple['total_recoil_hits'], 
                             150, 0, 150,
                             'Recoil Hit Multiplicity',
                             color=self.colors[0])
        index = 1
        for name, cut in cuts.iteritems():
            plt.create_root_hist('total_recoil_hits - %s' % name, 
                                 self.ntuple['total_recoil_hits'][cut], 
                                 150, 0, 150,
                                 'Recoil Hit Multiplicity', 
                                 color=self.colors[index])
            index += 1

        for layer_n in xrange(0, 10): 
            plt.plot_hists([
                self.ntuple['total_recoil_hits_l%s' % (layer_n + 1)],
                self.ntuple['total_recoil_hits_l%s' % (layer_n + 1)][single_track],
                self.ntuple['total_recoil_hits_l%s' % (layer_n + 1)][hcal_veto],
                self.ntuple['total_recoil_hits_l%s' % (layer_n + 1)][ecal_veto],
                self.ntuple['total_recoil_hits_l%s' % (layer_n + 1)][basic_veto]
                ],
                np.linspace(0, 150, 151),
                labels=['All', 'Single track', 'Hcal veto', 'Ecal veto', 'Basic veto'], 
                ylog=True,
                x_label='Recoil Hit Layer %s Multiplicity' % (layer_n + 1))

            plt.plot_hists([
                self.ntuple['total_charge_l%s' % (layer_n + 1)],
                self.ntuple['total_charge_l%s' % (layer_n + 1)][single_track],
                self.ntuple['total_charge_l%s' % (layer_n + 1)][hcal_veto],
                self.ntuple['total_charge_l%s' % (layer_n + 1)][ecal_veto],
                self.ntuple['total_charge_l%s' % (layer_n + 1)][basic_veto]
                ],
                np.linspace(0, 50, 151),
                labels=['All', 'Single track', 'Hcal veto', 'Ecal veto', 'Basic veto'], 
                ylog=True,
                x_label='Total Charge Layer %s' % (layer_n + 1))

            plt.create_root_hist('total_recoil_hits_l%s' % (layer_n + 1), 
                                 self.ntuple['total_recoil_hits_l%s' % (layer_n + 1)],
                                 150, 0, 150,
                                 'Recoil Hit Layer %s Multiplicity' % (layer_n + 1), 
                                 color=self.colors[0])
            plt.create_root_hist('total_charge_l%s' % (layer_n + 1), 
                                 self.ntuple['total_charge_l%s' % (layer_n + 1)],
                                 150, 0, 50,
                                 'Total Charge Layer %s' % (layer_n + 1),
                                 color=self.colors[0])
            index = 0
            for name, cut in cuts.iteritems():
                plt.create_root_hist('total_recoil_hits_l%s - %s' % ((layer_n + 1), name),
                                     self.ntuple['total_recoil_hits_l%s' % (layer_n + 1)][cut],
                                     150, 0, 150,
                                     'Recoil Hit Layer %s Multiplicity' % (layer_n + 1), 
                                     color=self.colors[index])
                plt.create_root_hist('total_charge_l%s - %s' % ((layer_n + 1), name),
                                     self.ntuple['total_charge_l%s' % (layer_n + 1)][cut],
                                     150, 0, 50,
                                     'Total Charge Layer %s' % (layer_n + 1),
                                     color=self.colors[index])
                index += 1

        plt.close()
       
        output = file(self.file_prefix + "_results.txt", 'w')
        cut = self.ntuple['track_count'] == 1
        pass_track_veto = len(self.ntuple['track_count'][cut])
        pass_track_veto_frac = pass_track_veto/len(self.ntuple['track_count'])
        output.write('Events that pass track veto: %s/%s = %s\n' % (
            pass_track_veto, 
            len(self.ntuple['track_count']), pass_track_veto_frac))

        cut = self.ntuple['stub_count'] == 1
        pass_stub_veto = len(self.ntuple['stub_count'][cut])
        pass_stub_veto_frac = pass_stub_veto/len(self.ntuple['stub_count'])
        output.write('Events passing stub veto: %s/%s = %s\n' % (
            pass_stub_veto, 
            len(self.ntuple['stub_count']), pass_stub_veto_frac))

        cut = self.ntuple['passes_ecal_veto'] == 1
        pass_ecal_veto = len(self.ntuple['passes_ecal_veto'][cut])
        pass_ecal_veto_frac = pass_ecal_veto/len(self.ntuple['passes_ecal_veto'])
        output.write('Events passing Ecal veto: %s/%s = %s\n' % (
            pass_ecal_veto, 
            len(self.ntuple['passes_ecal_veto']), pass_ecal_veto_frac))

        cut = self.ntuple['passes_hcal_veto'] == 1
        pass_hcal_veto = len(self.ntuple['passes_hcal_veto'][cut])
        pass_hcal_veto_frac = pass_hcal_veto/len(self.ntuple['passes_hcal_veto'])
        output.write('Events passing Hcal veto: %s/%s = %s\n' % (
            pass_hcal_veto, 
            len(self.ntuple['passes_hcal_veto']), pass_hcal_veto_frac))

        cut = ((self.ntuple['track_count'] == 1) 
                & (self.ntuple['passes_hcal_veto'] == 1)
                & (self.ntuple['passes_ecal_veto'] == 1))
        pass_veto = len(self.ntuple['passes_hcal_veto'][cut])
        pass_veto_frac = pass_veto/len(self.ntuple['passes_hcal_veto'])
        output.write('Events passing veto: %s/%s = %s\n' % (pass_veto, 
                len(self.ntuple['passes_hcal_veto']), pass_veto_frac))
       
        output.write('Events failing veto\n')
        for event_n in self.ntuple['events'][cut]:
            output.write('Event: %s\n' % event_n)

        output.close()
        '''
        output = file(self.file_prefix + "_results.txt", 'w')
        output.write( "Total events processed: %s\n" % len(self.ntuple['track_count']))
        
        pass_tp_veto = len(self.ntuple['down_tp_energy'][self.ntuple['down_tp_energy'] < 0.75])
        pass_tp_veto_per = pass_tp_veto/len(self.ntuple['down_tp_energy'])
        output.write( 'Pass down TP veto: %s (%s)\n' % (pass_tp_veto, pass_tp_veto_per*100))

        pass_tp_up_veto = len(self.ntuple['up_tp_energy'][self.ntuple['up_tp_energy'] < 0.5])
        pass_tp_up_veto_per = pass_tp_up_veto/len(self.ntuple['up_tp_energy'])
        output.write( 'Pass up TP veto: %s (%s)\n' % (pass_tp_up_veto, pass_tp_up_veto_per*100))

        output.write( 'Events that pass track veto: %s (%s)\n' % (pass_track_veto, pass_track_veto_per*100))


        pass_bs_veto = len(self.ntuple['axial_count'][self.ntuple['axial_count'] == 0])
        pass_bs_veto_per = pass_bs_veto/len(self.ntuple['axial_count'])
        output.write( 'Events passing bs veto: %s (%s)\n' % (pass_bs_veto, pass_bs_veto_per*100))

        pass_hcal_veto = len(self.ntuple['total_hcal_pe'][self.ntuple['total_hcal_pe'] <= 8])
        pass_hcal_veto_per = pass_hcal_veto/len(self.ntuple['total_hcal_pe'])
        output.write( 'Events passing hcal veto: %s (%s)\n' % (pass_hcal_veto, pass_hcal_veto_per*100))

        pass_ecal_sum_cut = len(self.ecal_summed_det[self.ecal_summed_det < 25])
        pass_ecal_sum_cut_per = pass_ecal_sum_cut/len(self.ecal_summed_det)
        output.write( 'Events passing ecal sum cut: %s (%s)\n' % (pass_ecal_sum_cut, pass_ecal_sum_cut_per*100))
        veto = ((self.ntuple['down_tp_energy'] < 0.75) & (self.ntuple['track_count'] == 1) & (self.ntuple['up_tp_energy'] < 0.5) 
                & (self.ntuple['stub_count'] == 0) & (self.ntuple['axial_count'] == 0))

        output.write( "Event passing upstream veto: %s\n" % len(self.ntuple['track_count'][veto]))
        
        veto = ((self.ntuple['down_tp_energy'] < 0.75) & (self.ntuple['track_count'] == 1) & (self.ntuple['up_tp_energy'] < 0.5) 
                & (self.ntuple['stub_count'] == 0) & (self.ntuple['axial_count'] == 0) & (self.ntuple['total_hcal_pe'] <= 8))
        
        output.write( "Event passing upstream veto + hcal: %s\n" % len(self.ntuple['track_count'][veto]))
        
        veto = ((self.ntuple['down_tp_energy'] < 0.75) & (self.ntuple['track_count'] == 1) & (self.ntuple['up_tp_energy'] < 0.5) 
                & (self.ntuple['stub_count'] == 0) & (self.ntuple['axial_count'] == 0) & (self.ntuple['total_hcal_pe'] <= 8)
                & (self.ntuple['passes_ecal_veto'] == 1))
        
        output.write( "Event passing upstream veto + hcal + ecal bdt: %s\n" % len(self.ntuple['track_count'][veto]))
        
        veto = ((self.ntuple['down_tp_energy'] < 0.75) & (self.ntuple['track_count'] == 1) & (self.ntuple['up_tp_energy'] < 0.5) 
                & (self.ntuple['stub_count'] == 0) & (self.ntuple['axial_count'] == 0) & (self.ntuple['total_hcal_pe'] <= 8)
                & (self.ntuple['passes_ecal_veto'] == 1) & (self.ecal_summed_det < 25))
        
        output.write( "Event passing upstream veto + hcal + ecal bdt + ecal sum: %s\n" % len(self.ntuple['track_count'][veto]))
        output.close()
        '''
