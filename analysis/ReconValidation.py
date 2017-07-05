from __future__ import division

import math
import ROOT as r 
import numpy as np
import Plotter

from numpy import linalg as la
from scipy.stats import norm

class ReconValidation(object): 

    def __init__(self): 
        self.initialize()
    
    def is_recoil(self, particle) : 
        return (particle.getPdgID() == 11) & (particle.getParentCount() == 0)

    def initialize(self) :

        self.ntuple = {}
        self.variables = [
                'events',
                'track_count', 'stub_count', 'axial_count', 
                'ecal_hit_energy', 'total_ecal_energy', 'passes_ecal_veto',
                'hcal_hit_energy', 'total_hcal_energy', 'total_hcal_pe',
                'passes_hcal_veto', 'up_tp_energy', 'down_tp_energy',
        ]
        self.variables = [
                'events',
                'track_count', 'stub_count', 'axial_count', 
                'ecal_hit_energy', 'total_ecal_energy', 'passes_ecal_veto',
                'hcal_hit_energy', 'total_hcal_energy', 'total_hcal_pe',
                'passes_hcal_veto', 'up_tp_energy', 'down_tp_energy',
                'total_recoil_hits'
        ]
        for layer_n in xrange(0, 10):
            self.variables.append('total_recoil_hits_l%s' % (layer_n + 1))
            self.variables.append('total_charge_l%s' % (layer_n + 1))

        for variable in self.variables: 
            self.ntuple[variable] = []
        
        self.event_count = -1
    
    def process(self, event) :
        
        self.event_count += 1
        self.ntuple['events'].append(event.get_event_number())
        
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

        plt = Plotter.Plotter('recon_validation')
       
        plt.plot_hists([
                        self.ntuple['track_count'],
                       ], 
                       np.linspace(0, 10, 11),
                       labels=['All'],
                       ylog=True,
                       x_label='Track Multiplicity')

        plt.create_root_hist('track_count', 
                             self.ntuple['track_count'], 
                             10, 0, 10,
                             'Track Multiplicity', 
                             color=r.kRed+2)
        
        plt.plot_hists([
                        self.ntuple['stub_count'],
                       ],  
                       np.linspace(0, 10, 11),
                       labels=['All'],
                       ylog=True,
                       x_label='Stub Multiplicity')

        plt.create_root_hist('stub_count', 
                             self.ntuple['stub_count'], 
                             10, 0, 10,
                             'Stub Multiplicity',
                             color=r.kRed+2)
        
        plt.plot_hists([
                        self.ntuple['axial_count'],
                       ],  
                       np.linspace(0, 10, 11),
                       labels=['All'],
                       ylog=True,
                       x_label='Axial Multiplicity')

        plt.create_root_hist('axial_count', 
                             self.ntuple['axial_count'], 
                             10, 0, 10,
                             'Axial Multiplicity',
                             color=r.kRed+2)
        
        plt.plot_hists([
                        self.ntuple['down_tp_energy'],
                       ], 
                       np.linspace(0, 100, 201),
                       labels=['All'],
                       ylog=True, 
                       x_label='Energy Deposited in Downstream Trigger Pad (MeV)')
        
        plt.create_root_hist('down_tp_energy', 
                             self.ntuple['down_tp_energy'], 
                             200, 0, 100,
                             'Energy Deposited in Downstream Trigger Pad (MeV)',
                             color=r.kRed+2)
        
        plt.plot_hists([
                        self.ntuple['up_tp_energy'],
                       ],
                       np.linspace(0, 100, 201),
                       labels=['All'],
                       ylog=True, 
                       x_label='Energy Deposited in Upstream Trigger Pad (MeV)')

        plt.create_root_hist('up_tp_energy', 
                             self.ntuple['up_tp_energy'], 
                             200, 0, 100,
                             'Energy Deposited in Upstream Trigger Pad (MeV)',
                             color=r.kRed+2)
        
        plt.plot_hists([
                        self.ntuple['total_ecal_energy'],
                       ], 
                       np.linspace(0, 140, 141),
                       labels=['All'],
                       ylog=True,
                       x_label='Total Energy Deposited in Ecal Si (MeV)')
        
        plt.create_root_hist('total_ecal_energy', 
                             self.ntuple['total_ecal_energy'], 
                             140, 0, 140,
                             'Energy Deposited in Upstream Trigger Pad (MeV)',
                             color=r.kRed+2)
        
        plt.plot_hists([
                        self.ntuple['ecal_hit_energy'],
                       ], 
                       np.linspace(0, 100, 201),
                       labels=['All'],
                       ylog=True,
                       x_label='Readout Hit Energy Deposited in Ecal Si (MeV)')
        
        plt.create_root_hist('ecal_hit_energy', 
                             self.ntuple['ecal_hit_energy'], 
                             200, 0, 100,
                             'Readout Hit Energy Deposited in Ecal Si (MeV)',
                             color=r.kRed+2)
       
        plt.plot_hists([
                        self.ntuple['hcal_hit_energy'],
                       ], 
                       np.linspace(0, 150, 301),
                       labels=['All'],
                       ylog=True,
                       x_label='Readout Hit Energy Deposited in Hcal Scint (MeV)')
        
        plt.create_root_hist('hcal_hit_energy', 
                             self.ntuple['hcal_hit_energy'], 
                             300, 0, 150,
                             'Readout Hit Energy Deposited in Ecal Si (MeV)',
                             color=r.kRed+2)
        
        plt.plot_hists([
                        self.ntuple['total_hcal_energy'],
                       ],
                       np.linspace(0, 400, 401),
                       labels=['All'],
                       ylog=True,
                       x_label='Total Energy Deposited in Hcal Scint (MeV)')

        plt.create_root_hist('total_hcal_energy', 
                             self.ntuple['total_hcal_energy'], 
                             400, 0, 400,
                             'Total Energy Deposited in Hcal Scint (MeV)',
                             color=r.kRed+2)

        plt.plot_hists([
                        self.ntuple['total_recoil_hits'],
                       ],
                       np.linspace(0, 150, 151),
                       labels=['All'],
                       ylog=True,
                       x_label='Recoil Hit Multiplicity')

        plt.create_root_hist('total_recoil_hits', 
                             self.ntuple['total_recoil_hits'], 
                             150, 0, 150,
                             'Recoil Hit Multiplicity',
                             color=r.kRed+2)

        for layer_n in xrange(0, 10): 
            plt.plot_hists([
                self.ntuple['total_recoil_hits_l%s' % (layer_n + 1)],
                ],
                np.linspace(0, 150, 151),
                labels=['All'],
                ylog=True,
                x_label='Recoil Hit Layer %s Multiplicity' % (layer_n + 1))

            plt.plot_hists([
                self.ntuple['total_charge_l%s' % (layer_n + 1)],
                ],
                np.linspace(0, 50, 151),
                labels=['All'],
                ylog=True,
                x_label='Total Charge Layer %s' % (layer_n + 1))

            plt.create_root_hist('total_recoil_hits_l%s' % (layer_n + 1), 
                                 self.ntuple['total_recoil_hits_l%s' % (layer_n + 1)],
                                 150, 0, 150,
                                 'Recoil Hit Layer %s Multiplicity' % (layer_n + 1), 
                                 color=r.kRed+2)
            plt.create_root_hist('total_charge_l%s' % (layer_n + 1), 
                                 self.ntuple['total_charge_l%s' % (layer_n + 1)],
                                 150, 0, 50,
                                 'Total Charge Layer %s' % (layer_n + 1),
                                 color=r.kRed+2)

        plt.close()
