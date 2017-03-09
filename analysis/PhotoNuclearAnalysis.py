
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
from numpy import linalg as la
from scipy.stats import norm

class PhotoNuclearAnalysis(object) : 

    def __init__(self) : 

        self.initialize()

        print Plotter.__file__

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
        #print "Vertex: %s" % particle.getVertex()[2]
        if abs(particle.getVertex()[2]) < 0.550 : return True 
        return False
    
    def initialize(self) :

        self.events              = []

        self.pn_gamma_energy     = []
        self.pn_particle_mult    = []
        self.pn_interaction_z    = []
            
        self.n_tracks = []
        self.n_snippets = []

        self.recoil_e_truth_p  = []
        self.recoil_e_truth_pt = []
        self.recoil_e_truth_px = []
        self.recoil_e_truth_py = []
        self.recoil_e_truth_pz = []

        self.lead_neutron_ke     = []
        self.lead_neutron_pdgid  = []
        self.lead_neutron_theta  = []
        self.lead_neutron_phi    = []
        self.lead_em_ke     = []
        self.lead_em_pdgid  = []
        self.lead_em_theta  = []
        self.lead_em_phi    = []
        
        self.n_charged = []

        self.ecal_hit_energy = []
        self.total_energy_ecal = []
        self.total_energy_hcal = []
        self.hcal_hit_energy = []
        self.hcal_z_pos = []

        self.up_tp_energy = []
        self.down_tp_energy = []
       
        self.target_total_energy = []
        
        self.event_count = -1

    def process(self, event) :
        
        self.event_count += 1

        # Get the collection of MC particles from the even
        particles = event.get_collection('SimParticles')
        
        recoil_e = None
        pn_gamma = None

        for particle in particles:

            #print 'PDG ID: %s' % particle.getPdgID()
            #print 'Parent count: %s' % particle.getParentCount()
            #print 'Daughtercount: %s' % partiyeah, it could be a local icle.getDaughterCount()
            #print 'Vertex: %s' % particle.getVertex()[2]
            
            if self.is_recoil(particle):
                recoil_e = particle
                break
        
        # Get the FindableTracks collection from the event
        findable_tracks = event.get_collection('FindableTracks')
        findable_dic = {}
        snippet_dic = {}

        # Create a map between a sim particle and a findable track. This particular
        # map will only contain tracks that can be found using the 4s and 3s1a
        # strategies.
        for findable_track in findable_tracks:
            if ((findable_track.is4sFindable()) 
                    or (findable_track.is3s1aFindable())):
                findable_dic[findable_track.getSimParticle()] = findable_track
            elif findable_track.is2sFindable() and findable_track.getSimParticle() != recoil_e:
                snippet_dic[findable_track.getSimParticle()] = findable_track

        # Use the recoil electron to find the PN gamma
        #print "Daughter count: %s" % recoil_e.getDaughterCount()
        for daughter_count in xrange(0, recoil_e.getDaughterCount()):
            daughter = recoil_e.getDaughter(daughter_count)
            #print "Event: %s, Energy %s, Count %s" % (event.get_event_number(), daughter.getEnergy(), daughter.getDaughterCount())
            if ((daughter.getPdgID() == 22) 
                & self.created_within_target(daughter)
                #& (daughter.getEndPoint()[2] < .55)):
                & (daughter.getEnergy() >= 2500)):
                #print "Event: %s, Energy %s, Count %s" % (event.get_event_number(), daughter.getEnergy(), daughter.getDaughterCount())
                pn_gamma = daughter
        
        if pn_gamma is None: 
            #print "Skipping event"
            #for daughter_count in xrange(0, recoil_e.getDaughterCount()):
            #    daughter = recoil_e.getDaughter(daughter_count)
            #    print "Event: %s, Energy %s, Count %s, Vertex %s" % (event.get_event_number(), daughter.getEnergy(), daughter.getDaughterCount(), daughter.getVertex()[2] )
            return

        self.n_tracks.append(len(findable_dic))
        self.n_snippets.append(len(snippet_dic))

        self.pn_gamma_energy.append(pn_gamma.getEnergy())
        self.pn_particle_mult.append(pn_gamma.getDaughterCount())
        self.pn_interaction_z.append(pn_gamma.getEndPoint()[2])

        # Calculate some kinematic variables for the recoil electron
        recoil_e_pvec = recoil_e.getMomentum()
        self.recoil_e_truth_p.append(la.norm(recoil_e_pvec))
        self.recoil_e_truth_pt.append(
                math.sqrt(recoil_e_pvec[0]*recoil_e_pvec[0] + recoil_e_pvec[1]*recoil_e_pvec[1]))
        self.recoil_e_truth_px.append(recoil_e_pvec[0]) 
        self.recoil_e_truth_py.append(recoil_e_pvec[1]) 
        self.recoil_e_truth_pz.append(recoil_e_pvec[2]) 

        trigger_pad_hits = event.get_collection('TriggerPadSimHits')
        total_up_tp_energy = 0
        total_down_tp_energy = 0
        for trigger_pad_hit in trigger_pad_hits:
            if (trigger_pad_hit.getPosition()[2] > 0): 
                total_down_tp_energy += trigger_pad_hit.getEdep()
            if (trigger_pad_hit.getPosition()[2] < 0) & (trigger_pad_hit.getContrib(0).pdgCode != 11): 
                total_up_tp_energy += trigger_pad_hit.getEdep()

        self.down_tp_energy.append(total_down_tp_energy)
        self.up_tp_energy.append(total_up_tp_energy)

        ecal_hits = event.get_collection('EcalHits')
        energy_counter = 0
        for ecal_hit in ecal_hits: 
            self.ecal_hit_energy.append(ecal_hit.getEnergy())
            energy_counter += ecal_hit.getEnergy()

        self.total_energy_ecal.append(energy_counter)

        hcal_hits = event.get_collection('HcalHits')
        energy_counter = 0
        for hcal_hit in hcal_hits: 
            self.hcal_hit_energy.append(hcal_hit.getEnergy())
            energy_counter += hcal_hit.getEnergy()

        self.total_energy_hcal.append(energy_counter)

    def finalize(self) :

        up_tp_cuts, up_tp_eff = self.calc_sig_eff(np.array(self.down_tp_energy), 0, 10, 0.1, False)
        #print 'TP cuts: %s' % up_tp_cuts
        #print 'TP up eff: %s' % up_tp_eff
        down_tp_cuts, down_tp_eff = self.calc_sig_eff(np.array(self.up_tp_energy), 0, 10, 0.1, False)
        #print 'TP down eff: %s' % down_tp_eff

        self.pn_gamma_energy = np.array(self.pn_gamma_energy)
        self.pn_particle_mult = np.array(self.pn_particle_mult)
        self.n_tracks = np.array(self.n_tracks)
        self.down_tp_energy = np.array(self.down_tp_energy)
        self.up_tp_energy = np.array(self.up_tp_energy)
        self.n_snippets = np.array(self.n_snippets)
        self.ecal_hit_energy = np.array(self.ecal_hit_energy)
        self.total_energy_ecal = np.array(self.total_energy_ecal)
        self.total_energy_hcal = np.array(self.total_energy_hcal)
        self.hcal_hit_energy = np.array(self.hcal_hit_energy)

        plt = Plotter.Plotter("photo_nuclear_analysis")

        veto = ((self.down_tp_energy < 0.75) & (self.n_tracks == 1) & (self.up_tp_energy < 0.5) 
                & (self.n_snippets == 0))

        plt.plot_hists([self.pn_gamma_energy, 
                        self.pn_gamma_energy[self.down_tp_energy < 0.75],
                        self.pn_gamma_energy[self.up_tp_energy < 0.5],
                        self.pn_gamma_energy[self.n_tracks == 1], 
                        self.pn_gamma_energy[self.n_snippets == 0], 
                        self.pn_gamma_energy[veto] 
                        ],
                      np.linspace(0, 4000, 41),
                      labels=['All', 
                              'Downstream TP Energy < 0.75 MeVt',
                              'Upstream TP Energy < 0.5 MeV',
                              'N tracks == 1',
                              'N 2s tracks == 0',
                              'Upstream Veto'
                              ],
                      ylog=True,
                      x_label='$E(\gamma)$ (MeV)')
                  

        plt.plot_hists([self.pn_particle_mult,
                        self.pn_particle_mult[self.down_tp_energy < 0.75],
                        self.pn_particle_mult[self.up_tp_energy < 0.5],
                        self.pn_particle_mult[self.n_tracks == 1], 
                        self.pn_particle_mult[self.n_snippets == 0], 
                        self.pn_particle_mult[veto] 
                        ],
                      np.linspace(0, 100, 101),
                      labels=['All', 
                              'Downstream TP Energy < 0.75 MeV',
                              'Upstream TP Energy < 0.5 MeV',
                              'N tracks == 1',
                              'N 2s tracks == 0',
                              'Upstream Veto'
                              ],
                      ylog=True,
                      x_label='PN Multiplicity')

        plt.plot_hists([self.n_tracks, 
                        self.n_tracks[self.down_tp_energy < 0.75],
                        self.n_tracks[self.up_tp_energy < 0.5],
                        self.n_tracks[self.n_tracks == 1], 
                        self.n_tracks[self.n_snippets == 0], 
                        self.n_tracks[veto] 
                        ],
                      np.linspace(0, 10, 11),
                      labels=['All', 
                              'Downstream TP Energy < 0.75 MeV',
                              'Upstream TP Energy < 0.5 MeV',
                              'N tracks == 1',
                              'N 2s tracks == 0',
                              'Upstream Veto'
                              ],
                      ylog=True,
                      x_label='Track Multiplicity')

        plt.plot_hists([self.n_snippets, 
                       self.n_snippets[self.down_tp_energy < 0.75], 
                       self.n_snippets[self.up_tp_energy < 0.5],
                       self.n_snippets[self.n_tracks == 1], 
                       self.n_snippets[self.n_snippets == 0], 
                       self.n_snippets[veto] 
                       ],
                      np.linspace(0, 10, 11),
                      labels=['All', 
                              'Downstream TP Energy < 0.75 MeV',
                              'Upstream TP Energy < 0.5 MeV',
                              'N tracks == 1',
                              'N 2s tracks == 0',
                              'Upstream Veto'
                              ],
                      ylog=True,
                      x_label='Snippet Multiplicity')
        
        
        plt.plot_hists([self.down_tp_energy, 
                       self.down_tp_energy[self.down_tp_energy < 0.75], 
                       self.down_tp_energy[self.up_tp_energy < 0.5],
                       self.down_tp_energy[self.n_tracks == 1], 
                       self.down_tp_energy[self.n_snippets == 0], 
                       self.down_tp_energy[veto] 
                       ],
                      np.linspace(0, 100, 101),
                      labels=['All', 
                              'Downstream TP Energy < 0.75 MeV',
                              'Upstream TP Energy < 0.5 MeV',
                              'N tracks == 1',
                              'N 2s tracks == 0',
                              'Upstream Veto'
                              ],
                      ylog=True, 
                      x_label='Energy Deposited in Downstream Trigger Pad (MeV)')

        plt.plot_graph(down_tp_cuts, down_tp_eff, 0, 0,
                       y_label='Cut Efficiency (%)', 
                       x_label='Energy Deposited in Downstream Trigger Pad (MeV)')

        plt.plot_hists([self.up_tp_energy, 
                       self.up_tp_energy[self.down_tp_energy < 0.75], 
                       self.up_tp_energy[self.up_tp_energy < 0.5],
                       self.up_tp_energy[self.n_tracks == 1], 
                       self.up_tp_energy[self.n_snippets == 0], 
                       self.up_tp_energy[veto] 
                       ],
                      np.linspace(0, 100, 101),
                      labels=['All', 
                              'Downstream TP Energy < 0.75 MeV',
                              'Upstream TP Energy < 0.5 MeV',
                              'N tracks == 1',
                              'N 2s tracks == 0',
                              'Upstream Veto'
                              ],
                      ylog=True, 
                      x_label='Energy Deposited in Upstream Trigger Pad (MeV)')

        plt.plot_graph(up_tp_cuts, up_tp_eff, 0, 0,
                       y_label='Cut Efficiency (%)', 
                       x_label='Energy Deposited in Upstream Trigger Pad (MeV)')
        
        plt.plot_hist(self.pn_interaction_z,
                      np.linspace(0, 350, 101),
                      ylog=True,
                      norm=True,
                      x_label='PN $\gamma$ Interaction Point z (mm)')
       
        '''
        plt.plot_hists([self.ecal_hit_energy, 
                       self.ecal_hit_energy[self.down_tp_energy < 0.75], 
                       self.ecal_hit_energy[self.up_tp_energy < 0.5],
                       self.ecal_hit_energy[self.n_tracks == 1], 
                       self.ecal_hit_energy[self.n_snippets == 0], 
                       self.ecal_hit_energy[veto] 
                       ],
                      np.linspace(0, 50, 51),
                      ylog=True,
                      x_label='Readout Hit Energy Deposited in Ecal Si (MeV)')
        '''

        plt.plot_hists([self.total_energy_ecal, 
                       self.total_energy_ecal[self.down_tp_energy < 0.75], 
                       self.total_energy_ecal[self.up_tp_energy < 0.5],
                       self.total_energy_ecal[self.n_tracks == 1], 
                       self.total_energy_ecal[self.n_snippets == 0], 
                       self.total_energy_ecal[veto] 
                       ],
                      np.linspace(0, 140, 35),
                      labels=['All', 
                              'Downstream TP Energy < 0.75 MeV',
                              'Upstream TP Energy < 0.5 MeV',
                              'N tracks == 1',
                              'N 2s tracks == 0',
                              'Upstream Veto'
                              ],
                      ylog=True,
                      x_label='Total Energy Deposited in Ecal Si (MeV)')


        '''
        plt.plot_hists([self.hcal_hit_energy, 
                       self.hcal_hit_energy[self.down_tp_energy < 0.75], 
                       self.hcal_hit_energy[self.up_tp_energy < 0.5],
                       self.hcal_hit_energy[self.n_tracks == 1], 
                       self.hcal_hit_energy[self.n_snippets == 0], 
                       self.hcal_hit_energy[veto] 
                       ],
                      np.linspace(0, 50, 51),
                      ylog=True,
                      x_label='Readout Hit Energy Deposited in Hcal Scint (MeV)')
        '''

        plt.plot_hists([self.total_energy_hcal, 
                       self.total_energy_hcal[self.down_tp_energy < 0.75], 
                       self.total_energy_hcal[self.up_tp_energy < 0.5],
                       self.total_energy_hcal[self.n_tracks == 1], 
                       self.total_energy_hcal[self.n_snippets == 0], 
                       self.total_energy_hcal[veto]
                       ],
                      np.linspace(0, 140, 35),
                      labels=['All', 
                              'Downstream TP Energy < 0.75 MeV',
                              'Upstream TP Energy < 0.5 MeV',
                              'N tracks == 1',
                              'N 2s tracks == 0',
                              'Upstream Veto'
                              ],
                      ylog=True,
                      x_label='Total Energy Deposited in Hcal Scint (MeV)')


        plt.close()
        
        print "Total events processed: %s" % len(self.n_tracks)
        
        pass_tp_veto = len(self.down_tp_energy[self.down_tp_energy < 0.75])
        pass_tp_veto_per = pass_tp_veto/len(self.down_tp_energy)
        print 'Pass down TP veto: %s (%s)' % (pass_tp_veto, pass_tp_veto_per*100)

        pass_tp_up_veto = len(self.up_tp_energy[self.up_tp_energy < 0.5])
        pass_tp_up_veto_per = pass_tp_up_veto/len(self.up_tp_energy)
        print 'Pass up TP veto: %s (%s)' % (pass_tp_up_veto, pass_tp_up_veto_per*100)

        pass_track_veto = len(self.n_tracks[self.n_tracks == 1])
        pass_track_veto_per = pass_track_veto/len(self.n_tracks)
        print 'Events that pass track veto: %s (%s)' % (pass_track_veto, pass_track_veto_per*100)

        pass_snippet_veto = len(self.n_snippets[self.n_snippets == 0])
        pass_snippet_veto_per = pass_snippet_veto/len(self.n_snippets)
        print 'Events passing snippet veto: %s (%s)' % (pass_snippet_veto, pass_track_veto_per*100)

        veto = ((self.down_tp_energy < 0.75) & (self.n_tracks == 1) & (self.up_tp_energy < 0.5) 
                & (self.n_snippets == 0))

        print "Event passing the veto: %s" % len(self.n_tracks[veto])
