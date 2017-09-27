from __future__ import division

import numpy as np
import math
import ROOT as r
import Plotter

from numpy import linalg as la

class PnReWeightingAnalysis:

    def __init__(self):
        self.initialize()

    def initialize(self): 
        
        self.ntuple = {}
        self.variables = [
            'events', 
            'pn_mult',
            'nucleon_ke', 'nucleon_theta', 'nucleon_w',
            'nucleon_ke_w', 'nucleon_theta_w', 'nucleon_w_w',
            'hnucleon_ke', 'hnucleon_theta', 'hnucleon_w',
            'hnucleon_ke_w', 'hnucleon_theta_w', 'hnucleon_w_w',
            'proton_ke', 'proton_theta', 'proton_w',
            'hproton_ke', 'hproton_theta', 'hproton_w',
            'neutron_ke', 'neutron_theta', 'neutron_w',
            'hneutron_ke', 'hneutron_theta', 'hneutron_w',
            'hw_nucleons', 'hw_nucleons_all', 'hw_nucleons_sum', 
            'hw_nucleons_theta', 'hw_nucleons_ke', 
        ]

        self.delta_variables = []
        delta = .1
        for idelta in xrange(25):
            self.delta_variables.extend(
                    ['nucleon_ke_delta_%s' % delta, 
                     'nucleon_theta_delta_%s' % delta, 
                     'nucleon_w_delta_%s' % delta])
            delta += .1

        self.theta_variables = { 
            'chnucleon_ke_theta':{}, 'chproton_ke_theta':{}, 'chneutron_ke_theta':{},  
            'chnucleon_w_theta':{}, 'chproton_w_theta':{}, 'chneutron_w_theta':{},  
            'chw_nucleon_theta':{}, 'chw_proton_theta':{}, 'chw_neutron_theta':{}
        }

        theta = 20
        for itheta in xrange(0, 8):
            for ptype in ['nucleon', 'proton', 'neutron']:
                self.variables.append('%s_ke_theta_%s_%s' % (ptype, theta, theta + 20))
                self.variables.append('%s_w_theta_%s_%s' % (ptype, theta, theta + 20))
                self.variables.append('h%s_ke_theta_%s_%s' % (ptype, theta, theta + 20))
                self.variables.append('h%s_w_theta_%s_%s' % (ptype, theta, theta + 20))
                self.variables.append('hw_%s_theta_%s_%s' % (ptype, theta, theta + 20))

            for key in self.theta_variables: 
                self.theta_variables[key]['theta_%s_%s' % (theta, theta + 20)] = -9999

            theta += 20

        for variable in self.variables: 
            self.ntuple[variable] = []
        
        for variable in self.delta_variables: 
            self.ntuple[variable] = []

        self.colors = [r.kAzure + 2, r.kGreen - 2, r.kRed + 2, r.kOrange + 8,
                       r.kMagenta - 4, r.kAzure + 10, r.kYellow, r.kBlack, r.kRed]
        
        self.event_count = 0
        self.file_prefix = None

    def reset_variables(self): 
        for var_key in self.theta_variables:
            for key in self.theta_variables[var_key]:
                self.theta_variables[var_key][key] = -9999

    def process(self, event):

        self.event_count += 1

        if self.event_count == 1: 
            self.file_prefix = event.get_file_name()[
                    event.get_file_name().rfind('/') + 1:-5]
            print self.file_prefix

        weights = event.get_collection('PNweight_recon')

        event_weight = 1.0
        for weight in weights:
            event_weight = weight.getWeight() 

        sim_particles = event.get_collection('SimParticles_sim')

        recoil_e = None
        for sim_particle in sim_particles:
            if ((sim_particle.getPdgID() == 11) 
                    & (sim_particle.getParentCount() == 0)):
                recoil_e = sim_particle
                #print '[ PnReWeightAnalysis ]: Recoil e- found.'
                break

        if not recoil_e: sys.exit('[ PnReWeightingAnalysis ]: Failed to find recoil e-.')

        pn_gamma = None
        for idaughter in xrange(recoil_e.getDaughterCount()): 
            daughter = recoil_e.getDaughter(idaughter)
            if daughter.getDaughterCount() == 0: continue
            if daughter.getDaughter(0).getProcessType() == 9: 
                pn_gamma = daughter
                break

        if not pn_gamma: sys.exit('[ PnReWeightingAnalysis ]: Failed to find PN gamma.')
        


        #print '[ PnReWeightAnalysis ]: Photo-nuclear multiplicity: %s' % pn_gamma.getDaughterCount()
        self.ntuple['pn_mult'].append(pn_gamma.getDaughterCount())
       
        self.reset_variables()        

        hnucleon_ke = -9999
        hnucleon_theta = -9999
        hw = -9999
        hnucleon = None

        hproton_ke = -9999
        hproton_theta = -9999
        hw_proton = -9999
        hproton = None

        hneutron_ke = -9999
        hneutron_theta = -9999
        hw_neutron = -9999
        hneutron = None

        for idaughter in xrange(pn_gamma.getDaughterCount()): 
            
            daughter = pn_gamma.getDaughter(idaughter)
            ke = daughter.getEnergy() - daughter.getMass()
            pvec = daughter.getMomentum()
            p = la.norm(pvec)
            theta = math.acos(pvec[2]/p)*180/3.14159
            pdg_id = abs(daughter.getPdgID())
            w = self.calculate_w(daughter, 0.5)
        
            if ((pdg_id == 2212) or (pdg_id == 2112)): 
                
                self.ntuple['nucleon_ke'].append(ke)
                self.ntuple['nucleon_theta'].append(theta)
                self.ntuple['nucleon_w'].append(w)
                self.ntuple['nucleon_ke_w'].append(ke*event_weight)
                self.ntuple['nucleon_theta_w'].append(theta*event_weight)
                self.ntuple['nucleon_w_w'].append(w*event_weight)
              
                
                delta = .1
                for idelta in xrange(25):
                    self.ntuple['nucleon_w_delta_%s' % delta].append(
                        self.calculate_w(daughter, delta))

                    delta += .1

                if ke > hnucleon_ke: 
                    hnucleon_ke = ke
                    hnucleon = daughter
                    hnucleon_theta = theta
                
                if pdg_id ==2212: 
                
                    self.ntuple['proton_ke'].append(ke)
                    self.ntuple['proton_theta'].append(theta)
                    self.ntuple['proton_w'].append(w)
                    
                    if ke > hproton_ke:
                        hproton_ke = ke
                        hproton = daughter
                        hproton_theta = theta
                    
                    if w > hw_proton: 
                        hw_proton = w

                if pdg_id ==2112: 
                    
                    self.ntuple['neutron_ke'].append(ke)
                    self.ntuple['neutron_theta'].append(theta)
                    self.ntuple['neutron_w'].append(w)
                    
                    if ke > hneutron_ke:
                        hneutron_ke = ke
                        hneutron = daughter
                        hneutron_theta = theta

                    if w > hw_neutron:
                        hw_neutron = w

                if w > hw: 
                    hw = w

               
                if (theta >= 20) & (theta < 40):
                    self.fill_angular_plots(pdg_id, ke, w, theta, 20, 40)
                elif (theta >= 40) & (theta < 60):
                    self.fill_angular_plots(pdg_id, ke, w, theta, 40, 60)
                elif (theta >= 60) & (theta < 80): 
                    self.fill_angular_plots(pdg_id, ke, w, theta, 60, 80)
                elif (theta >= 80) & (theta < 100): 
                    self.fill_angular_plots(pdg_id, ke, w, theta, 80, 100)
                elif (theta >= 100) & (theta < 120): 
                    self.fill_angular_plots(pdg_id, ke, w, theta, 100, 120)
                elif (theta >= 120) & (theta < 140): 
                    self.fill_angular_plots(pdg_id, ke, w, theta, 120, 140)
                elif (theta >= 140) & (theta < 160): 
                    self.fill_angular_plots(pdg_id, ke, w, theta, 140, 160)
                elif (theta >= 160) & (theta < 180): 
                    self.fill_angular_plots(pdg_id, ke, w, theta, 160, 180)
        
        self.ntuple['hnucleon_ke'].append(hnucleon_ke)
        self.ntuple['hnucleon_theta'].append(hnucleon_theta)
        if hnucleon: self.ntuple['hnucleon_w'].append(self.calculate_w(hnucleon, 0.5))
        else: self.ntuple['hnucleon_w'].append(-9999)
        self.ntuple['hproton_ke'].append(hproton_ke)
        self.ntuple['hproton_theta'].append(hproton_theta)
        if hproton: self.ntuple['hproton_w'].append(self.calculate_w(hproton, 0.5))
        else: self.ntuple['hproton_w'].append(-9999)
        self.ntuple['hneutron_ke'].append(hneutron_ke)
        self.ntuple['hneutron_theta'].append(hneutron_theta)
        if hneutron: self.ntuple['hneutron_w'].append(self.calculate_w(hneutron, 0.5))
        else: self.ntuple['hneutron_w'].append(-9999)

        self.ntuple['hw_nucleons'].append(hw)
        self.ntuple['hw_nucleons_all'].append(hw_proton)
        self.ntuple['hw_nucleons_all'].append(hw_neutron)
        if hw_proton == -9999: hw_proton = 0
        if hw_neutron == -9999: hw_neutron = 0
        self.ntuple['hw_nucleons_sum'].append(hw_proton + hw_neutron)

        theta = 20
        for itheta in xrange(0, 8): 
            for ptype in ['nucleon', 'proton', 'neutron']:
                self.ntuple['h%s_ke_theta_%s_%s' % (ptype, theta, theta + 20)].append(
                        self.theta_variables['ch%s_ke_theta' % ptype]['theta_%s_%s' % (theta, theta + 20)])
                self.ntuple['h%s_w_theta_%s_%s' % (ptype, theta, theta + 20)].append(
                        self.theta_variables['ch%s_w_theta' % ptype]['theta_%s_%s' % (theta, theta + 20)])
                self.ntuple['hw_%s_theta_%s_%s' % (ptype, theta, theta + 20)].append(
                        self.theta_variables['chw_%s_theta' % ptype]['theta_%s_%s' % (theta, theta + 20)])

            theta += 20

    def finalize(self):

        print 'Total: %s' % len(self.ntuple['hnucleon_ke'])


        for variable in self.variables: 
            self.ntuple[variable] = np.array(self.ntuple[variable])
        
        plt = Plotter.Plotter(self.file_prefix + '_plots')

        plt.plot_hist(self.ntuple['pn_mult'], 
                      np.linspace(0, 200, 201),
                      norm=True,
                      ylog=True,
                      x_label='Photonuclear Multiplicity')
        
        plt.create_root_hist('pn_mult', 
                             self.ntuple['pn_mult'], 
                             200, 0, 200,
                             'Photonuclear Multiplicity', 
                             color=self.colors[0])
        
        plt.plot_hists([
            self.ntuple['nucleon_ke'], 
            self.ntuple['proton_ke'], 
            self.ntuple['neutron_ke']],
            np.linspace(0, 5000, 501), 
            labels=['Nucleons', 'Protons', 'Neutrons'], 
            ylog=True, 
            x_label='Kinetic Energy, Inclusive (MeV)')

        for sparticle in ['nucleon', 'proton', 'neutron']:
            plt.create_root_hist('%s_ke' % sparticle, 
                                 self.ntuple['%s_ke' % sparticle], 
                                 500, 0, 5000,
                                 'Kinetic Energy, Inclusive (mev)', 
                                 color=self.colors[0])

        '''
        plt.plot_hists([
            self.ntuple['nucleon_ke_w']] 
            np.linspace(0, 5000, 501), 
            labels=['Nucleons', 'Protons', 'Neutrons'], 
            ylog=True, 
            x_label='Kinetic Energy, Inclusive, Weighted (MeV)')
        '''
        for sparticle in ['nucleon']:
            plt.create_root_hist('%s_ke_w' % sparticle, 
                                 self.ntuple['%s_ke_w' % sparticle], 
                                 500, 0, 5000,
                                 'Kinetic Energy, Inclusive, Weighted (mev)', 
                                 color=self.colors[0])


        theta_cut = 70
        while theta_cut != 180:
            plt.plot_hists([
                self.ntuple['nucleon_ke'][self.ntuple['nucleon_theta'] > theta_cut], 
                self.ntuple['proton_ke'][self.ntuple['proton_theta'] > theta_cut], 
                self.ntuple['neutron_ke'][self.ntuple['neutron_theta'] > theta_cut]],
                np.linspace(0, 5000, 501), 
                labels=['Nucleons, $\theta > %s$' % theta_cut,
                        'Protons $\theta > %s$' % theta_cut,
                        'Neutrons $\theta > %s$' % theta_cut], 
                ylog=True, 
                x_label='Kinetic Energy, Inclusive, $\theta > %s$ (MeV)' % theta_cut)

            for sparticle in ['nucleon', 'proton', 'neutron']:
                plt.create_root_hist('%s_ke_theta_cut_%s' % (sparticle, theta_cut),  
                                     self.ntuple['%s_ke' % sparticle][self.ntuple['%s_theta' % sparticle] > theta_cut], 
                                     500, 0, 5000,
                                     'Kinetic Energy, Inclusive, #theta > %s (MeV)' % theta_cut, 
                                     color=self.colors[0])
            theta_cut += 10

        plt.plot_hists([
            self.ntuple['hnucleon_ke'], 
            self.ntuple['hproton_ke'], 
            self.ntuple['hneutron_ke']],
            np.linspace(0, 5000, 501), 
            labels=['Hardest Nucleon', 'Hardest Proton', 'Hardest Neutron'], 
            ylog=True, 
            x_label='Kinetic Energy, Hardest (MeV)')

        for sparticle in ['hnucleon', 'hproton', 'hneutron']:
            plt.create_root_hist('%s_ke' % sparticle, 
                                 self.ntuple['%s_ke' % sparticle], 
                                 500, 0, 5000,
                                 'Kinetic Energy, Hardest (MeV)', 
                                 color=self.colors[0])

        plt.plot_hists([
            self.ntuple['hnucleon_ke'][self.ntuple['hnucleon_theta'] > 100], 
            self.ntuple['hproton_ke'][self.ntuple['hproton_theta'] > 100], 
            self.ntuple['hneutron_ke'][self.ntuple['hneutron_theta'] > 100]],
            np.linspace(0, 5000, 501), 
            labels=['Hardest Nucleon, $theta > 100$', 
                    'Hardest Proton, $theta > 100$',
                    'Hardest Neutron, $theta > 100$'], 
            ylog=True, 
            x_label='Kinetic Energy (MeV)')

        plt.create_root_hist('hnucleon_ke_theta_100', 
                             self.ntuple['hnucleon_ke'][self.ntuple['hnucleon_theta'] > 100], 
                             500, 0, 5000,
                             'Kinetic Energy (#theta > 100) (MeV)', 
                             color=self.colors[0])

        plt.create_root_hist('hproton_ke_theta_100', 
                             self.ntuple['hproton_ke'][self.ntuple['hnucleon_theta'] > 100], 
                             500, 0, 5000,
                             'Kinetic Energy (#theta > 100) (MeV)', 
                             color=self.colors[1])

        plt.create_root_hist('hneutron_ke_theta_100', 
                             self.ntuple['hneutron_ke'][self.ntuple['hnucleon_theta'] > 100], 
                             500, 0, 5000,
                             'Kinetic Energy (#theta > 100) (MeV)', 
                             color=self.colors[2])

        plt.plot_hists([
            self.ntuple['nucleon_w'], 
            self.ntuple['proton_w'], 
            self.ntuple['neutron_w']],
            np.linspace(0, 5000, 501), 
            labels=['Nucleons', 'Protons', 'Neutrons'], 
            ylog=True, 
            x_label='W, Inclusive (MeV)')

        for sparticle in ['nucleon', 'proton', 'neutron']:
            plt.create_root_hist('%s_w' % sparticle, 
                                 self.ntuple['%s_w' % sparticle], 
                                 500, 0, 5000,
                                 'W, inclusive (MeV)', 
                                 color=self.colors[0])
        

        '''
        plt.plot_hists([
            self.ntuple['nucleon_w_w']] 
            np.linspace(0, 5000, 501), 
            labels=['Nucleons', 'Protons', 'Neutrons'], 
            ylog=True, 
            x_label='W, Inclusive, Weighted (MeV)')
        '''

        for sparticle in ['nucleon']:
            plt.create_root_hist('%s_w_w' % sparticle, 
                                 self.ntuple['%s_w_w' % sparticle], 
                                 500, 0, 5000,
                                 'W, inclusive, Weighted (MeV)', 
                                 color=self.colors[0])

        theta_cut = 70
        while theta_cut != 180:
            plt.plot_hists([
                self.ntuple['nucleon_w'][self.ntuple['nucleon_theta'] > theta_cut], 
                self.ntuple['proton_w'][self.ntuple['proton_theta'] > theta_cut], 
                self.ntuple['neutron_w'][self.ntuple['neutron_theta'] > theta_cut]],
                np.linspace(0, 5000, 501), 
                labels=['Nucleons, $\theta > %s$' % theta_cut,
                        'Protons $\theta > %s$' % theta_cut,
                        'Neutrons $\theta > %s$' % theta_cut], 
                ylog=True, 
                x_label='W, Inclusive, $\theta > %s$ (MeV)' % theta_cut)

            for sparticle in ['nucleon', 'proton', 'neutron']:
                plt.create_root_hist('%s_w_theta_cut_%s' % (sparticle, theta_cut),  
                                     self.ntuple['%s_w' % sparticle][self.ntuple['%s_theta' % sparticle] > theta_cut], 
                                     500, 0, 5000,
                                     'W, Inclusive, #theta > %s (MeV)' % theta_cut, 
                                     color=self.colors[0])
           
            '''
            plt.plot_hists([
                self.ntuple['nucleon_w_w'][self.ntuple['nucleon_theta'] > theta_cut], 
                self.ntuple['proton_w_w'][self.ntuple['proton_theta'] > theta_cut], 
                self.ntuple['neutron_w_w'][self.ntuple['neutron_theta'] > theta_cut]],
                np.linspace(0, 5000, 501), 
                labels=['Nucleons, Weighted, $\theta > %s$' % theta_cut,
                        'Protons, Weighted, $\theta > %s$' % theta_cut,
                        'Neutrons, Weighted $\theta > %s$' % theta_cut], 
                ylog=True, 
                x_label='W, Inclusive, Weighted, $\theta > %s$ (MeV)' % theta_cut)
            '''

            for sparticle in ['nucleon']:
                plt.create_root_hist('%s_w_theta_cut_%s_w' % (sparticle, theta_cut),  
                                     self.ntuple['%s_w_w' % sparticle][self.ntuple['%s_theta' % sparticle] > theta_cut], 
                                     500, 0, 5000,
                                     'W, Inclusive, Weighted, #theta > %s (MeV)' % theta_cut, 
                                     color=self.colors[0])
            theta_cut += 10

        plt.plot_hists([
            self.ntuple['hnucleon_w'], 
            self.ntuple['hproton_w'], 
            self.ntuple['hneutron_w']],
            np.linspace(0, 5000, 501), 
            labels=['Hardest Nucleon', 'Hardest Proton', 'Hardest Neutron'], 
            ylog=True, 
            x_label='W (MeV)')

        plt.create_root_hist('hnucleon_w', 
                             self.ntuple['hnucleon_w'], 
                             500, 0, 5000,
                             'W (MeV)', 
                             color=self.colors[0])

        plt.create_root_hist('hproton_w', 
                             self.ntuple['hproton_w'], 
                             500, 0, 5000,
                             'W (MeV)', 
                             color=self.colors[1])

        plt.create_root_hist('hneutron_w', 
                             self.ntuple['hneutron_w'], 
                             500, 0, 5000,
                             'W (MeV)', 
                             color=self.colors[2])

        plt.plot_hists([
            self.ntuple['hnucleon_w'][self.ntuple['hnucleon_theta'] > 100], 
            self.ntuple['hproton_w'][self.ntuple['hproton_theta'] > 100], 
            self.ntuple['hneutron_w'][self.ntuple['hneutron_theta'] > 100]],
            np.linspace(0, 5000, 501), 
            labels=['Hardest Nucleon, $theta > 100$', 
                    'Hardest Proton, $theta > 100$',
                    'Hardest Neutron, $theta > 100$'], 
            ylog=True, 
            x_label='W (MeV)')

        plt.create_root_hist('hnucleon_w_theta_100', 
                             self.ntuple['hnucleon_w'][self.ntuple['hnucleon_theta'] > 100], 
                             500, 0, 5000,
                             'W (#theta > 100) (MeV)', 
                             color=self.colors[0])

        plt.create_root_hist('hproton_w_theta_100', 
                             self.ntuple['hproton_w'][self.ntuple['hnucleon_theta'] > 100], 
                             500, 0, 5000,
                             'W (#theta > 100) (MeV)', 
                             color=self.colors[1])

        plt.create_root_hist('hneutron_w_theta_100', 
                             self.ntuple['hneutron_w'][self.ntuple['hnucleon_theta'] > 100], 
                             500, 0, 5000,
                             'W (#theta > 100) (MeV)', 
                             color=self.colors[2])

        plt.plot_hist(self.ntuple['hw_nucleons'], 
                np.linspace(0, 5000, 501), 
                ylog=True,
                x_label='Highest W (MeV)')

        plt.create_root_hist('hw_nucleons', 
                             self.ntuple['hw_nucleons'], 
                             500, 0, 5000,
                             'Highest W (MeV)', 
                             color=self.colors[0])

        plt.plot_hist(self.ntuple['hw_nucleons_all'], 
                np.linspace(0, 5000, 501), 
                ylog=True,
                x_label='Highest W, All Nucleons (MeV)')

        plt.create_root_hist('hw_nucleons_all', 
                             self.ntuple['hw_nucleons_all'], 
                             500, 0, 5000,
                             'Highest W, All Nucleons (MeV)', 
                             color=self.colors[0])
        
        plt.plot_hist(self.ntuple['hw_nucleons_sum'], 
                np.linspace(0, 5000, 501), 
                ylog=True,
                x_label='Sum of Highest proton W + Highest neutron W (MeV)')

        plt.create_root_hist('hw_nucleons_sum', 
                             self.ntuple['hw_nucleons_sum'], 
                             500, 0, 5000,
                             'Sum of Highest W (MeV)', 
                             color=self.colors[0])

        theta = 20
        nucleon_ke_theta_array = []
        nucleon_w_theta_array = []
        labels = []
        for itheta in xrange(0, 8): 
            nucleon_ke_theta_array.append(self.ntuple['nucleon_ke_theta_%s_%s' % (theta, theta + 20)]) 
            nucleon_w_theta_array.append(self.ntuple['nucleon_w_theta_%s_%s' % (theta, theta + 20)]) 
            labels.append('$%s \leq \theta < %s$' % (theta, theta + 20))
            
            plt.create_root_hist('nucleon_ke_theta_%s_%s' % (theta, theta + 20), 
                             self.ntuple['nucleon_ke_theta_%s_%s' % (theta, theta + 20)], 
                             500, 0, 5000,
                             'Kinetic Energy, Inclusive, %s < #theta %s, (MeV)' % (theta, theta + 20), 
                             color=self.colors[itheta])

            plt.create_root_hist('nucleon_w_theta_%s_%s' % (theta, theta + 20), 
                             self.ntuple['nucleon_w_theta_%s_%s' % (theta, theta + 20)], 
                             500, 0, 5000,
                             'W, Inclusive, %s < #theta %s, (MeV)' % (theta, theta + 20), 
                             color=self.colors[itheta])
            theta += 20

        plt.plot_hists(nucleon_ke_theta_array, 
                      np.linspace(0, 3000, 301), 
                      ylog=True,
                      labels=labels,
                      norm=True,
                      x_label='Kinetic Energy, Inclusive (MeV)')

        plt.plot_hists(nucleon_w_theta_array, 
                      np.linspace(0, 3000, 301), 
                      ylog=True,
                      labels=labels, 
                      norm=True,
                      x_label='W, Inclusive (MeV)')

        theta = 20
        for itheta in xrange(0, 8):
           
            for pindex, ptype in enumerate(['nucleon', 'proton', 'neutron']):

                plt.plot_hist(self.ntuple['h%s_ke_theta_%s_%s' % (ptype, theta, theta + 20)], 
                    np.linspace(0, 5000, 501), 
                    ylog=True,
                    x_label='Kinetic Energy, Hardest (%s < $\theta$ < %s) (MeV)' % (theta, theta + 20))

                plt.create_root_hist('h%s_ke_theta_%s_%s' % (ptype, theta, theta + 20), 
                                 self.ntuple['h%s_ke_theta_%s_%s' % (ptype, theta, theta + 20)], 
                                 500, 0, 5000,
                                 'Kinetic Energy, Hardest (%s < #theta %s) (MeV)' % (theta, theta + 20), 
                                 color=self.colors[pindex])

                plt.plot_hist(self.ntuple['h%s_w_theta_%s_%s' % (ptype, theta, theta + 20)], 
                    np.linspace(0, 5000, 501), 
                    ylog=True,
                    x_label='W, Hardest (%s < $\theta$ < %s) (MeV)' % (theta, theta + 20))

                plt.create_root_hist('h%s_w_theta_%s_%s' % (ptype, theta, theta + 20), 
                                 self.ntuple['h%s_w_theta_%s_%s' % (ptype, theta, theta + 20)], 
                                 500, 0, 5000,
                                 'W, Hardest (%s < #theta %s) (MeV)' % (theta, theta + 20), 
                                 color=self.colors[pindex])

                plt.plot_hist(self.ntuple['hw_%s_theta_%s_%s' % (ptype, theta, theta + 20)], 
                    np.linspace(0, 5000, 501), 
                    ylog=True,
                    x_label='Highest W, %s < theta < %s (MeV)' % (theta, theta + 20))

                plt.create_root_hist('hw_%s_theta_%s_%s' % (ptype, theta, theta + 20), 
                                 self.ntuple['hw_%s_theta_%s_%s' % (ptype, theta, theta + 20)], 
                                 500, 0, 5000,
                                 'Highest W, %s < theta < %s (MeV)' % (theta, theta + 20), 
                                 color=self.colors[pindex])
            
            theta = theta + 20

        delta = .1
        for idelta in xrange(25):

            plt.plot_hist(self.ntuple['nucleon_w_delta_%s' % delta],
                          np.linspace(0, 5000, 501), 
                          ylog=True, 
                          x_label='W, Inclusive, delta = %s (MeV)' % delta)

            plt.create_root_hist('nucleon_w_delta_%s' % delta,
                                 self.ntuple['nucleon_w_delta_%s' % delta],
                                 500, 0, 5000,
                                 'W, Inclusive, delta = %s (MeV)' % delta, 
                                 color=self.colors[0])

            delta += .1

        plt.close()

    def calculate_w(self, particle, delta): 
        pvec = particle.getMomentum()
        p = la.norm(pvec)
        ke = particle.getEnergy() - particle.getMass()
        return 0.5*(p + ke)*(math.sqrt(1 + (delta*delta)) - delta*(pvec[2]/p))

    def fill_angular_plots(self, pdg_id, ke, w, theta, theta_min, theta_max): 
        
        self.ntuple['nucleon_ke_theta_%s_%s' % (theta_min, theta_max)].append(ke)
        self.ntuple['nucleon_w_theta_%s_%s' % (theta_min, theta_max)].append(w)

        if ke > self.theta_variables['chnucleon_ke_theta']['theta_%s_%s' % (theta_min, theta_max)]:
             self.theta_variables['chnucleon_ke_theta']['theta_%s_%s' % (theta_min, theta_max)] = ke
             self.theta_variables['chnucleon_w_theta']['theta_%s_%s' % (theta_min, theta_max)] = w

        if w > self.theta_variables['chw_nucleon_theta']['theta_%s_%s' % (theta_min, theta_max)]:
            self.theta_variables['chw_nucleon_theta']['theta_%s_%s' % (theta_min, theta_max)] = w

        if pdg_id == 2212: 
            if ke > self.theta_variables['chproton_ke_theta']['theta_%s_%s' % (theta_min, theta_max)]:
                self.theta_variables['chproton_ke_theta']['theta_%s_%s' % (theta_min, theta_max)] = ke
                self.theta_variables['chproton_w_theta']['theta_%s_%s' % (theta_min, theta_max)] = w

            if w > self.theta_variables['chw_proton_theta']['theta_%s_%s' % (theta_min, theta_max)]:
                self.theta_variables['chw_proton_theta']['theta_%s_%s' % (theta_min, theta_max)] = w
        
        if pdg_id == 2112: 
            if ke > self.theta_variables['chneutron_ke_theta']['theta_%s_%s' % (theta_min, theta_max)]:
                self.theta_variables['chneutron_ke_theta']['theta_%s_%s' % (theta_min, theta_max)] = ke
                self.theta_variables['chneutron_w_theta']['theta_%s_%s' % (theta_min, theta_max)] = w

            if w > self.theta_variables['chw_neutron_theta']['theta_%s_%s' % (theta_min, theta_max)]:
                self.theta_variables['chw_neutron_theta']['theta_%s_%s' % (theta_min, theta_max)] = w

