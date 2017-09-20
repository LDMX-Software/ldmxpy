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
            'hnucleon_ke', 'hnucleon_theta', 'hnucleon_w',
            'hproton_ke', 'hproton_theta', 'hproton_w',
            'hneutron_ke', 'hneutron_theta', 'hneutron_w',
            'hw_nucleons', 'hw_nucleons_all', 'hw_nucleons_sum', 
            'hw_nucleons_theta', 'hw_nucleons_ke', 
            'hnucleon_ke_theta_30_50',
            'hnucleon_ke_theta_50_70',
            'hnucleon_ke_theta_70_90',
            'hnucleon_ke_theta_90_110',
            'hnucleon_ke_theta_110_130',
            'hnucleon_ke_theta_130_150',
            'hnucleon_ke_theta_150',
            'hnucleon_w_theta_30_50',
            'hnucleon_w_theta_50_70',
            'hnucleon_w_theta_70_90',
            'hnucleon_w_theta_90_110',
            'hnucleon_w_theta_110_130',
            'hnucleon_w_theta_130_150',
            'hnucleon_w_theta_150',
            'hw_nucleon_theta_30_50',
            'hw_nucleon_theta_50_70',
            'hw_nucleon_theta_70_90',
            'hw_nucleon_theta_90_110',
            'hw_nucleon_theta_110_130',
            'hw_nucleon_theta_130_150',
            'hw_nucleon_theta_150',
        ]

        for variable in self.variables: 
            self.ntuple[variable] = []
        
        self.colors = ['#348ABD', '#A60628', '#7A68A6', '#467821', '#D55E00',
                       '#CC79A7', '#56B4E9', '#009E73', '#F0E442', '#0072B2']
        
        self.event_count = 0
        self.file_prefix = None
       
    def process(self, event):

        self.event_count += 1

        if self.event_count == 1: 
            self.file_prefix = event.get_file_name()[
                    event.get_file_name().rfind('/') + 1:-5]
            print self.file_prefix

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

        theta = 30
        hnucleon_ke_theta = {}
        hnucleon_particle_theta = {}
        hw_nucleon_theta = {}
        hproton_ke_theta = {}
        hw_proton_theta = {}
        hneutron_ke_theta = {}
        hw_neutron_theta = {}
        for itheta in xrange(0, 6): 
            hnucleon_ke_theta['theta_%s_%s'% (theta, theta + 20)] = -9999 
            hnucleon_particle_theta['theta_%s_%s'% (theta, theta + 20)] = None 
            hw_nucleon_theta['theta_%s_%s'% (theta, theta + 20)] = -9999  
            hproton_ke_theta['theta_%s_%s'% (theta, theta + 20)] = -9999  
            hw_proton_theta['theta_%s_%s'% (theta, theta + 20)] = -9999  
            hneutron_ke_theta['theta_%s_%s'% (theta, theta + 20)] = -9999  
            hw_neutron_theta['theta_%s_%s'% (theta, theta + 20)] = -9999  
            theta = theta + 20

        hnucleon_ke_theta['theta_150'] = -9999 
        hnucleon_particle_theta['theta_150'] = None
        hw_nucleon_theta['theta_150'] = -9999
        hproton_ke_theta['theta_150'] = -9999
        hw_proton_theta['theta_150'] = -9999
        hneutron_ke_theta['theta_150'] = -9999
        hw_neutron_theta['theta_150'] = -9999

        for idaughter in xrange(pn_gamma.getDaughterCount()): 
            daughter = pn_gamma.getDaughter(idaughter)
            ke = daughter.getEnergy() - daughter.getMass()
            pvec = daughter.getMomentum()
            p = la.norm(pvec)
            theta = math.acos(pvec[2]/p)*180/3.14159
            pdg_id = abs(daughter.getPdgID())
            w = self.calculate_w(daughter)
        
            if ((pdg_id == 2212) or (pdg_id == 2112)): 
                if ke > hnucleon_ke: 
                    hnucleon_ke = ke
                    hnucleon = daughter
                    hnucleon_theta = theta
                
                if pdg_id ==2212: 
                    if ke > hproton_ke:
                        hproton_ke = ke
                        hproton = daughter
                        hproton_theta = theta
                    
                    if w > hw_proton: 
                        hw_proton = w

                if pdg_id ==2112: 
                    if ke > hneutron_ke:
                        hneutron_ke = ke
                        hneutron = daughter
                        hneutron_theta = theta

                    if w > hw_neutron:
                        hw_neutron = w

                if w > hw: 
                    hw = w

                if (theta >= 30) & (theta < 50):

                    if ke > hnucleon_ke_theta['theta_30_50']: 
                        hnucleon_ke_theta['theta_30_50'] = ke
                        hnucleon_particle_theta['theta_30_50'] =  daughter
                    
                    if w > hw_nucleon_theta['theta_30_50']:
                        hw_nucleon_theta['theta_30_50'] = w

                    if pdg_id == 2212: 
                        if ke > hproton_ke_theta['theta_30_50']:
                            hproton_ke_theta['theta_30_50'] = ke

                        if w > hw_proton_theta['theta_30_50']:
                            hw_proton_theta['theta_30_50'] = w

                    if pdg_id == 2112: 
                        if ke > hneutron_ke_theta['theta_30_50']:
                            hneutron_ke_theta['theta_30_50'] = ke

                        if w > hw_neutron_theta['theta_30_50']:
                            hw_neutron_theta['theta_30_50'] = w

                elif (theta >= 50) & (theta < 70):
                    
                    if ke > hnucleon_ke_theta['theta_50_70']: 
                        hnucleon_ke_theta['theta_50_70'] = ke
                        hnucleon_particle_theta['theta_50_70'] =  daughter
                    
                    if w > hw_nucleon_theta['theta_50_70']:
                        hw_nucleon_theta['theta_50_70'] = w
                    
                    if pdg_id == 2212: 
                        if ke > hproton_ke_theta['theta_50_70']:
                            hproton_ke_theta['theta_50_70'] = ke

                        if w > hw_proton_theta['theta_50_70']:
                            hw_proton_theta['theta_50_70'] = w

                    if pdg_id == 2112: 
                        if ke > hneutron_ke_theta['theta_50_70']:
                            hneutron_ke_theta['theta_50_70'] = ke
                        
                        if w > hw_neutron_theta['theta_50_70']:
                            hw_neutron_theta['theta_50_70'] = w
                elif (theta >= 70) & (theta < 90):
                    
                    if ke > hnucleon_ke_theta['theta_70_90']: 
                        hnucleon_ke_theta['theta_70_90'] = ke
                        hnucleon_particle_theta['theta_70_90'] =  daughter
                    
                    if w > hw_nucleon_theta['theta_70_90']:
                        hw_nucleon_theta['theta_70_90'] = w
                    
                    if pdg_id == 2212: 
                        if ke > hproton_ke_theta['theta_70_90']:
                            hproton_ke_theta['theta_70_90'] = ke

                        if w > hw_proton_theta['theta_70_90']:
                            hw_proton_theta['theta_70_90'] = w
                    
                    if pdg_id == 2112: 
                        if ke > hneutron_ke_theta['theta_70_90']:
                            hneutron_ke_theta['theta_70_90'] = ke
                        
                        if w > hw_neutron_theta['theta_70_90']:
                            hw_neutron_theta['theta_70_90'] = w
                elif (theta >= 90) & (theta < 110):
                    
                    if ke > hnucleon_ke_theta['theta_90_110']: 
                        hnucleon_ke_theta['theta_90_110'] = ke
                        hnucleon_particle_theta['theta_90_110'] =  daughter
                    
                    if w > hw_nucleon_theta['theta_90_110']:
                        hw_nucleon_theta['theta_90_110'] = w
                    
                    if pdg_id == 2212: 
                        if ke > hproton_ke_theta['theta_90_110']:
                            hproton_ke_theta['theta_90_110'] = ke
                        
                        if w > hw_proton_theta['theta_90_110']:
                            hw_proton_theta['theta_90_110'] = w

                    if pdg_id == 2112: 
                        if ke > hneutron_ke_theta['theta_90_110']:
                            hneutron_ke_theta['theta_90_110'] = ke
                        
                        if w > hw_neutron_theta['theta_90_110']:
                            hw_neutron_theta['theta_90_110'] = w
                elif (theta >= 110) & (theta < 130):
                    
                    if ke > hnucleon_ke_theta['theta_110_130']: 
                        hnucleon_ke_theta['theta_110_130'] = ke
                        hnucleon_particle_theta['theta_110_130'] =  daughter
                    
                    if w > hw_nucleon_theta['theta_110_130']:
                        hw_nucleon_theta['theta_110_130'] = w
                    
                    if pdg_id == 2212: 
                        if ke > hproton_ke_theta['theta_110_130']:
                            hproton_ke_theta['theta_110_130'] = ke

                        if w > hw_proton_theta['theta_110_130']:
                            hw_proton_theta['theta_110_130'] = w

                    if pdg_id == 2112: 
                        if ke > hneutron_ke_theta['theta_110_130']:
                            hneutron_ke_theta['theta_110_130'] = ke
                        
                        if w > hw_neutron_theta['theta_110_130']:
                            hw_neutron_theta['theta_110_130'] = w
                elif (theta >= 130) & (theta < 150):
                    
                    if ke > hnucleon_ke_theta['theta_130_150']: 
                        hnucleon_ke_theta['theta_130_150'] = ke
                        hnucleon_particle_theta['theta_130_150'] =  daughter
                    
                    if w > hw_nucleon_theta['theta_130_150']:
                        hw_nucleon_theta['theta_130_150'] = w
                    
                    if pdg_id == 2212: 
                        if ke > hproton_ke_theta['theta_130_150']:
                            hproton_ke_theta['theta_130_150'] = ke

                        if w > hw_proton_theta['theta_130_150']:
                            hw_proton_theta['theta_130_150'] = w

                    if pdg_id == 2112: 
                        if ke > hneutron_ke_theta['theta_130_150']:
                            hneutron_ke_theta['theta_130_150'] = ke
                        
                        if w > hw_neutron_theta['theta_130_150']:
                            hw_neutron_theta['theta_130_150'] = w
                elif (theta >= 150):
                    
                    if ke > hnucleon_ke_theta['theta_150']: 
                        hnucleon_ke_theta['theta_150'] = ke
                        hnucleon_particle_theta['theta_150'] =  daughter
                    
                    if w > hw_nucleon_theta['theta_150']:
                        hw_nucleon_theta['theta_150'] = w

                    if pdg_id == 2212: 
                        if ke > hproton_ke_theta['theta_150']:
                            hproton_ke_theta['theta_150'] = ke

                        if w > hw_proton_theta['theta_150']:
                            hw_proton_theta['theta_150'] = w

                    if pdg_id == 2112: 
                        if ke > hneutron_ke_theta['theta_150']:
                            hneutron_ke_theta['theta_150'] = ke
                        
                        if w > hw_neutron_theta['theta_150']:
                            hw_neutron_theta['theta_150'] = w


        self.ntuple['hnucleon_ke'].append(hnucleon_ke)
        self.ntuple['hnucleon_theta'].append(hnucleon_theta)
        if hnucleon: self.ntuple['hnucleon_w'].append(self.calculate_w(hnucleon))
        else: self.ntuple['hnucleon_w'].append(-9999)
        self.ntuple['hproton_ke'].append(hproton_ke)
        self.ntuple['hproton_theta'].append(hproton_theta)
        if hproton: self.ntuple['hproton_w'].append(self.calculate_w(hproton))
        else: self.ntuple['hproton_w'].append(-9999)
        self.ntuple['hneutron_ke'].append(hneutron_ke)
        self.ntuple['hneutron_theta'].append(hneutron_theta)
        if hneutron: self.ntuple['hneutron_w'].append(self.calculate_w(hneutron))
        else: self.ntuple['hneutron_w'].append(-9999)

        self.ntuple['hw_nucleons'].append(hw)
        self.ntuple['hw_nucleons_all'].append(hw_proton)
        self.ntuple['hw_nucleons_all'].append(hw_neutron)
        if hw_proton == -9999: hw_proton = 0
        if hw_neutron == -9999: hw_neutron = 0
        self.ntuple['hw_nucleons_sum'].append(hw_proton + hw_neutron)

        self.ntuple['hnucleon_ke_theta_30_50'].append(hnucleon_ke_theta['theta_30_50'])
        self.ntuple['hnucleon_ke_theta_50_70'].append(hnucleon_ke_theta['theta_50_70'])
        self.ntuple['hnucleon_ke_theta_70_90'].append(hnucleon_ke_theta['theta_70_90'])
        self.ntuple['hnucleon_ke_theta_90_110'].append(hnucleon_ke_theta['theta_90_110'])
        self.ntuple['hnucleon_ke_theta_110_130'].append(hnucleon_ke_theta['theta_110_130'])
        self.ntuple['hnucleon_ke_theta_130_150'].append(hnucleon_ke_theta['theta_130_150'])
        self.ntuple['hnucleon_ke_theta_150'].append(hnucleon_ke_theta['theta_150'])

        if hnucleon_particle_theta['theta_30_50']: self.ntuple['hnucleon_w_theta_30_50'].append(self.calculate_w(hnucleon_particle_theta['theta_30_50']))
        else: self.ntuple['hnucleon_w_theta_30_50'].append(-9999)
        if hnucleon_particle_theta['theta_50_70']: self.ntuple['hnucleon_w_theta_50_70'].append(self.calculate_w(hnucleon_particle_theta['theta_50_70']))
        else: self.ntuple['hnucleon_w_theta_50_70'].append(-9999)
        if hnucleon_particle_theta['theta_70_90']: self.ntuple['hnucleon_w_theta_70_90'].append(self.calculate_w(hnucleon_particle_theta['theta_70_90']))
        else: self.ntuple['hnucleon_w_theta_70_90'].append(-9999)
        if hnucleon_particle_theta['theta_90_110']: self.ntuple['hnucleon_w_theta_90_110'].append(self.calculate_w(hnucleon_particle_theta['theta_90_110']))
        else: self.ntuple['hnucleon_w_theta_90_110'].append(-9999)
        if hnucleon_particle_theta['theta_110_130']: 
            self.ntuple['hnucleon_w_theta_110_130'].append(self.calculate_w(hnucleon_particle_theta['theta_110_130']))
        else: self.ntuple['hnucleon_w_theta_110_130'].append(-9999)
        if hnucleon_particle_theta['theta_130_150']: 
            self.ntuple['hnucleon_w_theta_130_150'].append(self.calculate_w(hnucleon_particle_theta['theta_130_150']))
        else: self.ntuple['hnucleon_w_theta_130_150'].append(-9999)
        if hnucleon_particle_theta['theta_150']: 
            self.ntuple['hnucleon_w_theta_150'].append(self.calculate_w(hnucleon_particle_theta['theta_150']))
        else: self.ntuple['hnucleon_w_theta_150'].append(-9999)

        self.ntuple['hw_nucleon_theta_30_50'].append(hw_nucleon_theta['theta_30_50'])
        self.ntuple['hw_nucleon_theta_50_70'].append(hw_nucleon_theta['theta_50_70'])
        self.ntuple['hw_nucleon_theta_70_90'].append(hw_nucleon_theta['theta_70_90'])
        self.ntuple['hw_nucleon_theta_90_110'].append(hw_nucleon_theta['theta_90_110'])
        self.ntuple['hw_nucleon_theta_110_130'].append(hw_nucleon_theta['theta_110_130'])
        self.ntuple['hw_nucleon_theta_130_150'].append(hw_nucleon_theta['theta_130_150'])
        self.ntuple['hw_nucleon_theta_150'].append(hw_nucleon_theta['theta_150'])

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
            self.ntuple['hnucleon_ke'], 
            self.ntuple['hproton_ke'], 
            self.ntuple['hneutron_ke']],
            np.linspace(0, 5000, 501), 
            labels=['Hardest Nucleon', 'Hardest Proton', 'Hardest Neutron'], 
            ylog=True, 
            x_label='Kinetic Energy (MeV)')

        plt.create_root_hist('hnucleon_ke', 
                             self.ntuple['hnucleon_ke'], 
                             500, 0, 5000,
                             'Kinetic Energy (MeV)', 
                             color=self.colors[0])

        plt.create_root_hist('hproton_ke', 
                             self.ntuple['hproton_ke'], 
                             500, 0, 5000,
                             'Kinetic Energy (MeV)', 
                             color=self.colors[1])

        plt.create_root_hist('hneutron_ke', 
                             self.ntuple['hneutron_ke'], 
                             500, 0, 5000,
                             'Kinetic Energy (MeV)', 
                             color=self.colors[2])

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

        theta = 30
        for itheta in xrange(0, 6):
            plt.plot_hist(self.ntuple['hnucleon_ke_theta_%s_%s' % (theta, theta + 20)], 
                np.linspace(0, 5000, 501), 
                ylog=True,
                x_label='Kinetic Energy (%s < $\theta$ < %s) (MeV)' % (theta, theta + 20))

            plt.create_root_hist('hnucleon_ke_theta_%s_%s' % (theta, theta + 20), 
                             self.ntuple['hnucleon_ke_theta_%s_%s' % (theta, theta + 20)], 
                             500, 0, 5000,
                             'Kinetic Energy (%s < #theta %s) (MeV)' % (theta, theta + 20), 
                             color=self.colors[0])

            plt.plot_hist(self.ntuple['hnucleon_w_theta_%s_%s' % (theta, theta + 20)], 
                np.linspace(0, 5000, 501), 
                ylog=True,
                x_label='W (%s < $\theta$ < %s) (MeV)' % (theta, theta + 20))

            plt.create_root_hist('hnucleon_w_theta_%s_%s' % (theta, theta + 20), 
                             self.ntuple['hnucleon_w_theta_%s_%s' % (theta, theta + 20)], 
                             500, 0, 5000,
                             'W (%s < #theta %s) (MeV)' % (theta, theta + 20), 
                             color=self.colors[0])


            plt.plot_hist(self.ntuple['hw_nucleon_theta_%s_%s' % (theta, theta + 20)], 
                np.linspace(0, 5000, 501), 
                ylog=True,
                x_label='Highest W, %s < theta < %s (MeV)' % (theta, theta + 20))

            plt.create_root_hist('hw_nucleon_theta_%s_%s' % (theta, theta + 20), 
                             self.ntuple['hw_nucleon_theta_%s_%s' % (theta, theta + 20)], 
                             500, 0, 5000,
                             'Highest W, %s < theta < %s (MeV)' % (theta, theta + 20), 
                             color=self.colors[0])
            
            theta = theta + 20

        plt.plot_hist(self.ntuple['hnucleon_ke_theta_150'], 
                np.linspace(0, 5000, 501), 
                ylog=True,
                x_label='Kinetic Energy, theta > 150 (MeV)')

        plt.create_root_hist('hnucleon_ke_theta_150', 
                             self.ntuple['hnucleon_ke_theta_150'], 
                             500, 0, 5000,
                             'Kinetic Energy, theta > 150 (MeV)', 
                             color=self.colors[0])

        plt.plot_hist(self.ntuple['hnucleon_w_theta_150'], 
                np.linspace(0, 5000, 501), 
                ylog=True,
                x_label='W ($\theta > 150$) (MeV)')

        plt.create_root_hist('hnucleon_w_theta_150', 
                             self.ntuple['hnucleon_ke_theta_150'], 
                             500, 0, 5000,
                             'W, (#theta > 150) (MeV)', 
                             color=self.colors[0])

        plt.plot_hist(self.ntuple['hw_nucleon_theta_150'], 
                np.linspace(0, 5000, 501), 
                ylog=True,
                x_label='Highest W, theta > 150 (MeV)')

        plt.create_root_hist('hw_nucleon_theta_150', 
                             self.ntuple['hw_nucleon_theta_150'], 
                             500, 0, 5000,
                             'Highest W, theta > 150 (MeV)', 
                             color=self.colors[0])


        plt.close()

    def calculate_w(self, particle): 
        pvec = particle.getMomentum()
        p = la.norm(pvec)
        ke = particle.getEnergy() - particle.getMass()
        return 0.5*(p + ke)*(1.12 - 0.5*(pvec[2]/p))
