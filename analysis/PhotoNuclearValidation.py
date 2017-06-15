from __future__ import division

import numpy as np
import Plotter

class PhotoNuclearValidation(object):

    def __init__(self):
        self.initialize()

    def initialize(self):
        self.events = []

        self.hardest_hadron_ke  = []
        self.hardest_hadron_mwp = []
        self.hardest_hadron_fwp = []
        self.hardest_hadron_theta = []
        
        self.weight = []

    def process(self, event):
        
        # Get the collection of MC particles from the event
        particles = event.get_collection('SimParticles_sim')

        weights = event.get_collection('pnWeight_recon')

        for weight in weights:
            if weight.getKineticEnergy() < 100: continue    
            self.hardest_hadron_ke.append(weight.getKineticEnergy())
            self.hardest_hadron_mwp.append(weight.getMeasuredWp())
            self.hardest_hadron_fwp.append(weight.getFitWp())
            self.hardest_hadron_theta.append(weight.getTheta())
            self.weight.append(weight.getWeight())

    def finalize(self):
        plt = Plotter.Plotter('pn_validation')

        self.hardest_hadron_ke = np.array(self.hardest_hadron_ke)
        self.hardest_hadron_mwp = np.array(self.hardest_hadron_mwp)
        self.hardest_hadron_fwp = np.array(self.hardest_hadron_fwp)
        self.hardest_hadron_theta = np.array(self.hardest_hadron_theta)
        self.weight = np.array(self.weight)
       
        plt.plot_hists([self.hardest_hadron_ke, 
                        self.hardest_hadron_ke[self.hardest_hadron_theta > 100]
                       ], 
                      np.linspace(0, 4500, 251), 
                      labels=['All', '$\theta > 100$'],
                      ylog=True, 
                      x_label='$T_{p}$ (MeV)')

        plt.plot_hists([np.sqrt(np.power(self.hardest_hadron_ke, 2) + 2*self.hardest_hadron_ke*938.272), 
                        np.sqrt(np.power(self.hardest_hadron_ke, 2) + 2*self.hardest_hadron_ke*938.272)[self.hardest_hadron_theta > 100]
                       ], 
                      np.linspace(0, 4500, 251), 
                      labels=['All', '$\theta > 100$'],
                      ylog=True, 
                      x_label='$p$ (MeV)')

        plt.plot_hists([self.hardest_hadron_ke[(self.hardest_hadron_theta > 25) &
                                               (self.hardest_hadron_theta < 35) 
                                              ],
                        self.hardest_hadron_ke[(self.hardest_hadron_theta > 55) &
                                               (self.hardest_hadron_theta < 65) 
                                              ],
                        self.hardest_hadron_ke[(self.hardest_hadron_theta > 85) &
                                               (self.hardest_hadron_theta < 95) 
                                              ],
                        self.hardest_hadron_ke[(self.hardest_hadron_theta > 115) &
                                               (self.hardest_hadron_theta < 125) 
                                              ],
                        self.hardest_hadron_ke[(self.hardest_hadron_theta > 155) &
                                               (self.hardest_hadron_theta < 165)
                                              ]
                       ], 
                      np.linspace(0, 4500, 251), 
                      labels=['30 degrees', '60 degrees', 
                              '90 degrees', '120 degrees', '160 degrees'],
                      ylog=True, 
                      x_label='$T_{p}$ (MeV)')

        x = np.linspace(0, 4500, 251)
        y = 17803.2*np.exp(-0.00807561*(x - 791.244))  
        plt.plot_hists([self.hardest_hadron_mwp,
                        self.hardest_hadron_mwp[self.hardest_hadron_theta > 100]
                       ],
                      np.linspace(0, 4500, 251), 
                      labels=['All', '$\theta > 100$'],
                      ylog=True, 
                      x_label='Measured $W_{p}$ (MeV)')

        plt.plot_graph(x, y, 0, 0,
                       ylog=True,
                       x_label='Fit $W_{p}$ (MeV)')

        plt.plot_hists([self.hardest_hadron_fwp,
                        self.hardest_hadron_fwp[self.hardest_hadron_theta > 100]
                       ],
                      np.linspace(0, 4500, 251), 
                      labels=['All', '$\theta > 100$'],
                      ylog=True, 
                      x_label='Fit $W_{p}$ (MeV)')
        
        plt.plot_hist(self.hardest_hadron_theta, 
                      np.linspace(0, 180, 360), 
                      ylog=True, 
                      x_label='$\theta$ (degrees)')

        plt.plot_hists([self.weight,
                        self.weight[self.hardest_hadron_theta > 100]
                       ],
                      np.linspace(0, 2, 100), 
                      labels=['All', '$\theta > 100$'],
                      ylog=True, 
                      x_label='PN Weight')

        plt.close()
