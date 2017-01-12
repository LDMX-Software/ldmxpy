
import matplotlib
import matplotlib.pyplot as plt

from matplotlib.backends.backend_pdf import PdfPages

class Plotter(object):

    def __init__(self, file_path): 
        
        plt.style.use('bmh')
        matplotlib.rcParams.update({'font.size': 12})
        matplotlib.rcParams['axes.facecolor'] = 'white'
        matplotlib.rcParams['legend.numpoints'] = 1
        
        self.pdf = PdfPages(file_path)

    def plot_hist(self, values, bins, **params):
        
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 10))
        
        label=None
        if 'label' in params:
            label=params['label']
            ax.legend()
        
        if 'ylog' in params:
            ax.set_yscale('symlog')

        if 'x_label' in params:
            ax.set_xlabel(params['x_label'])

        ax.hist(values, bins, histtype='step', lw=1.5, label=label)

        self.pdf.savefig(bbox_inches='tight')
        plt.close()

    def close(self):
        self.pdf.close()



        

