
from rootpy.tree import Tree

from EventModels import HcalEvent

class HcalAnalysis(object): 

    def __init__(self): 
        self.tree = None

    def initialize(self, params): 
        self.tree = Tree('hcal_ntuple', model=HcalEvent)
        self.event_count = 0
        self.veto_count = 0

    def process(self, event): 

        self.event_count += 1
        #print 'Event: %s' % self.event_count

        hcal_hits = event.get_collection('hcalDigis_recon')

        # Loop through all of the Hcal hits and calculate the total 
        # photoelectrons in the event
        total_pe = 0
        total_pe_fid = 0
        max_pe = 0
        max_pe_fid = 0
        max_pe_layer = 0
        max_pe_layer_fid = 0
        for hit in hcal_hits:
            #if hcal_hit.getSection() == 0: 
                #total_pe_fid += hcal_hit.getPE()
                #if max_pe_fid < hcal_hit.getPE(): 
                    #max_pe_fid = hcal_hit.getPE()
                    #max_pe_layer_fid = hcal_hit.getLayer()
            total_pe += hit.getPE()
            max_pe = max(max_pe, hit.getPE())

        self.tree.max_pe = max_pe
        self.tree.max_pe_fid = max_pe_fid
        self.tree.max_pe_layer = max_pe_layer
        self.tree.max_pe_layer_fid = max_pe_layer_fid
        self.tree.total_hits = hcal_hits.GetEntriesFast()
        self.tree.total_pe = total_pe
        self.tree.total_pe_fid = total_pe_fid

        self.tree.passes_hcal_veto =  1
        if max_pe >= 8: self.tree.passes_hcal_veto =  0
        else: self.veto_count += 1
        
        self.tree.fill(reset=True)

    def finalize(self):

        self.tree.write()
        
        print '[ HcalAnalysis ]: Total Events that have a max PE < 8: %s' % self.veto_count
            
