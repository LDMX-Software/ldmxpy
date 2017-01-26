
import ROOT as r

class Event(object):

    def __init__(self, event_lib_path):
    
        r.gSystem.Load(event_lib_path)
        
        self.rfile = None
        self.tree = None
        self.entry = 0
       
        self.event_header = r.ldmx.EventHeader()

        self.collections = {}
        self.collections['SimParticles'] = r.TClonesArray('ldmx::SimParticle')
        self.collections['EcalHits'] = r.TClonesArray('ldmx::EcalHit')
        self.collections['HcalHits'] = r.TClonesArray('ldmx::HcalHit')
    
    def load_file(self, rfile_path):
        self.rfile = r.TFile(rfile_path)
        
        self.tree = self.rfile.Get("LDMX_Events")
        self.tree.SetBranchAddress("EventHeader",  r.AddressOf(self.event_header))
        self.tree.SetBranchAddress("SimParticles_sim", r.AddressOf(self.collections['SimParticles']))
        self.tree.SetBranchAddress("ecalDigis_recon", r.AddressOf(self.collections['EcalHits']))
        self.tree.SetBranchAddress("hcalDigis_recon", r.AddressOf(self.collections['HcalHits']))

        self.entry = 0

    def close_file(self):
        if self.rfile: self.rfile.Close()

    def next_event(self):
        if self.entry >= self.tree.GetEntries(): return False
        
        if (self.entry)%100 == 0 : print "Event %s" % (self.entry + 1)
        
        self.tree.GetEntry(self.entry)
        self.entry += 1
        return True

    def get_collection(self, collection_name):
        return self.collections[collection_name]

    def get_event_number(self):
        return self.event_header.getEventNumber()
