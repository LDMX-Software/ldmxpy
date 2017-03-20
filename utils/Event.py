
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
        self.collections['TriggerPadSimHits'] = r.TClonesArray('ldmx::SimCalorimeterHit')
        self.collections['EcalSimHits'] = r.TClonesArray('ldmx::SimCalorimeterHit')
        self.collections['RecoilSimHits'] = r.TClonesArray('ldmx::SimTrackerHit')
        self.collections['EcalHits'] = r.TClonesArray('ldmx::EcalHit')
        self.collections['HcalHits'] = r.TClonesArray('ldmx::HcalHit')
        self.collections['EcalVeto'] = r.TClonesArray('ldmx::EcalVetoResult')
        self.collections['Trigger'] = r.TClonesArray('ldmx::TriggerResult')
        self.collections['FindableTracks'] = r.TClonesArray('ldmx::FindableTrackResult')

    def load_file(self, rfile_path):
        self.rfile = r.TFile(rfile_path)
        
        self.tree = self.rfile.Get("LDMX_Events")
        self.tree.SetBranchAddress("EventHeader",      r.AddressOf(self.event_header))
        self.tree.SetBranchAddress("SimParticles_sim", r.AddressOf(self.collections['SimParticles']))
        self.tree.SetBranchAddress('TriggerPadSimHits_sim', r.AddressOf(self.collections['TriggerPadSimHits']))
        self.tree.SetBranchAddress('EcalSimHits_sim', r.AddressOf(self.collections['EcalSimHits']))
        self.tree.SetBranchAddress('RecoilSimHits_sim', r.AddressOf(self.collections['RecoilSimHits']))
        self.tree.SetBranchAddress("ecalDigis_recon",  r.AddressOf(self.collections['EcalHits']))
        self.tree.SetBranchAddress("hcalDigis_recon",  r.AddressOf(self.collections['HcalHits']))
        self.tree.SetBranchAddress('EcalVeto_recon', r.AddressOf(self.collections['EcalVeto']))
        self.tree.SetBranchAddress("Trigger_recon",  r.AddressOf(self.collections['Trigger']))
        self.tree.SetBranchAddress("FindableTracks_recon", r.AddressOf(self.collections['FindableTracks']))

        self.entry = 0

    def close_file(self):
        if self.rfile: self.rfile.Close()

    def next_event(self):
        if self.entry >= self.tree.GetEntries(): return False
        
        if (self.entry)%1000 == 0 : print "Event %s" % (self.entry + 1)
        
        self.tree.GetEntry(self.entry)
        self.entry += 1
        return True

    def get_collection(self, collection_name):
        return self.collections[collection_name]

    def get_event_number(self):
        return self.event_header.getEventNumber()

    def get_tree(self):
        return self.tree

    def get_file_name(self): 
        return self.rfile.GetName()

