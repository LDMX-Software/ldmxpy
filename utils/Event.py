
import ROOT as r 

from rootpy.io import root_open
from rootpy.io import DoesNotExist

class Event(object):

    def __init__(self, config):
    
        # Get the path for the event lib
        r.gSystem.Load(config['EventLib'][0])
        
        self.rfile = None
        self.tree = None
        self.entry = 0
       
        self.event_header = r.ldmx.EventHeader()

        self.collections = {}
        for collection in config['Collections']:
            self.collections[collection.keys()[0]] = r.TClonesArray(collection.values()[0])

    def load_file(self, rfile_path, tree_name):
        
        self.rfile = root_open(rfile_path)
        try: 
            self.tree = self.rfile.Get(tree_name)
        except DoesNotExist: 
            print 'Tree does not exist.'
        
        for name, collection in self.collections.iteritems():
            self.tree.SetBranchAddress(name, collection)
        self.tree.SetBranchAddress("EventHeader", 
                r.AddressOf(self.event_header))

        self.entry = 0

    def close_file(self):
        if self.rfile: self.rfile.Close()

    def next_event(self):
        if self.entry >= self.tree.GetEntries(): return False
        
        if (self.entry)%1000 == 0 : print "Event %s" % (self.entry + 1)
        
        self.tree.GetEntry(self.entry)
        self.entry += 1
        return True

    def collection_exist(self, collection_name):
        if collection_name in self.collections: return True
        else: return False

    def get_collection(self, collection_name):
        return self.collections[collection_name]

    def get_event_number(self):
        return self.event_header.getEventNumber()

    def get_tree(self):
        return self.tree

    def get_file_name(self): 
        return self.rfile.GetName()

    def get_weight(self): 
        return self.event_header.getWeight()
