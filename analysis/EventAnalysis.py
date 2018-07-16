
from rootpy.tree import Tree

from EventModels import Event

class EventAnalysis(object): 

    def __init__(self): 
        self.tree = None

    def initialize(self): 
        self.tree = Tree('event_ntuple', model=Event)
        self.count = 0

    def process(self, event):
        
        self.tree.event_number = event.get_event_number()
        self.count += 1
        self.tree.event_count = self.count
        self.tree.fill(reset=True)

    def finalize(self): 
        self.tree.write()
        


        
