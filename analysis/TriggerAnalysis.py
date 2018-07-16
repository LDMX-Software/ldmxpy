
from rootpy.tree import Tree

from EventModels import TriggerEvent

class TriggerAnalysis(object): 

    def __init__(self): 
        self.tree = None

    def initialize(self): 
        self.tree = Tree('trigger_ntuple', model=TriggerEvent)

    def process(self, event): 

        # Get the collection of trigger results from the collection
        trigger_results = event.get_collection('Trigger_recon')

        # Set the trigger flag.
        if trigger_results[0].passed(): self.tree.triggered =  1
        else: self.tree.triggered =  0
       
        # Fill the tree
        self.tree.fill(reset=True)

    def finalize(self): 
      
        # Write the tree to the file
        self.tree.write()

