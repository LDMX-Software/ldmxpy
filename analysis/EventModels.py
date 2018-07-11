
from rootpy.tree import TreeModel
from rootpy.tree import IntCol, FloatCol

import rootpy.stl as stl, ROOT

class Event(TreeModel): 

    event_number = IntCol()
    event_count = IntCol()

class HcalEvent(TreeModel):

    passes_hcal_veto    = IntCol() 

    max_pe = FloatCol()
    max_pe_fid = FloatCol()
    max_pe_layer = FloatCol()
    max_pe_layer_fid = FloatCol()
    total_hits = FloatCol()
    total_pe = FloatCol()
    total_pe_fid = FloatCol()

class EcalEvent(TreeModel):

    ecal_dhit_count         = IntCol()
    
    average_ecal_layer_hit  = FloatCol()
    ecal_layer_std          = FloatCol()
    ecal_max_denergy_cell   = FloatCol()
    ecal_max_layer_hit      = FloatCol()
    ecal_x_pos_mean         = FloatCol()
    ecal_x_pos_std          = FloatCol()
    ecal_y_pos_mean         = FloatCol()
    ecal_y_pos_std          = FloatCol()
    ecal_shower_rms         = FloatCol()
    ecal_summed_tight_iso   = FloatCol()
    
    total_ecal_denergy      = FloatCol()

    ecal_layer1_hit_count  = FloatCol()
    ecal_layer1_energy_sum = FloatCol()
    ecal_layer1_hit_count_no_rec  = FloatCol()
    ecal_layer1_energy_sum_no_rec = FloatCol()
    

    #ecal_dhit_energy        = stl.vector("double")
    #ecal_dhit_layer         = stl.vector("double")
    
    #ecal_hit_energy         = stl.vector("double")
    #ecal_hit_x              = stl.vector("double")
    #ecal_hit_y              = stl.vector("double")
    #ecal_is_recoil_e_hit    = stl.vector('int') 
    
    recoil_e_ecal_sp_x = FloatCol()
    recoil_e_ecal_sp_y = FloatCol()
    recoil_e_ecal_sp_z = FloatCol()
    
    recoil_e_ecal_sp_p  = FloatCol()
    recoil_e_ecal_sp_px = FloatCol()
    recoil_e_ecal_sp_py = FloatCol()
    recoil_e_ecal_sp_pz = FloatCol()

    bdt_prob = FloatCol()
    passes_ecal_veto    = IntCol() 

class TriggerEvent(TreeModel): 

    triggered = IntCol()

class TriggerPadEvent(TreeModel): 

    total_hits = IntCol()

class TrackerEvent(TreeModel):

    recoil_track_count       = IntCol()
    recoil_loose_track_count = IntCol()
    recoil_axial_track_count = IntCol()

    recoil_hits_count  = IntCol()

    recoil_hits_count_l1 = IntCol()
    recoil_hits_count_l2 = IntCol()
    recoil_hits_count_l3 = IntCol()
    recoil_hits_count_l4 = IntCol()
    recoil_hits_count_l5 = IntCol()
    recoil_hits_count_l6 = IntCol()
    recoil_hits_count_l7 = IntCol()
    recoil_hits_count_l8 = IntCol()
    recoil_hits_count_l9 = IntCol()
    recoil_hits_count_l10 = IntCol()

    recoil_hits_count_l10_no_track = IntCol()

    recoil_charge_total_l1  = FloatCol()
    recoil_charge_total_l2  = FloatCol()
    recoil_charge_total_l3  = FloatCol()
    recoil_charge_total_l4  = FloatCol()
    recoil_charge_total_l5  = FloatCol()
    recoil_charge_total_l6  = FloatCol()
    recoil_charge_total_l7  = FloatCol()
    recoil_charge_total_l8  = FloatCol()
    recoil_charge_total_l9  = FloatCol()
    recoil_charge_total_l10 = FloatCol()

    tagger_hits_count  = IntCol()
    
    tagger_hits_count_l1 = IntCol()
    tagger_hits_count_l2 = IntCol()
    tagger_hits_count_l3 = IntCol()
    tagger_hits_count_l4 = IntCol()
    tagger_hits_count_l5 = IntCol()
    tagger_hits_count_l6 = IntCol()
    tagger_hits_count_l7 = IntCol()
    tagger_hits_count_l8 = IntCol()
    tagger_hits_count_l9 = IntCol()
    tagger_hits_count_l10 = IntCol()
    tagger_hits_count_l11 = IntCol()
    tagger_hits_count_l12 = IntCol()
    tagger_hits_count_l13 = IntCol()
    tagger_hits_count_l14 = IntCol()

    tagger_charge_total_l1  = FloatCol()
    tagger_charge_total_l2  = FloatCol()
    tagger_charge_total_l3  = FloatCol()
    tagger_charge_total_l4  = FloatCol()
    tagger_charge_total_l5  = FloatCol()
    tagger_charge_total_l6  = FloatCol()
    tagger_charge_total_l7  = FloatCol()
    tagger_charge_total_l8  = FloatCol()
    tagger_charge_total_l9  = FloatCol()
    tagger_charge_total_l10 = FloatCol()
    tagger_charge_total_l11 = FloatCol()
    tagger_charge_total_l12 = FloatCol()
    tagger_charge_total_l13 = FloatCol()
    tagger_charge_total_l14 = FloatCol()

class SignalEvent(TreeModel): 

    ap_mass           = FloatCol()
    
    recoil_e_truth_p  = FloatCol() 
    recoil_e_truth_pt = FloatCol() 
    recoil_e_truth_px = FloatCol() 
    recoil_e_truth_py = FloatCol() 
    recoil_e_truth_pz = FloatCol()

class PhotoNuclearEvent(TreeModel): 
    
    pn_particle_mult    = IntCol()
    pn_gamma_energy     = FloatCol()
    pn_gamma_vertex_z   = FloatCol()
    pn_gamma_int_z      = FloatCol()
    
    recoil_e_truth_p  = FloatCol() 
    recoil_e_truth_pt = FloatCol() 
    recoil_e_truth_px = FloatCol() 
    recoil_e_truth_py = FloatCol() 
    recoil_e_truth_pz = FloatCol()
    recoil_e_vertex_x = FloatCol()
    recoil_e_vertex_y = FloatCol()
    recoil_e_vertex_z = FloatCol()

    pn_weight = FloatCol()

    lead_hadron_ke       = FloatCol()
    lead_hadron_pdg_id   = FloatCol()
    lead_hadron_theta_z  = FloatCol()

    lead_proton_ke       = FloatCol(default=-9999)
    lead_proton_theta_z  = FloatCol(default=-9999)

    lead_neutron_ke      = FloatCol(default=-9999)
    lead_neutron_theta_z = FloatCol(default=-9999)

    lead_pion_ke         = FloatCol(default=-9999)
    lead_pion_theta_z    = FloatCol(default=-9999)
    
    single_trk_p   = FloatCol()
    single_trk_pt  = FloatCol()
    single_trk_pdg = IntCol()

    is_single_neutron = IntCol(default=-9999)
    is_dineutron      = IntCol(default=-9999)
    is_diproton       = IntCol(default=-9999)
    is_pn             = IntCol(default=-9999)
    is_ghost          = IntCol(default=0)

class ElectroNuclearEvent(TreeModel): 
    
    en_particle_mult = IntCol()
    en_reaction_z    = FloatCol()

    single_trk_p   = FloatCol()
    single_trk_pdg = IntCol()

    recoil_e_tp  = FloatCol(default=-9999) 
    recoil_e_tpt = FloatCol(default=-9999) 
    recoil_e_tpx = FloatCol(default=-9999) 
    recoil_e_tpy = FloatCol(default=-9999) 
    recoil_e_tpz = FloatCol(default=-9999)
    recoil_e_vx  = FloatCol(default=-9999)
    recoil_e_vy  = FloatCol(default=-9999)
    recoil_e_vz  = FloatCol(default=-9999)

    en_weight = FloatCol()

    lead_hadron_ke      = FloatCol()
    lead_hadron_pdg_id  = FloatCol()
    lead_hadron_theta_z = FloatCol()

    highest_w_nucleon_w = FloatCol()
