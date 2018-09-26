
from rootpy.tree import TreeModel
from rootpy.tree import IntCol, FloatCol, FloatArrayCol

import cppyy


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

    ecal_dhit_count         = IntCol(default=0)
    
    average_ecal_layer_hit  = FloatCol(default=-9999)
    ecal_max_denergy_cell   = FloatCol(default=-9999)
    ecal_max_layer_hit      = FloatCol(default=-9999)
  
    vecal_dhit_count         = FloatCol(default=-9999)
    vtotal_ecal_denergy      = FloatCol(default=0)
    vecal_summed_tight_iso   = FloatCol(default=-9999)
    vecal_max_denergy_cell   = FloatCol(default=-9999)
    vecal_shower_rms         = FloatCol(default=-9999)
    vecal_x_pos_std          = FloatCol(default=-9999)
    vecal_y_pos_std          = FloatCol(default=-9999)
    vaverage_ecal_layer_hit  = FloatCol(default=-9999)
    vecal_max_layer_hit      = FloatCol(default=-9999)
    vecal_layer_std          = FloatCol(default=-9999)

    cylinder_0_1_layer_0_0   = FloatCol(default=-9999)
    cylinder_0_1_layer_1_2   = FloatCol(default=-9999)
    cylinder_0_1_layer_3_6   = FloatCol(default=-9999)
    cylinder_0_1_layer_7_14  = FloatCol(default=-9999)
    cylinder_0_1_layer_15    = FloatCol(default=-9999)

    cylinder_1_3_layer_0_0   = FloatCol(default=-9999)
    cylinder_1_3_layer_1_2   = FloatCol(default=-9999)
    cylinder_1_3_layer_3_6   = FloatCol(default=-9999)
    cylinder_1_3_layer_7_14  = FloatCol(default=-9999)
    cylinder_1_3_layer_15    = FloatCol(default=-9999)

    cylinder_3_5_layer_0_0   = FloatCol(default=-9999)
    cylinder_3_5_layer_1_2   = FloatCol(default=-9999)
    cylinder_3_5_layer_3_6   = FloatCol(default=-9999)
    cylinder_3_5_layer_7_14  = FloatCol(default=-9999)
    cylinder_3_5_layer_15    = FloatCol(default=-9999)

    cylinder_5_layer_0_0   = FloatCol(default=-9999)
    cylinder_5_layer_1_2   = FloatCol(default=-9999)
    cylinder_5_layer_3_6   = FloatCol(default=-9999)
    cylinder_5_layer_7_14  = FloatCol(default=-9999)
    cylinder_5_layer_15    = FloatCol(default=-9999)

    ecal_dhit_energy = cppyy.gbl.std.vector('double') 
    ecal_dhit_layer  = cppyy.gbl.std.vector('int')

    total_ecal_denergy      = FloatCol(default=0)

    ecal_layer1_hit_count  = FloatCol(default=0)
    ecal_layer1_energy_sum = FloatCol(default=0)

    trigger_energy_sum = FloatCol(default=0)
    
    recoil_e_ecal_sp_x = FloatCol()
    recoil_e_ecal_sp_y = FloatCol()
    recoil_e_ecal_sp_z = FloatCol()
    
    recoil_e_ecal_sp_p  = FloatCol()
    recoil_e_ecal_sp_px = FloatCol()
    recoil_e_ecal_sp_py = FloatCol()
    recoil_e_ecal_sp_pz = FloatCol()

    bdt_prob         = FloatCol(default=-9999)
    passes_ecal_veto = IntCol(default=0) 

class TriggerEvent(TreeModel): 

    triggered = IntCol(default=-1)
    layer_sum = FloatCol(default=-9999)

class TriggerPadEvent(TreeModel): 

    total_hits = IntCol()

class TrackerEvent(TreeModel):

    recoil_track_count       = IntCol(default=-9999)
    recoil_loose_track_count = IntCol(default=-9999)
    recoil_axial_track_count = IntCol(default=-9999)

    recoil_hits_count  = IntCol(default=-9999)

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
    
    n_electrons       = IntCol()

    recoil_e_truth_p   = cppyy.gbl.std.vector('double') 
    recoil_e_truth_pt  = cppyy.gbl.std.vector('double')
    recoil_e_truth_px  = cppyy.gbl.std.vector('double')
    recoil_e_truth_py  = cppyy.gbl.std.vector('double')
    recoil_e_truth_pz  = cppyy.gbl.std.vector('double')
    recoil_e_vertex_x  = cppyy.gbl.std.vector('double')
    recoil_e_vertex_y  = cppyy.gbl.std.vector('double') 
    recoil_e_vertex_z  = cppyy.gbl.std.vector('double')

    recoil_e_sig_p  = FloatCol() 
    recoil_e_sig_pt = FloatCol() 
    recoil_e_sig_px = FloatCol() 
    recoil_e_sig_py = FloatCol() 
    recoil_e_sig_pz = FloatCol()

class PhotoNuclearEvent(TreeModel): 
    
    pn_particle_mult    = IntCol()
    pn_gamma_energy     = FloatCol()
    pn_gamma_vertex_z   = FloatCol()
    pn_gamma_int_z      = FloatCol()

    recoil_e_calc_p  = FloatCol() 
    recoil_e_calc_pt = FloatCol() 
    recoil_e_calc_px = FloatCol() 
    recoil_e_calc_py = FloatCol() 
    recoil_e_calc_pz = FloatCol()

    recoil_e_truth_p  = FloatCol() 
    recoil_e_truth_pt = FloatCol() 
    recoil_e_truth_px = FloatCol() 
    recoil_e_truth_py = FloatCol() 
    recoil_e_truth_pz = FloatCol()
    recoil_e_vertex_x  = cppyy.gbl.std.vector('double')
    recoil_e_vertex_y  = cppyy.gbl.std.vector('double') 
    recoil_e_vertex_z  = cppyy.gbl.std.vector('double')

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

    n_electrons       = IntCol(default=0)

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

    q2 = FloatCol(default=-9999)

    en_weight = FloatCol()

    lead_hadron_ke      = FloatCol()
    lead_hadron_pdg_id  = FloatCol()
    lead_hadron_theta_z = FloatCol()

    highest_w_nucleon_w = FloatCol()
    
    is_single_neutron = IntCol(default=-9999)
    is_dineutron      = IntCol(default=-9999)
    is_diproton       = IntCol(default=-9999)
    is_pn             = IntCol(default=-9999)
    is_ghost          = IntCol(default=0)
