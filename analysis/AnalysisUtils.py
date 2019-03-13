import math

from numpy import linalg as la

def get_kinetic_energy(particle): 
    return (particle.getEnergy() - particle.getMass())

def get_recoil_electrons(particles): 

    # Loop through all of the particles and search for the recoil electron i.e.
    # an electron which doesn't have any parents
    return [particle for particle in particles if ((particle.getPdgID() == 11) & (particle.getGenStatus() == 1))]

def get_pn_gamma(recoils): 
    
    # Search for the PN gamma and use it to get the PN daughters
    pn_gamma = None
    
    for recoil_e in recoils: 
        for daughter_count in xrange(0, recoil_e.getDaughterCount()):
            daughter = recoil_e.getDaughter(daughter_count)
        
            #daughter.Print()
            #if (daughter.getDaughterCount() > 0): 
                #print '###'
                #daughter.getDaughter(0).Print()

            if (daughter.getDaughterCount() > 0) and \
                    (daughter.getDaughter(0).getProcessType() == 9): 
                pn_gamma = daughter
                break
    
    # When looking at photonuclear events, a recoil electron should always have
    # a PN gamma associated with it.
    if not pn_gamma: 
        raise RuntimeError('PN gamma was not found!')

    return pn_gamma

def get_ap(particles): 

    # Loop through all of the particles and search for the A' i.e. a particle
    # with a PDG ID of 622
    aprime = None
    for particle in particles: 
        if particle.getPdgID() == 622: 
           aprime = particle
           break
    
    return aprime

def calculate_w(particle, delta): 

    pz = particle.getMomentum()[2]
    p = la.norm(particle.getMomentum())
    ke = get_kinetic_energy(particle)
    
    return 0.5*(p + ke)*(math.sqrt(1 + (delta*delta)) - delta*(pz/p))

def get_findable_tracks_map(findable_tracks):

    findable_dic = {}
    loose_dic = {}
    axial_dic = {}

    # Create a map between a sim particle and a findable track.
    for findable_track in findable_tracks:
        if (findable_track.is4sFindable() 
            or findable_track.is3s1aFindable()
                or findable_track.is2s2aFindable()):
            findable_dic[findable_track.getSimParticle()] = findable_track

        if findable_track.is2sFindable(): 
            loose_dic[findable_track.getSimParticle()] = findable_track   

        if findable_track.is2aFindable(): 
            axial_dic[findable_track.getSimParticle()] = findable_track

    return findable_dic, loose_dic, axial_dic

def classify_event(particles, threshold): 
    
    neutron_count = 0
    proton_count = 0
    pion_count = 0
    pi0_count = 0
    exotic_count = 0
    for particle in particles: 
        if get_kinetic_energy(particle) <= threshold: continue

        pdgid = abs(particle.getPdgID())
        if pdgid == 2112: neutron_count += 1
        elif pdgid == 2212: proton_count += 1
        elif pdgid == 211: pion_count += 1
        elif pdgid == 111: pi0_count += 1
        else: 
            print pdgid
            exotic_count += 1

    count = neutron_count + proton_count + pion_count + pi0_count + exotic_count
    count_a = proton_count + pion_count + pi0_count + exotic_count
    count_b = pion_count + pi0_count + exotic_count
    count_c = proton_count + neutron_count + pi0_count + exotic_count
    count_d = neutron_count + pi0_count + exotic_count
    count_e = proton_count + pi0_count + exotic_count
    count_f = neutron_count + proton_count + pion_count + exotic_count
    count_g = neutron_count + pion_count + pi0_count + exotic_count
    count_h = neutron_count + proton_count + pion_count + pi0_count
    if count == 0: return 0
    if neutron_count == 1: 
            if (count_a == 0): return 2
            elif (proton_count == 1) & (count_b == 0): return 9
    if (neutron_count == 2) & (count_a == 0): return 2
    if (neutron_count >= 3) & (count_a == 0): return 3
    if (pion_count == 1):
            if count_c == 0: return 4
            elif (proton_count == 1) & (count_d == 0): return 7   
            elif (neutron_count == 1) & (count_e == 0): return 7   
    if (pion_count == 2) & (count_c == 0): return 5
    if (pi0_count == 1) & (count_f == 0): return 6
    if (proton_count == 1) & (count_g == 0): return 8 
    if (proton_count == 2) & (count_g == 0): return 9
    if (exotic_count > 0) & (count_h == 0): return 10

    if ((neutron_count > 0) 
        & ((proton_count > 0) or (pion_count > 0) 
            or (pi0_count > 0) or (exotic_count > 0))): return 11
    elif ((proton_count > 0) 
        & ((neutron_count > 0) or (pion_count > 0) 
            or (pi0_count > 0) or (exotic_count > 0))): return 11
    elif ((pion_count > 0) 
        & ((neutron_count > 0) or (proton_count > 0) 
            or (pi0_count > 0) or (exotic_count > 0))): return 11
    elif ((pi0_count > 0) 
        & ((neutron_count > 0) or (proton_count > 0) 
            or (pion_count > 0) or (exotic_count > 0))): return 11
    elif ((exotic_count > 0) 
        & ((neutron_count > 0) or (proton_count > 0) 
            or (pion_count > 0) or (pi0_count > 0))): return 11

    return -9999
