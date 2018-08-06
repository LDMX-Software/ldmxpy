import math

from numpy import linalg as la

def get_kinetic_energy(particle): 
    return (particle.getEnergy() - particle.getMass())

def get_theta_z(particle): 
    pvec = particle.getMomentum()
    p = la.norm(pvec)
    return math.acos(pvec[2]/p)*(180/math.pi)

def get_recoil_electrons(particles): 

    # Loop through all of the particles and search for the recoil electron i.e.
    # an electron which doesn't have any parents
    recoils = []
    for particle in particles:
            
        # If the particle is an electron and has no parents, a recoil
        # has been found
        if (particle.getPdgID()) & (particle.getParentCount() == 0):
            recoils.append(particle)
            #recoil_e = particle
            #break

    # All events should contain a recoil electron
    if len(recoils) == 0: 
        raise RuntimeError('Recoil electron was not found!')

    return recoils

def get_pn_gamma(recoil_e): 
    
    # Search for the PN gamma and use it to get the PN daughters
    pn_gamma = None
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

def get_pt(p): 
    return math.sqrt(p[0]*p[0] + p[1]*p[1])
