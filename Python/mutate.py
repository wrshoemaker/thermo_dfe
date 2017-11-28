from __future__ import division
import random, math, os, imp, sys, argparse
import numpy as np
from pyrosetta import init, pose_from_pdb, get_fa_scorefxn, \
                    standard_packer_task, Pose, MoveMap
from rosetta import protocols
from pyrosetta.toolbox.cleaning import cleanATOM
#from pyrosetta.toolbox import mutate_residue
from rosetta.protocols.simple_moves import RotamerTrialsMover
from rosetta.protocols.simple_moves import MinMover
from rosetta.utility import vector1_bool
from rosetta.core.chemical import aa_from_oneletter_code

#protocols.simple_moves.MinMover()
#protocols.simple_moves.RotamerTrialsMover()

#toolbox.cleaning.cleanATOM()

#cleaning = imp.load_source('cleaning', cleaning_path)

mydir = os.path.expanduser("~/GitHub/thermo_dfe/")

def mutate(pdb_name, n_muts = 100):
    OUT = open(mydir + '/data/deltas.txt', 'w')
    #takes name of pdb file without the extention
    init(extra_options='-mute basic -mute core')
    # Constants
    PACK_RADIUS = 10.0
    #Amino acids, notice there is no C
    AAs = ("A","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
    #Number of mutations to accept
    max_accept_mut = 600
    #Population size
    N = 100
    #Beta (temp term)
    beta = 1
    pdb = pdb_name + ".pdb"
    cleanATOM(pdb)
    pdb_clean_name = pdb_name + '.clean.pdb'
    initial_pose = pose_from_pdb(pdb_clean_name)

    #Set up ScoreFunction
    sf = get_fa_scorefxn()

    #Set up MoveMap.
    mm = MoveMap()
    #change these for more or less flexability
    mm.set_bb(True)
    mm.set_chi(True)

    #Pack and minimize initial pose to remove clashes.
    pre_pre_packing_score = sf(initial_pose)

    task = standard_packer_task(initial_pose)
    task.restrict_to_repacking()
    task.or_include_current(True)
    pack_rotamers_mover = RotamerTrialsMover(sf, task)
    pack_rotamers_mover.apply(initial_pose)

    min_mover = MinMover()
    min_mover.movemap(mm)
    min_mover.score_function(sf)
    min_mover.min_type('dfpmin_armijo_nonmonotone')
    min_mover.apply(initial_pose)

    post_pre_packing_score = sf(initial_pose)

    #Set threshold for selection
    threshold = pre_pre_packing_score/2

    #number of residues to select from
    n_res = initial_pose.total_residue()

    max_accept_mut = 1500
    #start sim
    i=0
    gen=0
    columns = ['Number', 'Selection', 'P_fix', 'P_fix_MH']
    print>> OUT, '\t'.join(columns)
    while i < max_accept_mut:
        #update the number of generations that have pased
        gen+=1

        #pick a place to mutate
        mut_location = random.randint(1, n_res)

        #get the amino acid at that position
        res = initial_pose.residue(mut_location)
        #don't mess with C, just choose again
        if res.name1() == 'C':
            mut_location = random.randint(1, n_res)
            #get the amino acid at that position
            res = initial_pose.residue(mut_location)

        #choose the amino acid to mutate to
        new_mut_key = random.randint(0,len(AAs)-1)

        proposed_res = AAs[new_mut_key]

        #don't bother mutating to the same amino acid it just takes more time
        if proposed_res == res.name1():
            new_mut_key = random.randint(0,len(AAs)-1)
            proposed_res = AAs[new_mut_key]

        #make the mutation
        mutant_pose = mutate_residue(initial_pose, mut_location, proposed_res, PACK_RADIUS, sf)
        #score mutant
        variant_score = sf(mutant_pose)

        probability = calc_prob_fix(variant_score, post_pre_packing_score, N, beta, threshold)
        probability_mh = calc_prob_mh(variant_score, post_pre_packing_score, N, beta, threshold)
        f_i = calc_x_fix(post_pre_packing_score, beta, threshold)
        f_j = calc_x_fix(variant_score, beta, threshold)
        s = math.log(f_j) - math.log(f_i)
        #test to see if mutation is accepted
        if np.isnan(probability) == True:
            continue
        rndm = random.random()
        if (i < 100):
            if (rndm < probability):
                initial_pose = mutant_pose
                post_pre_packing_score = variant_score
            else:
                continue
        # Assuming 100 burn in phase, take this if out if you want to store everything
        else:
            data = [str(i), str(s), str(probability), str(probability_mh)]
            print>> OUT, '\t'.join(data)
            #save name and energy change
            #data.append(variant_name + "," + str(variant_score) + "," + str(variant_score - post_pre_packing_score) + "," + str(probability) + "," + str(gen) + "\n")
            #pdb_name=str(i)+".pdb"
            #mutant_pose.dump_pdb(pdb_name)

        print i, s, probability, rndm
        #update number of accepts
        i += 1

    OUT.close()





# replaces the residue at  <resid>  in  <pose>  with  <new_res>  with repacking
def mutate_residue(pose, mutant_position, mutant_aa,
        pack_radius = 0.0, pack_scorefxn = '' ):
    """
    Replaces the residue at  <mutant_position>  in  <pose>  with  <mutant_aa>
        and repack any residues within  <pack_radius>  Angstroms of the mutating
        residue's center (nbr_atom) using  <pack_scorefxn>
    note: <mutant_aa>  is the single letter name for the desired ResidueType

    example:
        mutate_residue(pose, 30, A)
    See also:
        Pose
        PackRotamersMover
        MutateResidue
        pose_from_sequence
    """
    #### a MutateResidue Mover exists similar to this except it does not pack
    ####    the area around the mutant residue (no pack_radius feature)
    #mutator = MutateResidue(mutant_position, mutant_aa)
    #mutator.apply(test_pose)

    if pose.is_fullatom() == False:
        IOError( 'mutate_residue only works with fullatom poses' )

    test_pose = Pose()
    test_pose.assign(pose)

    # create a standard scorefxn by default
    if not pack_scorefxn:
        pack_scorefxn = get_fa_scorefxn() #  create_score_function('standard')

    task = standard_packer_task(test_pose)

    # the Vector1 of booleans (a specific object) is needed for specifying the
    #    mutation, this demonstrates another more direct method of setting
    #    PackerTask options for design
    #aa_bool = rosetta.utility.vector1_bool()
    aa_bool = vector1_bool()
    # PyRosetta uses several ways of tracking amino acids (ResidueTypes)
    # the numbers 1-20 correspond individually to the 20 proteogenic amino acids
    # aa_from_oneletter returns the integer representation of an amino acid
    #    from its one letter code
    # convert mutant_aa to its integer representation
    #mutant_aa = core.chemical.aa_from_oneletter_code(mutant_aa)
    mutant_aa = aa_from_oneletter_code(mutant_aa)

    # mutation is performed by using a PackerTask with only the mutant
    #    amino acid available during design
    # to do this, construct a Vector1 of booleans indicating which amino acid
    #    (by its numerical designation, see above) to allow
    for i in range(1, 21):
        # in Python, logical expression are evaluated with priority, thus the
        #    line below appends to aa_bool the truth (True or False) of the
        #    statement i == mutant_aa
        aa_bool.append( i == int(mutant_aa) )

    # modify the mutating residue's assignment in the PackerTask using the
    #    Vector1 of booleans across the proteogenic amino acids
    task.nonconst_residue_task(mutant_position
        ).restrict_absent_canonical_aas(aa_bool)

    # prevent residues from packing by setting the per-residue "options" of
    #    the PackerTask
    center = pose.residue(mutant_position).nbr_atom_xyz()
    for i in range(1, pose.total_residue() + 1):
        # only pack the mutating residue and any within the pack_radius
        if not i == mutant_position or center.distance_squared(
                test_pose.residue(i).nbr_atom_xyz()) > pack_radius**2:
            task.nonconst_residue_task(i).prevent_repacking()

    # apply the mutation and pack nearby residues
    packer = protocols.simple_moves.PackRotamersMover(pack_scorefxn, task)
    packer.apply(test_pose)

    return test_pose


#score functions for met-hastings selection
def calc_prob_mh(stab_mut, stab_org, N, beta, thresholds):

  xi = calc_x_fix(stab_org, beta, thresholds)
  xj = calc_x_fix(stab_mut, beta, thresholds)

  if xj > xi:
    return((1.0))
  else:
    #change if you want to use exp expression, you also need to change how the score is calculated
    #exponent = -2 * float(N) * (xi - xj)
    ans=pow(float(xj/xi),2*float(N)-2)
    return(ans)
    #return(safe_calc(exponent))

#score functions for met-hastings selection
def calc_prob_fix(stab_mut, stab_org, N, beta, thresholds):

  xi = calc_x_fix(stab_org, beta, thresholds)
  xj = calc_x_fix(stab_mut, beta, thresholds)

  if xi == xj:
    return(1/float(N))

  if(xj==0.0):
    return(0.0)

  try:
    p =((1-pow( (xi/xj),2)) /(1-pow( (xi/xj), (2 * float(N) )) ) )
  except OverflowError as e:
    p = 0.0
  return (p)


#score considering the thresh value
def calc_x_fix(data, beta, threshold):
  total = 0
  exponent = float(beta) * (float(data) - float(threshold))
  total = 1/(safe_calc(exponent) + 1)
  return(total)


def safe_calc(exponent):
  if exponent > 700:
    #print("system maxed")
    return(sys.float_info.max)
  else:
    return(math.exp(exponent))

mutate('re_numb_1qhw')
