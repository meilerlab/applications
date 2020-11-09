""" This script takes as input a pdb filename, places the helix from a certain distance from the membrane
    and runs the helix through the membrane while sampling different z (depth) and x (alpha angle) rotation values.
    If multiple tilt angles are requested for the scan, these calculations are run for each different tilt angle (beta angle)
    and output files are reported for every tilt angle separately.
    
    Last edited 11/06/2020
    Alican Gulsevin"""

#load pyrosetta and other system libraries
from pyrosetta import *
from modules import AddSpanlessMembraneMover
from modules import HelixTools
from modules import HydrophobicMoment
import numpy as np
import os
import sys
import timeit
import pandas as pd
import argparse

#start the timer
start = timeit.default_timer()

#initiate pyrosetta and pymol
init('-ex1 -ex2aro -mute core.pack.interaction_graph.interaction_graph_factory core.pack.pack_rotamers core.pack.pack_rotamers core.packer.task protocols.grafting.util')
#init()
pymol = PyMOLMover()

#shorter forms of pyrosetta modules
rotate = pyrosetta.rosetta.protocols.rigid.WholeBodyRotationMover
translate = pyrosetta.rosetta.protocols.rigid.WholeBodyTranslationMover
Vector = pyrosetta.rosetta.numeric.xyzVector_double_t

#Arguments for the code 
parser = argparse.ArgumentParser()
parser.add_argument('-input_pdb', dest='input_pdb', type=str, help="Name of the input pdb.")
parser.add_argument('-silent', dest='silent', action='store_true', default=False, help="Activate the silent mode with no z-value outputs. default: False")
parser.add_argument('-x_final', dest='x_final', type=float, default=180, help="Final x-angle value. default: 180 degrees")
parser.add_argument('-y_final', dest='y_final', type=float, default=0, help="Final y-angle value. default: 0 (xz-only scan)")
parser.add_argument('-x_increment', dest='x_increment', type=int, default=5, help="x-angle increment. default: 5 degrees")
parser.add_argument('-y_increment', dest='y_increment', type=int, default=1, help="y-angle increment. default: 1 degrees")
parser.add_argument('-z_increment', dest='z_increment', type=float, default=0.1, help="Increment of movement along the z-axis. default: 0.1 Angstrom")
parser.add_argument('-show_pymol', dest='show_pymol', action='store_true', default=False, help="The option for real-time visualization with pymol. default: False")
parser.add_argument('-keep_history', dest='keep_history', action='store_true', default=False, help="The option to stack the scan steps at pymol. Useful to make movies. default: false")
parser.add_argument('-repack_helix', dest='repacking', action='store_true', default=False, help="The option to repack helices at every rotation step. default: False")
parser.add_argument('-relax_helix', dest='relax_helix', action='store_true', default=False, help="The option to relax helices at every x-rotation step. Don't use with repacking. default:False")
parser.add_argument('-final_relax', dest='final_relax', action='store_true', default=False, help="Relax the best structure before dumping.")
parser.add_argument('-keep_starting_tilt', dest='keep_tilt', action='store_true', default=False, help="Whether to keep the starting tilt angle from the input pdb file. Useful for xz only scans")
parser.add_argument('-thickness_from_file', dest='thickness_from_file', action='store_true', default=False, help="Whether the thickness should be read from a file.")
parser.add_argument('-helix_range', dest='helix_range', type=str, default=None, help="Define a helix region. Useful only if there are multiple helices in the system.")
args = parser.parse_args()

protein = args.input_pdb
protein_tag = protein.split(sep='.')[0]

#create the folders to output the results
if not os.path.exists('results'):
    os.makedirs('results')
if not os.path.exists('results/{}'.format(protein_tag)):
    os.makedirs('results/{}'.format(protein_tag))
if not os.path.exists('results/{}/txt'.format(protein_tag)):
    os.makedirs('results/{}/txt'.format(protein_tag))
if not os.path.exists('results/{}/csv'.format(protein_tag)):
    os.makedirs('results/{}/csv'.format(protein_tag))
if not os.path.exists('results/{}/output_pdbs'.format(protein_tag)):
    os.makedirs('results/{}/output_pdbs'.format(protein_tag))
        
#load the pose based on the input value
Pose = pose_from_pdb('{}'.format(protein))

#scan parameters  
z_final = 200
starting_z = -100

#if asked, load the membrane thickness from a file. otherwise, use the default thickness of 15A
if args.thickness_from_file is True:
    thickness_from_file = None
    #tinkering with the membrane thickness
    with open('metric_data/membrane_thickness.txt', 'r') as thickness_file:
         with open('metric_data/membrane_thickness.txt', 'r') as thickness_file:
            for thicc in thickness_file:
                if protein_tag in thicc:
                    thickness_from_file = round(float(thicc.split(',', 1)[1]))
                    print("Membrane thickness is set to: {}".format(thickness_from_file))
else:
    thickness_from_file = 15
        
#initiate the spanless membrane mover fom the modules file
fm = AddSpanlessMembraneMover()
fm.thickness = thickness_from_file
fm.membrane_core = thickness_from_file
fm.add_membrane_virtual(Pose)
fm.apply(Pose)

#Create a clone for the pose
rePose = Pose.clone()

#if a helix range is defined by the user, use that as the range, else, use the pose size
if args.helix_range is not None:
    range_start = int(args.helix_range.split(sep='-')[0])
    range_end = int(args.helix_range.split(sep='-')[1])
else:
    range_start = 1
    range_end = rePose.size() - 1
    
#calculate the center of mass of the protein
cmass = pyrosetta.rosetta.core.pose.center_of_mass(rePose, range_start, range_end)
print(range_start, range_end)
   
#move the helix to the center of mass (0,0,0)
realign1 = translate(cmass.negated())
realign1.apply(rePose)
if args.show_pymol is True:
    pymol.apply(rePose)
 
#move the helix center to the starting distance if one was set
if starting_z is not 0:
    realign2 = translate(Vector(0,0,starting_z))
    realign2.apply(rePose)
    if args.show_pymol is True:
        pymol.apply(rePose)

selected_helix_pose = pyrosetta.rosetta.protocols.grafting.return_region(rePose, range_start, range_end)
print(selected_helix_pose)

#initiate the helix tools from the modules file and calculate the angles between the helix
#screw axis and the xyz coordinates
ht = HelixTools()
helix_normal = ht.calculate_screw_axis(selected_helix_pose)
angle_with_x = ht.calc_angle(helix_normal,'x') 
angle_with_y = ht.calc_angle(helix_normal,'y') 
angle_with_z = ht.calc_angle(helix_normal,'z')
if args.silent is False:
    print(angle_with_x, angle_with_y, angle_with_z)

#Rotate the helix along the z-axis to align with the x-axis
if angle_with_y < 90:
    align_x = rotate(Vector(0,0,1),
                     pyrosetta.rosetta.core.pose.center_of_mass(rePose, range_start, range_end),
                     -angle_with_x)
else:
    align_x = rotate(Vector(0,0,1),
                     pyrosetta.rosetta.core.pose.center_of_mass(rePose, range_start, range_end),
                     angle_with_x)        

align_x.apply(rePose)
if args.show_pymol is True:
    pymol.apply(rePose)

#calculate angles after aligning the x-axis
selected_helix_pose = pyrosetta.rosetta.protocols.grafting.return_region(rePose, range_start, range_end)
helix_normal = ht.calculate_screw_axis(selected_helix_pose)
angle_with_x = ht.calc_angle(helix_normal, 'x')
angle_with_y = ht.calc_angle(helix_normal, 'y')
angle_with_z = ht.calc_angle(helix_normal, 'z')
if args.silent is False:
    print(angle_with_x, angle_with_y, angle_with_z)

#if keep_tilt is selected this part is skipped and the screw axis is not forced to be 
#parallel to the x-axis. Otherwise, additional rotation opertions are done to make the 
#helix parallel to the x-axis and orthogonal to the y-axis. 
if args.keep_tilt is False:
    if angle_with_y < 90:
        align_y = rotate(Vector(1,0,0),
                         pyrosetta.rosetta.core.pose.center_of_mass(rePose, range_start, range_end),
                         -(90-angle_with_y))
    else:
        align_y = rotate(Vector(1,0,0),
                         pyrosetta.rosetta.core.pose.center_of_mass(rePose, range_start, range_end),
                         (90-angle_with_y))

    align_y.apply(rePose)
    if args.show_pymol is True:
        pymol.apply(rePose)
    
    #calculate angles after aligning the y-axis
    selected_helix_pose = pyrosetta.rosetta.protocols.grafting.return_region(rePose, range_start, range_end)
    helix_normal = ht.calculate_screw_axis(selected_helix_pose)
    angle_with_x = ht.calc_angle(helix_normal, 'x')
    angle_with_y = ht.calc_angle(helix_normal, 'y')
    angle_with_z = ht.calc_angle(helix_normal, 'z')
    if args.silent is False:
        print(angle_with_x, angle_with_y, angle_with_z)    
        
    if angle_with_z < 90:
        align_z = rotate(Vector(0,1,0),
                         pyrosetta.rosetta.core.pose.center_of_mass(rePose, range_start, range_end),
                         90-angle_with_z)
    else:
        align_z = rotate(Vector(0,1,0),
                         pyrosetta.rosetta.core.pose.center_of_mass(rePose, range_start, range_end),
                         -(angle_with_z-90)) #-(90-angle_with_z)        
    align_z.apply(rePose)
    if args.show_pymol is True:
        pymol.apply(rePose)
    
    #calculate angles after aligning the z-axis
    selected_helix_pose = pyrosetta.rosetta.protocols.grafting.return_region(rePose, range_start, range_end)
    helix_normal = ht.calculate_screw_axis(selected_helix_pose)
    angle_with_x = ht.calc_angle(helix_normal, 'x')
    angle_with_y = ht.calc_angle(helix_normal, 'y')
    angle_with_z = ht.calc_angle(helix_normal, 'z')
    if args.silent is False:
        print(angle_with_x, angle_with_y, angle_with_z)
                    
#set the membrane scoring function
mem_sfxn = pyrosetta.get_fa_scorefxn()
mem_sfxn.add_weights_from_file("mpframework_smooth_fa_2012.wts") #this membrane mover won't work with franklin2019

#define the scan region
distance_from_membrane = 10
thickness = thickness_from_file
scan_start_region = -starting_z - thickness - distance_from_membrane
scan_stop_region = -starting_z + thickness + distance_from_membrane

#scan a range of tilt angles. If the default value is selected, then only y=0 will be scanned
yAngRange = np.arange(0, args.y_final + args.y_increment, args.y_increment, dtype=int)
for y_ang in yAngRange:
    if args.silent is False:
        print("Running tilt angle {}".format(y_ang))
    if os.path.exists("results/{}/txt/{}_{}_scores.txt".format(protein_tag, protein.split(".",1)[0], y_ang)):
        os.remove("results/{}/txt/{}_{}_scores.txt".format(protein_tag, protein.split(".",1)[0], y_ang))
    if os.path.exists("results/{}/csv/{}_{}.csv".format(protein_tag, protein.split(".",1)[0], y_ang)):
        os.remove("results/{}/csv/{}_{}.csv".format(protein_tag, protein.split(".",1)[0], y_ang))
        
    #set starting values
    best_z = 0
    best_xang = 0
    score_best = 999999
    
    #create a clone of the translated pose and rotate it to set the tilt angle
    rePose2 = rePose.clone()
    rot_y  = rotate(Vector(0,1,0),
                pyrosetta.rosetta.protocols.geometry.center_of_mass(rePose2, range_start, range_end),
                y_ang)
    rot_y.apply(rePose2)
    if args.show_pymol is True:
        pymol.apply(rePose2)

    #start the scan along the z-axis
    zRange = np.arange(args.z_increment, z_final + args.z_increment, args.z_increment)
    for dist in zRange:
        rePose3 = rePose2.clone()
        trans_x  = translate(Vector(0,0,dist)) 
        trans_x.apply(rePose3)
        if args.show_pymol is True: 
            pymol.apply(rePose3)
        
        #if not within 10A of the membrane surface
        if dist < scan_start_region or dist > scan_stop_region:
            continue
    
        #else if within 10A of the membrane at either end
        elif scan_start_region <= dist <= scan_stop_region:
            xAngRange = np.arange(0, args.x_final + args.x_increment, args.x_increment)
            for x_ang in xAngRange:
                rePose4 = rePose3.clone()                    
                selected_helix_pose = pyrosetta.rosetta.protocols.grafting.return_region(rePose3, range_start, range_end) 
                rot_x  = rotate(ht.calculate_screw_axis(selected_helix_pose),
                                pyrosetta.rosetta.protocols.geometry.center_of_mass(rePose4, range_start, range_end),
                                x_ang)
                rot_x.apply(rePose4)
                if args.show_pymol is True:
                    pymol.apply(rePose4)
                
                #repack the side-chains if requested
                if args.repacking is True:
                    task_pack = standard_packer_task(rePose4)
                    task_pack.restrict_to_repacking()
                    pack_mover = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover(mem_sfxn, task_pack)
                    pack_mover.apply(rePose4)
                
                elif args.relax_helix is True:
                    #define the movemap
                    move_map = pyrosetta.rosetta.core.kinematics.MoveMap()
                    #keep the N- and C-terminal residue backbones fixed. Originally conceived for capped helices.
                    move_map.set_bb_true_range(2,rePose4.size()-1)
                    move_map.set_chi_true_range(1,rePose4.size())
                    
                    #push the movemap to relax
                    relax = pyrosetta.rosetta.protocols.relax.FastRelax()
                    relax.set_scorefxn(mem_sfxn)
                    relax.set_movemap(move_map)
                    relax.apply(rePose4)
            
                #determine the best score over all x- and z-values
                score_new = mem_sfxn(rePose4)

                #check if better than the previous rounds. If yes, store the best parameters and save as the best pose
                if score_new < score_best:
                    score_best = score_new
                    best_z = dist
                    best_xang = x_ang
                    best_pose = rePose4.clone()
                
                #write the scores at the new y-angle
                with open("results/{}/txt/{}_{}_scores.txt".format(protein_tag, protein.split(".",1)[0], y_ang), "a+") as scores:
                    scores.write("The score corresponding to the z value: "  + str(dist) + \
                                    " and x-angle " + str(x_ang) \
                                    + " " + str(mem_sfxn(rePose4)) + "\n")
                with open("results/{}/csv/{}_{}_xz_scores.csv".format(protein_tag, protein.split(".",1)[0], y_ang), "a+") as xz_scores:
                    xz_scores.write("{}, {}, {}\n".format(dist, x_ang, mem_sfxn(rePose4)))
                                            
            #record the best x- and z-value information for each parameter
            if dist == scan_stop_region:
                with open("results/{}/txt/{}_{}_scores.txt".format(protein_tag, protein.split(".",1)[0], y_ang), "a+") as best:
                    best.write("\nThe best score of {} belongs to {} Angstroms and {} x-degrees.\n".format(round(score_best,4 ), best_z, best_xang))
        
        elif args.final_relax is True:
            #define the movemap
            move_map = pyrosetta.rosetta.core.kinematics.MoveMap()
            #keep the N- and C-terminal residue backbones fixed. Originally conceived for capped helices.
            move_map.set_bb_true_range(2,best_pose.size()-1)
            move_map.set_chi_true_range(1,best_pose.size())
            
            #push the movemap to relax
            relax = pyrosetta.rosetta.protocols.relax.FastRelax()
            relax.set_scorefxn(mem_sfxn)
            relax.set_movemap(move_map)
            relax.apply(best_pose)
    #save the best pose and write the best parameters    
    best_pose.dump_pdb('results/{}/output_pdbs/best_{}_{}.pdb'.format(protein_tag, protein.split(".",1)[0], y_ang))
    with open("results/{}/txt/{}_{}_scores.txt".format(protein_tag, protein.split(".",1)[0], y_ang), "a+") as best:
        best.write("Time elapsed for completion: {} seconds.\n".format(timeit.default_timer() - start))
        
