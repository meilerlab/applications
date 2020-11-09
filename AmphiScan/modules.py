""" This script includes movers and other functions that are necessary to run AmphiScan calculations or post-calculation
analyses. 

Last edited: 11/06/2020
Alican Gulsevin """

from pyrosetta import *
import numpy as np
#import matplotlib.pyplot as plt
import pandas as pd
#import seaborn as sns
import os
import string

#Defining the spanless mover
class AddSpanlessMembraneMover(rosetta.protocols.moves.Mover):
    """Creating a new membrane mover that works without a span file based on the default parameters of mpframework2012."""
    def __init__(self):
        rosetta.protocols.moves.Mover.__init__(self)
        self.include_lips = False
        self.topology = ""
        self.anchor_rsd = 1
        self.membrane_rsd = ""
        self.thickness = 15
        self.steepness = 10
        self.membrane_core = 15
        self.user_defined = True
        self.got_spans_from_xml = False
    

    """Apply the membrane mover to initiate the membrane object"""    
    def apply(self, Pose):
        print("=" * 80)
        print("WELCOME TO THE WORLD OF MEMBRANE PROTEINS...")
        print("=" * 80 + "\n")

        #Step 1: Initialize the membrane Info Object
        membrane_pos = Pose.size()         
        numjumps = Pose.fold_tree().num_jump()
        topology = pyrosetta.rosetta.core.conformation.membrane.SpanningTopology()
        #notify the user of the spanlessness
        mem_info = pyrosetta.rosetta.core.conformation.membrane.MembraneInfo(membrane_pos, 
                                                                             numjumps, 
                                                                             self.membrane_core, 
                                                                             self.thickness,
                                                                             self.steepness,
                                                                             topology)
        #Step 2: Add membrane info object to the pose conformation
        Pose.conformation().set_membrane_info(mem_info)
        #Step 3: Calculate the center and normal values as global variables to be used later
        if self.user_defined:
            print("Setting initial membrane center and normal to position used by the user-provided membrane residue")
            global center
            center = Pose.membrane_info().membrane_center(Pose.conformation())
            global normal
            normal = Pose.membrane_info().membrane_normal(Pose.conformation())


    
    def add_membrane_virtual(self, Pose):
        """ The function to add the virtual membrane residue """    
        #1) Grab the current residue typeset and create a new residue
        residue_set = Pose.residue_type_set_for_pose()        
        
        #2) Create a new residue of type MEM 
        membrane = residue_set.get_representative_type_name3("MEM")        
        rsd = rosetta.core.conformation.ResidueFactory.create_residue(membrane)

        #3) Append residue by jump, creating a new chain
        Pose.append_residue_by_jump(rsd, self.anchor_rsd, "", "", True)
        
        #4) Reorder the fold tree to make the membrane the root
        FoldTree = Pose.fold_tree()
        FoldTree.reorder(Pose.size())
        Pose.fold_tree(FoldTree)
        
        #5) Update chain record in PDB Info
        if Pose.pdb_info() is not None:
            #set the membrane chain id consequent to the protein chain id
            chain = list(string.ascii_uppercase)
            curr_chain = Pose.pdb_info().chain(Pose.size() - 1)
            index = chain.index(curr_chain)
            new_chain = chain[index + 1]

            Pose.pdb_info().number(Pose.size())
            Pose.pdb_info().chain(Pose.size(), new_chain)
            Pose.pdb_info().obsolete(False)
        
        return Pose.size()


class HelixTools:
    """ Some python functions to calculate helix parameters """    
    def __init__(self):
        self.diff_vec = None
        self.degree = 0

 
    def calculate_screw_axis(self, Pose):
        """ In order to calculate the normal, you need to first calculate the list of 
        alpha carbons. Because the last residue of a membrane-initiated pose is the
        membrane itself, the last protein residue is the n-1th residue. The screw axis
        is calculated as the average coordinates of the helix residues 2-4 and (n-1)-(n-4).
        Note tht this won't work for helices shorter than 7 amino acids. """
        #copy the pose and remove virtual residues
        pose_clone = Pose.clone()
        pyrosetta.rosetta.core.pose.remove_nonprotein_residues(pose_clone)
        #calculate alpha carbon centroids for the first and last 4 residues
        ca_list = []
        for resno in range(1,pose_clone.size()+1):
            ca_list.append(pose_clone.residue(resno).xyz('CA'))
        first = ca_list[1:5]
        last = ca_list[pose_clone.size()-5:pose_clone.size()-1]
        
        average_first = average_vectors(first)        
        average_last = average_vectors(last)
        
        self.screw_vector = average_last - average_first
        if self.screw_vector.is_zero() is False:
            self.screw_vector.normalize()
            
        return self.screw_vector
        
    
    def calc_angle(self, v1, axis):
        """ Calculate the angle between a vector and an axis unit vector """
        if axis is 'x':
            v2 = pyrosetta.rosetta.numeric.xyzVector_double_t(1,0,0)
        elif axis is 'y':
            v2 = pyrosetta.rosetta.numeric.xyzVector_double_t(0,1,0)
        elif axis is 'z':
            v2 = pyrosetta.rosetta.numeric.xyzVector_double_t(0,0,1)
        else:
            print("Enter 'x', 'y', or 'z'")
        
        v1_mag = np.linalg.norm(v1)    
        v2_mag = np.linalg.norm(v2)
        v1_v2_dot = np.dot(v1,v2)
        v1_v2_mag_dot = np.dot(v1_mag, v2_mag)
        
        costh = v1_v2_dot / v1_v2_mag_dot
        self.degree = np.arccos(costh) * 57.2958 #radian to degree conversion

        return self.degree
        
        
class HydrophobicMoment:
    """ Tools to calculate hydrophobic moment of a helix """
    def __init__(self):
        self.modifier = {}
        self.eisenberg = {'I': 0.73, 'F': 0.61, 'V': 0.54, 'L': 0.53, 'W': 0.37, 'M': 0.26,
                          'A': 0.25, 'G': 0.16, 'C': 0.04, 'Y': 0.02, 'P': -0.07, 'T': -0.18,
                          'S': -0.26, 'H': -0.40, 'E': -0.62, 'N': -0.64, 'Q': -0.69, 'D': -0.72,
                          'K': -1.10, 'R': -1.80}
        self.kytedoolittle = {'I': 4.5, 'F': 2.8, 'V': 4.2, 'L': 3.8, 'W': -0.9, 'M': 1.9,
                              'A': 1.8, 'G': -0.4, 'C': 2.5, 'Y': -1.3, 'P': -1.6, 'T': -0.7,
                              'S': -0.8, 'H': -3.2, 'E': -3.5, 'N': -3.5, 'Q': -3.5, 'D': -3.5,
                              'K': -3.9, 'R': -4.5}
        self.HM_vector = [0., 0., 0.]
        self.HM_vector_normalized = [0., 0., 0.]
        self.HM_vector_magnitude = 0.0
        self.center_list = []
        if mtype is "eisenberg":
            self.modifier = self.eisenberg
        elif mtype is "kytedoolittle":
            self.modifier = self.kytedoolittle
        else:
            print("The modifier type has to be 'eisenberg' or 'kytedoolittle'")
            sys.exit()
    
    def calculate_hydrophobic_moment(self, Pose, mtype='kytedoolittle'):
        """ Calculate the hydrophobic moment of a helix using the difference 
            vector between the center of mass of a side chain and the corresponding
            alpha carbon of the same residue (Eisenberg et. al, 1982). """
        if mtype is "eisenberg":
            self.modifier = self.eisenberg
        elif mtype is "kytedoolittle":
            self.modifier = self.kytedoolittle
        else:
            print("The modifier type has to be 'eisenberg' or 'kytedoolittle'")
            sys.exit()
            
        pose_clone = Pose.clone()
        pyrosetta.rosetta.core.pose.remove_nonprotein_residues(pose_clone)
        hm_magnitude = 0
        for resno in range(1, pose_clone.size()+1):
            Si = calculate_sc_center_of_mass(Pose, resno) - pose_clone.residue(resno).xyz('CA')
            hm_magnitude += np.dot(Si,pyrosetta.rosetta.numeric.xyzVector_double_t(self.modifier[pose_clone.residue(resno).name1()]))
            
        return hm_magnitude
    
    
    def HMFromStructure(self, pose, mtype='kytedoolittle'):
        """ Calculate the hydrophobic moment of a helix from the xyz coordinates of its alpha carbons """
        pose_clone = pose.clone()
        pyrosetta.rosetta.core.pose.remove_nonprotein_residues(pose_clone)
        protein_center = pyrosetta.rosetta.protocols.geometry.center_of_mass(pose_clone,1,pose_clone.size())
        number_of_residues = pose_clone.size()
        sequence = pose_clone.sequence()
        
        center_list = []
        for resno in range(1, number_of_residues + 1):
            center_list.append(pose_clone.residue(resno).xyz('CA'))
        
        
        for resno in range(1, number_of_residues + 1):
            residue_pose = pyrosetta.rosetta.core.pose.Pose(pose, resno, resno)
            residue_center = pyrosetta.rosetta.core.pose.get_center_of_mass(residue_pose)
            residue_distance_from_center = residue_center - protein_center
            residue_distance_from_center = np.array(residue_distance_from_center)
            self.center_list.append(residue_distance_from_center)
        
        self.center_list = np.array(self.center_list)
        
        if mtype is "eisenberg":
            self.modifier = self.eisenberg
        elif mtype is "kytedoolittle":
            self.modifier = self.kytedoolittle
        else:
            print("The modifier type has to be 'eisenberg' or 'kyte-doolittle'")
            sys.exit()
        
        for index, residue in enumerate(sequence):
            HM_per_residue = self.center_list[index] * self.modifier[residue]
            self.HM_vector = self.HM_vector + HM_per_residue
            self.HM_vector = np.array(self.HM_vector)
            self.HM_vector_magnitude = np.sqrt(np.square(self.HM_vector[0] - protein_center[0]) + np.square(self.HM_vector[1] - protein_center[1]) + np.square(self.HM_vector[2] - protein_center[2])) 

        return pyrosetta.rosetta.numeric.xyzVector_double_t(self.HM_vector[0], self.HM_vector[1], self.HM_vector[2])
        
    
#Additional utility scripts

def average_vectors(list_name):
    """ Takes a list of Rosetta vectors as input and returns the average 
        as a Rosetta vector object. """
    col1 = sum([item[0] for item in list_name])/len(list_name)
    col2 = sum([item[1] for item in list_name])/len(list_name)
    col3 = sum([item[2] for item in list_name])/len(list_name)    
    averaged_vector = pyrosetta.rosetta.numeric.xyzVector_double_t(col1, col2, col3)
    
    return averaged_vector


def calculate_sc_center_of_mass(Pose, residue_number):
    """ Takes pose name and residue number and calculates the coordinates of the 
        heavy atom centroid of any sidechain in the pose. For glycine residues,
        the position of the sidechain hydrogen is returned. """
    pose_clone = Pose.clone()
    pyrosetta.rosetta.core.pose.remove_nonprotein_residues(pose_clone)
    atom_positions = []    
    #glycine has no sidechain heavy atoms so an exception is defined
    if pose_clone.residue(residue_number).name1() == 'G':
        sc_center_of_mass = pose_clone.residue(residue_number).xyz('2HA')
    
    else:
        #The first 4 atoms consistently belong to the backbone heavy atoms
        for i in range(5, pose_clone.residue(residue_number).nheavyatoms()+1):
            atom_positions.append(pose_clone.residue(residue_number).xyz(i))
        sc_center_of_mass = average_vectors(atom_positions)
        
    return sc_center_of_mass
        

def sns_heatmap(protein, y_angle):
    """ Creates a heatmap from the helix rotation angle, membrane depth,
        and score belonging to a particular y-angle value. Bar positions are
        set based on 0.5 A z step size. """
    df = pd.read_csv('csv/{}_{}_xz_scores.csv'.format(protein, y_angle), names=['X', 'Y', 'Z'])
    df.iloc[:,0] = df.iloc[:,0] - 100
    
    data_pivoted = df.pivot_table(index='X', columns='Y', values='Z')
    
    ax = sns.heatmap(data_pivoted, cmap='jet')
    #membrane positive
    ax.hlines(y=20.5, xmin=0, xmax=180, linewidth=2, color='black', linestyle='solid')
    ax.hlines(y=80.5, xmin=0, xmax=180, linewidth=2, color='black', linestyle='solid')
    #membrane negative
    ax.hlines(y=26.5, xmin=0, xmax=180, linewidth=1, color='blue', linestyle='dashed')
    ax.hlines(y=74.5, xmin=0, xmax=180, linewidth=1, color='blue', linestyle='dashed')
    
    plt.xlabel('x-angle (degrees)')
    plt.ylabel('Distance from membrane (A)')
    plt.savefig('{}_{}.png'.format(protein, y_angle), dpi=800)


def analyze_y(protein):
    """ Plots the y-angles versus the best Rosetta score from each structure. 
        Takes the name of the protein as a string. """
    y_values = []
    scores = []
    
    path = "results/{}/txt".format(protein)
    
    os.system("grep best %s/*.txt | awk '{print $1, $5}' > %s/y_scores.txt" % (path, path)) 
    os.system("sed -i '$d' %s/y_scores.txt" % path) 
    
    with open('{}/y_scores.txt'.format(path), 'r') as file:
        lines = file.readlines()
        for line in lines:
            words = line.split()
            #find the numeric portion of the garbage file name
            elements = words[0].split("_")
            for element in elements:
                if element.isnumeric() is True:
                    y_angle = int(element)
            #keep the y_angle and score info
            y_values.append(float(y_angle))
            scores.append(float(words[1]))
    
    df1 = pd.DataFrame({'scores': scores})
    df2 = pd.DataFrame({'y_angle': y_values})
    df3 = df2.join(df1)
    
    df3 = df3.sort_values(['y_angle'], ascending=True)
    df3.to_csv('{}/y_angle_vs_scores.csv'.format(path), index=False, header=None)
    
    plt.scatter(df3['y_angle'], df3['scores'])
    plt.xlabel("Tilt angle (A)")
    plt.ylabel("Scores (REU)")
    plt.savefig("{}/y-angle_plot.png".format(path), dpi=600)
    
    best_df3 = min(df3['scores'])
    idx = df3.loc[df3['scores'] == best_df3].index[0]
    best_y_angle = int(df3.loc[idx][0])
    
    os.system("rm %s/y_scores.txt" % path)
    return best_y_angle 
