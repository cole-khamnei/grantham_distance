import os
import pathlib

########################################################################################################################
### System Constants ###
########################################################################################################################

GRANTHAM_DISTANCE_DIR_PATH = pathlib.Path(__file__).parent.resolve()

RESOURCES_DIR_PATH = os.path.join(GRANTHAM_DISTANCE_DIR_PATH, "resources")
COVER_PLOTS_PATH = os.path.join(GRANTHAM_DISTANCE_DIR_PATH, "cover_plots")

HLA_CLASS_I_ALLELE_AMINO_ACID_SEQUENCES_FULL_PATH = os.path.join(RESOURCES_DIR_PATH, "HLA_class_I_allele_AA_sequences.fas")
HLA_CLASS_I_ALLELE_AMINO_ACID_SEQUENCES_ABS_PATH = os.path.join(RESOURCES_DIR_PATH, "HLA_class_I_allele_AA_sequences_ABS.fas")
HLA_CLASS_II_ALLELE_AMINO_ACID_SEQUENCES_ABS_PATH = os.path.join(RESOURCES_DIR_PATH, "HLA_class_II_allele_AA_sequences_ABS.fas")

GRANTHAM_DISTANCE_MATRIX_PATH = os.path.join(RESOURCES_DIR_PATH, "grantham_distance_AA_matrix.cnv")

########################################################################################################################
### End ###
########################################################################################################################
