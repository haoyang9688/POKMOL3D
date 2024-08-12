import os
import re
import glob
import math
from rdkit import Chem
from rdkit.Chem import QED
from rdkit.Chem import rdMolDescriptors
from itertools import combinations
from rdkit.Chem import AllChem
from rdkit.DataStructs import TanimotoSimilarity
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit import DataStructs
import pickle
import gzip
import subprocess
import sys
import argparse
from argparse import Namespace
from pathlib import Path
import yaml
from easydict import EasyDict

# path
sys.path.append('/data/XXXX/POKMOL-3D')
class Target_failure_rate_Calculator:
    def __init__(self, config):
        self.config = config
class Sampling_success_rate_Calculator:
    def __init__(self, config):
        self.config = config        
class Validity_Calculator:
    def __init__(self, config):
        self.config = config
class Uniqueness_Calculator:
    def __init__(self, config):
        self.config = config
class Usability_Calculator:
    def __init__(self, config):
        self.config = config
class QED_Calculator:
    def __init__(self, config):
        self.config = config
class SA_Calculator:
    def __init__(self, config):
        self.config = config   
class Tdiv_Calculator:
    def __init__(self, config):
        self.config = config
class General_molecular_quality_metrics_Calculator:
    def __init__(self, config):
        self.config = config
class Structural_properties_2D_Calculator:
    def __init__(self, config):
        self.config = config
class Structural_properties_3D_Calculator:
    def __init__(self, config):
        self.config = config
class structural_summary_2D3D_Calculator:
     def __init__(self, config):
        self.config = config
class JSD_Calculator:
    def __init__(self, config):
        self.config = config        
class active_recovery_Calculator:
    def __init__(self, config):
        self.config = config  
class ReDocking_Calculator:
    def __init__(self, config):
        self.config = config
class In_situ_Docking_Calculator:
    def __init__(self, config):
        self.config = config
class RMSD_Calculator:
    def __init__(self, config):
        self.config = config
class Target_binding_metrics_Calculator:
    def __init__(self, config):
        self.config = config                                          
class MainProcessor:
    def __init__(self, config):
        self.config = config
##################Model-quelity########################################################
    def Target_failure_rate(self):
        from utils.Target_failure_rate import calculate_Target_failure_rate
        try:
            calculate_Target_failure_rate(self.config)
            print("Target_failure_rate completed")  
        except Exception as e:
            print(f"Target_failure_rate failed: {e}")

    def Sampling_success_rate(self):
        from utils.Sampling_success_rate import calculate_Sampling_success_rate
        try:
            calculate_Sampling_success_rate(self.config)
            print("Sampling_success_rate completed")  
        except Exception as e:
            print(f"Sampling_success_rate failed: {e}")        
##################Validity###################################################################
    def Validity(self):
        from utils.Validity import calculate_validity
        from utils.validity_summary import calculate_validity_summary
        from utils.Convert_smi import convert_smi
        try:
            calculate_validity(self.config)
            print("Validity calculation completed")
        except Exception as e:
            print(f"Validity calculation failed: {e}")
            raise  
        try:
            calculate_validity_summary(self.config)
            print("Validity summary calculation completed")
        except Exception as e:
            print(f"Validity summary calculation failed: {e}")
            raise  
        try:
            convert_smi(self.config)
            print("SMI conversion completed")
        except Exception as e:
            print(f"SMI conversion failed: {e}")
            raise  
##################Uniqueness###################################################################
    def Uniqueness(self):
        from utils.Uniqueness import calculate_uniqueness
        from utils.Convert_txt import calculate_convert_txt
        from utils.uniqueness_summary import calculate_uniqueness_summary
        from utils.Uniqueness_sdf_match import match_and_copy_sdf_files
        try:
            calculate_uniqueness(self.config)
            print("Uniqueness calculation completed")
        except Exception as e:
            print(f"Uniqueness calculation failed: {e}")
            raise  
        try:
            calculate_convert_txt(self.config)
            print("Convert txt calculation completed")
        except Exception as e:
            print(f"Convert txt calculation failed: {e}")
            raise  
        try:
            calculate_uniqueness_summary(self.config)
            print("Uniqueness summary calculation completed")
        except Exception as e:
            print(f"Uniqueness summary calculation failed: {e}")
            raise  
        try:
            match_and_copy_sdf_files(self.config)
            print("Match and copy SDF files completed")
        except Exception as e:
            print(f"Match and copy SDF files failed: {e}")
            raise 
##################Usability#####################################################################
    def Usability(self):
        from utils.Usability import calculate_usability
        from utils.usability_summary import calculate_usability_summary

        try:
            calculate_usability(self.config)
            print("Usability calculation completed")
        except Exception as e:
            print(f"Usability calculation failed: {e}")
            raise 
        try:
            calculate_usability_summary(self.config)
            print("Usability summary completed")
        except Exception as e:
            print(f"Usability summary failed: {e}")
            raise  
##################QED########################################################################
    def QED(self):
        from utils.QED import calculate_QED
        from utils.QED_merge import calculate_QED_merge
        try:
            calculate_QED(self.config)
            print("QED calculation completed")
        except Exception as e:
            print(f"QED calculation failed: {e}")
            raise
        try:
            calculate_QED_merge(self.config)
            print("QED merge completed")
        except Exception as e:
            print(f"QED merge failed: {e}")
            raise 
##################SA score###################################################################
    def SA(self):
        from utils.SA import calculate_SA
        from utils.SA_merge import calculate_SA_merge

        try:
            calculate_SA(self.config)
            print("SA_score completed")
        except Exception as e:
            print(f"SA_score failed: {e}")
            raise  

        try:
            calculate_SA_merge(self.config)
            print("SA_merge completed")
        except Exception as e:
            print(f"SA_merge failed: {e}")
            raise 
##################Tdiv##################################################################
    def Tdiv(self):
        from utils.Morgan_Tdiv import calculate_Morgan_fingerprint_Tdiv
        from utils.Scaffold_Tdiv import calculate_Scaffold_Tdiv
        from utils.Tdiv_merge import calculate_tdiv_merge
        try:
            calculate_Morgan_fingerprint_Tdiv(self.config)
            print("Morgan fingerprint Tdiv calculation completed")
        except Exception as e:
            print(f"Morgan fingerprint Tdiv calculation failed: {e}")
            raise  
        try:
            calculate_Scaffold_Tdiv(self.config)
            print("Scaffold Tdiv calculation completed")
        except Exception as e:
            print(f"Scaffold Tdiv calculation failed: {e}")
            raise  
        try:
            calculate_tdiv_merge(self.config)
            print("Tdiv merge calculation completed")
        except Exception as e:
            print(f"Tdiv merge calculation failed: {e}")
            raise 

    def General_molecular_quality_metrics(self):     
        from utils.GMQMS import calculate_Summary
        try:
            calculate_Summary(self.config)
            #print("General_molecular_quality_metrics completed")
        except Exception as e:
            #print(f"General_molecular_quality_metrics failed: {e}")
            raise  
##################structural_properties_2D#####################################################
    def structural_2D(self):
        try:
            from utils.Structural_properties_2D import calculate_structural_properties_2D
            from utils.Structural_properties_2D_merge import calculate_structural_properties_2D_merge
            try:
                calculate_structural_properties_2D(self.config)
                print("structural_properties_2D calculation completed")
            except Exception as e:
                print(f"structural_properties_2D calculation failed: {e}")
                raise  
            try:
                calculate_structural_properties_2D_merge(self.config)
                print("structural_properties_2D_merge calculation completed")
            except Exception as e:
                print(f"structural_properties_2D_merge calculation failed: {e}")
                raise  
        except Exception as e:
            print(f"An error occurred in the structural_2D process: {e}")
            raise  
##################structural_properties_3D########################
    def structural_3D(self):
        if self.config['evaluate']['structural_3D']:
            minimization_force_field = self.config['settings']['Structural_3D_settings']['minimization_force_field']
        try:
            # Force field optimization
            if minimization_force_field == "OPLS3":
                from utils.OPLS3_preparation import calculate_OPLS3_preparation
                from utils.delete_additional_files import delete_additional_files
                try:
                    calculate_OPLS3_preparation(self.config)
                    print("OPLS3_preparation completed")
                except Exception as e:
                    print(f"OPLS3_preparation failed: {e}")
                    raise

                try:
                    delete_additional_files(self.config)
                    print("delete_additional_files completed")
                except Exception as e:
                    print(f"delete_additional_files failed: {e}")
                    raise

            elif minimization_force_field == "MMFF":
                from utils.MMFF_preparation import calculate_MMFF_preparation
                try:
                    calculate_MMFF_preparation(self.config)
                    print("MMFF_force_field success")
                except Exception as e:
                    print(f"MMFF_force_field failed: {e}")
                    raise
            else:
                raise ValueError(f"Unknown force field: {minimization_force_field}")

            # Bond length
            from utils.Bond_length import calculate_bond_length
            try:
                calculate_bond_length(self.config)
                print("Bond length calculation completed")
            except Exception as e:
                print(f"Bond length failed: {e}")
                raise

            # Bond angles
            from utils.Bond_angles import calculate_bond_angles
            try:
                calculate_bond_angles(self.config)
                print("Bond angles calculation completed")
            except Exception as e:
                print(f"Bond angles failed: {e}")
                raise

            # Dihedral angles
            from utils.Dihedral_angles import calculate_dihedral_angles
            try:
                calculate_dihedral_angles(self.config)
                print("Dihedral angles calculation completed")
            except Exception as e:
                print(f"Dihedral angles failed: {e}")
                raise

            # JSD
            from utils.JSD_bond_length import calculate_bljs_divergence
            from utils.JSD_bond_angles import calculate_bajs_divergence
            from utils.JSD_dihedral_angles import calculate_dajs_divergence
            try:
                calculate_bljs_divergence(self.config)
                print("Bond_length_JSD completed")
            except Exception as e:
                print(f"Bond_length_JSD failed: {e}")
                raise

            try:
                calculate_bajs_divergence(self.config)
                print("Bond_angles_JSD completed")
            except Exception as e:
                print(f"Bond_angles_JSD failed: {e}")
                raise

            try:
                calculate_dajs_divergence(self.config)
                print("Dihedral_angles_JSD completed")
            except Exception as e:
                print(f"Dihedral_angles_JSD failed: {e}")
                raise

        except Exception as e:
            print(f"{minimization_force_field}_force_field_failed: {e}")
            raise
        else:
            print("structural properties 3D evaluation not enabled")
    
    def structural_summary_2D3D(self):     
        from utils.structural_summary import calculate_summary
        try:
            calculate_summary(self.config)
            #print("structural_summary completed")
        except Exception as e:
            #print(f"structural_summary failed: {e}")
            raise
##################Recovery###########################################################################      
    def active_recovery(self):
        from utils.Molecule_recovery import calculate_molecule_recovery
        from utils.Scaffold_recovery import calculate_scafold_recovery
        try:
            calculate_molecule_recovery(self.config)
            calculate_scafold_recovery(self.config)
            print("active_recovery calculation completed")
        except Exception as e:
            print(f"active_recovery failed: {e}")
        from utils.recovery_summary import calculate_recovery
        try:
            calculate_recovery(self.config)
            #print("recovery_summary calculation completed")
        except Exception as e:
            pass

##################Redocking######################################################################
    def Redocking_score(self):
        try:
            if self.config['evaluate']['Redocking_score']:
                docking_method = self.config['settings']['Redocking_settings']['docking_method']
                ligprep_needed = self.config['settings']['Redocking_settings']['ligprep']

                if docking_method == "glide":
                    from utils.Ligprep_glide import glide_ligprep_process
                    from utils.delete_additional_files import delete_additional_files
                    from utils.Glide_Docking import glide_redocking_main
                    from utils.Glide_docking_score import calculate_docking_score

                    if ligprep_needed:
                        try:
                            glide_ligprep_process(self.config)
                            print("glide_ligprep_process completed")
                        except Exception as e:
                            print(f"glide_ligprep_process failed: {e}")
                            return
                    try:
                        delete_additional_files(self.config)
                    except Exception as e:
                        print(f"delete_additional_files failed: {e}")
                        return
                    try:
                        glide_redocking_main(self.config)
                        print("glide_redocking_main completed")
                    except Exception as e:
                        print(f"glide_redocking_main failed: {e}")
                        return
                    try:
                        delete_additional_files(self.config)
                    except Exception as e:
                        print(f"delete_additional_files failed: {e}")
                        return
                    try:
                        calculate_docking_score(self.config)
                        print("docking_score completed")
                    except Exception as e:
                        print(f"docking_score failed: {e}")
                        return
                elif docking_method == "vina":
                    from utils.Ligprep_vina import vina_ligprep_process
                    from utils.Vina_Docking import vina_redocking_main   
                    from utils.Vina_docking_score import calculate_docking_score 

                    if ligprep_needed:
                        try:
                            vina_ligprep_process(self.config)
                            print("vina_ligprep_process completed")
                        except Exception as e:
                            print(f"vina_ligprep_process failed: {e}")
                            return
                    try:
                        vina_redocking_main(self.config)
                        print("vina_redocking_main completed")
                    except Exception as e:
                        print(f"vina_redocking_main failed: {e}")
                        return
                    try:
                        calculate_docking_score(self.config)
                        print("docking_score completed")
                    except Exception as e:
                        print(f"docking_score failed: {e}")
                        return
            else:
                print("Not meeting the execution conditions, do not perform any operation")
        except ImportError as e:
            print(f"Redocking_score failed: {e}")
        except Exception as e:
            print(f"Redocking_score failed: {e}")

##################In_situ_Docking######################################################################
#Attention-please:SDF files need to H atoms added in advance
    def In_situ_Docking(self):
        from utils.add_h import calculate_add_h_main
        if self.config['evaluate']['In_situ_Docking']:
            docking_method = self.config['settings']['In_situ_Docking_settings']['docking_method']
            try:
                 #Perform hydrogenation atom first
                calculate_add_h_main(self.config)
                print("add_h completed")
                #Select the method to execute based on Method
                if docking_method == "glide":
                    from utils.Glide_In_situ_docking import glide_in_situ_docking_main
                    from utils.Glide_In_situ_docking_score import calculate_in_situ_docking_score
                    try:
                        glide_in_situ_docking_main(self.config)
                        print("glide_in_situ_docking_main completed")
                    except Exception as e:
                        print(f"glide_in_situ_docking_main failed: {e}")
                        return  # Stop the program if glide_in_situ_docking_main fails
                    try:
                        calculate_in_situ_docking_score(self.config)
                        print("calculate_in_situ_docking_score completed")
                    except Exception as e:
                        print(f"calculate_in_situ_docking_score failed: {e}")
                        return  # Stop the program if calculate_in_situ_docking_score fails
                elif docking_method == "vina":
                    from utils.Ligprep_vina_in_situ import vina_insitu_ligprep
                    from utils.Vina_In_situ_docking import vina_in_situ_docking_main
                    from utils.Vina_In_situ_docking_score import calculate_score
                    try:
                        vina_insitu_ligprep(self.config)
                        print("ligprep_vina_in_situ completed")
                    except Exception as e:
                        print(f"ligprep_vina_in_situ failed: {e}")
                        return  # Stop the program if ligprep_vina_in_situ fails
                    try:
                        vina_in_situ_docking_main(self.config)
                        print("vina_in_situ_docking_main completed")
                    except Exception as e:
                        print(f"vina_in_situ_docking_main failed: {e}")
                        return  # Stop the program if vina_in_situ_docking_main fails
                    try:
                        calculate_score(self.config)
                        print("vina_in_situ_score_summary completed")
                    except Exception as e:
                        print(f"vina_in_situ_score_summary failed: {e}")
                else:
                    raise ValueError(f"Unsupported docking method: {docking_method}")
            except Exception as e:

                return  # Stop the program if any other step fails
##################RMSD######################################################################
    def RMSD(self):
        try:
            redocking_method = self.config['settings']['RMSD_settings']['redocking_method']
            rmsd_method = self.config['settings']['RMSD_settings']['rmsd_method']

            if redocking_method == "glide":
                if rmsd_method == "openeye":
                    
                    from utils.Glide_Openeye_RMSD import calculate_openeye_rmsd
                    try:
                        calculate_openeye_rmsd(self.config)
                        print("Glide-Openeye-RMSD completed")
                    except Exception as e:
                        print(f"Glide-Openeye-RMSD failed: {e}")
                        raise  
                elif rmsd_method == "rdkit":
                    from utils.Glide_RDKit_RMSD import calculate_rdkit_rmsd
                    try:
                        calculate_rdkit_rmsd(self.config)
                        print("Glide-RDKit-RMSD completed")
                    except Exception as e:
                        print(f"Glide-RDKit-RMSD failed: {e}")
                        raise  
                else:
                    raise ValueError(f"Unsupported RMSD method: {rmsd_method}")

            elif redocking_method == "vina":
                if rmsd_method == "openeye":
                    from utils.Pdbqt_Sdf import convert_pdbqt_sdf
                    try:
                        convert_pdbqt_sdf(self.config)
                        print("Pdbqt-Sdf completed")
                    except Exception as e:
                        print(f"Pdbqt-Sdf failed: {e}")
                        raise  

                    from utils.Vina_Openeye_RMSD import calculate_rmsd
                    try:
                        calculate_rmsd(self.config)
                        print("Vina-Openeye-RMSD completed")
                    except Exception as e:
                        print(f"Vina-Openeye-RMSD failed: {e}")
                        raise  

                elif rmsd_method == "rdkit":
                    from utils.Pdbqt_Sdf import convert_pdbqt_sdf
                    try:
                        convert_pdbqt_sdf(self.config)
                        print("Pdbqt-Sdf completed")
                    except Exception as e:
                        print(f"Pdbqt-Sdf failed: {e}")
                        raise  

                    from utils.Vina_RDKit_RMSD import calculate_rmsd
                    try:
                        calculate_rmsd(self.config)
                        print("Vina-RDKit-RMSD completed")
                    except Exception as e:
                        print(f"Vina-RDKit-RMSD failed: {e}")
                        raise  

                else:
                    raise ValueError(f"Unsupported RMSD method: {rmsd_method}")

            else:
                raise ValueError(f"Unsupported redocking method: {redocking_method}")

            from utils.rmsd_score import rmsd_score_summary
            try:
                rmsd_score_summary(self.config)
                print("RMSD-Score-Extract completed")
            except Exception as e:
                print(f"RMSD-Score-Extract failed: {e}")
                raise 
        except Exception as e:
            print(f"RMSD calculation process failed: {e}")
            raise  

    def Target_binding_metrics(self):    
        from utils.Target_binding_metrics import calculate_metrics
        try:
            calculate_metrics(self.config)
            #print("Target-binding-metrics calculation completed")
        except Exception as e:
            #print(f"Target-binding-metrics calculation failed: {e}")
            raise  
################################################################################################
    def process(self):

        # Path checking
        POKMOL3D_env = config['POKMOL-3D_env']
        if not os.path.exists(POKMOL3D_env):
           print("POKMOL-3D_env setting error!!!")
           return None

        input_path = config['data']['input_path']
        if not os.path.exists(input_path):
           print("input path does not exist!!!")
           return None

        # step1:Model-quality
        if config['evaluate']['Target_failure_rate']:
            self.Target_failure_rate()
        if config['evaluate']['Sampling_success_rate']:
            self.Sampling_success_rate()    

        # step2:Validity
        if config['evaluate']['Validity']:
            self.Validity()

        # step2:Uniqueness
        if config['evaluate']['Uniqueness']:
            self.Uniqueness()

        # step4:Usability
        if config['evaluate']['Usability']:
            self.Usability()
          
        # step5:QED
        if config['evaluate']['QED']:
            self.QED()

        # step6:SA
        if config['evaluate']['SA']:
            self.SA()

        # step7:Tdiv
        if config['evaluate']['Tdiv']:
            self.Tdiv()

        self.General_molecular_quality_metrics()

        # step8:structural_2D
        if config['evaluate']['structural_2D']:
            self.structural_2D()

        # step9:structural_3D
        if config['evaluate']['structural_3D']:
            self.structural_3D()

        self.structural_summary_2D3D()

        # step10:active_recovery
        if config['evaluate']['active_recovery']:
            self.active_recovery()
        
        # step11:Redocking_score
        if config['evaluate']['Redocking_score']:
             self.Redocking_score()

        # step12:In_situ_Docking
        if config['evaluate']['In_situ_Docking']:
             self.In_situ_Docking()

        # step13:RMSD
        if config['evaluate']['RMSD']:
            self.RMSD()
        self.Target_binding_metrics()
        
if __name__ == "__main__":
    # Main program entrance
    # Create a command-line parameter parser
    p = argparse.ArgumentParser()
    # Add command line parameter --config
    p.add_argument('--config', type=str, required=True)
    # Parse command line parameters
    args = p.parse_args()
    # Open and read the configuration file
    with open(args.config, 'r') as f:
        config = yaml.safe_load(f)
    config = EasyDict(config)
    processor = MainProcessor(config)
    try:
        processor.process()
    except FileNotFoundError as e:
        print(e)

