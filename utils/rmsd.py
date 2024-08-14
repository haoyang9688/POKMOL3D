import sys
import os
import yaml

def main(argv=[__name__]):
    with open('config.yml', 'r') as file:
        config = yaml.safe_load(file)
  
    licence_path = config['settings']['RMSD_settings']['OE_LICENSE']
  
    os.environ['OE_LICENSE'] = licence_path

    from openeye import oechem

    itf = oechem.OEInterface(InterfaceData, argv)

    if not itf.GetBool("-verbose"):
        oechem.OEThrow.SetLevel(oechem.OEErrorLevel_Warning)

    ref_file = itf.GetString("-ref")
    input_file = itf.GetString("-in")

    automorph = itf.GetBool("-automorph")
    heavy = itf.GetBool("-heavyonly")

    overlay = False

    ifs = oechem.oemolistream()
    if not ifs.open(ref_file):
        oechem.OEThrow.Fatal("error" % ref_file)

    rmol = oechem.OEGraphMol()
    if not oechem.OEReadMolecule(ifs, rmol):
        oechem.OEThrow.Fatal("Unable to read reference molecule")

    ifs = oechem.oemolistream()
    if not ifs.open(input_file):
        oechem.OEThrow.Fatal("error" % input_file)

    
    results = []  

    for mol in ifs.GetOEMols():
        rmsds = oechem.OEDoubleArray(mol.GetMaxConfIdx())

     
        oechem.OERMSD(rmol, mol, rmsds, automorph, heavy, overlay)

        result_str = f"results\tfor\t{mol.GetTitle()}: " + "\t".join([f"RMSD-{i}:\t{rmsds[i]:.2f}" for i in range(len(rmsds))])
        results.append(result_str)

    
    if itf.HasString("-out"):
        output_file = itf.GetString("-out")
        with open(output_file, 'w') as txt_file:
            for result in results:
                txt_file.write(result + "\n")

    if itf.HasString("-sdfout"):
        sdf_output_file = itf.GetString("-sdfout")
        ofs = oechem.oemolostream()
        if not ofs.open(sdf_output_file):
            oechem.OEThrow.Fatal("error" % sdf_output_file)
        
        for mol in ifs.GetOEMols():
            oechem.OEWriteMolecule(ofs, mol)

    return 0

# InterfaceData 
InterfaceData = """\
!BRIEF [options] [-ref <mol file>] [-in <mol file>] [-sdfout <sdf file>] [-out <txt file>]

!CATEGORY "input/output options"

  !PARAMETER -ref
    !TYPE string
    !REQUIRED true
    !BRIEF 
    !KEYLESS 1
  !END

  !PARAMETER -in
    !ALIAS -i
    !TYPE string
    !REQUIRED true
    !BRIEF 
    !KEYLESS 2
  !END

  !PARAMETER -sdfout
    !ALIAS -so
    !TYPE string
    !REQUIRED false
    !BRIEF 
    !KEYLESS 3
  !END

  !PARAMETER -out
    !ALIAS -o
    !TYPE string
    !REQUIRED false
    !BRIEF 
    !KEYLESS 4
  !END

!END

!CATEGORY "options"

  !PARAMETER -automorph
    !TYPE bool
    !DEFAULT true
    !BRIEF 
    !DETAIL
  !END

  !PARAMETER -heavyonly
    !TYPE bool
    !DEFAULT true
    !BRIEF 
  !END

  !PARAMETER -verbose
    !ALIAS -v
    !TYPE bool
    !DEFAULT false
    !BRIEF 
  !END

!END
"""

if __name__ == "__main__":
    sys.exit(main(sys.argv))

