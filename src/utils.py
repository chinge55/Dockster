import prody
import os
from fastapi import HTTPException
import pdbfixer
from openmm import app
import openmm
from openmm import unit
import os
from Bio.PDB import PDBParser
import numpy as np

import qrotate
import subprocess

temporary_file_dir = "temp"
if not os.path.exists(temporary_file_dir):
    os.makedirs(temporary_file_dir)


def convert_bytes_prody(bytes_pdb):
    """This function takes a bytes, i.e uploaded pdb file and converts it to prody AtomGroup.
    Currently, doing this by creating temporary file, Will do better in the future
    TODO:
    """
    with open(f"{temporary_file_dir}/protein.pdb", "wb") as f:
        f.write(bytes_pdb)
    pdb = prody.parsePDB(f"{temporary_file_dir}/protein.pdb")
    os.remove(f"{temporary_file_dir}/protein.pdb")
    if pdb is None:
        raise HTTPException(status_code=404, detail="problem in PDB File")
    return pdb


def set_coordinates(protein, ligand):
    # Get the coordinates of protein and ligand
    prot_xyz = protein.getCoords()
    lig_xyz = ligand.getCoords()

    # Centre both on the COM of the ligand:
    prot_xyz -= lig_xyz.mean(0)
    lig_xyz -= lig_xyz.mean(0)

    al = qrotate.Align()
    angles = al.align_pcl(lig_xyz, get_angles=True)

    # Apply the rotation matrices

    prot_xyz = prot_xyz.dot(angles[0]).dot(angles[1])
    lig_xyz = lig_xyz.dot(angles[0]).dot(angles[1])

    protein.setCoords(prot_xyz)
    ligand.setCoords(lig_xyz)

    prody.writePDB(f"{temporary_file_dir}/ligand.pdb", ligand)
    prody.writePDB(f"{temporary_file_dir}/protein.pdb", protein)

    fixer = pdbfixer.PDBFixer(f"{temporary_file_dir}/protein.pdb")

    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)
    return fixer


def minimize_energy(fixer):
    forcefield = app.ForceField("amber14-all.xml", "amber14/tip3pfb.xml")
    system = forcefield.createSystem(
        fixer.topology,
        nonbondedMethod=app.CutoffNonPeriodic,
        nonbondedCutoff=0.9 * unit.nanometer,
    )
    atom_elements = [atom.element.name for atom in fixer.topology.atoms()]
    for i in range(system.getNumParticles()):
        if atom_elements[i] != "hydrogen":
            system.setParticleMass(i, 0.0)

    integrator = openmm.LangevinIntegrator(
        298 * unit.kelvin, 1 / unit.picosecond, 1 * unit.femtosecond
    )
    platform = openmm.Platform.getPlatformByName("CPU")

    simulation = app.Simulation(fixer.topology, system, integrator, platform)
    simulation.context.setPositions(fixer.positions)
    simulation.minimizeEnergy()

    positions = simulation.context.getState(getPositions=True).getPositions()

    app.PDBFile.writeFile(
        fixer.topology, positions, open(f"{temporary_file_dir}/proteinH.pdb", "w")
    )

    return True


def get_center_of_mass(pdb_file, chain_id):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)
    model = structure[0]  # assuming only one model in PDB file
    chain = model[chain_id]

    coords = np.array([atom.get_coord() for atom in chain.get_atoms()])
    center = np.mean(coords, axis=0)

    return center


def get_box_size(pdb_file, chain_id, buffer=10.0):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)
    model = structure[0]
    chain = model[chain_id]

    coords = np.array([atom.get_coord() for atom in chain.get_atoms()])
    center = np.mean(coords, axis=0)
    distances = np.linalg.norm(coords - center, axis=1)
    max_distance = np.max(distances)

    return max_distance + buffer


def convert_molecules():
    # Can't risk having variables here
    ligand_command = "docker run -v `pwd`/temp:/home/obabel --rm myobabel obabel ligand.pdb -O ligand.pdbqt"
    protein_command = "docker run -v `pwd`/temp:/home/obabel --rm myobabel obabel proteinH.pdb -xr -O protein.pdbqt"
    os.system(ligand_command)
    os.system(protein_command)


def perform_docking(exhaustiveness=8, num_modes=10):
    center = get_center_of_mass(f"{temporary_file_dir}/proteinH.pdb", "A")
    box_size = get_box_size(f"{temporary_file_dir}/proteinH.pdb", "A")
    print(center, box_size)
    with open(f"{temporary_file_dir}/config.txt", "w") as f:
        f.write(f"receptor = protein.pdbqt")
        f.write("\n")
        f.write(f"ligand = ligand.pdbqt")
        f.write("\n")
        f.write("\n")
        f.write(f"center_x = {center[0]}")
        f.write("\n")
        f.write(f"center_y = {center[1]}")
        f.write("\n")
        f.write(f"center_z = {center[2]}")
        f.write("\n")
        f.write("\n")
        f.write("\n")
        f.write(f"size_x = {int(box_size)}")
        f.write("\n")
        f.write(f"size_y = {int(box_size)}")
        f.write("\n")
        f.write(f"size_z = {int(box_size)}")
        f.write("\n")

        f.write(f"exhaustiveness = {exhaustiveness}")
        f.write("\n")
        f.write(f"num_modes = {num_modes}")
        f.write("\n")
    pwd = os.getcwd()
    result = subprocess.run(
        [
            "docker",
            "run",
            "-v",
            f"{pwd}/temp:/data",
            "--rm",
            "ghcr.io/metaphorme/vina:release",
            "vina",
            "--config",
            "config.txt",
        ],
        stdout=subprocess.PIPE,
    )
    return result.stdout
