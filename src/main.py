from fastapi import FastAPI, HTTPException, File, UploadFile
import prody
from typing_extensions import Annotated
from typing import Union
import tempfile
import os
from utils import (
    convert_bytes_prody,
    set_coordinates,
    minimize_energy,
    convert_molecules,
    perform_docking,
)


app = FastAPI()


@app.get("/")
async def root():
    return {"message": "Hello World"}


@app.post("/enterPDBID")
async def enter_pdbid(pdb_id: str):
    file_name = prody.fetchPDB(pdb_id, compressed=False)
    if file_name is None:
        raise HTTPException(status_code=404, detail="PDB ID not found")
    return {"fileName": file_name}


@app.post("/uploadPDBFiles")
async def upload_ligand_pdb(
    protein_pdb_bytes: Annotated[Union[bytes, None], File()],
    ligand_pdb_bytes: Annotated[Union[bytes, None], File()],
):
    protein_pdb = convert_bytes_prody(protein_pdb_bytes)
    protein = protein_pdb.select("protein")
    if protein is None:
        protein = protein_pdb
    ligand = convert_bytes_prody(ligand_pdb_bytes)
    fixer = set_coordinates(protein, ligand)
    sim_status = minimize_energy(fixer)
    if not sim_status:
        raise HTTPException(status_code=400, detail="Some error occured")
    convert_molecules()
    result = perform_docking()
    return result
