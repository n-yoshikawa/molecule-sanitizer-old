import os
import json
import uuid
import threading
import pandas as pd
import subprocess
from subprocess import PIPE

from flask import Flask, render_template, request, send_file, jsonify

from rdkit import Chem
from rdkit.Chem import Draw, AllChem, rdCoordGen
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')


app = Flask(__name__)
app.secret_key = 'hogefuga'

@app.route('/')
def index():
    return render_template("index.html")

@app.route('/viewer/<list_id>')
def viewer(list_id):
    lst = []
    if not os.path.isfile(f'{list_id}.sub'):
        return render_template("viewer.html", is_ready=False, list_id=list_id)
    with open(f'{list_id}.sub') as f:
        for line in f:
            ID,description,nodes,edges,s_abs,s_rel,c_abs,c_rel = line.split(',')
            mol = Chem.MolFromSmiles(description)
            if mol is not None:
                smiles = Chem.MolToSmiles(mol)
                rdCoordGen.AddCoords(mol)
                d2d = rdMolDraw2D.MolDraw2DSVG(180, 180)
                d2d.DrawMolecule(mol)
                d2d.FinishDrawing()
                svg = d2d.GetDrawingText()
                lst.append([int(ID),description,int(nodes),int(edges),int(s_abs),float(s_rel),int(c_abs),float(c_rel), smiles, svg])
    df = pd.DataFrame(lst, columns = ['ID', 'description', 'nodes', 'edges', 's_abs', 's_rel', 'c_abs', 'c_rel','smiles', 'svg'])
    df = df.set_index('ID')
    df_s = df.sort_values('s_abs', ascending=False)
    substructures = []
    for index, row in df_s.iterrows():
        substructures.append({
            "ID": index,
            "SMILES": row['smiles'],
            "SVG": row['svg'] })
        if len(substructures) > 100:
            break
    return render_template("viewer.html", is_ready=True, substructures=substructures, list_id=list_id)

@app.route('/select', methods=['POST'])
def select():
    data = request.get_json()
    if data['ID'] == 'all':
        send_all = True
    else:
        send_all = False
        ssID = int(data['ID'])
    list_id = data['list_id']
    if not os.path.isfile(f'{list_id}.in'):
        return jsonify(values=json.dumps({"molecules": [], "exists": False}))
    smis = []
    with open(f'{list_id}.in') as f:
        for line in f:
            _, _, smi = line.rstrip('\n').split(',')
            smis.append(smi)
    molecules = []
    if send_all or not (os.path.isfile(f'{list_id}.sub') and os.path.isfile(f'{list_id}.ids')):
        for idx, smi in enumerate(smis):
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            AllChem.Compute2DCoords(mol)
            d2d = rdMolDraw2D.MolDraw2DSVG(280, 210)
            d2d.DrawMolecule(mol)
            d2d.FinishDrawing()
            svg = d2d.GetDrawingText()
            molecules.append({
                "smiles": Chem.MolToSmiles(mol),
                "svg_html": svg,
                "name": f"checkbox{idx}" })
    else:
        with open(f'{list_id}.sub') as f:
            next(f)
            for line in f:
                ID,description,nodes,edges,s_abs,s_rel,c_abs,c_rel = line.split(',')
                if int(ID) == ssID:
                    template_smi = description
                    break
        with open(f'{list_id}.ids') as f:
            next(f)
            for line in f:
                ID, id_list = line.split(':')
                if int(ID) == ssID:
                    ids = id_list.split(',')
                    ss_mols = [Chem.MolFromSmiles(smis[int(i)]) for i in ids]
                    break
        template = Chem.MolFromSmiles(template_smi)
        AllChem.Compute2DCoords(template)
        for idx, mol in enumerate(ss_mols):
            AllChem.Compute2DCoords(mol)
            #AllChem.GenerateDepictionMatching2DStructure(mol, template)
            highlights = mol.GetSubstructMatch(template)
            d2d = rdMolDraw2D.MolDraw2DSVG(280, 210)
            d2d.DrawMolecule(mol, highlightAtoms=highlights)
            d2d.FinishDrawing()
            svg = d2d.GetDrawingText()
            molecules.append({
                "smiles": Chem.MolToSmiles(mol),
                "svg_html": svg,
                "name": f"checkbox{idx}" })
    return jsonify(values=json.dumps({"molecules": molecules, "exists": True}))

@app.route('/send', methods=['POST'])
def send():
    data = request.get_json()
    list_id = uuid.uuid4().hex
    with open(f'{list_id}.smi', 'w') as f:
        f.write(data['smiles_list'])
    thread = threading.Thread(target=analyze_ss, args=(list_id,))
    thread.start()
    return jsonify(values=json.dumps({"url": f"/viewer/{list_id}"}))

def analyze_ss(list_id):
    smis = []
    with open(f'{list_id}.smi') as f:
        for line in f:
            smi = line.split()[0]
            if Chem.MolFromSmiles(smi) is not None:
                smis.append(smi)
    with open(f'{list_id}.in', 'w') as f:
        for i, smi in enumerate(smis):
            f.write(f'{i},0,{smi}\n')
    proc = subprocess.run(f"java -cp moss.jar moss.Miner -s0.1 -S0.5 -m12 -R -B {list_id}.in {list_id}.sub {list_id}.ids >> log.txt 2>&1", shell=True, text=True)

if __name__ == '__main__':
    app.run(debug=True)
    # For deployment
    # app.run(debug=False, host='0.0.0.0', port=12345)
