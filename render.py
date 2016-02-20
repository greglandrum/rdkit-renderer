#
#  Created by Greg Landrum (greg.landrum@t5informatics.com), Feb 2016
#
from flask import Flask,make_response,request
from rdkit import Chem
from rdkit import rdBase
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
import json

app = Flask(__name__)

@app.route('/')
def health():
    res = dict(rdkitVersion=rdBase.rdkitVersion,boostVersion=rdBase.boostVersion)
    return json.dumps(res)

@app.route('/canon_smiles/<smiles>')
def canon_smiles(smiles):
    " returns canonical SMILES for an input SMILES "
    return Chem.CanonSmiles(smiles)

def _moltosvg(mol,molSize=(450,200),kekulize=True,drawer=None,**kwargs):
    mc = rdMolDraw2D.PrepareMolForDrawing(mol,kekulize=kekulize)
    if drawer is None:
        drawer = rdMolDraw2D.MolDraw2DSVG(molSize[0],molSize[1])
    drawer.DrawMolecule(mc,**kwargs)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    #return svg
    return svg.replace('svg:','')
def _moltopng(mol,molSize=(450,200),kekulize=True,drawer=None,**kwargs):
    mc = rdMolDraw2D.PrepareMolForDrawing(mol,kekulize=kekulize)
    if drawer is None:
        drawer = rdMolDraw2D.MolDraw2DCairo(molSize[0],molSize[1])
    drawer.DrawMolecule(mc,**kwargs)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()

def _render(mol,renderer,size=(150,100),**kwargs):
    sz = int(request.values.get('w',size[0])),int(request.values.get('h',size[1]))
    return renderer(mol,molSize=sz,**kwargs)

@app.route('/mol_to_png/mol.png', methods=['GET', 'POST'])
def mol_to_png():
    " returns a PNG for a mol block "
    molblock = request.values.get('mol')
    m = Chem.MolFromMolBlock(molblock)
    response = make_response(_render(m,_moltopng))
    response.headers['Content-Type'] = 'image/png'
    return response
@app.route('/mol_to_svg/mol.svg', methods=['GET', 'POST'])
def mol_to_svg():
    " returns an SVG for a mol block "
    molblock = request.values.get('mol')
    m = Chem.MolFromMolBlock(molblock)
    response = make_response(_render(m,_moltosvg))
    return response

@app.route('/smiles_to_png/<smiles>.png', methods=['GET', 'POST'])
def smiles_to_png(smiles):
    " returns a PNG for a SMILES "
    m = Chem.MolFromSmiles(smiles)
    response = make_response(_render(m,_moltopng))
    response.headers['Content-Type'] = 'image/png'
    return response

@app.route('/smiles_to_svg/<smiles>.svg', methods=['GET', 'POST'])
def smiles_to_svg(smiles):
    " returns an SVG for a SMILES "
    m = Chem.MolFromSmiles(smiles)
    response = make_response(_render(m,_moltosvg))
    #response.headers['Content-Type'] = 'image/svg+xml' # SVG doesn't render properly in chrome if we send this
    return response


if __name__ == '__main__':
    # FIX: turn this off pre-deployment
    app.debug  = True
    app.run()
