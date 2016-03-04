#
#  Created by Greg Landrum (greg.landrum@t5informatics.com), Feb 2016
#
from flask import Flask,make_response,request,jsonify
from rdkit import Chem
from rdkit import rdBase
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from io import StringIO
import json,sys

Chem.WrapLogs()
app = Flask(__name__)

# error handline example from the Flask docs
# (http://flask.pocoo.org/docs/0.10/patterns/apierrors/)
class InvalidUsage(Exception):
    status_code = 400

    def __init__(self, message, status_code=None, payload=None):
        Exception.__init__(self)
        self.message = message
        if status_code is not None:
            self.status_code = status_code
        self.payload = payload

    def to_dict(self):
        rv = dict(self.payload or ())
        rv['message'] = self.message
        return rv
@app.errorhandler(InvalidUsage)
def handle_invalid_usage(error):
    response = jsonify(error.to_dict())
    response.status_code = error.status_code
    return response


@app.route('/')
def health():
    res = dict(rdkitVersion=rdBase.rdkitVersion,boostVersion=rdBase.boostVersion)
    return json.dumps(res)

def _molfromrequest():
    # get errors on stderr:
    sio = sys.stderr = StringIO()
    if 'smiles' in request.values:
        mol = Chem.MolFromSmiles(request.values.get('smiles'))
    elif 'mol' in request.values:
        mol = Chem.MolFromMolBlock(request.values.get('mol'))
    else:
        raise InvalidUsage("Neither 'smi' nor 'mol' present.", status_code=410)
    if mol is None:
        errm = sio.getvalue()
        errm = errm.replace('RDKit ERROR: \n','') # some errors leave blank lines
        raise InvalidUsage("Molecule could not be processed. Error message was:\n%s"%errm, status_code=411)
    return mol

@app.route('/canon_smiles', methods=['GET', 'POST'])
def canon_smiles():
    " returns canonical SMILES for input data "
    mol = _molfromrequest()
    return Chem.MolToSmiles(mol,True)

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

@app.route('/to_img/mol.png', methods=['GET', 'POST'])
def to_png():
    " returns a PNG for input data "
    mol = _molfromrequest()
    response = make_response(_render(mol,_moltopng))
    response.headers['Content-Type'] = 'image/png'
    return response
@app.route('/to_img/mol.svg', methods=['GET', 'POST'])
def to_svg():
    " returns an SVG for input data "
    mol = _molfromrequest()
    response = make_response(_render(mol,_moltosvg))
    return response


if __name__ == '__main__':
    # FIX: turn this off pre-deployment
    app.debug  = True
    app.run()
