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
    res = dict(rdkitVersion=rdBase.rdkitVersion,boostVersion=rdBase.boostVersion,pythonVersion=sys.version)
    return json.dumps(res)

def _molfromrequest():
    # get errors on stderr:
    sio = sys.stderr = StringIO()
    sanitize = int(request.values.get('sanitize',1))
    removeHs = int(request.values.get('removeHs',1))
    if 'smiles' in request.values:
        mol = Chem.MolFromSmiles(request.values.get('smiles'),sanitize=sanitize)
    elif 'mol' in request.values:
        mol = Chem.MolFromMolBlock(request.values.get('mol'),sanitize=sanitize,removeHs=removeHs)
    else:
        raise InvalidUsage("Neither 'smiles' nor 'mol' present.", status_code=410)
    if mol is None:
        errm = sio.getvalue()
        errm = errm.replace('RDKit ERROR: \n','') # some errors leave blank lines
        raise InvalidUsage("Molecule could not be processed. Error message was:\n%s"%errm, status_code=411)
    if not sanitize:
        mol.UpdatePropertyCache(False)
        Chem.FastFindRings(mol)
        Chem.SetConjugation(mol)
        Chem.SetHybridization(mol)
    return mol

@app.route('/canon_smiles', methods=['GET', 'POST'])
def canon_smiles():
    " returns canonical SMILES for input data "
    mol = _molfromrequest()
    return Chem.MolToSmiles(mol,True)

def _paramToList(text):
    if not text:
        return None
    if text[-1] in ')]':
        text = text[:-1]
    if text[0] in '[(':
        text = text[0:]
    text = text.split(',')
    return [int(x) for x in text]


def _drawHelper(mol,kekulize,drawer,**kwargs):
    sanit = bool(int(request.values.get('sanitize',1)))
    if 'highlightAtoms' not in kwargs:
        tmp = _paramToList(request.values.get('highlightAtoms',None))
        if tmp:
            kwargs['highlightAtoms'] = tmp
    if 'highlightBonds' not in kwargs:
        tmp = _paramToList(request.values.get('highlightBonds',None))
        if tmp:
            kwargs['highlightBonds'] = tmp
    if 'legend' not in kwargs:
        kwargs['legend'] = request.values.get('legend','')
    mc = rdMolDraw2D.PrepareMolForDrawing(mol,kekulize=kekulize&sanit,
                                          addChiralHs=sanit)
    drawer.DrawMolecule(mc,**kwargs)
    drawer.FinishDrawing()

def _moltosvg(mol,molSize=(450,200),kekulize=True,drawer=None,**kwargs):
    if drawer is None:
        drawer = rdMolDraw2D.MolDraw2DSVG(molSize[0],molSize[1])
    _drawHelper(mol,kekulize,drawer,**kwargs)
    svg = drawer.GetDrawingText()
    return svg.replace('svg:','')

def _moltopng(mol,molSize=(450,200),kekulize=True,drawer=None,**kwargs):
    if drawer is None:
        drawer = rdMolDraw2D.MolDraw2DCairo(molSize[0],molSize[1])
    _drawHelper(mol,kekulize,drawer,**kwargs)
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

def _gen3d_sdf(mol,**kwargs):
    from rdkit.Chem import AllChem
    minimize = int(request.values.get('minimize',0))
    mh = Chem.AddHs(mol)
    mh.SetProp("_Name","2D to 3D output")
    cid=AllChem.EmbedMolecule(mh,enforceChirality=True,useExpTorsionAnglePrefs=True,
                     useBasicKnowledge=True)
    if cid<0:
        raise InvalidUsage("Molecule could not be embedded.", status_code=418)
    if minimize:
        AllChem.MMFFOptimizeMolecule(mh)
    res = Chem.MolToMolBlock(mh)
    return res

@app.route('/to_3d/mol.sdf', methods=['GET', 'POST'])
def to_3d_sdf():
    " returns a mol block for input data "
    mol = _molfromrequest()
    response = make_response(_gen3d_sdf(mol))
    return response


if __name__ == '__main__':
    # FIX: turn this off pre-deployment
    app.debug  = True
    app.run()
