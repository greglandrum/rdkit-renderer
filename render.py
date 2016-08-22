#
#  Created by Greg Landrum (greg.landrum@t5informatics.com), Feb 2016
#
from flask import Flask,make_response,request,jsonify
try:
    from flasgger import Swagger
except ImportError:
    Swagger = None

from rdkit import Chem
from rdkit import rdBase
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from io import StringIO
import json,sys

Chem.WrapLogs()
app = Flask(__name__)
if Swagger is not None:
    Swagger(app)

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
    '''
    Simple health check call
    ---
    responses:
      500:
        description: Error!
      200:
        description: Everything is fine
        schema:
          id: what?
          properties:
            rdkitVersion:
              type: string
              description: version of the RDKit being used
            boostVersion:
              type: string
              description: The version of boost being used
            pythonVersion:
              type: string
              description: The version of python being used
    '''
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
        print(request.values.get('mol'))
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
def _queryfromrequest(suffix='_query'):
    # get errors on stderr:
    sio = sys.stderr = StringIO()
    if 'smiles'+suffix in request.values:
        mol = Chem.MolFromSmiles(request.values.get('smiles'+suffix),sanitize=False)
        if mol is not None:
            try:
                Chem.SanitizeMol(mol)
            except:
                mol = None
    elif 'smarts'+suffix in request.values:
        mol = Chem.MolFromSmarts(request.values.get('smarts'+suffix))
    elif 'mol'+suffix in request.values:
        mol = Chem.MolFromMolBlock(request.values.get('mol'+suffix),removeHs=False)
        mol = Chem.AdjustQueryProperties(mol)
    else:
        return None
    if mol is None:
        errm = sio.getvalue()
        errm = errm.replace('RDKit ERROR: \n','') # some errors leave blank lines
        raise InvalidUsage("Molecule could not be processed. Error message was:\n%s"%errm, status_code=411)
    return mol

@app.route('/canon_smiles', methods=['GET', 'POST'])
def canon_smiles():
    '''
    Returns the canonical SMILES for a molecule
    ---
    parameters:
      - name: smiles
        in: query
        type: string
        required: false
        description: input SMILES. Provide either this or mol
      - name: mol
        in: query
        type: string
        required: false
        description: input CTAB. Provide either this or smiles

    responses:
      500:
        description: Error!
      410:
        description: no molecule data provided
      411:
        description: input molecule could not be processed
      200:
        description: Everything is fine
        schema:
          id: what?
          properties:
            smiles:
              type: string
              description: canonical SMILES
    '''
    mol = _molfromrequest()
    return json.dumps(dict(smiles=Chem.MolToSmiles(mol,True)))

import json
def _loadJSONParam(text):
    if not text:
        return None
    return json.loads(text)

def _drawHelper(mol,kekulize,drawer,**kwargs):
    sanit = bool(int(request.values.get('sanitize',1)))
    for arg in ('highlightAtomColors', 'highlightBondColors'):
        if arg not in kwargs:
            tmp = _loadJSONParam(request.values.get(arg,None))
            if tmp:
                for k in list(tmp):
                    tmp[int(k)] = tuple(tmp[k])
                    del tmp[k]
                kwargs[arg] = tmp
    for arg in ('highlightAtoms','highlightBonds'):

        if arg not in kwargs:
            tmp = _loadJSONParam(request.values.get(arg,None))
            if tmp:
                kwargs[arg] = tmp
    if 'legend' not in kwargs:
        kwargs['legend'] = request.values.get('legend','')
    if 'highlightSubstruct' not in kwargs:
        if 'highlightColor' not in kwargs:
            hcolor = _loadJSONParam(request.values.get('highlightColor',None))
            if hcolor:
                highlightColor = tuple(hcolor)
            else:
                highlightColor = None
        else:
            highlightColor = kwargs['higlightColor']
            del kwargs['higlightColor']
        if bool(int(request.values.get('highlightSubstruct',0))):
            qmol = _queryfromrequest(suffix='_highlight')
            if qmol is not None:
                qMatches = mol.GetSubstructMatches(qmol)
                highlights=kwargs.get('highlightAtoms',[])
                if highlightColor is not None:
                    highlightBonds=kwargs.get('highlightBonds',[])
                else:
                    highlightBonds = None

                for entry in qMatches:
                    if highlightColor is not None:
                        had = kwargs.get('highlightAtomColors',{})
                        for v in entry:
                            had[v] = highlightColor
                        kwargs['highlightAtomColors'] = had

                        hbd = kwargs.get('highlightBondColors',{})
                        emap = dict(enumerate(entry))
                        for bond in qmol.GetBonds():
                            mbond = mol.GetBondBetweenAtoms(emap[bond.GetBeginAtomIdx()],emap[bond.GetEndAtomIdx()])
                            highlightBonds.append(mbond.GetIdx())
                            hbd[mbond.GetIdx()] = highlightColor
                        kwargs['highlightBondColors'] = hbd
                    highlights.extend(list(entry))
                if highlights:
                    kwargs['highlightAtoms']=highlights
                    if highlightBonds is not None:
                        kwargs['highlightBonds']=highlightBonds
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
    """
    Returns a PNG rendering of a molecule
    ---
    parameters:
      - name: smiles
        in: query
        type: string
        required: false
        description: input SMILES. Provide either this or mol
      - name: mol
        in: query
        type: string
        required: false
        description: input CTAB. Provide either this or smiles
      - name: w
        in: query
        type: integer
        required: false
        default: 150
        description: image width
      - name: h
        in: query
        type: integer
        required: false
        default: 100
        description: image height
      - name: legend
        in: query
        type: string
        required: false
        description: text to be displayed beneath the molecule
      - name: highlightAtoms
        in: query
        type: array
        items:
            type: integer
            description: indices (0-based) of the atoms to highlight
        required: false
        description: atoms to be highlighted (as a JSON array)
      - name: highlightBonds
        in: query
        type: array
        items:
            type: integer
            description: indices (0-based) of the bonds to highlight
        required: false
        description: bonds to be highlighted (as a JSON array)
    produces:
        image/png
    responses:
      500:
        description: Error!
      410:
        description: no molecule data provided
      411:
        description: input molecule could not be processed
      200:
        description: Everything is fine
    """
    mol = _molfromrequest()
    response = make_response(_render(mol,_moltopng))
    response.headers['Content-Type'] = 'image/png'
    return response
@app.route('/to_img/mol.svg', methods=['GET', 'POST'])
def to_svg():
    """
    Returns a SVG rendering of a molecule
    ---
    parameters:
      - name: smiles
        in: query
        type: string
        required: false
        description: input SMILES. Provide either this or mol
      - name: mol
        in: query
        type: string
        required: false
        description: input CTAB. Provide either this or smiles
      - name: w
        in: query
        type: integer
        required: false
        default: 150
        description: image width
      - name: h
        in: query
        type: integer
        required: false
        default: 100
        description: image height
      - name: legend
        in: query
        type: string
        required: false
        description: text to be displayed beneath the molecule
      - name: highlightAtoms
        in: query
        type: array
        items:
            type: integer
            description: indices (0-based) of the atoms to highlight
        required: false
        description: atoms to be highlighted (as a JSON array)
      - name: highlightBonds
        in: query
        type: array
        items:
            type: integer
            description: indices (0-based) of the bonds to highlight
        required: false
        description: bonds to be highlighted (as a JSON array)
    produces:
        image/svg+xml
    responses:
      500:
        description: Error!
      410:
        description: no molecule data provided
      411:
        description: input molecule could not be processed
      200:
        description: Everything is fine
    """
    mol = _molfromrequest()
    response = make_response(_render(mol,_moltosvg))
    #response.headers['Content-Type'] = 'image/svg+xml'
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
