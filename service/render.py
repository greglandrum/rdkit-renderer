#
#  Created by Greg Landrum (greg.landrum@t5informatics.com), Feb 2016
#
from flask import Flask, make_response, request, jsonify
try:
    from flasgger import Swagger
except ImportError:
    Swagger = None

from rdkit import Chem
from rdkit import rdBase
from rdkit.Chem import AllChem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from io import StringIO
import json
import sys

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
    res = dict(rdkitVersion=rdBase.rdkitVersion, boostVersion=rdBase.boostVersion,
               pythonVersion=sys.version)
    return json.dumps(res)


def _stringtobool(text):
    if text in ('0', 'false', 'False'):
        return False
    return bool(text)


def _molfromrequest():
    # get errors on stderr:
    sio = sys.stderr = StringIO()
    tgt = request.get_json()
    if tgt is None:
        tgt = request.values
    sanitize = _stringtobool(tgt.get('sanitize', True))
    removeHs = _stringtobool(tgt.get('removeHs', True))
    asSmarts = _stringtobool(tgt.get('asSmarts', False))
    asRxn = _stringtobool(tgt.get('asRxn', False))
    isMol = True
    if 'smiles' in tgt:
        if asRxn:
            mol = AllChem.ReactionFromSmarts(
                tgt.get('smiles'), useSmiles=not asSmarts)
        elif asSmarts:
            mol = Chem.MolFromSmarts(tgt.get('smiles'))
            sanitize = False
        else:
            mol = Chem.MolFromSmiles(tgt.get('smiles'), sanitize=sanitize)
    elif 'mol' in tgt:
        mol = Chem.MolFromMolBlock(
            tgt.get('mol'), sanitize=sanitize, removeHs=removeHs)
    else:
        raise InvalidUsage(
            "Neither 'smiles' nor 'mol' present.", status_code=410)

    if mol is None:
        errm = sio.getvalue()
        # some errors leave blank lines
        errm = errm.replace('RDKit ERROR: \n', '')
        raise InvalidUsage("Molecule could not be processed. Error message was:\n%s" % errm,
                           status_code=411)
    if not sanitize and not asRxn:
        mol.UpdatePropertyCache(False)
        Chem.FastFindRings(mol)
        Chem.SetConjugation(mol)
        Chem.SetHybridization(mol)
    return mol


def _queryfromrequest(suffix='_query'):
    # get errors on stderr:
    tgt = request.get_json()
    if tgt is None:
        tgt = request.values

    sio = sys.stderr = StringIO()
    if 'smiles' + suffix in tgt:
        mol = Chem.MolFromSmiles(tgt.get('smiles' + suffix), sanitize=False)
        if mol is not None:
            try:
                Chem.SanitizeMol(mol)
            except:
                mol = None
    elif 'smarts' + suffix in tgt:
        mol = Chem.MolFromSmarts(tgt.get('smarts' + suffix))
    elif 'mol' + suffix in tgt:
        mol = Chem.MolFromMolBlock(tgt.get('mol' + suffix), removeHs=False)
        mol = Chem.AdjustQueryProperties(mol)
    else:
        return None
    if mol is None:
        errm = sio.getvalue()
        # some errors leave blank lines
        errm = errm.replace('RDKit ERROR: \n', '')
        raise InvalidUsage("Molecule could not be processed. Error message was:\n%s" % errm,
                           status_code=411)
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
    return json.dumps(dict(smiles=Chem.MolToSmiles(mol, True)))


import json


def _loadJSONParam(text):
    if not text:
        return None
    return json.loads(text)

def _drawHelper(mol, drawer, **kwargs):
    tgt = request.get_json()
    if tgt is None:
        tgt = request.values
    sanit = _stringtobool(tgt.get('sanitize', True))
    kekulize = _stringtobool(tgt.get('kekulize', True))

    for arg in ('highlightAtomColors', 'highlightBondColors'):
        if arg not in kwargs:
            tmp = _loadJSONParam(tgt.get(arg, None))
            if tmp:
                for k in list(tmp):
                    tmp[int(k)] = tuple(tmp[k])
                    del tmp[k]
                kwargs[arg] = tmp
    for arg in ('highlightAtoms', 'highlightBonds'):
        if arg not in kwargs:
            tmp = _loadJSONParam(tgt.get(arg, None))
            if tmp:
                kwargs[arg] = tmp
    if 'legend' not in kwargs:
        kwargs['legend'] = tgt.get('legend', '')
    # if 'lw' not in kwargs and 'lw' in request.values:
    #     kwargs['lw'] = int(request.values['lw'])
    if 'highlightSubstruct' not in kwargs:
        if 'highlightColor' not in kwargs:
            hcolor = _loadJSONParam(tgt.get('highlightColor', None))
            if hcolor:
                highlightColor = tuple(hcolor)
            else:
                highlightColor = None
        else:
            highlightColor = kwargs['highlightColor']
            del kwargs['higlightColor']
        if _stringtobool(tgt.get('highlightSubstruct', False)) or \
          'smiles_highlight' in tgt or \
          'smarts_highlight' in tgt or \
          'mol_highlight' in tgt:
            qmol = _queryfromrequest(suffix='_highlight')
            if qmol is not None:
                qMatches = mol.GetSubstructMatches(qmol)
                highlights = kwargs.get('highlightAtoms', [])
                if highlightColor is not None:
                    highlightBonds = kwargs.get('highlightBonds', [])
                else:
                    highlightBonds = None

                for entry in qMatches:
                    if highlightColor is not None:
                        had = kwargs.get('highlightAtomColors', {})
                        for v in entry:
                            had[v] = highlightColor
                        kwargs['highlightAtomColors'] = had

                        hbd = kwargs.get('highlightBondColors', {})
                        emap = dict(enumerate(entry))
                        for bond in qmol.GetBonds():
                            mbond = mol.GetBondBetweenAtoms(emap[bond.GetBeginAtomIdx()],
                                                            emap[bond.GetEndAtomIdx()])
                            highlightBonds.append(mbond.GetIdx())
                            hbd[mbond.GetIdx()] = highlightColor
                        kwargs['highlightBondColors'] = hbd
                    highlights.extend(list(entry))
                if highlights:
                    kwargs['highlightAtoms'] = highlights
                    if highlightBonds is not None:
                        kwargs['highlightBonds'] = highlightBonds
    from rdkit.Chem.Draw import rdDepictor
    rdDepictor.SetPreferCoordGen(True)
    # print(Chem.MolToMolBlock(mc))
    drawo = drawer.drawOptions()
    if _stringtobool(tgt.get('dummiesAreAttachments', False)):
        drawo.dummiesAreAttachments = True

    if 'backgroundColor' in tgt:
        tmp = _loadJSONParam(tgt.get('backgroundColor', None))

        drawo.SetBackgroundColor(tuple)

    if _stringtobool(tgt.get('atomIndices',False)):
      drawo.addAtomIndices = True
    if _stringtobool(tgt.get('bondIndices',False)):
      drawo.addBondIndices = True


    # if 'lw' in kwargs:
    #     drawer.SetLineWidth(int(lw))
    if _stringtobool(tgt.get('useGrid', '0')):
        frags = []
        for frag in Chem.GetMolFrags(mol, asMols=True):
            frags.append(rdMolDraw2D.PrepareMolForDrawing(frag,
                                                          kekulize=kekulize & sanit,
                                                          addChiralHs=sanit))

        drawer.DrawMolecules(frags)
    elif _stringtobool(tgt.get('asRxn', False)):
        drawer.DrawReaction(mol)
    else:
        mc = rdMolDraw2D.PrepareMolForDrawing(
            mol, kekulize=kekulize & sanit, addChiralHs=sanit)
        mc.SetProp("_Name", "temp")
        drawer.DrawMolecule(mc, **kwargs)
    drawer.FinishDrawing()


def _moltosvg(mol, molSize=(450, 200), drawer=None, **kwargs):
    if drawer is None:
        tgt = request.get_json()
        if tgt is None:
            tgt = request.values
        if _stringtobool(tgt.get('useGrid', '0')):
            print("grid")
            nMols = len(Chem.GetMolFrags(mol))
            nCols = int(tgt.get('molsPerRow', 3))
            nRows = nMols // nCols
            if nMols % nCols:
                nRows += 1
            drawer = rdMolDraw2D.MolDraw2DSVG(
                molSize[0]*nCols, molSize[1]*nRows, molSize[0], molSize[1])
        else:
            drawer = rdMolDraw2D.MolDraw2DSVG(molSize[0], molSize[1])
    _drawHelper(mol, drawer, **kwargs)
    svg = drawer.GetDrawingText()
    return svg


def _moltopng(mol, molSize=(450, 200), drawer=None, **kwargs):
    if drawer is None:
        tgt = request.get_json()
        if tgt is None:
            tgt = request.values
        if _stringtobool(tgt.get('useGrid', '0')):
            print("grid")
            nMols = len(Chem.GetMolFrags(mol))
            nCols = int(tgt.get('molsPerRow', 3))
            nRows = nMols // nCols
            if nMols % nCols:
                nRows += 1
            drawer = rdMolDraw2D.MolDraw2DCairo(
                molSize[0]*nCols, molSize[1]*nRows, molSize[0], molSize[1])
        else:
            drawer = rdMolDraw2D.MolDraw2DCairo(molSize[0], molSize[1])
    _drawHelper(mol, drawer, **kwargs)
    return drawer.GetDrawingText()


def _render(mol, renderer, size=(150, 100), **kwargs):
    tgt = request.get_json()
    if tgt is None:
        tgt = request.values
    asRxn = _stringtobool(tgt.get('asRxn', False))
    if 'w' not in tgt and asRxn:
      s0 = size[0]*(mol.GetNumReactantTemplates()+mol.GetNumProductTemplates()+1)
      size = (s0,size[1])
    sz = int(tgt.get('w', size[0])), int(tgt.get('h', size[1]))
    return renderer(mol, molSize=sz, **kwargs)


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
        - name: kekulize
          in: query
          type: boolean
          required: false
          default: true
          description: determines whether or not the molecule is kekulized before being rendered
        - name: sanitize
          in: query
          type: boolean
          required: false
          default: true
          description: determines whether or not the molecule is sanitized before being rendered
        - name: removeHs
          in: query
          type: boolean
          required: false
          default: true
          description: determines whether or not Hs are removed from the molecule before being rendered
        - name: legend
          in: query
          type: string
          required: false
          description: text to be displayed beneath the molecule
        - name: highlightSubstruct
          in: query
          type: boolean
          required: false
          default: false
          description: indicates that substructure highlighting should be done
        - name: atomIndicess
          in: query
          type: boolean
          required: false
          default: false
          description: indicates that atom indices should be shown      
        - name: bondIndices
          in: query
          type: boolean
          required: false
          default: false
          description: indicates that bond indices should be shown
        - name: smiles_highlight
          in: query
          type: string
          required: false
          description: SMILES describing the query to highlight. If highlightSubstruct is set you should provide this, smarts_highlight, or mol_highlight
        - name: smarts_highlight
          in: query
          type: string
          required: false
          description: SMARTS describing the query to highlight. If highlightSubstruct is set you should provide this, smiles_highlight, or mol_highlight
        - name: mol_highlight
          in: query
          type: string
          required: false
          description: CTAB describing the query to highlight. If highlightSubstruct is set you should provide this, smarts_highlight, or smiles_highlight
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
        - name: dummiesAreAttachments
          in: query
          type: boolean
          required: false
          default: false
          description: if true, dummy atoms are drawn as attachment points (with a squiggle line instead of an atom symbol)
        - name: asRxn
          in: query
          type: boolean
          required: false
          default: false
          description: if true, the input will be treated as a reaction, not a molecule
        - name: useGrid
          in: query
          type: boolean
          required: false
          default: false
          description: if true, draws multiple molecules using a grid
        - name: molsPerRow
          in: query
          type: integer
          required: false
          default: 3
          description: sets the number of columns in the grid
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
    response = make_response(_render(mol, _moltopng))
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
        - name: kekulize
          in: query
          type: boolean
          required: false
          default: true
          description: determines whether or not the molecule is kekulized before being rendered
        - name: sanitize
          in: query
          type: boolean
          required: false
          default: true
          description: determines whether or not the molecule is sanitized before being rendered
        - name: removeHs
          in: query
          type: boolean
          required: false
          default: true
          description: determines whether or not Hs are removed from the molecule before being rendered
        - name: legend
          in: query
          type: string
          required: false
          description: text to be displayed beneath the molecule
        - name: highlightSubstruct
          in: query
          type: boolean
          required: false
          default: false
          description: indicates that substructure highlighting should be done
        - name: atomIndicess
          in: query
          type: boolean
          required: false
          default: false
          description: indicates that atom indices should be shown      
        - name: bondIndices
          in: query
          type: boolean
          required: false
          default: false
          description: indicates that bond indices should be shown
        - name: smiles_highlight
          in: query
          type: string
          required: false
          description: SMILES describing the query to highlight. If highlightSubstruct is set you should provide this, smarts_highlight, or mol_highlight
        - name: smarts_highlight
          in: query
          type: string
          required: false
          description: SMARTS describing the query to highlight. If highlightSubstruct is set you should provide this, smiles_highlight, or mol_highlight
        - name: mol_highlight
          in: query
          type: string
          required: false
          description: CTAB describing the query to highlight. If highlightSubstruct is set you should provide this, smarts_highlight, or smiles_highlight
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
        - name: dummiesAreAttachments
          in: query
          type: boolean
          required: false
          default: false
          description: if true, dummy atoms are drawn as attachment points (with a squiggle line instead of an atom symbol)
        - name: asRxn
          in: query
          type: boolean
          required: false
          default: false
          description: if true, the input will be treated as a reaction, not a molecule
        - name: useGrid
          in: query
          type: boolean
          required: false
          default: false
          description: if true, draws multiple molecules using a grid
        - name: molsPerRow
          in: query
          type: integer
          required: false
          default: 3
          description: sets the number of columns in the grid
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
    response = make_response(_render(mol, _moltosvg))
    #response.headers['Content-Type'] = 'image/svg+xml'
    return response


def _gen3d_sdf(mol, **kwargs):
    from rdkit.Chem import AllChem
    tgt = request.get_json()
    if tgt is None:
        tgt = request.values
    minimize = _stringtobool(tgt.get('minimize', False))
    seed = int(tgt.get('randomSeed', -1))

    mh = Chem.AddHs(mol)
    mh.SetProp("_Name", "2D to 3D output")
    ps = AllChem.ETKDG()
    ps.randomSeed = seed
    cid = AllChem.EmbedMolecule(mh, ps)
    if cid < 0:
        raise InvalidUsage("Molecule could not be embedded.", status_code=418)
    if minimize:
        try:
            AllChem.MMFFOptimizeMolecule(mh)
        except:
            raise InvalidUsage(
                "Molecule could not be minimized.", status_code=419)

    if not _stringtobool(tgt.get('returnHs', True)):
        mh = Chem.RemoveHs(mh)
    res = Chem.MolToMolBlock(mh)
    return res


@app.route('/to_3d/mol.sdf', methods=['GET', 'POST'])
def to_3d_sdf():
    """
      Returns a 3D conformation for a molecule
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
        - name: sanitize
          in: query
          type: boolean
          required: false
          default: true
          description: determines whether or not the molecule is sanitized before being processed
        - name: removeHs
          in: query
          type: boolean
          required: false
          default: true
          description: determines whether or not Hs are removed from the input molecule. Hs will be added to the full molecule before generating the conformer.
        - name: returnHs
          in: query
          type: boolean
          required: false
          default: true
          description: determines whether or not Hs are present in the returned structure
        - name: minimize
          in: query
          type: boolean
          required: false
          default: false
          description: determines whether or not the generated conformer is minimized (with the MMFF94 force field) before being returned
        - name: randomSeed
          in: query
          type: integer
          required: false
          default: -1
          description: provides a random number seed to be used in the embedding
      produces:
          chemical/x-mdl-molfile
      responses:
        500:
          description: Error!
        410:
          description: no molecule data provided
        411:
          description: input molecule could not be processed
        418:
          description: embedding failed
        419:
          description: minimization failed
        200:
          description: Everything is fine
      """
    mol = _molfromrequest()
    response = make_response(_gen3d_sdf(mol))
    response.headers['Content-Type'] = 'chemical/x-mdl-molfile'
    return response


if __name__ == '__main__':
    # FIX: turn this off pre-deployment
    app.debug = True
    app.run()
