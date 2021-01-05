from rdkit import Chem  
def DrawMolAsSVG(m): 
    Chem.rdDepictor.Compute2DCoords(m)
    drawer = Draw.rdMolDraw2D.MolDraw2DSVG(200,200)
    drawer.DrawMolecule(m)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText().replace('svg:','')
    display(SVG(svg))