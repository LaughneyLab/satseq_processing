#!python3

import pymol
from pymol import cmd
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa
#import numpy as np
import matplotlib
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import argparse
import json

parser = argparse.ArgumentParser(prog='color_sting.py', usage='%(prog)s [options]')
parser.add_argument('--png', nargs='?', help='output png file, default=sting_color.png', default='sting_color.png', type=str)
parser.add_argument('--pdb', nargs='?', default='Mouse STING homodimer(ranked_0).pdb', help='protein file (PDB)', type=str)
parser.add_argument('--colors', nargs='?', help='JSON string mapping residues to colors directly, e.g. \'{"112": "#ff4949", "141": "#ff9393", "145": "#ffffff"}\'')
parser.add_argument('--scores', nargs='?', help='JSON string with scores, use with cmap/vmin/vmax to map to colors. Alternative to --colors. e.g, \'{"159": -1.26955099963628, "22": -0.8021170040184462}\'. Will color chains A and B the same', type=str)
parser.add_argument('--cmap', nargs='?', default='bwr', help='Matplotlib colormap to use, default=bwr', type=str)
parser.add_argument('--vlim', nargs=2, default=[-2,2], help='vmin/vmax for colormap, default=[-2,2]', type=float)

args = parser.parse_args()

#print(f'args.png={args.png}')
#print(f'args.pdb={args.pdb}')

residue_colors=None
residue_scores=None
if args.colors is not None:
    if args.scores is not None:
        print('Warning; should provide either colors or scores, not both. ignoring scores')
    residue_colors = json.loads(args.colors)
    residues = residue_colors.keys()
elif args.scores is not None:
    residue_scores = json.loads(args.scores)
    residues = residue_scores.keys()
else:
    print('Error: need to provide args.colors or args.scores')
    quit()

# Start PyMOL
pymol.finish_launching(['pymol', '-cq'])

# Load PDB
#print(f'pdb file={args.pdb}')
cmd.load(args.pdb, "prot")
cmd.color("gray70", "prot")
cmd.hide("everything", "all")
cmd.show("cartoon", "prot")



# Apply coloring
#print(f'cmap={args.cmap} vlim=({args.vlim[0]}, {args.vlim[1]})')
if residue_scores is not None:
    norm = mcolors.Normalize(vmin=args.vlim[0], vmax=args.vlim[1])
    cmap = matplotlib.colormaps[args.cmap]

for i in residues:
    if residue_scores is not None:
        color_hex = mcolors.to_hex(cmap(norm(residue_scores[i])))
    else:
        color_hex = residue_colors[i]
    pymol_color=f'score_color_{i}'
    cmd.set_color(pymol_color, mcolors.to_rgb(color_hex))
    for chain in ['A', 'B']:
        selection=f'chain {chain} and resi {i}'
        cmd.color(pymol_color, selection)

# Render PNG
cmd.bg_color("white")
cmd.set("ray_opaque_background", 0)
#cmd.set_view([-0.593658030,   -0.054702803,    0.802855909,
#    -0.625859380,   -0.595750630,   -0.503369391,
#     0.505839169,   -0.801303744,    0.319434971,
#     0.000000000,    0.000000000, -306.655059814,
#     0.225017548,   -0.535675049,   -0.549560547,
#   241.769195557,  371.540924072,  -20.000000000])
cmd.set_view([0.479985863,   -0.019433912,   -0.877059817,
     0.655858099,   -0.656033754,    0.373463154,
    -0.582640409,   -0.754482031,   -0.302139550,
     0.000000000,    0.000000000, -306.655059814,
     0.225017548,   -0.535675049,   -0.549560547,
   261.217132568,  352.092956543,  -20.000000000 ])
#cmd.turn('y', 180)
#cmd.turn('x', 180)
#cmd.turn('z', 180)
cmd.ray(1200, 900)
cmd.png(args.png, dpi=300)
#print(f'Wrote {args.png}')
cmd.quit()
