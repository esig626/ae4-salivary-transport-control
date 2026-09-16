#!/usr/bin/env python3
"""Reproduce pixel-to-axis digitisation; no model execution or fitted constants.
Run with the provided runtime Python (Pillow required). Coordinates refer to original
PMC JPEG pixels (origin upper left). Manual picks are persisted, never entered as
physical measurements; affine calibration maps create every reported value.
"""
from pathlib import Path
import csv, json, hashlib
from PIL import Image, ImageDraw
import numpy as np
P=Path(__file__).resolve().parent
axes={
 '4A':{'file':'figure4.jpg','x':[[128.5,0],[386.5,10]],'y':[[418.5,0],[40.5,3]],'x_unit':'min','y_unit':'uL/min'},
 '4B':{'file':'figure4.jpg','y':[[425.5,0],[41.5,15]],'y_unit':'uL/10min'},
 '5A':{'file':'figure5.jpg','y':[[300,0],[43,45]],'y_unit':'uL/10min'},
 '5B':{'file':'figure5.jpg','x':[[301.5,0],[473,10]],'y':[[294.5,0],[45.5,6]],'x_unit':'min','y_unit':'uL/min'},
 '5C':{'file':'figure5.jpg','x':[[611,0],[774,500]],'y':[[302,.8],[45,1.2]],'x_unit':'s','y_unit':'F0/F (label also V/V0)'},
 '5E':{'file':'figure5.jpg','y':[[681,0],[425,400]],'y_unit':'pA absolute change'}
}
# Raw manually localised marker centres/error-cap pixels, from zoomed original JPEG.
# Each row: x, y_mean, y_upper_SEM, y_lower_SEM. None means not resolved.
series={
 '4A_control':{'panel':'4A','n':6,'picks':[
 [128.5,418.5,None,None],[154,308,299,318],[180,280.5,270,291.5],[206,280.5,272,289],
 [232,291,281,300.5],[258,287.5,278.5,297.5],[284,296.5,285,308.5],[309,298,284,312],
 [335,290.5,274.5,308.5],[361,282.5,270,294.5],[387,301.5,288.5,314]]},
 '4A_Tmem16A_KO':{'panel':'4A','n':6,'picks':[
 [128.5,418.5,None,None],[154,335,321,350.5],[180,324,302.5,345.5],[206,318,297.5,339],
 [232,330,315,345],[258,328.5,310.5,346.5],[284,328.5,312.5,344.5],[309,322.5,304,340.5],
 [335,327,310,344.5],[361,335.5,319,352.5],[387,325.5,308.5,342.5]]},
 '4B_control':{'panel':'4B','n':6,'picks':[[650,169,150,None]]},
 '4B_Tmem16A_KO':{'panel':'4B','n':6,'picks':[[699,242,207,None]]},
 '5A_control':{'panel':'5A','n':8,'picks':[[113,191,170,None]]},
 '5A_bicarbonate_free':{'panel':'5A','n':8,'picks':[[135,286,280,None]]},
 '5A_low_chloride':{'panel':'5A','n':8,'picks':[[157,298,None,None]]},
 '5B_control':{'panel':'5B','n':8,'picks':[
 [301.5,294.5,None,None],[319,213,None,None],[336,227.5,219,None],[353,220.5,205,None],
 [370,221,205,None],[387,214.5,195,None],[405,216.5,199.5,None],[422,216.5,195,None],
 [439,219.5,202,None],[456,218.5,199,None],[473,218.5,195.5,None]]},
 '5B_Clcn2_KO':{'panel':'5B','n':8,'picks':[
 [301.5,294.5,None,None],[319,225,None,232],[336,239,None,247],[353,237,228,246],
 [370,235.5,None,243],[387,234.5,228,None],[405,240,None,247.5],[422,241.5,None,248.5],
 [439,242,233,250],[456,247.5,241,254],[473,243,None,249.5]]},
 '5B_DCPIB':{'panel':'5B','n':8,'picks':[
 [301.5,294.5,None,None],[319,275.5,None,None],[336,277.5,None,None],[353,280.5,None,None],
 [370,282.5,None,None],[387,284,None,None],[405,286,None,None],[422,287.5,None,None],
 [439,288,None,None],[456,288,None,None],[473,288,None,None]]},
 '5B_NPPB':{'panel':'5B','n':6,'picks':[[x,294.5,None,None] for x in [301.5,319,336,353,370,387,405,422,439,456,473]]},
 '5E_IPR_induced':{'panel':'5E','n':4,'picks':[[541,547.5,525,None]]},
 '5E_IPR_DCPIB_block_magnitude':{'panel':'5E','n':4,'picks':[[556,612,601,None]]},
 '5E_hypotonic_induced':{'panel':'5E','n':3,'picks':[[606,566.5,510,None]]},
 '5E_hypotonic_DCPIB_block_magnitude':{'panel':'5E','n':3,'picks':[[622,621.5,618,None]]}
}
def affine(z,pairs):
 (p0,v0),(p1,v1)=pairs
 return v0+(z-p0)*(v1-v0)/(p1-p0)
def yy(z,panel):return affine(z,axes[panel]['y'])
rows=[]
for name,s in series.items():
 panel=s['panel'];a=axes[panel]
 for j,(x,y,up,down) in enumerate(s['picks']):
  r={'series':name,'panel':panel,'index':j,'n':s['n'],'x_pixel':x,'y_pixel':y,'sem_upper_pixel':up,'sem_lower_pixel':down,
   'time':affine(x,a['x']) if 'x' in a else None,'time_unit':a.get('x_unit'),
   'mean':yy(y,panel),'sem_upper':yy(up,panel)-yy(y,panel) if up is not None else None,
   'sem_lower':yy(y,panel)-yy(down,panel) if down is not None else None,'unit':a['y_unit'],
   'digitisation_halfwidth':abs(yy(y+1.5,panel)-yy(y,panel)),
   'sem_note':'Unresolved/overlapped SEM is empty, never zero. One-sided caps retained without assuming symmetry.'}
  rows.append(r)
# Dense Fig5C: fixed vertical pixel windows capture the IPR trace envelope.
# This is a bounded raster trace observation, NOT a reconstruction of biological SEM.
im=np.array(Image.open(P/'figure5.jpg').convert('L'))
volume=[]
for t in range(25,451,25):
 x=round(611+163*t/500)
 # Broad image regions isolate the full upper trace, not a +/-10px centre corridor.
 # At the right edge the printed IPR label starts at y=92; the trace ends at y<=89.
 # Early IPR/CCh symbols overlap: retain an envelope and flag ambiguity.
 high=58; low=90 if x>=740 else 200
 mask=(im[high:low+1,x-1:x+2]<150).any(axis=1)
 ys=np.flatnonzero(mask)+high
 if not len(ys):raise RuntimeError((t,x,'empty trace'))
 ylo,yhi=float(ys.min()),float(ys.max());mid=(ylo+yhi)/2
 volume.append({'time_s_nominal':t,'x_pixel':x,'time_s_calibrated':affine(x,axes['5C']['x']),
 'trace_y_min_pixel':ylo,'trace_y_max_pixel':yhi,'trace_mid':yy(mid,'5C'),
 'trace_lower':yy(yhi+1.5,'5C'),'trace_upper':yy(ylo-1.5,'5C'),
 'uncertainty_type':'raster marker/trace envelope +1.5px; not SEM',
 'early_symbol_overlap':t<=75})
meta={'article_doi':'10.1073/pnas.1415739112','pmcid':'PMC4343136','extraction_date':'2026-09-16',
 'coordinate_origin':'upper-left original JPEG; x right, y down; not resized crops',
 'method':'Manual pixel picks persisted before affine conversion; Fig5C deterministic vertical-envelope extraction.',
 'version':2,
 'correction':'Version1 +/-10px corridors clipped portions of the Fig5C symbol envelope. Version2 uses broad image regions with explicit bar/text exclusion; all original Version1 files are preserved in digitisation_v1/.',
 'figure5C_regions':{'x_below_740_y_px':[58,200],'x_at_least_740_y_px':[58,90],'dark_threshold':150,'strip_width_px':3,'early_symbol_overlap_time_le_s':75},
 'axis_anchor_halfwidth_px':1.0,'point_pick_halfwidth_px':1.5,
 'reported_digitisation_error':'1.5px mapped to y units only; axis uncertainty is separately declared; not biological SEM.',
 'source_hashes':{f:hashlib.sha256((P/f).read_bytes()).hexdigest() for f in ['figure4.jpg','figure5.jpg']},
 'axes':axes,'series':series,
 'cautions':['Fig5E grey bars are magnitudes BLOCKED by DCPIB, not residual currents after blockade.',
 'Fig5B NPPB markers conceal tiny or absent SEM. Do not declare exact zero.',
 'Whole-gland secretion is not normalised per acinar cell. Different control cohorts have materially different outputs.',
 'Fig5C trace endpoint differs from narrative 12.5±0.2% swelling; both are preserved and discrepancy not resolved.',
 'Current and swelling were measured in different experiments; no paired gain is inferred.',
 'Digitisation does not licence an Ohmic law for a reported outward-rectifying conductance.']}
(P/'digitisation_metadata.json').write_text(json.dumps(meta,indent=2)+'\n')
for file,data in [('digitised_points.csv',rows),('digitised_volume_trace.csv',volume)]:
 with (P/file).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
for filename in ['figure4.jpg','figure5.jpg']:
 im0=Image.open(P/filename).convert('RGB');d=ImageDraw.Draw(im0)
 for name,s in series.items():
  if axes[s['panel']]['file']!=filename:continue
  for x,y,up,dn in s['picks']:
   d.line((x-2,y,x+2,y),fill='red',width=1);d.line((x,y-2,x,y+2),fill='red',width=1)
   for q in [up,dn]:
    if q is not None:d.line((x-3,q,x+3,q),fill='blue',width=1)
 if filename=='figure5.jpg':
  for r in volume:d.line((r['x_pixel'],r['trace_y_min_pixel'],r['x_pixel'],r['trace_y_max_pixel']),fill='red',width=1)
 im0.save(P/(filename.replace('.jpg','_picks.png')))
summary={r['series']:{k:r[k] for k in ['mean','sem_upper','unit','digitisation_halfwidth']} for r in rows if r['panel'] in ['4B','5A','5E']}
summary['Figure5C_endpoint']=volume[-1]
(P/'digitised_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
