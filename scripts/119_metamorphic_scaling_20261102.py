"""Check certificates under exact scaling (v,f)->(c*v,c*f), c in {-2,-1,1,2}."""
from pathlib import Path
from fractions import Fraction as F
from copy import deepcopy
import argparse,hashlib,json,os,sys
from arb_certification import exact_fraction
from arb_cubic_reconstruction import CubicVerifier
from arb_shaped_remainder import ShapedVerifier
ROOT=Path(__file__).resolve().parents[1];SCALES=(1,-1,2,-2);MESH=(32,64);DEGREES=(2,3)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def atom(p,obj):
 b=(json.dumps(obj,indent=2,allow_nan=False)+'\n').encode();tmp=p.with_suffix('.tmp')
 with tmp.open('xb') as f:f.write(b);f.flush();os.fsync(f.fileno())
 os.replace(tmp,p)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--fixtures',nargs='+',required=True,type=Path);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=False)
 source=(ROOT/'scripts/74_polynomial_forcing_prototype.py').read_text().split("fixture=ROOT/'configs/")[0];needle="OUT=ROOT.parent/'run_20261030/polynomial_forcing';OUT.mkdir(parents=True,exist_ok=False)"
 if source.count(needle)!=1:raise ValueError('prototype layout changed')
 ns={'__file__':str(ROOT/'scripts/74_polynomial_forcing_prototype.py')};exec(compile(source.replace(needle,'OUT=None'),'defs','exec'),ns)
 rows=[];fixtures=[]
 for index,path in enumerate(a.fixtures):
  digest=sha(path);nets=json.loads(path.read_bytes())['networks'];fixtures.append({'sha256':digest,'networks':len(nets)})
  derived=[]
  for net in nets:
   for scale in SCALES:
    n=deepcopy(net);n['label']=net['label']+'__scale_'+str(scale);layer=n['layers'][-1]
    if layer['activation']!='identity':raise ValueError('last layer not linear')
    def scale_exact(v):
     w=float(v)*scale
     if F.from_float(w)!=F.from_float(float(v))*scale:raise ValueError('scaled weight not exact binary64')
     return w
    layer['weight']=[[scale_exact(v) for v in row] for row in layer['weight']];layer['bias']=[scale_exact(v) for v in layer['bias']];n['source_coefficients']=[v*scale for v in net['source_coefficients']];derived.append(n)
  derived_path=a.out/f'derived_fixture_{index}.json';atom(derived_path,{'networks':derived});derived_sha=sha(derived_path)
  for net in nets:
   for mesh in MESH:
    for degree in DEGREES:
     baseline=None
     for scale in SCALES:
      label=net['label']+'__scale_'+str(scale);coeff=[c*scale for c in net['source_coefficients']];p=ns['PolynomialForcing'](derived_path,label,derived_sha,coeff);v=ShapedVerifier(p,'moment') if degree==2 else CubicVerifier(p,'moment');v.uniform(mesh);r=v.certify();lo,hi=exact_fraction(r['lower']),exact_fraction(r['upper']);normalized=(lo/abs(scale),hi/abs(scale))
      if baseline is None:baseline=normalized
      rows.append({'fixture_sha256':digest,'network':net['label'],'mesh':mesh,'degree':degree,'scale':scale,'lower':str(lo),'upper':str(hi),'normalized_lower':str(normalized[0]),'normalized_upper':str(normalized[1]),'baseline_lower':str(baseline[0]),'baseline_upper':str(baseline[1]),'equal_to_baseline':normalized==baseline,'overlap_baseline':normalized[0]<=baseline[1] and baseline[0]<=normalized[1]})
   print(net['label'],'done',flush=True)
 atom(a.out/'records.json',rows);n=sum(f['networks'] for f in fixtures);s={'fixtures':fixtures,'scales':SCALES,'meshes':MESH,'degrees':DEGREES,'expected':n*len(SCALES)*len(MESH)*len(DEGREES),'records':len(rows),'equal':sum(r['equal_to_baseline'] for r in rows),'overlap':sum(r['overlap_baseline'] for r in rows),'records_sha256':sha(a.out/'records.json'),'script_sha256':sha(Path(__file__)),'derived_fixture_sha256':[sha(a.out/f'derived_fixture_{i}.json') for i in range(len(fixtures))],'failures':[r for r in rows if not r['overlap_baseline']],'scope':'exact scale metamorphic check, same implementation; not independent validation or trained models'}
 if len(rows)!=s['expected']:raise ValueError('count mismatch')
 atom(a.out/'summary.json',s);print(json.dumps(s),flush=True);return int(bool(s['failures']))
if __name__=='__main__':sys.exit(main())
