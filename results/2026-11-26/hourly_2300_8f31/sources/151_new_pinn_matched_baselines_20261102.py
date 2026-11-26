"""Compare global, spatial and signed bounds at matched fixed meshes.

Global/spatial/signed quadratic share precisely the same quadratic cell data;
cubic separately uses its own interpolation nodes and fourth derivative.
No timing or generic superiority claim. Acceptance is one-sided U <= tau.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,json,os,sys
from flint import arb
from arb_certification import finite,exact_fraction
from arb_reconstruction import norm_enclosure
from arb_shaped_remainder import ShapedVerifier
from arb_cubic_reconstruction import CubicVerifier
ROOT=Path(__file__).resolve().parents[1];MESH=(16,32,64);METHODS=('global','spatial','signed_quadratic','signed_cubic');TOL=(F(1,1000),F(1,100000))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def atom(p,obj):
 b=(json.dumps(obj,indent=2,allow_nan=False)+'\n').encode();tmp=p.with_suffix('.tmp')
 with tmp.open('xb') as f:f.write(b);f.flush();os.fsync(f.fileno())
 os.replace(tmp,p)
def baseline(cells,depth):
 envelopes=[]
 for c in cells:
  r0,rm,r1=c['r0'],c['rm'],c['r1'];_,upper=norm_enclosure([r0,2*rm-(r0+r1)/2,r1],depth);envelopes.append(finite(upper+c['delta'],'rho').upper())
 glob=finite(max(envelopes)/8,'global').upper();slope=arb(0)
 for c,rho in zip(cells,envelopes):slope=finite(slope+rho*c['h']*(1-(c['a']+c['b'])/2),'slope')
 value=upper=arb(0)
 for c,rho in zip(cells,envelopes):
  h=c['h'];_,b=norm_enclosure([value,value+slope*h/2,value+slope*h-rho*h*h/2],depth);upper=max(upper,b);value=finite(value+slope*h-rho*h*h/2,'value');slope=finite(slope-rho*h,'slope update')
 if not value.contains(0):raise ArithmeticError('spatial envelope misses zero boundary')
 return exact_fraction(glob),exact_fraction(finite(upper,'spatial').upper())
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--training',type=Path,required=True);ap.add_argument('--certificates',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=False)
 fixture=a.training/'frozen_networks.json';h=sha(fixture);ts=json.loads((a.training/'summary.json').read_text());cs=json.loads((a.certificates/'summary.json').read_text());refp=a.certificates/'references.json';recp=a.certificates/'records.json'
 if h!=ts['fixture_sha256'] or h!=cs['fixture_sha256'] or sha(refp)!=cs['references_sha256'] or sha(recp)!=cs['records_sha256']:raise ValueError('input hash mismatch')
 nets=json.loads(fixture.read_bytes())['networks'];refs={r['network']:r for r in json.loads(refp.read_text())};old={(r['network'],r['degree'],r['intervals']):r for r in json.loads(recp.read_text())};labels={n['label'] for n in nets}
 if len(labels)!=len(nets) or len(refs)!=len(nets) or set(refs)!=labels:raise ValueError('input manifest')
 src=(ROOT/'scripts/74_polynomial_forcing_prototype.py').read_text().split("fixture=ROOT/'configs/")[0];needle="OUT=ROOT.parent/'run_20261030/polynomial_forcing';OUT.mkdir(parents=True,exist_ok=False)"
 if src.count(needle)!=1:raise ValueError('prototype changed')
 ns={'__file__':str(ROOT/'scripts/74_polynomial_forcing_prototype.py')};exec(compile(src.replace(needle,'OUT=None'),'defs','exec'),ns);rows=[]
 for net in nets:
  q=refs[net['label']];rl,ru=F(q['lower']),F(q['upper'])
  for mesh in MESH:
   p=ns['PolynomialForcing'](fixture,net['label'],h,net['source_coefficients']);v=ShapedVerifier(p,'moment');v.uniform(mesh);r=v.certify();counts=(v.residual_evaluations,v.curvature_evaluations);cells=v.prepare(v.nodes);g,s=baseline(cells,v.depth);unchanged=counts==(v.residual_evaluations,v.curvature_evaluations);ql,qu=exact_fraction(r['lower']),exact_fraction(r['upper'])
   p3=ns['PolynomialForcing'](fixture,net['label'],h,net['source_coefficients']);v3=CubicVerifier(p3,'moment');v3.uniform(mesh);r3=v3.certify();cl,cu=exact_fraction(r3['lower']),exact_fraction(r3['upper']);models={'global':(F(0),g),'spatial':(F(0),s),'signed_quadratic':(ql,qu),'signed_cubic':(cl,cu)}
   for method,(lo,hi) in models.items():
    prior=None
    if method.startswith('signed_'):
     z=old[(net['label'],2 if method=='signed_quadratic' else 3,mesh)];prior=(str(lo)==z['lower'] and str(hi)==z['upper'])
    rows.append({'network':net['label'],'source':net['source_name'],'depth':net['depth'],'seed':net['seed'],'intervals':mesh,'method':method,'lower':str(lo),'upper':str(hi),'reference_lower':str(rl),'reference_upper':str(ru),'reference_converged':q['converged'],'reference_overlap':lo<=ru and rl<=hi,'reference_contained':lo<=rl and ru<=hi,'accepted':{str(t):hi<=t for t in TOL},'quadratic_cache_reused':unchanged if method!='signed_cubic' else None,'matches_prior':prior,'residual_calls':counts[0] if method!='signed_cubic' else v3.residual_evaluations,'derivative_calls':counts[1] if method!='signed_cubic' else v3.curvature_evaluations})
  print(net['label'],'done',flush=True)
 atom(a.out/'records.json',rows);summary={'fixture_sha256':h,'input_records_sha256':cs['records_sha256'],'references_sha256':cs['references_sha256'],'meshes':MESH,'methods':METHODS,'tolerances':[str(t) for t in TOL],'expected':len(nets)*len(MESH)*len(METHODS),'records':len(rows),'reference_overlap':sum(r['reference_overlap'] for r in rows),'reference_contained':sum(r['reference_contained'] for r in rows),'shared_cache_passed':sum(r['quadratic_cache_reused'] is True for r in rows),'prior_matched':sum(r['matches_prior'] is True for r in rows),'counts':{},'records_sha256':sha(a.out/'records.json'),'script_sha256':sha(Path(__file__)),'failures':[{'network':r['network'],'method':r['method'],'intervals':r['intervals']} for r in rows if not r['reference_converged'] or not r['reference_overlap'] or r['quadratic_cache_reused'] is False or r['matches_prior'] is False],'scope':'one-sided acceptance at matched mesh; baselines and quadratic share cell data; no runtime comparison'}
 for mesh in MESH:
  for method in METHODS:
   group=[r for r in rows if r['intervals']==mesh and r['method']==method];summary['counts'][f'{mesh}/{method}']={'accepted':{str(t):sum(r['accepted'][str(t)] for r in group) for t in TOL},'reference_contained':sum(r['reference_contained'] for r in group)}
 if len(rows)!=summary['expected']:raise ValueError('record count')
 atom(a.out/'summary.json',summary);print(json.dumps(summary),flush=True);return int(bool(summary['failures']))
if __name__=='__main__':sys.exit(main())
