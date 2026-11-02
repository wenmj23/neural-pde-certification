"""Randomized matched-grid in-process timing for four upper-bound methods.

The global branch returns immediately after a quadratic residual envelope;
the spatial branch additionally integrates the spatial envelope. No complete
startup time or equal-accuracy time is measured.
"""
from pathlib import Path
from statistics import median
import argparse,hashlib,importlib.util,json,os,random,sys,time
from flint import arb
from arb_certification import finite,exact_fraction
from arb_reconstruction import norm_enclosure
from arb_shaped_remainder import ShapedVerifier
from arb_cubic_reconstruction import CubicVerifier
ROOT=Path(__file__).resolve().parents[1];MESH=(16,64);METHODS=('global','spatial','signed_quadratic','signed_cubic');REPEATS=3;SEED=20261104
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def atom(p,obj):
 b=(json.dumps(obj,indent=2,allow_nan=False)+'\n').encode();tmp=p.with_suffix('.tmp')
 with tmp.open('xb') as f:f.write(b);f.flush();os.fsync(f.fileno())
 os.replace(tmp,p)
def global_upper(cells,depth):
 env=[]
 for c in cells:
  r0,rm,r1=c['r0'],c['rm'],c['r1'];_,u=norm_enclosure([r0,2*rm-(r0+r1)/2,r1],depth);env.append(finite(u+c['delta'],'rho').upper())
 return exact_fraction(finite(max(env)/8,'global').upper())
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--training',required=True,type=Path);ap.add_argument('--baseline',required=True,type=Path);ap.add_argument('--out',required=True,type=Path);a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=False)
 fixture=a.training/'frozen_networks.json';h=sha(fixture);ts=json.loads((a.training/'summary.json').read_text());bs=json.loads((a.baseline/'summary.json').read_text());rp=a.baseline/'records.json'
 if h!=ts['fixture_sha256'] or h!=bs['fixture_sha256'] or sha(rp)!=bs['records_sha256']:raise ValueError('input hash mismatch')
 nets=json.loads(fixture.read_bytes())['networks'];old={(r['network'],r['intervals'],r['method']):r for r in json.loads(rp.read_text())};expected={(n['label'],k,m) for n in nets for k in MESH for m in METHODS}
 if not expected.issubset(old) or len(expected)!=len(nets)*len(MESH)*len(METHODS):raise ValueError('input manifest')
 path=ROOT/'scripts/151_new_pinn_matched_baselines_20261102.py';spec=importlib.util.spec_from_file_location('comparison151',path);bm=importlib.util.module_from_spec(spec);spec.loader.exec_module(bm)
 if sha(path)!=bs['script_sha256']:raise ValueError('baseline script hash mismatch')
 src=(ROOT/'scripts/74_polynomial_forcing_prototype.py').read_text().split("fixture=ROOT/'configs/")[0];needle="OUT=ROOT.parent/'run_20261030/polynomial_forcing';OUT.mkdir(parents=True,exist_ok=False)"
 if src.count(needle)!=1:raise ValueError('prototype changed')
 ns={'__file__':str(ROOT/'scripts/74_polynomial_forcing_prototype.py')};exec(compile(src.replace(needle,'OUT=None'),'defs','exec'),ns)
 def run(net,mesh,method):
  start=time.perf_counter_ns();cpu=time.process_time_ns();p=ns['PolynomialForcing'](fixture,net['label'],h,net['source_coefficients']);v=CubicVerifier(p,'moment') if method=='signed_cubic' else ShapedVerifier(p,'moment');v.uniform(mesh)
  if method=='global':u=global_upper(v.prepare(v.nodes),v.depth)
  elif method=='spatial':_,u=bm.baseline(v.prepare(v.nodes),v.depth)
  else:u=exact_fraction(v.certify()['upper'])
  return time.perf_counter_ns()-start,time.process_time_ns()-cpu,str(u)
 for net in nets:
  for mesh in MESH:
   for method in METHODS:run(net,mesh,method)
 tasks=[(n,mesh,method,rep) for n in nets for mesh in MESH for method in METHODS for rep in range(REPEATS)];random.Random(SEED).shuffle(tasks);rows=[]
 for seq,(net,mesh,method,rep) in enumerate(tasks):
  wall,cpu,u=run(net,mesh,method);rows.append({'sequence':seq,'network':net['label'],'intervals':mesh,'method':method,'repeat':rep,'wall_ns':wall,'cpu_ns':cpu,'upper':u,'matches_baseline':u==old[(net['label'],mesh,method)]['upper']})
 atom(a.out/'records.json',rows);pairs=[]
 for net in nets:
  for mesh in MESH:
   med={m:median(r['wall_ns'] for r in rows if r['network']==net['label'] and r['intervals']==mesh and r['method']==m) for m in METHODS};pairs.append({'network':net['label'],'intervals':mesh,'median_wall_ns':med,'quadratic_over_global':med['signed_quadratic']/med['global'],'quadratic_over_spatial':med['signed_quadratic']/med['spatial'],'cubic_over_quadratic':med['signed_cubic']/med['signed_quadratic']})
 atom(a.out/'pairs.json',pairs);summary={'fixture_sha256':h,'baseline_sha256':bs['records_sha256'],'baseline_script_sha256':bs['script_sha256'],'meshes':MESH,'methods':METHODS,'repeats':REPEATS,'shuffle_seed':SEED,'expected':len(nets)*len(MESH)*len(METHODS)*REPEATS,'records':len(rows),'matches_baseline':sum(r['matches_baseline'] for r in rows),'records_sha256':sha(a.out/'records.json'),'pairs_sha256':sha(a.out/'pairs.json'),'script_sha256':sha(Path(__file__)),'failures':[{'network':r['network'],'intervals':r['intervals'],'method':r['method']} for r in rows if not r['matches_baseline'] or r['wall_ns']<=0 or r['cpu_ns']<=0],'scope':'one machine matched fixed-mesh timing, 1 warmup plus 3 repeats; fresh verifiers; excludes startup, training, references; spatial includes calculation of global bound as well'}
 if len(rows)!=summary['expected']:raise ValueError('count mismatch')
 atom(a.out/'summary.json',summary);print(json.dumps(summary),flush=True);return int(bool(summary['failures']))
if __name__=='__main__':sys.exit(main())
