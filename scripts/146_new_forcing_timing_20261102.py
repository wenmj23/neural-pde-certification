"""Small, randomized, matched-mesh timing benchmark for new frozen networks.

Times problem construction, fixture reading, verifier construction and one
fixed mesh certificate; excludes Python startup, training and reference work.
One warmup per network/degree/mesh; three measured fresh instances, shuffled.
"""
from pathlib import Path
from fractions import Fraction as F
from statistics import median
import argparse,hashlib,json,os,random,sys,time
from arb_certification import exact_fraction
from arb_cubic_reconstruction import CubicVerifier
from arb_shaped_remainder import ShapedVerifier
ROOT=Path(__file__).resolve().parents[1];MESH=(16,64);DEGREES=(2,3);REPEATS=3;SEED=20261102
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def atom(p,obj):
 b=(json.dumps(obj,indent=2,allow_nan=False)+'\n').encode();t=p.with_suffix('.tmp')
 with t.open('xb') as f:f.write(b);f.flush();os.fsync(f.fileno())
 os.replace(t,p)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--training',required=True,type=Path);ap.add_argument('--certificates',required=True,type=Path);ap.add_argument('--out',required=True,type=Path);a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=False)
 fixture=a.training/'frozen_networks.json';h=sha(fixture);ts=json.loads((a.training/'summary.json').read_text());cs=json.loads((a.certificates/'summary.json').read_text());cp=a.certificates/'records.json'
 if h!=ts['fixture_sha256'] or h!=cs['fixture_sha256'] or sha(cp)!=cs['records_sha256']:raise ValueError('input hash mismatch')
 nets=json.loads(fixture.read_bytes())['networks'];prior={(r['network'],r['degree'],r['intervals']):r for r in json.loads(cp.read_text())};expected={(n['label'],d,k) for n in nets for d in DEGREES for k in MESH}
 if not expected.issubset(prior) or len(expected)!=len(nets)*4:raise ValueError('input manifest')
 src=(ROOT/'scripts/74_polynomial_forcing_prototype.py').read_text().split("fixture=ROOT/'configs/")[0];needle="OUT=ROOT.parent/'run_20261030/polynomial_forcing';OUT.mkdir(parents=True,exist_ok=False)"
 if src.count(needle)!=1:raise ValueError('prototype changed')
 ns={'__file__':str(ROOT/'scripts/74_polynomial_forcing_prototype.py')};exec(compile(src.replace(needle,'OUT=None'),'defs','exec'),ns)
 def run(net,degree,mesh):
  start=time.perf_counter_ns();cpu=time.process_time_ns();p=ns['PolynomialForcing'](fixture,net['label'],h,net['source_coefficients']);v=ShapedVerifier(p,'moment') if degree==2 else CubicVerifier(p,'moment');v.uniform(mesh);r=v.certify();cpu=time.process_time_ns()-cpu;elapsed=time.perf_counter_ns()-start;return elapsed,cpu,str(exact_fraction(r['lower'])),str(exact_fraction(r['upper']))
 for net in nets:
  for d in DEGREES:
   for mesh in MESH:run(net,d,mesh)
 tasks=[(net,d,mesh,rep) for net in nets for d in DEGREES for mesh in MESH for rep in range(REPEATS)];random.Random(SEED).shuffle(tasks);rows=[]
 for i,(net,d,mesh,rep) in enumerate(tasks):
  elapsed,cpu,l,u=run(net,d,mesh);old=prior[(net['label'],d,mesh)];rows.append({'sequence':i,'network':net['label'],'degree':d,'intervals':mesh,'repeat':rep,'wall_ns':elapsed,'cpu_ns':cpu,'lower':l,'upper':u,'matches_previous':l==old['lower'] and u==old['upper']})
 atom(a.out/'records.json',rows);pairs=[]
 for net in nets:
  for mesh in MESH:
   times={d:median([r['wall_ns'] for r in rows if r['network']==net['label'] and r['degree']==d and r['intervals']==mesh]) for d in DEGREES};pairs.append({'network':net['label'],'intervals':mesh,'degree2_median_ns':times[2],'degree3_median_ns':times[3],'ratio_cubic_over_quadratic':times[3]/times[2]})
 atom(a.out/'pairs.json',pairs);summary={'fixture_sha256':h,'certificates_sha256':cs['records_sha256'],'seed':SEED,'meshes':MESH,'degrees':DEGREES,'repeats':REPEATS,'expected':len(nets)*len(MESH)*len(DEGREES)*REPEATS,'records':len(rows),'endpoint_matches':sum(r['matches_previous'] for r in rows),'records_sha256':sha(a.out/'records.json'),'pairs_sha256':sha(a.out/'pairs.json'),'script_sha256':sha(Path(__file__)),'scope':'matched fixed mesh, fresh verifier, single machine, one warm-up and three timing repeats; not time-to-goal or universal speed result','failures':[{'network':r['network'],'degree':r['degree'],'intervals':r['intervals']} for r in rows if not r['matches_previous'] or r['wall_ns']<=0 or r['cpu_ns']<=0]}
 if len(rows)!=summary['expected']:raise ValueError('count mismatch')
 atom(a.out/'summary.json',summary);print(json.dumps(summary),flush=True);return int(bool(summary['failures']))
if __name__=='__main__':sys.exit(main())
