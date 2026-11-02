"""Independent readback of 153 timing keys, bounds, hashes, and medians."""
from pathlib import Path
from itertools import product
from statistics import median
import argparse,hashlib,json,math,sys
MESH=(16,64);METHODS=('global','spatial','signed_quadratic','signed_cubic');REPEATS=3
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def need(v,msg):
 if not v:raise ValueError(msg)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--training',type=Path,required=True);ap.add_argument('--baseline',type=Path,required=True);ap.add_argument('--timing',type=Path,required=True);a=ap.parse_args();fixture=a.training/'frozen_networks.json';h=sha(fixture);ts=json.loads((a.training/'summary.json').read_text());bs=json.loads((a.baseline/'summary.json').read_text());s=json.loads((a.timing/'summary.json').read_text());need(h==ts['fixture_sha256']==bs['fixture_sha256']==s['fixture_sha256'],'fixture');need(sha(a.baseline/'records.json')==bs['records_sha256']==s['baseline_sha256'],'baseline');need(sha(Path(__file__).with_name('151_new_pinn_matched_baselines_20261102.py'))==bs['script_sha256']==s['baseline_script_sha256'],'baseline script');need(sha(a.timing/'records.json')==s['records_sha256'] and sha(a.timing/'pairs.json')==s['pairs_sha256'] and sha(Path(__file__).with_name('153_matched_baseline_timing_20261102.py'))==s['script_sha256'],'result/script');need(s['meshes']==list(MESH) and s['methods']==list(METHODS) and s['repeats']==REPEATS,'config')
 nets=json.loads(fixture.read_text())['networks'];labels={n['label'] for n in nets};need(len(labels)==len(nets)==12,'labels');old={(r['network'],r['intervals'],r['method']):r for r in json.loads((a.baseline/'records.json').read_text())};rows=json.loads((a.timing/'records.json').read_text());keys=[(r['network'],r['intervals'],r['method'],r['repeat']) for r in rows];expected=set(product(labels,MESH,METHODS,range(REPEATS)));need(len(rows)==len(set(keys))==len(expected)==s['expected']==s['records']==288 and set(keys)==expected,'record manifest');need({r['sequence'] for r in rows}==set(range(288)),'sequence')
 for r in rows:
  prev=old[(r['network'],r['intervals'],r['method'])];need(type(r['wall_ns']) is int and r['wall_ns']>0 and type(r['cpu_ns']) is int and r['cpu_ns']>0,'duration');need(r['matches_baseline']==(r['upper']==prev['upper']) and r['matches_baseline'],'upper mismatch')
 pairs=json.loads((a.timing/'pairs.json').read_text());pk=[(r['network'],r['intervals']) for r in pairs];need(len(pairs)==len(set(pk))==24 and set(pk)==set(product(labels,MESH)),'pairs')
 for r in pairs:
  m={method:median(x['wall_ns'] for x in rows if x['network']==r['network'] and x['intervals']==r['intervals'] and x['method']==method) for method in METHODS};need(r['median_wall_ns']==m,'median');need(math.isclose(r['quadratic_over_global'],m['signed_quadratic']/m['global'],rel_tol=1e-15) and math.isclose(r['quadratic_over_spatial'],m['signed_quadratic']/m['spatial'],rel_tol=1e-15) and math.isclose(r['cubic_over_quadratic'],m['signed_cubic']/m['signed_quadratic'],rel_tol=1e-15),'ratios')
 need(s['matches_baseline']==sum(r['matches_baseline'] for r in rows)==288 and not s['failures'],'summary');out={}
 for n in MESH:
  g=[r for r in pairs if r['intervals']==n];out[str(n)]={name:median(r[name] for r in g) for name in ('quadratic_over_global','quadratic_over_spatial','cubic_over_quadratic')};out[str(n)]['median_ms']={m:median(r['median_wall_ns'][m] for r in g)/1e6 for m in METHODS}
 print(json.dumps({'ok':True,'records':len(rows),'pairs':len(pairs),'medians':out,'scope':'readback arithmetic, not universal performance'}));return 0
if __name__=='__main__':
 try:sys.exit(main())
 except (OSError,ValueError,KeyError,TypeError,ZeroDivisionError) as e:print(json.dumps({'ok':False,'error':str(e)}));sys.exit(1)
