"""Read back benchmark records and recompute fixed-mesh timing medians."""
from pathlib import Path
from itertools import product
from statistics import median
import argparse,hashlib,json,math,sys

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def need(b,msg):
 if not b:raise ValueError(msg)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--training',required=True,type=Path);ap.add_argument('--certificates',required=True,type=Path);ap.add_argument('--timing',required=True,type=Path);a=ap.parse_args();t,c,p=a.training,a.certificates,a.timing;s=json.loads((p/'summary.json').read_text());ts=json.loads((t/'summary.json').read_text());cs=json.loads((c/'summary.json').read_text());fixture=t/'frozen_networks.json';need(sha(fixture)==ts['fixture_sha256']==cs['fixture_sha256']==s['fixture_sha256'],'fixture mismatch');need(sha(c/'records.json')==cs['records_sha256']==s['certificates_sha256'],'certificate hash');need(sha(p/'records.json')==s['records_sha256'] and sha(p/'pairs.json')==s['pairs_sha256'],'timing hash');need(sha(Path(__file__).with_name('146_new_forcing_timing_20261102.py'))==s['script_sha256'],'script hash')
 nets=json.loads(fixture.read_text())['networks'];labels={n['label'] for n in nets};need(len(labels)==len(nets)==12,'network labels');prior={(r['network'],r['degree'],r['intervals']):r for r in json.loads((c/'records.json').read_text())};rows=json.loads((p/'records.json').read_text());keys=[(r['network'],r['degree'],r['intervals'],r['repeat']) for r in rows];expected=set(product(labels,(2,3),(16,64),range(3)));need(len(rows)==len(set(keys))==len(expected)==s['records']==s['expected']==144 and set(keys)==expected,'timing key coverage');need({r['sequence'] for r in rows}==set(range(144)),'sequence')
 for r in rows:
  q=prior[(r['network'],r['degree'],r['intervals'])];need(type(r['wall_ns']) is int and r['wall_ns']>0 and type(r['cpu_ns']) is int and r['cpu_ns']>0,'timing value');match=r['lower']==q['lower'] and r['upper']==q['upper'];need(r['matches_previous']==match and match,'endpoint mismatch')
 pairs=json.loads((p/'pairs.json').read_text());pairkeys=[(r['network'],r['intervals']) for r in pairs];need(len(pairs)==len(set(pairkeys))==24 and set(pairkeys)==set(product(labels,(16,64))),'pair coverage')
 for r in pairs:
  med={d:median(x['wall_ns'] for x in rows if x['network']==r['network'] and x['degree']==d and x['intervals']==r['intervals']) for d in (2,3)};need(r['degree2_median_ns']==med[2] and r['degree3_median_ns']==med[3] and math.isclose(r['ratio_cubic_over_quadratic'],med[3]/med[2],rel_tol=1e-15),'median mismatch')
 need(s['endpoint_matches']==144 and not s['failures'],'summary');out={}
 for mesh in (16,64):
  ratios=[r['ratio_cubic_over_quadratic'] for r in pairs if r['intervals']==mesh];out[str(mesh)]={'median':median(ratios),'minimum':min(ratios),'maximum':max(ratios),'cubic_faster':sum(x<1 for x in ratios),'pairs':len(ratios)}
 print(json.dumps({'ok':True,'records':len(rows),'pairs':len(pairs),'ratios':out,'scope':'readback arithmetic only, no performance generalization'}));return 0
if __name__=='__main__':
 try:sys.exit(main())
 except (OSError,ValueError,KeyError,TypeError,ZeroDivisionError) as e:print(json.dumps({'ok':False,'error':str(e)}));sys.exit(1)
