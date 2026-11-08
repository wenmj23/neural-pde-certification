"""Post-exit manifest and rational consistency of saved sampled errors."""
from fractions import Fraction as F
from pathlib import Path
from statistics import median
import hashlib,json
ROOT=Path(__file__).resolve().parents[1];B=ROOT/'results/2026-11-08';D=B/'midterm_independent_sample';OUT=B/'midterm_independent_sample_readback.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 if OUT.exists():raise FileExistsError(OUT)
 c=json.loads((D/'locked_config.json').read_text());s=json.loads((D/'summary.json').read_text());rows=json.loads((D/'records.json').read_text());train=ROOT/c['training_dir'];certdir=ROOT/c['certificates_dir'];fixture=train/'frozen_networks.json';refs=certdir/'references.json';cert=certdir/'records.json'
 if sha(D/'records.json')!=s['records_sha256'] or sha(D/'locked_config.json')!=sha(ROOT/'configs/midterm_independent_sample_2026-11-08.json') or {'fixture':sha(fixture),'references':sha(refs),'certificates':sha(cert)}!=s['inputs_sha256']:raise ValueError('digest')
 if any(sha(ROOT/p)!=h for p,h in s['source_hashes_before_run'].items()):raise ValueError('source hash')
 networks={n['label']:n for n in json.loads(fixture.read_text())['networks']};references={q['network']:q for q in json.loads(refs.read_text())};certs={(r['network'],r['degree'],r['intervals']):r for r in json.loads(cert.read_text())};expected={(n,g) for n in networks for g in c['grids']};seen=set();idx={}
 for r in rows:
  key=r['network'],r['grid_intervals']
  if key not in expected or key in seen:raise ValueError('manifest')
  seen.add(key);n=networks[key[0]];q=references[key[0]];x=F(r['sample_x']);a,b=F(r['sample_error_80']),F(r['sample_error_110']);rl,ru=F(q['lower']),F(q['upper']);stable=abs(a-b)<=F(c['relative_stability'])*max(F(1),abs(b))
  if r['source']!=n['source_name'] or not q['converged'] or (F(r['reference_lower']),F(r['reference_upper']))!=(rl,ru) or x<0 or x>1 or (x*key[1]).denominator!=1 or not 0<=a or not 0<=b or r['stable']!=stable or not stable or r['sample_below_reference_upper']!=(b<=ru) or b>ru:raise ValueError('sample/reference')
  for d in c['degrees']:
   item=r['bounds'][str(d)];old=certs[(key[0],d,c['certification_mesh'])];l,u=F(item['lower']),F(item['upper'])
   if (l,u)!=(F(old['lower']),F(old['upper'])) or item['sample_below_upper']!=(b<=u) or item['sample_above_lower']!=(b>=l) or b>u:raise ValueError('certificate')
  idx[key]=r
 if seen!=expected or len(networks)!=12 or s['records']!=24 or s['stable']!=24 or s['below_reference_upper']!=24 or s['below_certificate_upper']!=48:raise ValueError('summary')
 ratios=[];improved=0
 for n in networks:
  coarse=F(idx[(n,64)]['sample_error_80']);fine=F(idx[(n,256)]['sample_error_80']);q=references[n];mid=(F(q['lower'])+F(q['upper']))/2
  if fine<coarse or mid<=0:raise ValueError('nested grids')
  improved+=fine>coarse;ratios.append(float(fine/mid))
 out={'networks':len(networks),'grid_records':len(rows),'stable':len(rows),'below_reference_upper':len(rows),'below_certificate_upper':48,'fine_grid_increased_sample_max_networks':improved,'sample256_over_reference_midpoint':{'min':min(ratios),'median':median(ratios),'max':max(ratios)},'records_sha256':sha(D/'records.json'),'scope':'readback of saved finite numerical samples and existing bounds; not independent recomputation of sampled forward passes or a continuous proof'}
 with OUT.open('x') as f:json.dump(out,f,indent=2);f.write('\n')
 print(json.dumps(out,indent=2))
if __name__=='__main__':main()
