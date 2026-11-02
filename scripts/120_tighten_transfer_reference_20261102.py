"""Post-hoc refinement of reference intervals escaping certificate bounds."""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,json,os,sys
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def atom(p,obj):
 b=(json.dumps(obj,indent=2)+'\n').encode();t=p.with_suffix('.tmp')
 with t.open('xb') as f:f.write(b);f.flush();os.fsync(f.fileno())
 os.replace(t,p)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--fixtures',type=Path,nargs='+',required=True);ap.add_argument('--input',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=False)
 s=json.loads((a.input/'summary.json').read_text());records=json.loads((a.input/'records.json').read_text());refs=json.loads((a.input/'references.json').read_text())
 if sha(a.input/'records.json')!=s['records_sha256'] or sha(a.input/'references.json')!=s['references_sha256']:raise ValueError('input hash mismatch')
 fixtures={sha(p):(p,{n['label']:n for n in json.loads(p.read_bytes())['networks']}) for p in a.fixtures}
 if set(fixtures)!={x['sha256'] for x in s['fixtures']}:raise ValueError('fixture hash mismatch')
 src=(ROOT/'scripts/74_polynomial_forcing_prototype.py').read_text().split("fixture=ROOT/'configs/")[0];needle="OUT=ROOT.parent/'run_20261030/polynomial_forcing';OUT.mkdir(parents=True,exist_ok=False)"
 if src.count(needle)!=1:raise ValueError('prototype layout changed')
 ns={'__file__':str(ROOT/'scripts/74_polynomial_forcing_prototype.py')};exec(compile(src.replace(needle,'OUT=None'),'polydefs','exec'),ns)
 old={(x['fixture_sha256'],x['network'],x['source']):x for x in refs};selected={(x['fixture_sha256'],x['network'],x['source']) for x in records if not x['reference_contained']};out=[]
 for h,label,source in sorted(selected):
  raw=fixtures[h][1][label];ref=ns['PolynomialReference'](raw,s['sources'][source]).maximum(relative=F(1,10**10),absolute=F(1,10**15),cap=5000);lo,hi=ns['fr'](ref['lower']),ns['fr'](ref['upper']);prev=old[h,label,source];pair=[x for x in records if (x['fixture_sha256'],x['network'],x['source'])==(h,label,source)]
  for row in pair:
   L,U=F(row['lower']),F(row['upper']);out.append({'fixture_sha256':h,'network':label,'source':source,'degree':row['degree'],'old_reference_lower':prev['lower'],'old_reference_upper':prev['upper'],'new_reference_lower':str(lo),'new_reference_upper':str(hi),'certificate_lower':row['lower'],'certificate_upper':row['upper'],'reference_converged':bool(ref['converged']),'reference_contained':L<=lo and hi<=U,'reference_overlap':L<=hi and lo<=U,'new_overlaps_old':lo<=F(prev['upper']) and F(prev['lower'])<=hi,'splits':ref['splits']})
  print(label,source,'refined',flush=True)
 atom(a.out/'records.json',out);summary={'selected_reference_cases':len(selected),'expected':len(selected)*2,'records':len(out),'converged':sum(r['reference_converged'] for r in out),'contained':sum(r['reference_contained'] for r in out),'overlap':sum(r['reference_overlap'] for r in out),'old_new_overlap':sum(r['new_overlaps_old'] for r in out),'records_sha256':sha(a.out/'records.json'),'input_records_sha256':s['records_sha256'],'input_references_sha256':s['references_sha256'],'fixtures':[{'sha256':h,'networks':len(z[1])} for h,z in fixtures.items()],'script_sha256':sha(Path(__file__)),'failures':[r for r in out if not r['reference_converged'] or not r['reference_overlap'] or not r['new_overlaps_old']],'scope':'post-hoc tighter conditional continuous reference; original wider reference retained; not formal verification'}
 if len(out)!=summary['expected']:raise ValueError('unexpected count')
 atom(a.out/'summary.json',summary);print(json.dumps(summary),flush=True);return int(bool(summary['failures']))
if __name__=='__main__':sys.exit(main())
