"""Research prototype: resumable, validated decision-result persistence.

Journals one request per line and rejects incomplete/malformed/duplicate records.
Writes final results atomically; readers must check summary's record-file SHA.
--simulate-stop-after is a controlled interruption after a complete flushed row,
not a simulation of every kind of process, OS or power failure.
"""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
import argparse,hashlib,json,os,sys,time
from arb_cubic_reconstruction import CubicVerifier
from arb_shaped_remainder import ShapedVerifier
from arb_certification import exact_fraction
ROOT=Path(__file__).resolve().parents[1]
def sha(data):return hashlib.sha256(data).hexdigest()
def atomic(path,data):
 temp=path.with_name(path.name+'.tmp-'+str(os.getpid()))
 with temp.open('xb') as f:
  f.write(data);f.flush();os.fsync(f.fileno())
 os.replace(temp,path)
 try:
  fd=os.open(str(path.parent),os.O_RDONLY)
  try:os.fsync(fd)
  finally:os.close(fd)
 except OSError:pass  # directory fsync may not be available on every platform

def canonical(value):return (json.dumps(value,ensure_ascii=False,sort_keys=True,allow_nan=False)+'\n').encode()
def load_journal(path,plan):
 if not path.exists():return []
 data=path.read_bytes()
 if data and not data.endswith(b'\n'):raise ValueError('Journal has an incomplete last line')
 rows=[];seen=set()
 for lineno,line in enumerate(data.splitlines(),1):
  try:r=json.loads(line)
  except (ValueError,UnicodeDecodeError) as e:raise ValueError(f'Malformed journal line {lineno}') from e
  key=(r.get('network'),r.get('degree'),r.get('tolerance'))
  if key in seen:raise ValueError(f'Duplicate request at line {lineno}')
  if len(rows)>=len(plan) or key!=plan[len(rows)]:raise ValueError(f'Unexpected request sequence at line {lineno}')
  if not isinstance(r.get('trace'),list) or not r['trace']:raise ValueError(f'Empty trace at line {lineno}')
  tau=F(key[2])
  for entry in r['trace']:
   lo,hi=F(entry['lower']),F(entry['upper']);d='within' if hi<=tau else 'exceeds' if lo>tau else 'unknown'
   if lo<0 or hi<lo or entry['decision']!=d:raise ValueError(f'Invalid decision bounds at line {lineno}')
  if r.get('decision')!=r['trace'][-1]['decision']:raise ValueError(f'Final decision differs at line {lineno}')
  seen.add(key);rows.append(r)
 return rows

def verify_export(path):
 summ=json.loads((path/'summary.json').read_text());data=(path/'records.json').read_bytes()
 if sha(data)!=summ['records_sha256']:raise ValueError('records.json checksum differs from summary')
 rows=json.loads(data)
 if len(rows)!=summ['records_count'] or len(rows)!=summ['expected_count']:raise ValueError('record count differs from summary')
 return rows,summ

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--out',type=Path,required=True);parser.add_argument('--resume',action='store_true');parser.add_argument('--simulate-stop-after',type=int);args=parser.parse_args()
 if args.simulate_stop_after is not None and args.simulate_stop_after<1:raise ValueError('Stop count must be positive')
 out=args.out.resolve()
 if args.resume:
  if not out.is_dir():raise ValueError('Resume directory missing')
  if (out/'summary.json').exists():raise ValueError('Already finalized; do not overwrite')
 else:out.mkdir(parents=True,exist_ok=False)
 s=(ROOT/'scripts/74_polynomial_forcing_prototype.py').read_text().split("fixture=ROOT/'configs/")[0];s=s.replace("OUT=ROOT.parent/'run_20261030/polynomial_forcing';OUT.mkdir(parents=True,exist_ok=False)","OUT=None");ns={'__file__':str(ROOT/'scripts/74_polynomial_forcing_prototype.py')};exec(compile(s,'poly_forcing_definitions','exec'),ns)
 fixture=ROOT.parent/'run_20261101/poly_new_training/frozen_polynomial_pinns.json';data=fixture.read_bytes();fixture_sha=sha(data)
 if fixture_sha!=(fixture.parent/'fixture_sha256.txt').read_text().strip():raise ValueError('Fixture checksum mismatch')
 refpath=ROOT.parent/'run_20261101/poly_new_certification/references.json';refdata=refpath.read_bytes();refs=json.loads(refdata);nets=json.loads(data)['networks'];tols=[F(1,1000),F(1,100000),F(1,1000000)];plan=[(n['label'],degree,str(t)) for n in nets for degree in (2,3) for t in tols];mapping={n['label']:n for n in nets}
 if len(mapping)!=len(nets):raise ValueError('Duplicate network labels')
 config={'fixture_sha256':fixture_sha,'references_sha256':sha(refdata),'source_sha256':sha(Path(__file__).read_bytes()),'prototype_sha256':sha((ROOT/'scripts/74_polynomial_forcing_prototype.py').read_bytes()),'network_count':len(nets),'expected_count':len(plan),'tolerances':[str(t) for t in tols],'degrees':[2,3],'mesh':[16,32,64,128],'reference_usage':'loaded before loop, compared after decision only; not a blind test'}
 cfg=out/'config.json'
 if args.resume:
  if json.loads(cfg.read_text())!=config:raise ValueError('Resume inputs/config differ')
 else:atomic(cfg,canonical(config))
 journal=out/'records.jsonl';rows=load_journal(journal,plan);count_added=0
 with journal.open('ab') as jf:
  for idx in range(len(rows),len(plan)):
   name,degree,tau=plan[idx];net=mapping[name];tauF=F(tau);tick=time.perf_counter();trace=[];decision='unknown'
   p=ns['PolynomialForcing'](fixture,name,fixture_sha,net['source_coefficients']);v=ShapedVerifier(p,'moment') if degree==2 else CubicVerifier(p,'moment')
   for n in config['mesh']:
    v.uniform(n);r=v.certify();lo,hi=exact_fraction(r['lower']),exact_fraction(r['upper']);decision='within' if hi<=tauF else 'exceeds' if lo>tauF else 'unknown';trace.append({'intervals':n,'lower':str(lo),'upper':str(hi),'decision':decision})
    if decision!='unknown':break
   ref=refs[name];rl,ru=ns['fr'](ref['lower']),ns['fr'](ref['upper']);refdec='within' if ru<=tauF else 'exceeds' if rl>tauF else 'reference_unknown';ok=decision=='unknown' or (decision==refdec and ref.get('converged',False))
   row={'network':name,'source':net['source_name'],'seed':net['seed'],'depth':net['depth'],'degree':degree,'tolerance':tau,'decision':decision,'reference_decision':refdec,'reference_converged':ref.get('converged',False),'compatible':ok,'trace':trace,'elapsed_seconds':time.perf_counter()-tick,'error':None}
   jf.write(canonical(row));jf.flush();os.fsync(jf.fileno());rows.append(row);count_added+=1
   if args.simulate_stop_after is not None and count_added>=args.simulate_stop_after:
    print('SIMULATED_STOP',len(rows),flush=True);return 3
 # Read back journal and check against the full planned sequence before finalizing.
 rows=load_journal(journal,plan)
 if len(rows)!=len(plan):raise ValueError('Incomplete journal')
 records=(json.dumps(rows,ensure_ascii=False,indent=2,allow_nan=False)+'\n').encode();atomic(out/'records.json',records)
 summary={'records_count':len(rows),'expected_count':len(plan),'records_sha256':sha(records),'journal_sha256':sha(journal.read_bytes()),'outcomes':dict(Counter(x['decision'] for x in rows)),'compatible_count':sum(x['compatible'] for x in rows),'config_sha256':sha(cfg.read_bytes()),'durability_note':'fsync journal per row and atomically replace exports; not proof against all storage failure modes'};atomic(out/'summary.json',canonical(summary));verify_export(out);print('FINALIZED',json.dumps(summary),flush=True);return 0
if __name__=='__main__':raise SystemExit(main())
