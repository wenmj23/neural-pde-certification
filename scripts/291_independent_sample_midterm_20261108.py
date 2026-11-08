"""Sample analytic-solution errors using an independent mpmath forward pass."""
from fractions import Fraction as F
from pathlib import Path
from mpmath.libmp import to_rational
import hashlib,json,os,tempfile
import mpmath as mp
ROOT=Path(__file__).resolve().parents[1];CFG=ROOT/'configs/midterm_independent_sample_2026-11-08.json';B=ROOT/'results/2026-11-08';OUT=B/'midterm_independent_sample'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def as_fraction(x):return F(*to_rational(x._mpf_))
def mp_error(n,i,den,dps):
 with mp.workdps(dps):
  x=mp.mpf(i)/den;a=[x]
  for layer in n['layers']:
   weights=[[mp.mpf(v) for v in row] for row in layer['weight']];bias=[mp.mpf(v) for v in layer['bias']];a=[sum((v*z for v,z in zip(row,a)),b) for row,b in zip(weights,bias)]
   if layer['activation']=='tanh':a=[mp.tanh(z) for z in a]
   elif layer['activation']!='identity':raise ValueError('activation')
  if len(a)!=1:raise ValueError('output width')
  v=x*(1-x)*a[0];u=sum((mp.mpf(v)*(x-x**(j+2))/((j+1)*(j+2)) for j,v in enumerate(n['source_coefficients'])),mp.mpf(0))
  return +abs(v-u)
def main():
 if OUT.exists():raise FileExistsError(OUT)
 c=json.loads(CFG.read_text());t=ROOT/c['training_dir'];d=ROOT/c['certificates_dir'];fixture=t/'frozen_networks.json';refs=d/'references.json';cert=d/'records.json';ts=json.loads((t/'summary.json').read_text());cs=json.loads((d/'summary.json').read_text())
 if sha(fixture)!=c['fixture_sha256'] or sha(fixture)!=ts['fixture_sha256'] or sha(refs)!=cs['references_sha256'] or sha(cert)!=cs['records_sha256']:raise ValueError('input hash')
 networks=json.loads(fixture.read_text())['networks'];ref={q['network']:q for q in json.loads(refs.read_text())};certificates={(r['network'],r['degree'],r['intervals']):r for r in json.loads(cert.read_text())}
 if len(networks)!=12 or {n['label'] for n in networks}!=set(ref):raise ValueError('manifest')
 src={str(p.relative_to(ROOT)):sha(p) for p in (CFG,Path(__file__))};rows=[]
 for n in networks:
  q=ref[n['label']];rl,ru=F(q['lower']),F(q['upper']);den=max(c['grids']);values=[as_fraction(mp_error(n,i,den,c['candidate_digits'])) for i in range(den+1)]
  for grid in c['grids']:
   if den%grid:raise ValueError('non-nested grids')
   indexes=range(0,den+1,den//grid);imax=max(indexes,key=lambda i:values[i]);v80=values[imax];v110=as_fraction(mp_error(n,imax,den,c['check_digits']));stable=abs(v80-v110)<=F(c['relative_stability'])*max(F(1),abs(v110));bounds={}
   for degree in c['degrees']:
    certrow=certificates[(n['label'],degree,c['certification_mesh'])];L,U=F(certrow['lower']),F(certrow['upper']);bounds[str(degree)]={'lower':str(L),'upper':str(U),'sample_below_upper':v110<=U,'sample_above_lower':v110>=L}
   rows.append({'network':n['label'],'source':n['source_name'],'grid_intervals':grid,'sample_x':str(F(imax,den)),'sample_error_80':str(v80),'sample_error_110':str(v110),'stable':stable,'reference_lower':str(rl),'reference_upper':str(ru),'sample_below_reference_upper':v110<=ru,'bounds':bounds})
  print(n['label'],'completed',flush=True)
 B.mkdir(parents=True,exist_ok=True);stage=Path(tempfile.mkdtemp(prefix='.sample_midterm_',dir=B));(stage/'locked_config.json').write_bytes(CFG.read_bytes());(stage/'records.json').write_text(json.dumps(rows,indent=2)+'\n');s={'expected':24,'records':len(rows),'stable':sum(r['stable'] for r in rows),'below_reference_upper':sum(r['sample_below_reference_upper'] for r in rows),'below_certificate_upper':sum(r['bounds'][str(d)]['sample_below_upper'] for r in rows for d in c['degrees']),'inputs_sha256':{'fixture':sha(fixture),'references':sha(refs),'certificates':sha(cert)},'source_hashes_before_run':src,'records_sha256':sha(stage/'records.json'),'scope':c['scope']};(stage/'summary.json').write_text(json.dumps(s,indent=2)+'\n');os.rename(stage,OUT);print(json.dumps(s));return int(s['records']!=24 or s['stable']!=24 or s['below_reference_upper']!=24 or s['below_certificate_upper']!=48)
if __name__=='__main__':raise SystemExit(main())
