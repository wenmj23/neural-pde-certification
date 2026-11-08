"""Paired same-network, mesh, tolerance decisions, using checked runs."""
from collections import Counter
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parents[1];B=ROOT/'results/2026-11-07';OUT=B/'midterm_decision_pairwise_readback.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 if OUT.exists():raise FileExistsError(OUT)
 out={}
 for name in ('poly','sine'):
  d=B/f'midterm_{name}_decisions';p=d/'records.json';c=json.loads((d/'locked_config.json').read_text());audit=json.loads((B/f'midterm_{name}_decisions_readback.json').read_text());rows=json.loads(p.read_text())
  if sha(p)!=audit['records_sha256']:raise ValueError('hash')
  idx={(r['network'],r['mesh'],r['degree']):r for r in rows};counts=Counter();by_mesh=Counter();examples=[]
  for r in rows:
   if r['degree']!=2:continue
   z=idx[(r['network'],r['mesh'],3)]
   for t in c['tolerances']:
    a,b=r['decisions'][t],z['decisions'][t]
    if a!='unknown' and b!='unknown':kind='both_agree' if a==b else 'contradiction'
    elif a=='unknown' and b=='unknown':kind='both_unknown'
    else:kind='only_quadratic' if b=='unknown' else 'only_cubic'
    counts[kind]+=1;by_mesh[(r['mesh'],kind)]+=1
    if kind in ('only_cubic','only_quadratic') and len(examples)<4:examples.append({'network':r['network'],'mesh':r['mesh'],'tolerance':t,'quadratic':a,'cubic':b})
  out[name]={'paired_settings':sum(counts.values()),'counts':dict(counts),'by_mesh':[{'mesh':k[0],'outcome':k[1],'count':v} for k,v in sorted(by_mesh.items())],'examples':examples,'records_sha256':sha(p)}
  if counts['contradiction'] or sum(counts.values())!=(288 if name=='poly' else 144):raise ValueError('pair mismatch')
 out['scope']='same previously studied networks, fixed meshes and thresholds; not independent networks or matched timing'
 with OUT.open('x') as f:json.dump(out,f,indent=2);f.write('\n')
 print(json.dumps(out,indent=2))
if __name__=='__main__':main()
