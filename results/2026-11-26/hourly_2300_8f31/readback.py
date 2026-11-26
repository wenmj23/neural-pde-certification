from pathlib import Path
from fractions import Fraction as F
from collections import Counter
import json,hashlib,itertools
D=Path(__file__).resolve().parent;sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();s=json.loads((D/'summary.json').read_text());c=json.loads((D/'protocol.json').read_text());rows=json.loads((D/'records.json').read_text());nets=json.loads((D/'inputs/frozen_networks.json').read_text())['networks'];refs={r['network']:r for r in json.loads((D/'inputs/references.json').read_text())}
assert sha(D/'inputs/frozen_networks.json')==s['fixture_sha256'] and sha(D/'protocol.json')==s['protocol_sha256'] and sha(D/'records.json')==s['records_sha256'] and sha(D/'run.py')==s['run_sha256']
for n,h in s['sources'].items():assert sha(D/'sources'/n)==h
for n,h in s['cases'].items():assert sha(D/n)==h
names=[n['label'] for n in nets];assert len(rows)==96 and {(r['network'],r['degree'],r['mode'],r['budget']) for r in rows}==set(itertools.product(names,c['degrees'],c['mesh_modes'],c['budgets']))
counts=Counter();by={}
for r in rows:
 assert r['status']=='finite';l,u,zl,zu,du=map(F,[r['lower'],r['upper'],r['z_lower'],r['z_upper'],r['remainder_upper']]);assert 0<=l<=u and l<=max(F(0),zl-du) and u>=zu+du;nodes=list(map(F,r['nodes']));assert len(nodes)==r['budget']+1 and nodes[0]==0 and nodes[-1]==1 and all(a<b for a,b in zip(nodes,nodes[1:]));q=refs[r['network']];rl=F(int(q['lower']['numerator']),int(q['lower']['denominator']));ru=F(int(q['upper']['numerator']),int(q['upper']['denominator']));counts['reference_contained']+=l<=rl<=ru<=u
 for t in ['1/10000','1/1000000','1/10000000']:
  tol=F(t);counts[t+':'+('within' if u<=tol else 'exceeds' if l>tol else 'unknown')]+=1
 by[r['network'],r['degree'],r['mode'],r['budget']]=r
pairs=[]
for n,d,b in itertools.product(names,c['degrees'],c['budgets']):
 u=by[n,d,'uniform',b];a=by[n,d,'adaptive',b];uw=F(u['upper'])-F(u['lower']);aw=F(a['upper'])-F(a['lower']);pairs.append(dict(network=n,degree=d,budget=b,winner='adaptive' if aw<uw else 'uniform' if uw<aw else 'equal',nodes_different=a['nodes']!=u['nodes'],adaptive_derivative_calls=a['derivative_interval_requests'],uniform_derivative_calls=u['derivative_interval_requests']))
out={'records':len(rows),'networks':len(names),'counts':dict(counts),'pairs':len(pairs),'winners':dict(Counter(x['winner'] for x in pairs)),'nodes_different':sum(x['nodes_different'] for x in pairs),'adaptive_more_derivative_calls':sum(x['adaptive_derivative_calls']>x['uniform_derivative_calls'] for x in pairs),'groups':{f'{d}_{b}':dict(Counter(x['winner'] for x in pairs if x['degree']==d and x['budget']==b)) for d,b in itertools.product(c['degrees'],c['budgets'])}}
(D/'pairs.json').write_text(json.dumps(pairs,indent=2)+'\n');(D/'readback.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
