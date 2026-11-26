import json
from fractions import Fraction as F
from pathlib import Path
p=Path(__file__).resolve().parent;rows=json.loads((p/'records.json').read_text());refs={x['network']:x for x in json.loads((p/'inputs/references.json').read_text())};out=[]
for x in rows:
 r=refs[x['network']];rl=F(int(r['lower']['numerator']),int(r['lower']['denominator']));ru=F(int(r['upper']['numerator']),int(r['upper']['denominator']));l,u=F(x['lower']),F(x['upper']);out.append({'network':x['network'],'degree':x['degree'],'mode':x['mode'],'budget':x['budget'],'reference_fully_contained':l<=rl and ru<=u,'overlap':max(l,rl)<=min(u,ru),'certificate_contained_in_reference':rl<=l and u<=ru})
(p/'overlap_records.json').write_text(json.dumps(out,indent=2)+'\n');print('rows',len(out),'reference-contained',sum(x['reference_fully_contained'] for x in out),'overlap',sum(x['overlap'] for x in out),'non-contained',sum(not x['reference_fully_contained'] for x in out))
