"""Independently read back script 111's point and full-interval result bundle."""
from fractions import Fraction as F
from pathlib import Path
import argparse,hashlib,importlib.util,json,sys

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 ap=argparse.ArgumentParser();ap.add_argument('directory',type=Path);ap.add_argument('--fixtures',nargs='+',required=True,type=Path);a=ap.parse_args();p=a.directory
 spec=importlib.util.spec_from_file_location('guard',Path(__file__).with_name('108_verify_result_bundle_20261101.py'));guard=importlib.util.module_from_spec(spec);spec.loader.exec_module(guard)
 base=guard.verify(p,a.fixtures)
 if not base['ok']:raise ValueError('point result guard: '+str(base['errors']))
 s=json.loads((p/'summary.json').read_text());rows=json.loads((p/'intervals.json').read_text());points=json.loads((p/'records.json').read_text())
 if sha(p/'intervals.json')!=s['intervals_sha256']:raise ValueError('interval hash mismatch')
 fixtures={(sha(path),net['label']) for path in a.fixtures for net in json.loads(path.read_bytes())['networks']}
 cells=[('0','1/32'),('7/32','1/4'),('15/32','1/2'),('31/32','1')];orders=[0,2,3,4]
 expected={(h,n,c,j) for h,n in fixtures for c in cells for j in orders}
 keys=[(r['fixture_sha256'],r['network'],tuple(r['cell']),r['order']) for r in rows]
 if len(keys)!=len(expected) or set(keys)!=expected or len(keys)!=s['interval_expected'] or len(keys)!=s['interval_records']:raise ValueError('interval key coverage mismatch')
 sign_detection=0
 for r in rows:
  lo,hi,il,iu=map(F,(r['arb_lower'],r['arb_upper'],r['mpiv_lower'],r['mpiv_upper']))
  overlap=lo<=iu and il<=hi
  if lo>hi or il>iu or overlap!=r['overlap'] or not overlap:raise ValueError('interval bounds or overlap mismatch')
  sign_detection+=not (lo<=-il and -iu<=hi)
 if s['interval_overlap']!=sum(r['overlap'] for r in rows):raise ValueError('interval count mismatch')
 if s['mpiv_contained']!=sum(r['mpiv90_contained'] is True for r in points) or s['mpiv_contained']!=len(points):raise ValueError('point mpiv containment mismatch')
 print(json.dumps({'ok':True,'point_records':len(points),'interval_records':len(rows),'interval_overlap':s['interval_overlap'],'record_sha256':s['records_sha256'],'interval_sha256':s['intervals_sha256'],'sign_flip_detected':sign_detection,'scope':'file integrity and consistency, not proof of numerical correctness'}));return 0
if __name__=='__main__':sys.exit(main())
