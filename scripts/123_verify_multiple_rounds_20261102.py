"""Post-exit integrity/consistency checker for rounds 118-122.

Rechecks stored arithmetic and request manifests. It does not rerun numerical
or interval calculations or authenticate self-consistent modified results.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,json,sys
SOURCES={'quadratic_12x2':[0,0,12],'cubic_20x3':[0,0,0,20],'quintic_42x5':[0,0,0,0,0,42]}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def need(v,msg):
 if not v:raise ValueError(msg)
def data(root,name):
 p=root/name;s=json.loads((p/'summary.json').read_text());need(sha(p/'records.json')==s['records_sha256'],name+' hash mismatch');rows=json.loads((p/'records.json').read_text());need(len(rows)==s['records'],name+' count mismatch');need(not s['failures'],name+' reported failures');return p,s,rows
def key(r):return (r['fixture_sha256'],r['network'],r['source'])
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--root',required=True,type=Path);ap.add_argument('--fixtures',nargs='+',required=True,type=Path);a=ap.parse_args();root=a.root
 labels={(sha(p),n['label']) for p in a.fixtures for n in json.loads(p.read_bytes())['networks']};need(len(labels)==24,'unexpected fixture labels')
 p,s,rows=data(root,'transferred_forcing');need(s['sources']==SOURCES,'source manifest');need(sha(p/'references.json')==s['references_sha256'],'ref hash mismatch');refs=json.loads((p/'references.json').read_text());rmap={key(x):x for x in refs};basekeys={(h,n,src) for h,n in labels for src in SOURCES};need(len(refs)==len(rmap)==72 and set(rmap)==basekeys,'reference keys')
 allkeys={(h,n,src,d) for h,n,src in basekeys for d in (2,3)};rkeys=[(*key(r),r['degree']) for r in rows];need(len(rows)==len(allkeys)==144 and set(rkeys)==allkeys,'certificate keys');need(s['expected_records']==144 and s['expected_references']==72,'certificate expected')
 for r in rows:
  q=rmap[key(r)];lo,hi,rl,ru=map(F,(r['lower'],r['upper'],q['lower'],q['upper']));need(0<=lo<=hi and 0<=rl<=ru and r['intervals']==64,'bounds');need(r['reference_converged']==q['converged'] is True,'ref convergence');need(r['reference_overlap']==(lo<=ru and rl<=hi),'overlap');need(r['reference_contained']==(lo<=rl and ru<=hi),'containment');need(r['gap_1_02']==(lo>0 and hi<=F(51,50)*lo),'gap')
 need(s['reference_converged']==sum(r['converged'] for r in refs),'reference total');need(s['reference_overlap']==sum(r['reference_overlap'] for r in rows),'overlap total');need(s['reference_contained']==sum(r['reference_contained'] for r in rows),'containment total')
 for d in (2,3):need(s['gap_1_02'][str(d)]==sum(r['gap_1_02'] for r in rows if r['degree']==d),'gap total')
 original_rows=rows;original_s=s
 p,s,scale=data(root,'metamorphic_scaling');need(s['expected']==384,'scale expected');expected={(h,n,m,d,c) for h,n in labels for m in (32,64) for d in (2,3) for c in (1,-1,2,-2)};ks=[(r['fixture_sha256'],r['network'],r['mesh'],r['degree'],r['scale']) for r in scale];need(len(ks)==len(expected)==384 and set(ks)==expected,'scale keys')
 if len(a.fixtures)!=len(s['derived_fixture_sha256']):raise ValueError('derived fixture count')
 for i,(path,expected_sha) in enumerate(zip(a.fixtures,s['derived_fixture_sha256'])):
  d=p/f'derived_fixture_{i}.json';need(sha(d)==expected_sha,'derived fixture hash');derived={x['label']:x for x in json.loads(d.read_bytes())['networks']};orig=json.loads(path.read_bytes())['networks'];need(len(derived)==len(orig)*4,'derived count')
  for n in orig:
   for c in (1,-1,2,-2):
    t=derived[n['label']+'__scale_'+str(c)];need(t['source_coefficients']==[v*c for v in n['source_coefficients']],'derived source')
    for x,y in zip(n['layers'][-1]['weight'][0],t['layers'][-1]['weight'][0]):need(F.from_float(float(y))==F.from_float(float(x))*c,'derived weight')
    need(F.from_float(float(t['layers'][-1]['bias'][0]))==F.from_float(float(n['layers'][-1]['bias'][0]))*c,'derived bias');need(t['layers'][:-1]==n['layers'][:-1],'hidden layers changed')
 baseline={(r['fixture_sha256'],r['network'],r['mesh'],r['degree']):(F(r['normalized_lower']),F(r['normalized_upper'])) for r in scale if r['scale']==1}
 for r in scale:
  l,u,nl,nu,bl,bu=map(F,(r['lower'],r['upper'],r['normalized_lower'],r['normalized_upper'],r['baseline_lower'],r['baseline_upper']));need((nl,nu)==(l/abs(r['scale']),u/abs(r['scale'])),'scale normalization');need((bl,bu)==baseline[(r['fixture_sha256'],r['network'],r['mesh'],r['degree'])],'scale baseline');need(r['equal_to_baseline']==((nl,nu)==(bl,bu)),'scale equality flag');need(r['overlap_baseline']==(nl<=bu and bl<=nu),'scale overlap flag')
 need(s['equal']==sum(r['equal_to_baseline'] for r in scale) and s['overlap']==sum(r['overlap_baseline'] for r in scale),'scale totals')
 p,s,tight=data(root,'tightened_reference');need(s['input_records_sha256']==original_s['records_sha256'] and s['input_references_sha256']==original_s['references_sha256'],'tight inputs');selected={key(r) for r in original_rows if not r['reference_contained']};expected={(*k,d) for k in selected for d in (2,3)};keys=[(*key(r),r['degree']) for r in tight];need(len(keys)==len(expected)==s['expected'] and set(keys)==expected,'tight keys');origmap={(*key(r),r['degree']):r for r in original_rows}
 for r in tight:
  q=origmap[(*key(r),r['degree'])];old=rmap[key(r)];need(r['certificate_lower']==q['lower'] and r['certificate_upper']==q['upper'] and r['old_reference_lower']==old['lower'] and r['old_reference_upper']==old['upper'],'tight parent mismatch');lo,hi,rl,ru,ol,ou=map(F,(r['certificate_lower'],r['certificate_upper'],r['new_reference_lower'],r['new_reference_upper'],r['old_reference_lower'],r['old_reference_upper']));need(r['reference_converged'] is True and r['reference_contained']==(lo<=rl and ru<=hi) and r['reference_overlap']==(lo<=ru and rl<=hi) and r['new_overlaps_old']==(rl<=ou and ol<=ru),'tight flags')
 for x,y in [('converged','reference_converged'),('contained','reference_contained'),('overlap','reference_overlap'),('old_new_overlap','new_overlaps_old')]:need(s[x]==sum(r[y] for r in tight),'tight totals')
 p,s,dec=data(root,'decision_boundaries');need(s['input_sha256']==original_s['records_sha256'] and s['source_certificates']==144,'decision input');need(len(dec)==723,'decision count');seen=set()
 for r in dec:
  lo,hi,t=map(F,(r['lower'],r['upper'],r['tolerance']));expected='within_tolerance' if hi<=t else ('exceeds_tolerance' if lo>t else 'undetermined');need(r['expected']==r['actual']==expected and r['goal_met']==(expected!='undetermined') and r['passed'] is True,'decision flag')
  if r['input_index']!='synthetic':
   i=r['input_index'];need(type(i) is int and 0<=i<len(original_rows),'decision index');q=original_rows[i];need(lo==F(q['lower']) and hi==F(q['upper']),'decision bound');kind=r['kind'];targets={'at_lower':lo,'midpoint':(lo+hi)/2,'at_upper':hi,'above_upper':hi+(hi-lo)/4,'below_lower':lo/2};need(kind in targets and t==targets[kind], 'decision target');need((i,kind) not in seen,'duplicate decision');seen.add((i,kind))
 need(len(seen)==720 and s['passed']==723,'decision coverage')
 p,s,der=data(root,'transferred_derivatives');expected={(h,n,src,i,d) for h,n in labels for src in SOURCES for i in range(4) for d in (0,2,3,4)};keys=[(*key(r),r['cell_index'],r['order']) for r in der];need(len(keys)==len(expected)==s['expected'] and set(keys)==expected,'derivative keys')
 for r in der:
  lo,hi,il,iu,low,high=map(F,(r['arb_lower'],r['arb_upper'],r['mpiv_lower'],r['mpiv_upper'],r['mp50'],r['mp70']));need(lo<=hi and il<=iu,'derivative endpoints');checks={'overlap':lo<=iu and il<=hi,'stable':abs(low-high)<=F(1,10**40)*max(F(1),abs(high)),'point_in_arb':lo<=high<=hi,'point_in_mpiv':il<=high<=iu}
  for k,v in checks.items():need(r[k] is v and v,'derivative '+k)
 for k in ('overlap','stable','point_in_arb','point_in_mpiv'):need(s[k]==sum(r[k] for r in der),'derivative totals')
 print(json.dumps({'ok':True,'transfer_records':len(original_rows),'references':len(refs),'scaled_records':len(scale),'tightened_records':len(tight),'boundary_records':len(dec),'derivative_records':len(der),'scope':'readback integrity and stored arithmetic; not mathematical proof or authentication'}));return 0
if __name__=='__main__':
 try:sys.exit(main())
 except (ValueError,KeyError,TypeError,OSError,ZeroDivisionError) as e:print(json.dumps({'ok':False,'error':str(e)}));sys.exit(1)
