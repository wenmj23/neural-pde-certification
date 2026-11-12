"""Post-exit hash/manifest/readback checks for 24-network diagnostics."""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parents[1]
D=ROOT/'results/2026-11-12/extra_polynomial_green_audit'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def need(c,m):
    if not c:raise AssertionError(m)
def main():
    c=json.loads((D/'locked_config.json').read_text());s=json.loads((D/'summary.json').read_text())
    need(sha(D/'locked_config.json')==s['config_sha256'],'config')
    need(sha(ROOT/'scripts/304_extra_polynomial_green_audit_20261112.py')==s['script_sha256'],'producer')
    need(sha(ROOT/'scripts/296_green_quadrature_audit_20261112.py')==s['helper_sha256'],'helper')
    need(sha(D/'records.json')==s['records_sha256'],'records')
    nets=[];bounds={}
    for g in c['groups']:
        fp,op=ROOT/g['fixture_path'],ROOT/g['records_path']
        need(sha(fp)==g['fixture_sha256'] and sha(op)==g['records_sha256'],'input')
        ns=json.loads(fp.read_text())['networks'];nets.extend(ns)
        bounds.update({r['network']:r for r in json.loads(op.read_text()) if r['degree']==c['degree'] and r['intervals']==c['mesh']})
    labels=[n['label'] for n in nets];rows=json.loads((D/'records.json').read_text())
    keys=[(r['network'],r['point']) for r in rows]
    need(len(nets)==24 and len(rows)==72 and len(set(keys))==72 and set(keys)=={(n,p) for n in labels for p in c['points']},'manifest')
    need(len(s['case_sha256'])==72,'case count')
    for ni,n in enumerate(labels):
        for pi,p in enumerate(c['points']):
            name=f'case_{ni:02d}_{pi:02d}.json'
            need(sha(D/name)==s['case_sha256'][name],'case hash')
            need(json.loads((D/name).read_text())==next(r for r in rows if (r['network'],r['point'])==(n,p)),'case data')
    green=[];residual=[]
    for row in rows:
        need([z['digits'] for z in row['results']]==c['digits'],'precision')
        for z in row['results']:
            d=abs(F(z['green'])-F(z['direct_error']));need(d<F(1,10**30),'Green mismatch');green.append(d)
        ctrl=row['source_control']
        delta=abs(F(ctrl['correct_residual'])-F(ctrl['independent_minus_error_second_derivative']))
        need(delta<F(1,10**65),'residual mismatch');residual.append(delta)
        injected=F(ctrl['injected_omission_residual']);correct=F(ctrl['independent_minus_error_second_derivative'])
        need(ctrl['injected_fault_nonzero']==(abs(injected-correct)>F(1,10**50)),'fault flag')
        u=F(bounds[row['network']]['upper']);need(F(row['saved_cubic_upper'])==u,'old upper')
        need(row['direct_below_saved_cubic_upper']==(abs(F(row['results'][-1]['direct_error']))<=u),'upper flag')
    result={'records':len(rows),'case_file_checks':len(s['case_sha256']),
      'point_precision_checks':len(green),'point_residual_checks':len(residual),
      'below_saved_upper':sum(r['direct_below_saved_cubic_upper'] for r in rows),
      'injected_omission_detected':sum(r['source_control']['injected_fault_nonzero'] for r in rows),
      'max_green_direct_difference_as_saved':float(max(green)),
      'max_residual_difference_as_saved':float(max(residual)),
      'readback_script_sha256':sha(Path(__file__)),
      'scope':'post-exit hashes and stored-decimal arithmetic; injected fault is a controlled test, not a production defect; finite samples do not prove continuous coverage'}
    p=D/'independent_readback.json'
    if p.exists():raise FileExistsError(p)
    p.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
if __name__=='__main__':main()
