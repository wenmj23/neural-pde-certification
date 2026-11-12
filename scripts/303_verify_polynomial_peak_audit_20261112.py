"""Post-exit manifests, hash and exact stored-decimal peak checks."""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parents[1]
D=ROOT/'results/2026-11-12/polynomial_peak_audit'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def need(q,m):
    if not q:raise AssertionError(m)
def main():
    c=json.loads((D/'locked_config.json').read_text());s=json.loads((D/'summary.json').read_text())
    paths={k:ROOT/c[k+'_path'] for k in ['fixture','references','spots']}
    for k,p in paths.items():need(sha(p)==c[k+'_sha256']==s['inputs'][k],k+' hash')
    need(sha(D/'locked_config.json')==s['config_sha256'],'config hash')
    need(sha(ROOT/'scripts/302_polynomial_peak_audit_20261112.py')==s['script_sha256'],'producer hash')
    need(sha(ROOT/'scripts/296_green_quadrature_audit_20261112.py')==s['helper_sha256'],'helper hash')
    need(sha(D/'records.json')==s['records_sha256'],'records hash')
    nets=json.loads(paths['fixture'].read_text())['networks'];refs={r['network']:r for r in json.loads(paths['references'].read_text())}
    rows=json.loads((D/'records.json').read_text())
    need(len(rows)==12 and {r['network'] for r in rows}=={n['label'] for n in nets},'row manifest')
    need(len(s['cases'])==12,'case manifest')
    rootdiff=[];green=[]
    for i,row in enumerate(rows):
        name=f'case_{i:02d}.json';need(name in s['cases'] and sha(D/name)==s['cases'][name],'case hash')
        need(json.loads((D/name).read_text())==row,'case contents')
        need([r['digits'] for r in row['roots']]==c['root_digits'],'root digits')
        need([q['digits'] for q in row['quadrature']]==c['quadrature_digits'],'quadrature digits')
        a,b=row['roots'];r=refs[row['network']]
        need(F(row['reference_lower'])==F(r['lower']) and F(row['reference_upper'])==F(r['upper']),'reference bounds')
        within=F(r['lower'])<=abs(F(b['error']))<=F(r['upper'])
        need(row['numerical_peak_within_reference']==within,'reference flag')
        need(b['local_abs_max_curvature_test']==(F(b['error'])*F(b['curvature'])<0),'curvature flag')
        rootdiff.append(abs(F(a['root'])-F(b['root'])))
        for q in row['quadrature']:
            d=abs(F(q['green'])-F(q['direct_error']));need(d<F(1,10**30),'Green/direct discrepancy');green.append(d)
        need(F(b['abs_derivative'])<F(1,10**65),'stationary derivative')
    result={'records':len(rows),'case_file_checks':len(s['cases']),
      'numerical_peaks_within_reference':sum(r['numerical_peak_within_reference'] for r in rows),
      'local_max_curvature_checks':sum(r['roots'][-1]['local_abs_max_curvature_test'] for r in rows),
      'max_root_cross_precision_difference':float(max(rootdiff)),
      'max_green_direct_discrepancy_as_saved':float(max(green)),
      'readback_script_sha256':sha(Path(__file__)),
      'scope':'post-exit hash and stored arithmetic for selected numerical roots; no proof of global maxima'}
    p=D/'independent_readback.json'
    if p.exists():raise FileExistsError(p)
    p.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
if __name__=='__main__':main()
