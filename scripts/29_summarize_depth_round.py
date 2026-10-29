"""Verify completed run manifests and generate a compact depth-round digest."""
from datetime import datetime,timezone
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import statistics
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parents[1]


def latest(prefix):
    completed=[p.parent for p in (ROOT/'results/generated'/prefix).glob('*/config.json') if json.loads(p.read_text()).get('status')=='complete']
    if not completed:raise RuntimeError('Missing completed study: '+prefix)
    return max(completed)


def main():
    names=['23_depth_derivatives','24_depth_reference','25_depth_certification','26_depth_interior_audit','27_depth_refinements','28_cubic_depth','30_cubic_adaptivity','31_budget_sensitivity']
    runs={name:latest(name) for name in names};configs={name:json.loads((p/'config.json').read_text()) for name,p in runs.items()}
    validation=[]
    for name,p in runs.items():
        for file,sha in configs[name]['sources'].items():
            if hashlib.sha256((ROOT/file).read_bytes()).hexdigest()!=sha:raise AssertionError('Source changed: '+file)
        if (p/'checks.json').exists():
            checks=json.loads((p/'checks.json').read_text())
            if len(checks)!=configs[name]['checks_passed'] or not all(x['passed'] for x in checks):raise AssertionError('Failed checks in '+name)
        validation.append(dict(study=name,source_hashes_match=True,completed=True))
    out=ROOT/'results/generated/29_depth_summary'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ');out.mkdir(parents=True,exist_ok=False)
    primary=json.loads((runs['25_depth_certification']/'summary.json').read_text())
    follow=json.loads((runs['27_depth_refinements']/'summary.json').read_text())
    cubic=json.loads((runs['28_cubic_depth']/'summary.json').read_text())
    fixture=json.loads((ROOT/'configs/frozen_networks_depth_2026-10-29.json').read_text());labels=[r['label'] for r in fixture['networks']]
    def pfind(label,derivative,reconstruction):return next(r for r in primary if r['label']==label and r['derivatives']==derivative and r['reconstruction']==reconstruction)
    def cfind(label,method,target='51/50'):return next(r for r in cubic if r['label']==label and r['method']==method and r['target']==target)
    quad_ratios=[pfind(label,'chain','linear')['median_wall_ms']/pfind(label,'chain','quadratic_moment')['median_wall_ms'] for label in labels]
    cubic_ratios=[cfind(label,'quadratic_moment')['median_wall_ms']/cfind(label,'cubic_moment')['median_wall_ms'] for label in labels]
    alternatives={}
    for method in ['centered','value_centered','centered_lite','split']:
        ratios=[];wins=[]
        for label in labels:
            base=next(r for r in follow if r['label']==label and r['derivatives']=='chain')
            candidate=next(r for r in follow if r['label']==label and r['derivatives']==method)
            ratio=base['median_wall_ms']/candidate['median_wall_ms'];ratios.append(ratio)
            if ratio>1.05:wins.append(dict(label=label,speed_ratio=ratio))
        alternatives[method]=dict(median_speed_ratio=statistics.median(ratios),wins_over_5percent=wins)
    totals=dict(networks=len(labels),checks=sum(c.get('checks_passed',0) for c in configs.values()),
                fixed_grid_configurations=configs['25_depth_certification']['fixed_records']+configs['27_depth_refinements']['fixed_configurations']+configs['28_cubic_depth']['fixed_configurations'],
                timed_requests=sum(c.get('timed_requests',0) for c in configs.values()),successful_requests=sum(c.get('successful_requests',0) for c in configs.values()))
    def stats(ratios):return dict(min=min(ratios),median=statistics.median(ratios),max=max(ratios),faster_cases=sum(x>1 for x in ratios),cases=len(ratios))
    digest=dict(totals=totals,quadratic_vs_linear=stats(quad_ratios),cubic_vs_quadratic=stats(cubic_ratios),derivative_alternatives=alternatives,
                studies={k:str(p.relative_to(ROOT)) for k,p in runs.items()},manifest_validation=validation,
                selected_cubic=[cfind(label,method) for label in ['depth1_k1_seed7','depth2_k1_seed7','depth3_k1_seed7','depth3_k5_seed7'] for method in ['quadratic_moment','cubic_constant','cubic_moment']],
                elapsed_since_first_workspace_check_min=(datetime.now(timezone.utc)-datetime.fromisoformat('2026-10-29T04:41:22+00:00')).total_seconds()/60,
                caution='All ratios compare observed three-repeat medians, not population speedups; exact scope and trust conditions remain in source study configs.')
    (out/'digest.json').write_text(json.dumps(digest,indent=2)+'\n')
    fig,axes=plt.subplots(1,2,figsize=(12,4.7))
    for ax,k in zip(axes,(1,5)):
        selected=[r for r in fixture['networks'] if r['frequency']==k];x=np.arange(len(selected))
        for method,label in [('linear','Linear'),('quadratic_constant','Quadratic / constant'),('quadratic_moment','Quadratic / shaped')]:
            ys=[pfind(r['label'],'chain',method)['median_wall_ms'] for r in selected]
            ax.plot(x,ys,'-o',label=label,markersize=4)
        ax.set_yscale('log');ax.set_xticks(x,[f'D{r["depth"]}/S{r["seed"]}' for r in selected],rotation=45,ha='right',fontsize=8)
        ax.set(title=f'Poisson sine frequency {k}',ylabel='Complete request median (ms, log scale)',xlabel='Hidden layers / seed');ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.suptitle('Same frozen output, same gap target U/L <= 1.02 (study25)',fontsize=11);fig.tight_layout();fig.savefig(out/'depth_reconstruction_cost.png',dpi=160);plt.close(fig)
    fig,ax=plt.subplots(figsize=(12,4.6));x=np.arange(len(labels));colors=['#3b82f6' if r['depth']==1 else '#059669' if r['depth']==2 else '#d97706' for r in fixture['networks']]
    low=[];high=[]
    for label,ratio in zip(labels,cubic_ratios):
        q=cfind(label,'quadratic_moment');c=cfind(label,'cubic_moment')
        low.append(ratio-q['min_wall_ms']/c['max_wall_ms']);high.append(q['max_wall_ms']/c['min_wall_ms']-ratio)
    ax.bar(x,cubic_ratios,color=colors,alpha=.85);ax.errorbar(x,cubic_ratios,yerr=[low,high],fmt='none',ecolor='#374151',capsize=2,linewidth=1)
    ax.axhline(1,color='black',linewidth=1);ax.set_xticks(x,[f'D{r["depth"]}/k{r["frequency"]}/S{r["seed"]}' for r in fixture['networks']],rotation=60,ha='right',fontsize=8)
    ax.set(ylabel='Quadratic time / cubic time',title='Cubic versus quadratic shaped reconstruction (study28)',xlabel='Same U/L <= 1.02; bars use 3-run medians; whiskers use observed min/max combinations')
    ax.grid(axis='y',alpha=.2);fig.tight_layout();fig.savefig(out/'cubic_speed_comparison.png',dpi=160);plt.close(fig)
    (out/'config.json').write_text(json.dumps(dict(status='complete',source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),studies=digest['studies']),indent=2)+'\n')
    print('DIGEST',json.dumps(digest,indent=2));print('COMPLETED',out)


if __name__=='__main__':main()
