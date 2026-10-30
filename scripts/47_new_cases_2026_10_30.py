"""Pre-lock existing selectors, then generate new seed49 test networks.

Training only creates fixed verification inputs; keep every final network,
including low-quality ones. Continuous references are validation only.
Never overwrite an existing canonical fixture or prior results.
"""
from datetime import datetime, timezone
from fractions import Fraction as F
from pathlib import Path
import hashlib, importlib.util, json, platform, time
import torch

ROOT=Path(__file__).resolve().parents[1]
CANONICAL=ROOT/'configs/frozen_networks_holdout_2026-10-30.json'

def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m

def main():
    if CANONICAL.exists(): raise FileExistsError('Do not overwrite frozen inputs')
    trainer=load('training22',ROOT/'scripts/22_freeze_depth_networks.py')
    reference=load('reference24',ROOT/'scripts/24_depth_continuous_reference.py')
    out=ROOT/'results/generated/47_new_cases'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ');out.mkdir(parents=True,exist_ok=False)
    names=['arb_mesh_check_forecast.py','arb_mesh_check_stride.py','arb_mesh_selection.py',
           'arb_cubic_reconstruction.py','arb_reconstruction.py','arb_shaped_remainder.py',
           'arb_deep_network.py','arb_certification.py','arb_problem_families.py','arb_degree_selection.py',
           'interval_certification.py','18_continuous_reference_crosscheck.py','24_depth_continuous_reference.py',
           '22_freeze_depth_networks.py','mesh_problem_families.py','48_holdout_study.py']
    sources=[Path(__file__),*[ROOT/'scripts'/n for n in names]]
    snapshot=out/'source_snapshot';snapshot.mkdir()
    for p in sources:(snapshot/p.name).write_bytes(p.read_bytes())
    config=dict(status='running',seed=49,depths=[1,3],frequencies=[2,5],width=16,steps=2000,
                optimizer='Adam',learning_rate=.01,collocation_points=128,dtype='float64',threads=1,
                torch=torch.__version__,python=platform.python_version(),
                shapes={'shift_034_sigma003':[dict(center='0.34',width='0.003',weight='1')],
                        'shift_031_sigma040':[dict(center='0.31',width='0.04',weight='1')],
                        'double_021_068':[dict(center='0.21',width='0.009',weight='2/5'),dict(center='0.68',width='0.014',weight='3/5')]},
                policies=['cubic_uniform16','auto_heap16','nested_left125_16','nested_left_stride2_125_16',
                          'nested_left_stride3_125_16','nested_left_fcguard_125_16','nested_left_fcfast_125_16'],
                targets=['51/50','201/200'],repeats=3,max_intervals=4096,seconds=30,
                source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
                scope='Frozen seed49/k2,k5 networks generated after source/parameters fixed, same training and PDE family. All outputs retained. New analytic peaks are manufactured, not PINNs.')
    def save(name,obj):
        p=out/name;p.write_text(json.dumps(obj,indent=2)+'\n');assert json.loads(p.read_text())==obj
    save('config.json',config)
    torch.set_default_dtype(torch.float64);torch.set_num_threads(1)
    rows=[];refs={}
    try:
        for depth in (1,3):
            for k in (2,5):
                torch.manual_seed(49);model=trainer.PINN(depth);opt=torch.optim.Adam(model.parameters(),lr=.01)
                base=torch.linspace(0,1,130)[1:-1].reshape(-1,1);history=[];start=time.perf_counter()
                for step in range(1,2001):
                    x=base.detach().clone().requires_grad_(True);opt.zero_grad(set_to_none=True)
                    loss=trainer.residual(model,x,k).square().mean()
                    if not torch.isfinite(loss):raise ArithmeticError('Nonfinite training loss')
                    loss.backward();opt.step()
                    if step in (1,100,500,1000,2000):history.append(dict(step=step,loss_before_update=float(loss.detach())))
                layers=[dict(weight=l.weight.detach().tolist(),bias=l.bias.detach().tolist(),activation='tanh' if i<depth else 'identity') for i,l in enumerate(model.layers)]
                row=dict(label=f'depth{depth}_k{k}_seed49',depth=depth,width=16,frequency=k,seed=49,
                         parameter_count=sum(p.numel() for p in model.parameters()),layers=layers,
                         training_history=history,training_wall_s=time.perf_counter()-start)
                rows.append(row);save('networks_partial.json',{'networks':rows})
                print('FROZEN',row['label'],round(row['training_wall_s'],2),flush=True)
        fixture={'config':{'seed':49,'depths':[1,3],'frequencies':[2,5],'steps':2000,
                           'selection':'all final networks retained','generator_run':str(out.relative_to(ROOT))},'networks':rows}
        save('networks.json',fixture)
        CANONICAL.write_bytes((out/'networks.json').read_bytes())
        config['fixture_sha256']=hashlib.sha256(CANONICAL.read_bytes()).hexdigest()
        for row in rows:
            result=reference.MPDeepReference(row).maximum(relative=F(1,10**8),absolute=F(1,10**14),cap=5000)
            refs[row['label']]=result;save('network_references.json',refs)
            print('REFERENCE',row['label'],'converged',result['converged'],flush=True)
            if not result['converged']:raise ArithmeticError('Reference cap exhausted')
        for name,sha in config['source_sha256'].items():
            if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=sha:raise RuntimeError('Prelocked source changed')
        config.update(status='complete',frozen_cases=len(rows),continuous_references=len(refs));save('config.json',config)
    except Exception as exc:
        config.update(status='failed',error=repr(exc));save('config.json',config);raise
    print('COMPLETED',out,flush=True)

if __name__=='__main__':main()
