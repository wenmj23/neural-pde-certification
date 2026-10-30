"""Predefined training-only run for six polynomial Poisson PINNs.
No exact solution used for training, selection, or stopping. Keep all outcomes.
"""
from pathlib import Path
import hashlib,json,random,time,platform
import torch
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT.parent/'run_20261030/polynomial_pinn_training';OUT.mkdir(parents=True,exist_ok=False)
torch.set_num_threads(1);torch.set_default_dtype(torch.float64)
SOURCES={'constant_2':[2],'linear_6x':[0,6],'quartic_30x4':[0,0,0,0,30]};SEEDS=[101,102];STEPS=2000;LBFGS=80;N=64
config={'sources':SOURCES,'seeds':SEEDS,'architecture':'1-16-1 tanh, hard x*(1-x) boundary','dtype':'float64','device':'cpu','collocation':'64 equally spaced interior nodes i/65','adam_steps':STEPS,'adam_lr':'0.01','lbfgs_max_iter':LBFGS,'lbfgs_history':30,'lbfgs_line_search':'strong_wolfe','max_iter_stopping':'torch LBFGS internal criteria may stop early','sampling':'same deterministic collocation on every step','selection':'retain every planned network; no reference accuracy used','torch_version':torch.__version__,'python':platform.python_version(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()};(OUT/'locked_config.json').write_text(json.dumps(config,indent=2)+'\n')
rows=[];logs=[];xs0=torch.arange(1,N+1).reshape(-1,1)/(N+1)
for name,coeff in SOURCES.items():
 for seed in SEEDS:
  t=time.perf_counter();random.seed(seed);torch.manual_seed(seed);net=torch.nn.Sequential(torch.nn.Linear(1,16),torch.nn.Tanh(),torch.nn.Linear(16,1))
  def loss():
   x=xs0.detach().clone().requires_grad_(True);v=x*(1-x)*net(x);dv=torch.autograd.grad(v,x,torch.ones_like(v),create_graph=True)[0];ddv=torch.autograd.grad(dv,x,torch.ones_like(dv),create_graph=True)[0];f=torch.zeros_like(x)
   for c in reversed(coeff):f=f*x+c
   return ((-ddv-f)**2).mean()
  adam=torch.optim.Adam(net.parameters(),lr=0.01);history=[]
  for step in range(1,STEPS+1):
   adam.zero_grad();l=loss();l.backward();adam.step()
   if step in (1,500,1000,1500,2000):history.append({'step':step,'loss_before_update':float(l.detach())})
  bfgs=torch.optim.LBFGS(net.parameters(),lr=1,max_iter=LBFGS,history_size=30,tolerance_grad=1e-12,tolerance_change=1e-16,line_search_fn='strong_wolfe');calls=[0]
  def closure():
   bfgs.zero_grad();l=loss();l.backward();calls[0]+=1;return l
  bfgs.step(closure);last=float(loss().detach());layers=[]
  for layer,activation in ((net[0],'tanh'),(net[2],'identity')):
   layers.append({'weight':layer.weight.detach().tolist(),'bias':layer.bias.detach().tolist(),'activation':activation})
  label=f'poly_{name}_seed{seed}';row={'label':label,'frequency':1,'frequency_note':'legacy adapter requirement only; polynomial forcing specified separately','depth':1,'width':16,'seed':seed,'source_name':name,'source_coefficients':coeff,'layers':layers};rows.append(row);logs.append({'label':label,'adam_history':history,'lbfgs_closure_calls':calls[0],'final_training_collocation_loss':last,'training_wall_seconds':time.perf_counter()-t})
  (OUT/'partial_fixture.json').write_text(json.dumps({'config':config,'networks':rows},indent=2)+'\n');(OUT/'training_log.json').write_text(json.dumps(logs,indent=2)+'\n');print(label,'loss',last,'LBFGS_calls',calls[0],'seconds',round(time.perf_counter()-t,2),flush=True)
fixture=OUT/'frozen_polynomial_pinns.json';fixture.write_text(json.dumps({'config':config,'networks':rows},indent=2)+'\n');(OUT/'fixture_sha256.txt').write_text(hashlib.sha256(fixture.read_bytes()).hexdigest()+'\n');print('COMPLETED',len(rows),'fixture_sha256',hashlib.sha256(fixture.read_bytes()).hexdigest(),flush=True)
