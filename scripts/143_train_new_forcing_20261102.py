"""Train preconfigured PINNs on three sources; keep all outcomes locally.

This research output is not automatically published. Two seeds are reused
across sources and architectures; the 12 outcomes are not iid samples.
"""
from pathlib import Path
import argparse,hashlib,json,math,os,platform,random,sys,time,traceback
import torch
SOURCES={'quadratic_12x2':[0,0,12],'cubic_20x3':[0,0,0,20],'quintic_42x5':[0,0,0,0,0,42]};SEEDS=(401,402);DEPTHS=(1,2);N=64;STEPS=2000;LBFGS=80
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def atom(p,obj):
 b=(json.dumps(obj,indent=2,allow_nan=False)+'\n').encode();t=p.with_suffix('.tmp')
 with t.open('xb') as f:f.write(b);f.flush();os.fsync(f.fileno())
 os.replace(t,p)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--out',required=True,type=Path);a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=False);torch.set_num_threads(1);torch.set_default_dtype(torch.float64)
 config={'status':'planned','sources':SOURCES,'seeds':SEEDS,'depths':DEPTHS,'width':16,'dtype':'float64','device':'cpu','collocation':'64 interior points i/65','adam_steps':STEPS,'adam_lr':'0.01','lbfgs_max_iter':LBFGS,'lbfgs_history':30,'lbfgs_line_search':'strong_wolfe','selection':'all planned outcomes; no use of known solution in training or selection','torch_version':torch.__version__,'python':platform.python_version(),'script_sha256':sha(Path(__file__))};atom(a.out/'locked_config.json',config)
 x0=torch.arange(1,N+1).reshape(-1,1)/(N+1);nets=[];logs=[]
 for source,coeff in SOURCES.items():
  for depth in DEPTHS:
   for seed in SEEDS:
    label=f'poly_{source}_depth{depth}_seed{seed}';tick=time.perf_counter();random.seed(seed);torch.manual_seed(seed)
    try:
     dims=[1]+[16]*depth+[1];modules=[]
     for i,(n,m) in enumerate(zip(dims[:-1],dims[1:])):
      modules.append(torch.nn.Linear(n,m))
      if i<depth:modules.append(torch.nn.Tanh())
     net=torch.nn.Sequential(*modules)
     def loss():
      x=x0.detach().clone().requires_grad_(True);v=x*(1-x)*net(x);dv=torch.autograd.grad(v,x,torch.ones_like(v),create_graph=True)[0];ddv=torch.autograd.grad(dv,x,torch.ones_like(dv),create_graph=True)[0];f=torch.zeros_like(x)
      for c in reversed(coeff):f=f*x+c
      return ((-ddv-f)**2).mean()
     opt=torch.optim.Adam(net.parameters(),lr=.01);history=[]
     for step in range(1,STEPS+1):
      opt.zero_grad();l=loss();l.backward();opt.step()
      if step in (1,500,1000,1500,2000):history.append({'step':step,'loss_before_update':float(l.detach())})
     opt2=torch.optim.LBFGS(net.parameters(),lr=1,max_iter=LBFGS,history_size=30,tolerance_grad=1e-12,tolerance_change=1e-16,line_search_fn='strong_wolfe');calls=[0]
     def closure():
      opt2.zero_grad();ll=loss();ll.backward();calls[0]+=1;return ll
     opt2.step(closure);final=float(loss().detach());linear=[x for x in net if isinstance(x,torch.nn.Linear)];layers=[{'weight':l.weight.detach().tolist(),'bias':l.bias.detach().tolist(),'activation':'identity' if i==len(linear)-1 else 'tanh'} for i,l in enumerate(linear)]
     if not math.isfinite(final):raise ValueError('nonfinite training loss')
     row={'label':label,'frequency':1,'frequency_note':'legacy placeholder; actual source in source_coefficients','depth':depth,'width':16,'seed':seed,'source_name':source,'source_coefficients':coeff,'layers':layers};nets.append(row);entry={'label':label,'status':'trained','adam_history':history,'lbfgs_closure_calls':calls[0],'final_collocation_loss':final,'training_wall_seconds':time.perf_counter()-tick}
    except Exception as e:entry={'label':label,'status':'failed','error':str(e),'traceback':traceback.format_exc(),'training_wall_seconds':time.perf_counter()-tick}
    logs.append(entry);print(label,entry['status'],entry.get('final_collocation_loss'),flush=True)
 atom(a.out/'training_log.json',logs);expected=len(SOURCES)*len(SEEDS)*len(DEPTHS);failures=[x for x in logs if x['status']!='trained'];summary={'expected':expected,'trained':len(nets),'log_records':len(logs),'failed':failures,'training_log_sha256':sha(a.out/'training_log.json'),'config_sha256':sha(a.out/'locked_config.json')}
 if not failures and len(nets)==expected:
  atom(a.out/'frozen_networks.json',{'config':config,'networks':nets});summary['fixture_sha256']=sha(a.out/'frozen_networks.json')
 atom(a.out/'summary.json',summary);print(json.dumps(summary),flush=True);return int(bool(failures) or len(nets)!=expected)
if __name__=='__main__':sys.exit(main())
