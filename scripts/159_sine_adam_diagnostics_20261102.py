"""Replay k=5 Adam phase with fixed configuration and instrumented checkpoints.

This diagnostic is post-hoc after observing training failures. It samples
collocation points; it does not establish a root cause or repair the solver.
"""
from pathlib import Path
import argparse,hashlib,json,math,os,random,sys
import torch
STEPS=2000;N=64;CHECKPOINTS=(1,10,50,100,200,500,1000,1500,2000);COMPARE=(1,500,1000,1500,2000)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def atom(p,obj):
 b=(json.dumps(obj,indent=2,allow_nan=False)+'\n').encode();tmp=p.with_suffix('.tmp')
 with tmp.open('xb') as f:f.write(b);f.flush();os.fsync(f.fileno())
 os.replace(tmp,p)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--training',required=True,type=Path);ap.add_argument('--out',required=True,type=Path);a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=False);s=json.loads((a.training/'summary.json').read_text());config=json.loads((a.training/'locked_config.json').read_text());lp=a.training/'training_log.json'
 if sha(lp)!=s['training_log_sha256'] or sha(a.training/'locked_config.json')!=s['config_sha256']:raise ValueError('training hash mismatch')
 if config['adam_steps']!=STEPS or config['collocation']!='64 interior nodes i/65' or config['torch_version']!=torch.__version__:raise ValueError('training environment differs')
 old={x['label']:x for x in json.loads(lp.read_text())};torch.set_num_threads(1);torch.set_default_dtype(torch.float64);x0=torch.arange(1,N+1).reshape(-1,1)/(N+1);rows=[]
 for depth in (1,2):
  for seed in (501,502):
   label=f'sine_k5_depth{depth}_seed{seed}';random.seed(seed);torch.manual_seed(seed);dims=[1]+[16]*depth+[1];mods=[]
   for i,(u,v) in enumerate(zip(dims[:-1],dims[1:])):
    mods.append(torch.nn.Linear(u,v))
    if i<depth:mods.append(torch.nn.Tanh())
   net=torch.nn.Sequential(*mods)
   def loss():
    x=x0.detach().clone().requires_grad_(True);v=x*(1-x)*net(x);dv=torch.autograd.grad(v,x,torch.ones_like(v),create_graph=True)[0];ddv=torch.autograd.grad(dv,x,torch.ones_like(dv),create_graph=True)[0];f=(5*torch.pi)**2*torch.sin(5*torch.pi*x);return ((-ddv-f)**2).mean()
   opt=torch.optim.Adam(net.parameters(),lr=.01);original={x['step']:x['loss_before_update'] for x in old[label]['adam_history']}
   for step in range(1,STEPS+1):
    opt.zero_grad();l=loss();l.backward()
    if step in CHECKPOINTS:
     with torch.no_grad():
      hidden=[];state=x0
      for layer in net:
       state=layer(state)
       if isinstance(layer,torch.nn.Linear) and state.shape[1]!=1:
        z=state;hidden.append({'fraction_abs_z_ge_5':float((z.abs()>=5).double().mean()),'fraction_abs_z_ge_10':float((z.abs()>=10).double().mean()),'mean_tanh_slope':float((1-torch.tanh(z)**2).mean()),'max_abs_z':float(z.abs().max())})
      gradients=[]
      for layer in net:
       if isinstance(layer,torch.nn.Linear):
        gradients.append(float(torch.sqrt(sum((p.grad*p.grad).sum() for p in layer.parameters()))))
      observed=float(l.detach());n=net(x0);match=(math.isclose(observed,original[step],rel_tol=1e-12,abs_tol=1e-12) if step in COMPARE else None)
      rows.append({'network':label,'depth':depth,'seed':seed,'step_before_update':step,'loss_before_update':observed,'prior_loss_matches':match,'network_N_min':float(n.min()),'network_N_max':float(n.max()),'network_N_range':float(n.max()-n.min()),'linear_layer_gradient_norms':gradients,'hidden_layer_samples':hidden})
    opt.step()
   print(label,'Adam replay done',flush=True)
 atom(a.out/'records.json',rows);comparisons=[r for r in rows if r['prior_loss_matches'] is not None];summary={'training_log_sha256':s['training_log_sha256'],'training_config_sha256':s['config_sha256'],'expected':4*len(CHECKPOINTS),'records':len(rows),'prior_loss_comparisons':len(comparisons),'prior_loss_matches':sum(r['prior_loss_matches'] for r in comparisons),'records_sha256':sha(a.out/'records.json'),'script_sha256':sha(Path(__file__)),'scope':'post-hoc replay on same seed and CPU library, finite sampled collocation points; no root-cause proof'}
 if len(rows)!=summary['expected']:raise ValueError('count mismatch')
 atom(a.out/'summary.json',summary);print(json.dumps(summary),flush=True);return int(summary['prior_loss_matches']!=summary['prior_loss_comparisons'])
if __name__=='__main__':sys.exit(main())
