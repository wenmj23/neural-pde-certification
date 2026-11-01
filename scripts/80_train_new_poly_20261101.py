"""Prespecified 12-network polynomial Poisson study; retain all outcomes."""
from pathlib import Path
import json,hashlib,random,time,platform,traceback
import torch
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT.parent/'run_20261101'/'poly_new_training';OUT.mkdir(parents=True,exist_ok=False)
torch.set_num_threads(1);torch.set_default_dtype(torch.float64)
SOURCES={'constant_2':[2],'linear_6x':[0,6],'quartic_30x4':[0,0,0,0,30]};SEEDS=[203,204];DEPTHS=[1,2];N=64;STEPS=2000;LBFGS=80
config={'status':'running','sources':SOURCES,'seeds':SEEDS,'depths':DEPTHS,'hidden_width':16,'architecture':'fully connected tanh; hard x*(1-x) boundary; depth = number of hidden layers','dtype':'float64','device':'cpu','collocation':'64 equally spaced interior nodes i/65','adam_steps':STEPS,'adam_lr':'0.01','lbfgs_max_iter':LBFGS,'lbfgs_history':30,'lbfgs_line_search':'strong_wolfe','sampling':'same deterministic collocation on every step','selection':'all planned outcomes; known solution not used during training or selecting','torch_version':torch.__version__,'python':platform.python_version(),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'baseline_commit':'8e09c023419ef92dfff9b12f29867c915ae0a4fd'}
(OUT/'locked_config.json').write_text(json.dumps(config,indent=2)+'\n');rows=[];logs=[];xs0=torch.arange(1,N+1).reshape(-1,1)/(N+1)
for name,coeff in SOURCES.items():
 for depth in DEPTHS:
  for seed in SEEDS:
   label=f'poly_{name}_depth{depth}_seed{seed}';tick=time.perf_counter();random.seed(seed);torch.manual_seed(seed)
   try:
    layerdims=[1]+[16]*depth+[1];mods=[]
    for i,(a,b) in enumerate(zip(layerdims[:-1],layerdims[1:])):
     mods.append(torch.nn.Linear(a,b))
     if i<depth:mods.append(torch.nn.Tanh())
    net=torch.nn.Sequential(*mods)
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
     bfgs.zero_grad();ll=loss();ll.backward();calls[0]+=1;return ll
    bfgs.step(closure);final=float(loss().detach());linears=[x for x in net if isinstance(x,torch.nn.Linear)];layers=[{'weight':l.weight.detach().tolist(),'bias':l.bias.detach().tolist(),'activation':'identity' if i==len(linears)-1 else 'tanh'} for i,l in enumerate(linears)]
    row={'label':label,'frequency':1,'frequency_note':'legacy parser placeholder; actual source in source_coefficients','depth':depth,'width':16,'seed':seed,'source_name':name,'source_coefficients':coeff,'layers':layers};rows.append(row);log={'label':label,'status':'trained','adam_history':history,'lbfgs_closure_calls':calls[0],'final_training_collocation_loss':final,'training_wall_seconds':time.perf_counter()-tick}
   except Exception as e:log={'label':label,'status':'training_failed','error':str(e),'traceback':traceback.format_exc(),'training_wall_seconds':time.perf_counter()-tick}
   logs.append(log);(OUT/'training_log.json').write_text(json.dumps(logs,indent=2)+'\n');(OUT/'partial_fixture.json').write_text(json.dumps({'config':config,'networks':rows},indent=2)+'\n');print(label,log['status'],log.get('final_training_collocation_loss'),round(log['training_wall_seconds'],2),flush=True)
fixture=OUT/'frozen_polynomial_pinns.json';fixture.write_text(json.dumps({'config':config,'networks':rows},indent=2)+'\n');sha=hashlib.sha256(fixture.read_bytes()).hexdigest();(OUT/'fixture_sha256.txt').write_text(sha+'\n');summary={'planned':len(SOURCES)*len(SEEDS)*len(DEPTHS),'trained':len(rows),'failed':[x for x in logs if x['status']!='trained'],'fixture_sha256':sha};(OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print('SUMMARY',json.dumps(summary),flush=True)
