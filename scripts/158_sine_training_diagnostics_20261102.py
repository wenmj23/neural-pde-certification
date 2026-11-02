"""Grid diagnostics and best constant-network residual model, not certificates.

For v=c*x*(1-x), residual r=2c-f, and least squares on the training
collocation points is minimized at c=mean(f)/2. Grid samples never establish
continuous error bounds. The training loss itself is read from its saved log.
"""
from pathlib import Path
import argparse,hashlib,json,math,os,sys
import numpy as np

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def atom(p,obj):
 b=(json.dumps(obj,indent=2,allow_nan=False)+'\n').encode();tmp=p.with_suffix('.tmp')
 with tmp.open('xb') as f:f.write(b);f.flush();os.fsync(f.fileno())
 os.replace(tmp,p)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--training',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=False);fixture=a.training/'frozen_networks.json';ts=json.loads((a.training/'summary.json').read_text());lp=a.training/'training_log.json'
 if sha(fixture)!=ts['fixture_sha256'] or sha(lp)!=ts['training_log_sha256']:raise ValueError('training hash')
 nets=json.loads(fixture.read_text())['networks'];logs={x['label']:x for x in json.loads(lp.read_text())};x=np.linspace(0,1,10001,dtype=np.float64)[:,None];coll=np.arange(1,65,dtype=np.float64)/65;rows=[]
 for net in nets:
  val=x
  for layer in net['layers']:
   val=val@np.asarray(layer['weight'],dtype=np.float64).T+np.asarray(layer['bias'],dtype=np.float64)
   if layer['activation']=='tanh':val=np.tanh(val)
   elif layer['activation']!='identity':raise ValueError('activation')
  v=(x*(1-x)*val).ravel();truth=np.sin(net['frequency']*np.pi*x.ravel());f=(net['frequency']*np.pi)**2*np.sin(net['frequency']*np.pi*coll);c=float(np.mean(f)/2);constant_mse=float(np.mean((2*c-f)**2));observed_loss=logs[net['label']]['final_collocation_loss'];vals=val.ravel()
  rows.append({'network':net['label'],'frequency':net['frequency'],'depth':net['depth'],'seed':net['seed'],'grid_points':len(x),'network_N_grid_min':float(vals.min()),'network_N_grid_max':float(vals.max()),'network_N_grid_range':float(np.ptp(vals)),'v_grid_absmax':float(np.max(np.abs(v))),'error_grid_absmax':float(np.max(np.abs(v-truth))),'best_constant_c':c,'best_constant_collocation_mse':constant_mse,'final_collocation_loss':observed_loss,'difference_loss_to_constant':float(observed_loss-constant_mse)})
 atom(a.out/'records.json',rows);summary={'fixture_sha256':ts['fixture_sha256'],'training_log_sha256':ts['training_log_sha256'],'expected':len(nets),'records':len(rows),'grid_points':len(x),'numpy_version':np.__version__,'records_sha256':sha(a.out/'records.json'),'script_sha256':sha(Path(__file__)),'scope':'float64 grid observations, not continuous maximum or formal validated computation'}
 if len(rows)!=12 or len({r['network'] for r in rows})!=12:raise ValueError('manifest')
 atom(a.out/'summary.json',summary);print(json.dumps(summary),flush=True);return 0
if __name__=='__main__':sys.exit(main())
