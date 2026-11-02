"""Apply a deliberate missing-source test to the separate 163 audit."""
from pathlib import Path
import argparse,importlib.util,sys
p=argparse.ArgumentParser();p.add_argument('--fixtures',required=True,nargs='+',type=Path);p.add_argument('--out',required=True,type=Path);a=p.parse_args();f=Path(__file__).with_name('163_sine_derivative_point_audit_20261102.py');spec=importlib.util.spec_from_file_location('audit163',f);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);sys.exit(m.run(a.fixtures,a.out,omit_for_test=True))
