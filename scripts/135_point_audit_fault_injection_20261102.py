"""Run 134's point audit against an isolated deliberately faulty subclass.

The production implementation is untouched. A mismatch is expected exactly
where the polynomial source derivative is nonzero at the tested point.
"""
from pathlib import Path
import argparse,importlib.util,sys
p=argparse.ArgumentParser();p.add_argument('--fixtures',nargs='+',required=True,type=Path);p.add_argument('--out',required=True,type=Path);p.add_argument('--mode',choices=['original','transfer','all'],default='all');a=p.parse_args()
file=Path(__file__).with_name('134_point_audit_polynomial_20261102.py');spec=importlib.util.spec_from_file_location('point_audit134',file);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
sys.exit(m.run(a.fixtures,a.out,a.mode,omit_source_for_test=True))
