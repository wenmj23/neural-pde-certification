"""One entry point for the scoped frozen-network error certifier."""
import argparse
from pathlib import Path
import sys


def main():
    parser = argparse.ArgumentParser(description='冻结tanh网络的一维Poisson连续最大误差验证。')
    parser.add_argument('--case', required=True, type=Path, help='案例JSON；网络路径相对案例文件解析')
    parser.add_argument('--out', type=Path, help='新的结果目录；已有目录不会覆盖')
    parser.add_argument('--degree', type=int, choices=(2, 3), help='残差插值次数')
    parser.add_argument('--mesh', choices=('uniform', 'adaptive'))
    parser.add_argument('--goal', choices=('enclosure', 'decision'))
    parser.add_argument('--error-tolerance', help='最大绝对误差容差，如1e-3')
    parser.add_argument('--gap-ratio', help='误差上/下界比例要求，如1.02')
    parser.add_argument('--max-intervals', type=int)
    parser.add_argument('--max-attempts', type=int)
    parser.add_argument('--time-limit-seconds')
    args = parser.parse_args()
    overrides = {name: getattr(args, name) for name in
                 ('degree', 'mesh', 'goal', 'error_tolerance', 'gap_ratio', 'max_intervals', 'max_attempts', 'time_limit_seconds')
                 if getattr(args, name) is not None}
    try:
        from certification_pipeline import run_case, human_summary
        result, output = run_case(args.case, args.out, overrides)
    except ImportError as error:
        print('缺少验证依赖：请安装 requirements-certifier.txt。' + str(error), file=sys.stderr)
        return 1
    except (ValueError, KeyError, TypeError, OSError, RuntimeError) as error:
        print('未完成验证：' + str(error), file=sys.stderr)
        return 1
    for line in human_summary(result):
        print(line)
    print('结果目录：' + str(output))
    print('可读报告：' + str(output / 'report.md'))
    return 0 if result['goal_met'] else 2


if __name__ == '__main__':
    raise SystemExit(main())
