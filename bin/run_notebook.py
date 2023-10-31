import sys,os
import argparse

CONFIG_PATH = '.ipynb.config'

os.chdir(os.path.dirname(__file__))

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--notebook', required=True, help='Name of the Jupyter notebook to run')
    parser.add_argument('-o', '--outname', help='Name of the output HTML report')
    parser.add_argument('-f', '--format', default='html', help='Name of the output HTML report')
    parser.add_argument('remainder', nargs=argparse.REMAINDER, help='All arguments to be passed to the Jupyter notebook')
    
    return parser.parse_args()

def main():

    args = get_args()

    with open(CONFIG_PATH,'w') as f:
        f.write(' '.join([__file__] + args.remainder)+'\n')
    output = f'--output {args.outname}' if args.outname else ''
    os.system(f'jupyter nbconvert --execute --no-input --to {args.format} {output} {args.notebook}')
    return None

if __name__ == '__main__':
    main()

    