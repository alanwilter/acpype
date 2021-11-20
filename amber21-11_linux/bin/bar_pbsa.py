#!/home/conda/feedstock_root/build_artifacts/ambertools_1635195478443/_h_env_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_/bin/python

import yaml
import numpy as np
import pandas as pd
import argparse
import logging
import parmed.amber
import parmed.tools
from pathlib import Path
import pytraj as pt
import shutil
import multiprocessing as mp
import subprocess
import shlex
import os
import itertools
#import pymbar


def init_logger(log_path='mbar_pbsa.log'):

    logging.basicConfig(level=logging.INFO,
            format='%(asctime)s : %(levelname)s : %(message)s',
            filename=log_path)
    logging.getLogger().addHandler(logging.StreamHandler())


def parse_args():

    '''
    CLI to 1) Strip trajectories
           2) Prepare BAR/PBSA input parameters
           3) Run BAR/PBSA in parallel
           4) Post-process BAR data
    '''

    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest='command')

    parser_a = subparsers.add_parser('strip')
    parser_a.add_argument('strip_input', type=str, help='Path to input strip yaml file')

    parser_b = subparsers.add_parser('prep')
    parser_b.add_argument('prep_input', type=str, help='Path in input prep yaml file')

    parser_c = subparsers.add_parser('run')
    parser_c.add_argument('run_input', type=str, help='Path in input run yaml file')
    parser_c.add_argument('-n', '--num_processes', default=4, type=int,
                          help='Number of processors to run SANDER jobs in parallel')

    parser_d = subparsers.add_parser('calc')
    parser_d.add_argument('calc_input', type=str, help='Path in input calc yaml file (same as run input)')

    return parser.parse_args()


class strip_traj:

    '''
    Prepare explicit solvent trajectories for Sander PBSA by stripping water molecules 
    and ions, concatenating replicate trajectories, and autoimaging/RMSD aligning snapshots
    to first frame. Keep dummy counter-ions for decharging and ligand atoms through AMBER
    mask selection. Setup can ignore the first half frames from each lambda window to 
    avoid unequilibrated data.
    '''

    def __init__(self, yaml_path):

        self.yaml_path = yaml_path
        self.parse_input()


    def parse_input(self):

        '''Parse input yaml'''

        with open(self.yaml_path) as f:
            inp = yaml.safe_load(f)
            logging.info(inp)

        self.dest_path = inp['dest_path']
        self.complex_paths = inp['complex_paths']
        self.ligand_paths = inp['ligand_paths']
        self.ligand_mask = inp['ligand_mask']
        self.ion_decharge = inp['ion_decharge']
        self.last_half_frames = inp['last_half_frames']
        self.stride = inp['stride']


    def get_decharged_ion(self, run_path):

        '''Find index of decharged ion to keep'''

        assert self.ion_decharge

        # load last parm at 1.000
        last_parm = list(Path(run_path, '1.000').glob('*.parm7'))[0]
        parm_in = parmed.amber.AmberParm(str(last_parm))

        # find Na+ or Cl- with 0 charge
        ions = parmed.tools.printDetails(parm_in, ':Cl-,:Na+')
        ions.execute()
        ions_str = str(ions)

        # get all ions
        ions_idx = []
        for line in ions_str.split('\n'):
            if len(line.split()) > 9 and \
               line.split()[9] not in ['-1.0000', '1.0000'] and \
               'ATOM' not in line:

                ion_line = line

                # get index
                atom_idx = ion_line.split()[0]

                ions_idx.append(atom_idx)

        return '|@'.join(ions_idx)


    def make_folds(self):

        '''Make folders to save stripped trajectories'''

        self.dest_dir = Path(Path.cwd(), self.dest_path)

        # make for stripped complex and ligands
        self.com_folds = [f't{i}' for i in range(1, len(self.complex_paths) + 1)]
        self.lig_folds = [f't{i}' for i in range(1, len(self.ligand_paths) + 1)]

        self.lig_lamdas = sorted([x.name for x in Path(self.ligand_paths[0]).iterdir()])
        self.com_lamdas = sorted([x.name for x in Path(self.complex_paths[0]).iterdir()])

        for com_fold in self.com_folds:
            for com_lamda in self.com_lamdas:
                Path(self.dest_dir, 'complex', com_fold, com_lamda).mkdir(parents=True, exist_ok=True)
                logging.info(f"mkdir {Path(self.dest_dir, 'complex', com_fold, com_lamda)}")

        for lig_fold in self.lig_folds:
            for lig_lamda in self.lig_lamdas:
                Path(self.dest_dir, 'ligands', lig_fold, lig_lamda).mkdir(parents=True, exist_ok=True)
                logging.info(f"mkdir {Path(self.dest_dir, 'ligands', lig_fold, lig_lamda)}")

        # for concat traj
        for com_lamda in self.com_lamdas:
            Path(self.dest_dir, 'concat_complex', com_lamda).mkdir(parents=True, exist_ok=True)
            logging.info(f"mkdir {Path(self.dest_dir, 'concat_complex', com_lamda)}")

        for lig_lamda in self.lig_lamdas:
            Path(self.dest_dir, 'concat_ligands', lig_lamda).mkdir(parents=True, exist_ok=True)
            logging.info(f"mkdir {Path(self.dest_dir, 'concat_ligands', com_lamda)}")


    def strip_single_traj(self, ligcom, run_path, lamda):

        '''Strip waters and ions, use amber mask to select which atoms to keep'''

        assert ligcom in ['ligands', 'complex']

        parm = str(list(Path(run_path, lamda).glob('*.parm7'))[0])
        trajin = str(list(Path(run_path, lamda).glob('*.nc'))[0])

        traj_name = trajin.split('/')[-1]
        parm_name = parm.split('/')[-1]

        traj = pt.load(trajin, parm)
        n_frames = traj.n_frames

        # strip 
        if ligcom == 'ligands':
            # ligand only
            mask = f'{self.ligand_mask}'
            paths = self.ligand_paths

        elif ligcom == 'complex':
            # ligand and protein
            # residues from number of ca
            n_residues = len(traj.top.select('@CA'))
            mask = f':1-{n_residues}|{self.ligand_mask}'
            paths = self.complex_paths
            
        # get ion index
        if self.ion_decharge:
            atom_idx = self.get_decharged_ion(run_path)
            mask = mask + f'|@{atom_idx}'

        # select region to keep
        traj = traj[mask]
        logging.info(f'{run_path}, lambda: {lamda}, mask: {mask}')

        # last half
        if self.last_half_frames:
            traj = traj[n_frames//2:]

        # subsample
        if self.stride != 1:
            traj = traj[::self.stride]

        # save nc and restart
        path_idx = paths.index(run_path)

        if ligcom == 'ligands':
            fold = self.lig_folds[path_idx]
            out_traj = str(Path(self.dest_dir, 'ligands', fold, lamda, traj_name))
            out_parm = str(Path(self.dest_dir, 'ligands', fold, lamda, parm_name))

        elif ligcom == 'complex':
            fold = self.com_folds[path_idx]
            out_traj = str(Path(self.dest_dir, 'complex', fold, lamda, traj_name))
            out_parm = str(Path(self.dest_dir, 'complex', fold, lamda, parm_name))

        pt.write_traj(out_traj, traj)
        pt.write_parm(out_parm, traj.top)


    def align(self, ligcom, lamda):

        '''Concat all replicate trajectories, RMSD fit to first frame and autoimage'''

        # glob all trajs
        # keep corresponding parm at each lambda to maintain partial charges
        if ligcom == 'ligands':
            paths = sorted(list(Path(self.dest_dir, 'ligands').glob(f'*/{lamda}/*.nc')))
            parm = sorted(list(Path(self.dest_dir, 'ligands').glob(f'*/{lamda}/*.parm7')))
        elif ligcom == 'complex':
            paths = sorted(list(Path(self.dest_dir, 'complex').glob(f'*/{lamda}/*.nc')))
            parm = sorted(list(Path(self.dest_dir, 'complex').glob(f'*/{lamda}/*.parm7')))

        paths = list(map(str, paths))

        # read first parm
        parm = str(parm[0])

        # iterload all
        traj = pt.iterload(paths, parm)

        # rmsd align to first
        if ligcom == 'ligands':
            mask = self.ligand_mask
        elif ligcom == 'complex':
            n_residues = len(traj.top.select('@CA'))
            mask = f':1-{n_residues}@CA'

        rmsd = pt.rmsd(traj, mask=mask, ref=0)

        # autoimage
        traj = traj.autoimage()

        # save concat traj, restart, parm
        # bug? only works with ncrst restarts, not rst7
        traj_save = str(Path(self.dest_dir, f'concat_{ligcom}', lamda, 'ti001.nc'))
        restart_save = Path(self.dest_dir, f'concat_{ligcom}', lamda, 'ti.ncrst')
        parm_save = str(Path(self.dest_dir, f'concat_{ligcom}', lamda, 'ti.parm'))

        pt.write_traj(traj_save, traj)
        pt.write_traj(str(restart_save), traj, frame_indices=[0])
        pt.save(parm_save, traj.top)

        # fix restart name, remove suffix
        Path(f'{restart_save}.1').rename(restart_save)

        logging.info(f'align {self.dest_path} {ligcom} {lamda}')


    def strip_all_traj(self):

        for lig_path in self.ligand_paths:
            for lamda in self.lig_lamdas:
                self.strip_single_traj('ligands', lig_path, lamda)

        for com_path in self.complex_paths:
            for lamda in self.com_lamdas:
                self.strip_single_traj('complex', com_path, lamda)


    def align_all_traj(self):

        for lig_lambda in self.lig_lamdas:
            self.align('ligands', lig_lambda)

        for com_lambda in self.com_lamdas:
            self.align('complex', com_lambda)


    def run_all(self):

        self.make_folds()
        self.strip_all_traj()
        self.align_all_traj()


class prep_bar_pbsa:

    '''
    Write Sander PBSA input files with target surface area (radiscale, 
    protscale) and interior dielectric (epsin) parameters for post-processing 
    each trajectory. This is carried out for both the ligand and complex paths. 
    '''

    def __init__(self, yaml_path):

        self.yaml_path = yaml_path
        self.parse_input()


    def parse_input(self):

        '''Parse input yaml'''

        with open(self.yaml_path) as f:
            inp = yaml.safe_load(f)
            logging.info(inp)

        self.dest_path = inp['dest_path']
        self.ligand_res = inp['ligand_res']
        self.istrng = inp['istrng']
        self.epsin = inp['epsin']
        self.radiscale = inp['radiscale']
        self.protscale = inp['protscale']


    def get_n_frames(self, ligcom):

        '''Frame counts from first traj'''

        assert ligcom in ['ligands', 'complex']

        trajin = sorted(list(Path(self.dest_path, f'concat_{ligcom}').rglob('*.nc')))[0]
        parm = sorted(list(Path(self.dest_path, f'concat_{ligcom}').rglob('*.parm')))[0]

        traj = pt.iterload(str(trajin), str(parm))

        return traj.n_frames


    def write_sander(self, ligcom):

        # set fillratio based on ligcom
        if ligcom == 'ligands':
            fillratio = 4
            dest = Path(Path.cwd(), self.lig_dest, 'pb_input.txt')
        elif ligcom == 'complex':
            fillratio = 2
            dest = Path(Path.cwd(), self.com_dest, 'pb_input.txt')

        n_frames = self.get_n_frames(ligcom)

        template = f'''PB energy calculation
 &cntrl
   imin = 5, nstlim = {n_frames}, dt = 0.002,
   ntx = 2, irest = 0,
   ntb = 0,
   ntt = 1,
   temp0 = 298.0, ig = -1,
   ntp = 0,
   ntc = 2, ntf = 2,
   ioutfm = 1, iwrap = 0,
   ntwe = 1000, ntwx = 10000, ntpr = 1, ntwr = 10000,
   cut = 9999.0,
   maxcyc = 1,
   ipb=2, inp=0

/
&pb
   npbverb=0, istrng={self.istrng}, epsout=80.0, epsin={self.epsin}, space=.5,
   accept=0.001, radiopt=0, fillratio={fillratio},
   npbopt=0, bcopt=1, solvopt=1, maxitn=10000,
   frcopt=0, nfocus=2, radiscale={self.radiscale}, radires="{self.ligand_res}", protscale={self.protscale},
   dprob=1.4
/
'''

        Path(dest).write_text(template)
        logging.info(f'{self.dest_path} {ligcom}\n{template}')


    def make_folds(self):

        self.lig_path = Path(Path.cwd(), self.dest_path, 'param_sweep_ligands')
        self.com_path = Path(Path.cwd(), self.dest_path, 'param_sweep_complex')

        folder = f'e_{self.epsin}_r_{self.radiscale}_p_{self.protscale}'

        self.lig_dest = Path(self.lig_path, folder)
        self.com_dest = Path(self.com_path, folder)

        # copy clean dir over
        com_source = Path(self.dest_path, 'concat_complex')
        lig_source = Path(self.dest_path, 'concat_ligands')

        shutil.copytree(com_source, str(self.com_dest))
        shutil.copytree(lig_source, str(self.lig_dest))

        logging.info(f'setup {self.lig_dest}')
        logging.info(f'setup {self.com_dest}')


    def run_all(self):

        self.make_folds()
        self.write_sander('ligands')
        self.write_sander('complex')


class run_bar_pbsa:

    '''
    Run Sander BAR/PBSA calculations for neighboring lambdas with multiprocessing. 
    This stage is carried out separately for the complex and ligand paths due 
    to memory limitations, this way both sets of trajectories do not have to be 
    stored together.
    '''

    def __init__(self, yaml_path, np):

        self.yaml_path = yaml_path
        self.np = np
        self.parse_input()


    def parse_input(self):

        '''Parse input yaml'''

        with open(self.yaml_path) as f:
            inp = yaml.safe_load(f)
            logging.info(inp)

        self.dest_path = inp['dest_path']
        self.epsin = inp['epsin']
        self.radiscale = inp['radiscale']
        self.protscale = inp['protscale']
        self.ligcom = inp['ligcom']
        self.del_traj = inp['del_traj']

        self.run_path = Path(Path.cwd(), self.dest_path, f'param_sweep_{self.ligcom}', f'e_{self.epsin}_r_{self.radiscale}_p_{self.protscale}')


    def make_fold(self):

        '''Make folder to store bar output'''

        self.bar_out = Path(self.run_path, 'bar_out')

        if not self.bar_out.exists():
            self.bar_out.mkdir()


    def get_cross_terms(self):

        '''Get bar neighbor lambda window combinations'''

        pairs = []

        # ignore bar_out folder
        lamdas = sorted([x.name for x in self.run_path.iterdir() if x.is_dir() if x.name != 'bar_out'])

        for i, lamda in enumerate(lamdas):
            if i == 0:
                pairs.append([lamdas[i], lamdas[i]])
                pairs.append([lamdas[i], lamdas[i+1]])

            elif i == len(lamdas) - 1:
                pairs.append([lamdas[i], lamdas[i]])
                pairs.append([lamdas[i], lamdas[i-1]])

            else:
                pairs.append([lamdas[i], lamdas[i-1]])
                pairs.append([lamdas[i], lamdas[i]])
                pairs.append([lamdas[i], lamdas[i+1]])

        return pairs


    def run_traj_parm(self, traj, parm):

        '''Run cpu sander PBSA calculation for single traj-parm cross-term'''

        amberhome = os.environ['AMBERHOME']

        cmd = f'{amberhome}/bin/sander -i {self.run_path}/pb_input.txt ' \
              f'-c {self.run_path}/{traj}/ti.ncrst ' \
              f'-p {self.run_path}/{parm}/ti.parm ' \
              f'-O -o {self.bar_out}/{traj}_{parm}.out ' \
              f'-y {self.run_path}/{traj}/ti001.nc ' \
              f'-r {self.bar_out}/{traj}_{parm}.ncrst ' \
              f'-x {self.bar_out}/{traj}_{parm}_mdcrd ' \
              f'-e {self.bar_out}/{traj}_{parm}_mden '

        logging.info(cmd)
        subprocess.run(shlex.split(cmd))


    def parallel_sander(self):

        '''Run all sander pbsa traj-parm cross terms in parallel'''

        # make cross-terms for all traj-parm pairs
        pairs = self.get_cross_terms()

        # run parallel
        pool = mp.Pool(processes=self.np)
        out = [pool.apply_async(self.run_traj_parm, args=(pair[0], pair[1])) for pair in pairs]

        pool.close()
        pool.join()


    def clean_up(self):

        '''Delete replicate traj files, ncrst, mdcrd, mden'''

        restarts = self.bar_out.glob('*.ncrst')
        mdcrds = self.bar_out.glob('*_mdcrd')
        mdens = self.bar_out.glob('*_mden')
        trajs = self.run_path.rglob('ti001.nc')

        [x.unlink() for x in restarts]
        [x.unlink() for x in mdcrds]
        [x.unlink() for x in mdens]

        if self.del_traj:
            [x.unlink() for x in trajs]


    def run_all(self):

        self.make_fold()
        self.parallel_sander()
        self.clean_up()


class parse_mbar:

    '''
    Run BAR to calculate decharging energies at each lambda window and aggregate 
    data to obtain the decharging energy for the full process. Perform 
    individually for the complex and ligand paths.
    '''

    def __init__(self, yaml_path):

        self.yaml_path = yaml_path
        self.parse_input()


    def parse_input(self):

        '''Parse input yaml'''

        with open(self.yaml_path) as f:
            inp = yaml.safe_load(f)
            logging.info(inp)

        self.dest_path = inp['dest_path']
        self.epsin = inp['epsin']
        self.radiscale = inp['radiscale']
        self.protscale = inp['protscale']
        self.ligcom = inp['ligcom']

        self.run_path = Path(Path.cwd(), self.dest_path, f'param_sweep_{self.ligcom}', f'e_{self.epsin}_r_{self.radiscale}_p_{self.protscale}')


    def read_single(self, path):

        '''Parse energies from single AMBER output'''

        traj = path.parts[-1].split('_')[0]
        parm = path.parts[-1].split('_')[-1][:-4]

        content = path.read_text().splitlines()

        energies = []
        for line in content:
            if 'minimization completed, ENE=' in line:
                energy = float(line.split()[3])
                energies.append(energy)

        df = pd.DataFrame(energies, columns=['energy'])
        df['sample'] = self.dest_path
        df['radiscale'] = self.radiscale
        df['protscale'] = self.protscale
        df['traj'] = traj
        df['parm'] = parm
        df['ligcom'] = self.ligcom
        df['frame'] = list(range(1, len(energies) + 1))

        return df


    def get_outputs(self):

        '''Glob paths to AMBER outputs for every lambda window'''

        outputs = sorted(list(self.run_path.rglob('*.out')))
        return outputs


    def get_all_energies(self):

        '''
        Parse energies from AMBER outputs at all lambda windows and save 
        to single dataframe
        '''

        outputs = self.get_outputs()
        df_list = []

        for out in outputs:
            logging.info(f'parsing {out}')
            df = self.read_single(out)
            df_list.append(df)

        self.df = pd.concat(df_list, ignore_index=True)
        self.df.to_csv(Path(self.run_path, 'energies.csv'), index=False)


#    def run_bar(self, traj_lambda):
#
#        # wF = u01 - u00
#        # wR = u10 - u11
#
#        # get index of lambda, select next neighbor
#        all_lambdas = sorted(self.df['traj'].unique())
#        idx = all_lambdas.index(traj_lambda)
#        after = all_lambdas[idx+1]
#
#        u00 = self.df.query("traj == @traj_lambda and parm == @traj_lambda")['energy'].values
#        u01 = self.df.query("traj == @traj_lambda and parm == @after")['energy'].values
#        u10 = self.df.query("traj == @after and parm == @traj_lambda")['energy'].values
#        u11 = self.df.query("traj == @after and parm == @after")['energy'].values
#
#        wF = u01 - u00
#        wR = u10 - u11
#
#        # run bar
#        bar_en, stdev = pymbar.bar.BAR(wF, wR)
#
#        return bar_en, stdev


    def bar_single(self, const, u10, u11, u01, u00):

        numer = np.sum(1 / (1 + np.exp((u10 - u11 + const) / RT)))
        denom = np.sum(1 / (1 + np.exp((u01 - u00 - const) / RT)))

        dA = RT * np.log(numer/denom) + const

        return dA


    def bar_iter(self, const, u10, u11, u01, u00, thresh=0.001, max_iters=1000):

        old_dA = const
        new_dA = self.bar_single(old_dA, u10, u11, u01, u00)
        diff = np.absolute(new_dA - old_dA)

        bailout = 0
        while diff >= thresh:

            old_dA = new_dA
            new_dA = self.bar_single(old_dA, u10, u11, u01, u00)

            bailout += 1
            if bailout == max_iters:
                break

        return new_dA


    def run_bar(self, traj_lambda):

        # get index of lambda, select next neighbor
        all_lambdas = sorted(self.df['traj'].unique())
        idx = all_lambdas.index(traj_lambda)
        after = all_lambdas[idx+1]

        u00 = self.df.query("traj == @traj_lambda and parm == @traj_lambda")['energy'].values
        u01 = self.df.query("traj == @traj_lambda and parm == @after")['energy'].values
        u10 = self.df.query("traj == @after and parm == @traj_lambda")['energy'].values
        u11 = self.df.query("traj == @after and parm == @after")['energy'].values

        bar = self.bar_iter(0, u10, u11, u01, u00)

        return bar


    def run_all(self):

        #def sd_error_prop(arr):
        #    return np.sqrt(np.sum(np.square(arr)))

        self.get_all_energies()

        # can't run bar on last lambda
        lamdas = sorted(self.df['traj'].unique())[:-1]

        energies = []
        for lamda in lamdas:
            bar = self.run_bar(lamda)
            energies.append([lamda, bar])

        # raw bar energies at each lambda
        results = pd.DataFrame(energies, columns=['lambda', 'bar'])
        results.to_csv(Path(self.run_path, 'bar_energies.csv'), index=False)

        # add together for final bar energies
        bar_total = results['bar'].sum()
        #bar_total_std = sd_error_prop(results['bar_std'].values)

        with open(Path(self.run_path, 'bar_final.txt'), 'w') as fo:
            fo.write(f'bar_total: {bar_total}\n')
            #fo.write(f'bar_total_std: {bar_total_std}')

        logging.info(f'{self.dest_path} {self.ligcom}')
        logging.info(f'bar_total: {bar_total}')
        #logging.info(f'bar_total_std: {bar_total_std}')


if __name__ == '__main__':

    # constants in kcal/mol
    R = 0.001987
    T = 298.0
    RT = R*T

    init_logger()
    args = parse_args()
    logging.info(args)

    # >>> python bar_pbsa.py strip strip_input.yaml
    if args.command == 'strip':
        cleaner = strip_traj(args.strip_input)
        cleaner.run_all()

    # >>> python bar_pbsa.py prep prep_input.yaml
    if args.command == 'prep':
        prepper = prep_bar_pbsa(args.prep_input)
        prepper.run_all()

    # >>> python bar_pbsa.py run complex_run_input.yaml
    # >>> python bar_pbsa.py run ligand_run_input.yaml
    if args.command == 'run':

        # set number of cpus for parallel sander
        avail_cpus = mp.cpu_count()
        if args.num_processes > avail_cpus:
            n_cpus = avail_cpus
        else:
            n_cpus = args.num_processes

        runner = run_bar_pbsa(args.run_input, n_cpus)
        runner.make_fold()
        runner.run_all()

    # >>> python bar_pbsa.py calc complex_run_input.yaml
    # >>> python bar_pbsa.py calc ligand_run_input.yaml
    if args.command == 'calc':
        calculator = parse_mbar(args.calc_input)
        calculator.run_all()

