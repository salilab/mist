import unittest
import os
import sys
import utils
import contextlib
import subprocess
import shutil

TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(TOPDIR)
os.environ['PATH'] = TOPDIR + ':' + os.environ['PATH']

@contextlib.contextmanager
def mock_mdp_package():
    with utils.temporary_directory() as tmpdir:
        with open(os.path.join(tmpdir, 'mdp.py'), 'w') as fh:
            fh.write('pass\n')
        sys.path.insert(0, tmpdir)
        yield
        del sys.path[0]

class Tests(unittest.TestCase):
    def test_bad(self):
        """Test wrong arguments"""
        for args in (['x'],['x']*6):
            out = utils.check_output(['MiST.py'] + args,
                                     stderr=subprocess.STDOUT, retcode=2)

    def test_import(self):
        """Test simple import"""
        with mock_mdp_package():
            import MiST

    def test_complete_no_training(self):
        """Test complete run, no training"""
        with utils.temporary_working_directory() as tmpdir:
            shutil.copy(os.path.join(TOPDIR, 'test', 'input', 'input.tsv'), '.')
            subprocess.check_call(['MiST.py', 'input.tsv', 'test', '0', '0'])
            with open('test.log') as fh:
                contents = fh.readlines()
            self.assertTrue('Number of Preys: 5\n' in contents)
            self.assertTrue('Number of Experiments: 6\n' in contents)
            self.assertTrue('Number of Baits: 3\n' in contents)
            os.unlink("test.log")
            os.unlink("test_metrics.out")
            os.unlink("test_mist.out")

    def test_complete_training(self):
        """Test complete run, with training"""
        with utils.temporary_working_directory() as tmpdir:
            shutil.copy(os.path.join(TOPDIR, 'test', 'input', 'input.tsv'), '.')
            subprocess.check_call(['MiST.py', 'input.tsv', 'test', '0', '1'])
            with open('test.log') as fh:
                contents = fh.readlines()
            self.assertTrue('Number of Preys: 5\n' in contents)
            self.assertTrue('Number of Experiments: 6\n' in contents)
            self.assertTrue('Number of Baits: 3\n' in contents)
            os.unlink("test.log")
            os.unlink("test_metrics.out")
            os.unlink("test_mist.out")

if __name__ == '__main__':
    unittest.main()
