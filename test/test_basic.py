import unittest
import os
import sys
import utils
import contextlib
import subprocess
import shutil
import numpy

TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(TOPDIR)
MIST = os.path.join(TOPDIR, 'MiST.py')

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
            out = utils.check_output([sys.executable, MIST] + args,
                                     stderr=subprocess.STDOUT, retcode=2)

    def test_import(self):
        """Test simple import"""
        with mock_mdp_package():
            import MiST

    def test_complete_no_training(self):
        """Test complete run, no training"""
        with utils.temporary_working_directory() as tmpdir:
            shutil.copy(os.path.join(TOPDIR, 'test', 'input', 'input.tsv'), '.')
            subprocess.check_call([sys.executable, MIST, 'input.tsv',
                                   'test', '0', '0'])
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
            subprocess.check_call([sys.executable, MIST, 'input.tsv',
                                   'test', '0', '1'])
            with open('test.log') as fh:
                contents = fh.readlines()
            self.assertTrue('Number of Preys: 5\n' in contents)
            self.assertTrue('Number of Experiments: 6\n' in contents)
            self.assertTrue('Number of Baits: 3\n' in contents)
            os.unlink("test.log")
            os.unlink("test_metrics.out")
            os.unlink("test_mist.out")

    def test_complete_filtering(self):
        """Test complete run, with filtering"""
        with utils.temporary_working_directory() as tmpdir:
            shutil.copy(os.path.join(TOPDIR, 'test', 'input', 'input.tsv'), '.')
            subprocess.check_call([sys.executable, MIST, 'input.tsv',
                                   'test', '1', '0'])
            with open('test.log') as fh:
                contents = fh.readlines()
            self.assertTrue('Number of Preys: 5\n' in contents)
            self.assertTrue('Number of Experiments: 6\n' in contents)
            self.assertTrue('Number of Baits: 3\n' in contents)
            os.unlink("test.log")
            os.unlink("test_metrics.out")
            os.unlink("test_mist.out")

    def test_read_input_duplicate_experiment(self):
        """Test read_input() with duplicate experiment"""
        import MiST
        fname = os.path.join(TOPDIR, 'test', 'input',
                             'duplicate-experiment.tsv')
        self.assertRaises(MiST.MatrixFormatError, MiST.ReadInput, fname)

    def test_output_metrics_noout(self):
        """Test OutputMetrics() with no file output"""
        import MiST
        R = numpy.array(((1,2), (3,4)))
        A = numpy.array(((5,6), (7,8)))
        S = numpy.array(((9,10), (11,12)))
        B = ['bait1', 'bait2']
        P = ['prey1', 'prey2']
        matrix, pairs = MiST.OutputMetrics(R, A, S, B, P, out=0)
        self.assertAlmostEqual(matrix[0][0], 1., places=2)
        self.assertEqual(pairs,
                         [('bait1', 'prey1'), ('bait2', 'prey1'),
                          ('bait1', 'prey2'), ('bait2', 'prey2')])

if __name__ == '__main__':
    unittest.main()
