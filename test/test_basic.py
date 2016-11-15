import unittest
import os
import sys
import utils
import contextlib

TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(TOPDIR)

@contextlib.contextmanager
def mock_mdp_package():
    with utils.temporary_directory() as tmpdir:
        with open(os.path.join(tmpdir, 'mdp.py'), 'w') as fh:
            fh.write('pass\n')
        sys.path.insert(0, tmpdir)
        yield
        del sys.path[0]

class Tests(unittest.TestCase):
    def test_import(self):
        """Test simple import"""
        with mock_mdp_package():
            import MiST

if __name__ == '__main__':
    unittest.main()
