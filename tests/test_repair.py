"""Synthetic regression tests ONLY. No biological metrics are reported."""
import sys
from pathlib import Path
import unittest
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from repair_model import BiProfileEncoder, checked_splits, nested_cv


class RepairTests(unittest.TestCase):
    def test_invalid_residues_not_deleted(self):
        for bad in ['AAAXSKL', 'aaaaskl', 'AAAA*', 'AAA-AAAA', 'AAA']:
            with self.assertRaises(ValueError):
                BiProfileEncoder(window=4).fit([bad, 'AAAASKL'], [0, 1])

    def test_profiles_use_only_fit_input(self):
        encoder = BiProfileEncoder(window=3).fit(['AAASKL', 'AAASRL', 'AAAAAA', 'AAACCC'], [1, 1, 0, 0])
        before = [p.copy() for p in encoder.profiles_]
        encoder.transform(['TTTKKK', 'YYYSKL'])
        for p, original in zip(encoder.profiles_, before):
            np.testing.assert_array_equal(p, original)
        self.assertEqual(encoder.n_fit_records_, 4)
        # Two positive records end L: pseudocount 1, denominator 20+2.
        from repair_model import AA_INDEX
        self.assertAlmostEqual(encoder.profiles_[0][-1, AA_INDEX['L']], 3/22)

    def test_group_overlap_prevented(self):
        y = np.tile([0, 1], 30)
        groups = np.repeat(np.arange(30), 2)
        for tr, te in checked_splits(np.arange(60), y, groups, 3, 42):
            self.assertFalse(set(groups[tr]) & set(groups[te]))

    def test_full_nested_scope_and_determinism(self):
        rng = np.random.default_rng(21)
        aa = np.array(list('ACDEFGHIKLMNPQRSTVWY'))
        y = np.tile([0, 1], 45)
        seqs = [''.join(rng.choice(aa, 27)) + ('SKL' if v else 'AAA') for v in y]
        groups = np.repeat(np.arange(45), 2)
        kwargs = dict(window=14, outer_k=3, inner_k=2, calibration_k=2,
                      seed=42, grid={'svm__C':[1], 'svm__gamma':[.01]})
        result = nested_cv(seqs, y, groups, **kwargs)
        repeated = nested_cv(seqs, y, groups, **kwargs)
        np.testing.assert_array_equal(result['pred'], repeated['pred'])
        np.testing.assert_allclose(result['prob'], repeated['prob'], rtol=0, atol=0)
        self.assertTrue(np.all(np.isfinite(result['prob'])))
        self.assertTrue(np.all((result['prob'] >= 0) & (result['prob'] <= 1)))
        self.assertEqual(set(result['fold']), {0, 1, 2})
        for t in result['trace']:
            fit, cal, test = [set(t[k]) for k in ['fit_indices','calibration_indices','test_indices']]
            self.assertFalse(fit & cal or fit & test or cal & test)
            self.assertFalse(set(groups[list(fit)]) & set(groups[list(cal)]))
            for inner in t['inner_splits']:
                self.assertTrue(set(inner['train']).issubset(fit))
                self.assertTrue(set(inner['validation']).issubset(fit))
                self.assertFalse(set(inner['train']) & set(inner['validation']))


if __name__ == '__main__':
    unittest.main(verbosity=2)
