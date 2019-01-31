# TODO (Dan) rename file to match pwm.py
import numpy as np
import pandas as pd
import pkg_resources

from pathlib import Path
from itertools import product
from collections import defaultdict

from seekr.pwm import CountsWeighter


class TestCountsWeighter:

    def minimal_pwm(self):
        pwm = {'A': {0: .5, 1: .5, 2: .95},
               'G': {0: .1, 1: .2, 2: .05},
               'T': {0: .3, 1: 0., 2: 0},
               'C': {0: .1, 1: .3, 2: 0}}
        return pwm

    def test_get_adj_ndarray(self):
        cw = CountsWeighter(k=1)
        array = np.array([[1, 2, 3, 4], [5, 6, 7, 8]])
        counts = cw.get_counts(array)
        expected = pd.DataFrame(array, columns=['A', 'G', 'T', 'C'])
        assert expected.equals(counts)

    def test_get_adj_df(self):
        cw = CountsWeighter(k=1)
        array = np.array([[1, 2, 3, 4], [5, 6, 7, 8]])
        expected = pd.DataFrame(array, columns=['A', 'G', 'T', 'C'])
        counts = cw.get_counts(expected)
        assert expected.equals(counts)

    def test_get_adj_str_ndarray(self, tmpdir):
        cw = CountsWeighter(k=1)
        array = np.array([[1, 2, 3, 4], [5, 6, 7, 8]])
        array_path = str(Path(tmpdir, 'out.npy'))
        np.save(array_path, array)
        counts = cw.get_counts(array_path)
        expected = pd.DataFrame(array, columns=['A', 'G', 'T', 'C'])
        assert expected.equals(counts)

    def test_get_pwm_dicts(self):
        pwm_dir = pkg_resources.resource_filename('seekr', 'tests/data/pwms/')
        cw = CountsWeighter(pwm_dir, k=1)
        pwm_path, pwm = next(cw.gen_pwm_dicts())
        expected = Path(pkg_resources.resource_filename('seekr', 'tests/data/pwms/M001_0.6.txt'))
        assert pwm_path == expected
        assert len(pwm) == 4
        assert len(pwm['A']) == 7
        assert pwm['A'][0] == 0.39532879396435

    def test_set_kmer2weight(self):
        kmer2weight = defaultdict(int)
        pwm = self.minimal_pwm()
        cw = CountsWeighter(k=2)
        for kmer in (''.join(p) for p in product('AGTC', repeat=2)):
            cw.set_kmer2weight(kmer2weight, pwm, kmer, kmer, 2)
        assert kmer2weight['AA'] == (.5 * .5) + (.5 * .95)
        assert kmer2weight['GG'] == (.1 * .2) + (.2 * .05)
        assert kmer2weight['CC'] == (.1 * .3) + (.3 * 0)
        assert kmer2weight['AG'] == (.5 * .2) + (.5 * 0.05)
        assert len(kmer2weight) == 16

    def test_build_weights_dict_minimal(self):
        pwm = self.minimal_pwm()
        cw = CountsWeighter(k=2)
        kmer2weight = cw.build_weights_dict(pwm)
        assert kmer2weight['AA'] == (.5 * .5) + (.5 * .95)
        assert kmer2weight['GG'] == (.1 * .2) + (.2 * .05)
        assert kmer2weight['CC'] == (.1 * .3) + (.3 * 0)
        assert kmer2weight['AG'] == (.5 * .2) + (.5 * 0.05)
        assert len(kmer2weight) == 16

    def test_build_weights_dict_full(self):
        pwm = pkg_resources.resource_filename('seekr', 'tests/data/pwms/M001_0.6.txt')
        pwm = pd.read_csv(pwm, sep='\t').rename(columns={'U': 'T'}).to_dict()
        cw = CountsWeighter(k=2)
        kmer2weight = cw.build_weights_dict(pwm)
        assert kmer2weight['AA'] == 0.9749391864711447
        assert kmer2weight['CG'] == 0.00473602191097646
        assert kmer2weight['CT'] == 0.1186835711075973
        assert kmer2weight['AG'] == 0.02190317370024123
        assert len(kmer2weight) == 16

    def test_weight_counts(self):
        kmers = 'AGTC'
        kmer2weight = dict(zip(kmers, range(4)))
        cw = CountsWeighter(k=1)
        counts = [[1, 2, 3, 4],
                  [1, 1, 1, 1],
                  [1, 2, 1, 2]]
        counts = pd.DataFrame(np.array(counts), columns=list(kmers))
        cw.counts = counts
        sum_scores = cw.weight_counts(kmer2weight)
        assert np.all(sum_scores == np.array([20,  6, 10]))

    def test_run(self):
        pwm_dir = pkg_resources.resource_filename('seekr', 'tests/data/pwms/')
        counts = pkg_resources.resource_filename('seekr', 'tests/data/example_2mers.npy')
        cw = CountsWeighter(pwm_dir, counts, k=2)
        cw.run()
        expected = np.array([[4.56194883, 5.65991574, 6.5880013, 6.47169944, 8.14747599]])
        assert np.allclose(cw.df.values, expected)
