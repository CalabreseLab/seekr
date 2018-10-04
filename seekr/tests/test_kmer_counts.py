import numpy as np
import pkg_resources

from seekr import kmer_counts

class TestBasicCounter:

    def _create_basic_counter_with_data(self, **kwargs):
        infasta = 'tests/data/example.fa'
        infasta = pkg_resources.resource_filename('seekr', infasta)
        counter = kmer_counts.BasicCounter(infasta=infasta,
                                           silent=True,
                                           **kwargs)
        return counter

    def test_counter_init(self):
        counter = self._create_basic_counter_with_data()
        assert len(counter.seqs) == 5
        assert counter.seqs[0] == 'AAAAAA'

    def test_occurrences_k1(self):
        counter = self._create_basic_counter_with_data(k=1)
        row = np.zeros(4)
        expected = row.copy()
        expected[0] = 1000
        row = counter.occurrences(row, counter.seqs[0])
        assert np.allclose(row, expected)

        row = np.zeros(4)
        expected=row.copy()
        expected[1] = 500
        expected[2] = 500
        row = counter.occurrences(row, counter.seqs[1])
        assert np.allclose(row, expected)

    def test_occurrences_k2(self):
        counter = self._create_basic_counter_with_data(k=2)
        row = np.zeros(16)
        expected=row.copy()
        expected[5] = 416.666666666666
        expected[9] = 83.333333333333
        expected[10] = 416.666666666666
        row = counter.occurrences(row, counter.seqs[1])
        assert np.allclose(row, expected)

    def test_center_true(self):
        counter = self._create_basic_counter_with_data(k=1)
        counts = np.array([[1,2,3,4], [1, -2, 5, 10]], dtype=np.float32)
        counter.counts = counts
        expected = np.array([[0, 2, -1, -3], [0, -2, 1, 3]], dtype=np.float32)
        counter.center()
        assert np.allclose(counter.counts, expected)

    def test_center_vector(self):
        counter = self._create_basic_counter_with_data(k=1)
        counts = np.array([[1, 2, 3, 4], [1, -2, 5, 10]], dtype=np.float32)
        counter.counts = counts
        mean = np.ones(4)
        mean[3] = -1
        counter.mean = mean
        expected = np.array([[0, 1, 2, 5], [0, -3, 4, 11]], dtype=np.float32)
        counter.center()
        assert np.allclose(counter.counts, expected)

    def test_standardize_true(self):
        counter = self._create_basic_counter_with_data(k=1)
        counts = np.array([[1, 2, 3, 4], [0, -2, 5, 10]], dtype=np.float32)
        counter.counts = counts
        expected = np.array([[2, 1, 3, 4/3], [0, -1, 5, 10/3]], dtype=np.float32)
        counter.standardize()
        assert np.allclose(counter.counts, expected)

    def test_standardize_vector(self):
        counter = self._create_basic_counter_with_data(k=1)
        counts = np.array([[1, 2, 3, 4], [0, -2, 5, 10]], dtype=np.float32)
        counter.counts = counts
        std = np.arange(1, 5)
        counter.std = std
        expected = np.array([[1, 1, 1, 1], [0, -1, 5/3, 2.5]], dtype=np.float32)
        counter.standardize()
        assert np.allclose(counter.counts, expected)

    def test_get_counts(self):
        counter = self._create_basic_counter_with_data(k=1)
        counter.get_counts()
        expected = np.array([[ 1.8932464 , -1.4252993 , -1.6378726 , -1.2105556 ],
                             [-1.0827353 ,  1.6793344 ,  1.4974962 , -1.2105556 ],
                             [-0.3387398 ,  0.12701759, -0.07018813,  0.8703152 ],
                             [-0.32881987, -0.01786529, -0.04928571,  1.0090398 ],
                             [-0.14295156, -0.36318776,  0.25985077,  0.5417567 ]],
                             dtype=np.float32)
        assert np.allclose(counter.counts, expected)

    def test_get_counts_raw(self):
        counter = self._create_basic_counter_with_data(
            k=2, mean=False, std=False)
        counter.get_counts()
        expected = np.zeros((5, 16))
        expected[0] = counter.occurrences(counter.counts[0], counter.seqs[0])
        expected[1] = counter.occurrences(counter.counts[1], counter.seqs[1])
        expected[2] = counter.occurrences(counter.counts[2], counter.seqs[2])
        expected[3] = counter.occurrences(counter.counts[3], counter.seqs[3])
        expected[4] = counter.occurrences(counter.counts[4], counter.seqs[4])
        assert np.allclose(counter.counts, expected)