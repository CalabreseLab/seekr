import numpy as np

from seekr import pearson


class TestPearson:

    def test_pearson(self):
        counts1 = np.array([[8, 5, 6, 9, 2],
                            [8, 3, 6, 6, 7],
                            [7, 7, 3, 3, 7]])
        counts2 = np.array([[ 2,  8, -9, -1, -8],
                            [-4,  1,  2, -1,  2],
                            [ 5, -3, -7,  2, -9]])
        dist = pearson.pearson(counts1, counts2)
        expected = np.array([[0.3217847, -0.71611487, 0.85110363],
                             [-0.52756992, -0.47172818, 0.22652512],
                             [0.43762719, -0.17902872, 0.01547461]])
        assert np.allclose(dist, expected)

    def test_pearson_one(self):
        counts1 = np.array([[1, 2, 3, 4], [2, 4, 6, 8]])
        dist = pearson.pearson(counts1, counts1)
        expected = np.ones((2, 2))
        assert np.allclose(dist, expected)


class TestDomainPearson:
    # TODO (jessime) Add tests!
    pass
