# Testing

Before distributing, the test suite should pass all tests. Run:

```
$ pytest -p no:warnings -v
```

to generate a report. 

# Create new PYPI version

You must have a pypi account, and be a collaborator. 
Sign up [here](https://pypi.org/account/register/) if necessary.

Install `twine` with:

```
$ pip install twine
```

From within the `seekr` home directory, do:

1. Bump the version number in `seekr/__version__.py` appropriately.
2. `$ rm -rf build/ dist/`
3. `$ python setup.py sdist bdist_wheel`
4. `twine upload dist/*`
