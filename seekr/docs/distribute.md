# Saving changes

Before distributing, you need to save your changes.
Git can be a lot to learn, but there are good resources on the internet.
Here's one beginner tutorial:

https://git-scm.com/docs/gittutorial

You should read the tutorial above before proceeding,
but here's a one liner for saving any changes and adding them to Github

```
 $ cd ~/seekr  # Navigate to the cloned seekr repository
 $ git add . && git commit -m "Replace this with a clear description of changes." && git push
```
# Testing

Before distributing, the test suite should pass all tests. Run:

```
$ cd /path/to/seekr/
$ pytest -p no:warnings -v
```

to generate a report.

After you've added the new version to PyPI, run:

```
$./seekr/tests/integration.sh
```

to make sure the main pipeline runs without errors.

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
