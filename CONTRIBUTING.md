# Want to make changes to seekr?

A great place to start is by adding more documentation to this file. 
It's a little underwhelming at the moment.

## Running tests

Install pytest with `pip install pytest`. 
Run from the repo home directory:

```
pytest
```

Also make sure you can successfully execute `seekr/tests/integration.sh`.
This is always worth checking, but doubly so if you're messing with integration code.
The unit tests don't fully cover the interfaces.

## Pushing changes

If you're ready to push changes, please ensure you've commited them to a non-master branch.
Push the branch to the remote repository, and submit a [PR](https://github.com/CalabreseLab/seekr/pulls).
Reviews are good for everyone. :)

## Publishing changes

You'll need `twine` installed via `pip install twine`.
 
Once changes have been approved and are in master

1. Bump the version in `__version__.py`
2. Build the package: `python setup.py sdist bdist_wheel`
3. Upload to pypi: `twine upload dist/*`
4. ... Pass the username and password... Right now this is tide to my private account. If someone makes it this far, email me at jessime.kirk@gmail.com and we'll figure out what we need to do. :)