# uorf-predictor

Predict uORF traduction

## Setup

Setup for using the package is easy. Just install with pip.

```shell
cd uorf-predictor
python -m pip install .
```


### For developers

Install in editable mode to see the updates without having to reinstall. Just restart the kernel.

```shell
cd uorf-predictor
python -m pip install --editable .
```

### Run tests

To make sure the installation went OK, run the tests:

```shell
cd uorf-predictor
python -m pip install .[test]

# run tests
pytest
```