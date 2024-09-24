## Documentation building: Read the docs
`sphinx-quickstart` ## in /docs

`sphinx-apidoc -o docs/source src/sciduck/`

`make clean` ## in /docs

`make html`  ## in /docs

## Package building: hatchling
### Bump version of package then:
`hatch build .`

`hatch publish -u __token__ -a $SCIDUCK_PACKAGE_TOKEN`
