## Resubmission
This is a resubmission. In this version I have:

* Fixed the NOTE "Found the following (possibly) invalid file URI: URI: www.robbievanaert.com From: README.md" by using http://www.robbievanaert.com. Building via win-builder did not result in an ERROR, WARNING or NOTE (although this did also not happen for the initial submission).

## Test environments
* Windows 10, R 3.6.1 (local)
* ubuntu 16.04 (on travis-ci), R 3.6.1
* osx 10.13.3 (on travis-ci), R 3.6.1
* win-builder (devel and release)

## R CMD check results
There were no ERRORs, WARNINGs or NOTEs.

## Downstream dependencies
I ran reverse dependency checks and all packages passed the checks.