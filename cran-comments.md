## Test environments
* ubuntu 16.04 (on travis-ci), R 3.3.1
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 2 NOTE:

* checking CRAN incoming feasibility ... NOTE
Maintainer: ‘Hu Sheng <shenghu@nju.edu.cn>’

  I have no idea with this note, but I think it doesn't matter.

* checking top-level files ... NOTE
Files ‘README.md’ or ‘NEWS.md’ cannot be checked without ‘pandoc’ being installed.

  If I add the README.md file to .Rbuildignore, this note will be disappear, but I don't think it's a good idea.
