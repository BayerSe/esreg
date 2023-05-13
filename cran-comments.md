## Test environments
* local macOS install, R 4.3.0
* win-builder (devel and release)
* r-hub

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE:

Compilation used the following non-portable flag(s):
  ‘-Werror=format-security’ ‘-Wformat’ ‘-Wp,-D_FORTIFY_SOURCE=2’
  ‘-Wp,-D_GLIBCXX_ASSERTIONS’ ‘-march=x86-64’

## Downstream dependencies
There is one downstream dependency of which I am the author: esback. 
It is compatible with the current version of esreg.
