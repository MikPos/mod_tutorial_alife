SHELL=/bin/bash
summary/all.tex: gen
	echo ''
	echo 'Summary compilation disabled by user. Run'
	echo '   mod_post --mode compileSummary'
	echo 'to manually invoke it.'
summary.pdf: summary/all.tex
	mod_post --mode compileSummary
include summary/MakefileImpl
include summary/MakefileDep
include summary/MakefileAllTarget
include summary/MakefileClean
