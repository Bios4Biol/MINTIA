#/bin/bash

if [ -e $CONDA_PREFIX/.perl5lib ]
then
	export PERL5LIB=`cat $CONDA_PREFIX/.perl5lib`
	rm $CONDA_PREFIX/.perl5lib
fi
if [ -e $CONDA_PREFIX/.per5lib ]
then
	export PER5LIB=`cat $CONDA_PREFIX/.per5lib`
	rm $CONDA_PREFIX/.per5lib
fi
if [ -e $CONDA_PREFIX/.perl_local_lib_root ]
then
	export PERL_LOCAL_LIB_ROOT=`cat $CONDA_PREFIX/.perl_local_lib_root`
	rm $CONDA_PREFIX/.perl_local_lib_root
fi
if [ -e $CONDA_PREFIX/.perl_mb_opt ]
then
	export PERL_MB_OPT=`cat $CONDA_PREFIX/.perl_mb_opt`
	rm $CONDA_PREFIX/.perl_mb_opt
fi
if [ -e $CONDA_PREFIX/.perl_mm_opt ]
then
	export PERL_MM_OPT=`cat $CONDA_PREFIX/.perl_mm_opt`
	rm $CONDA_PREFIX/.perl_mm_opt
fi

if [ -e $CONDA_PREFIX/.path ]
then
	export PATH=`cat $CONDA_PREFIX/.path`
	rm $CONDA_PREFIX/.path
fi
