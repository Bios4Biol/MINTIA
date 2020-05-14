#/bin/bash

if [ -n $PERL5LIB ]
then
	echo $PERL5LIB > $CONDA_PREFIX/.perl5lib
	export PERL5LIB=""
fi
if [ -n $PER5LIB ]
then
	echo $PER5LIB > $CONDA_PREFIX/.per5lib
	export PER5LIB=""
fi
if [ -n $PERL_LOCAL_LIB_ROOT ]
then
	echo $PERL_LOCAL_LIB_ROOT > $CONDA_PREFIX/.perl_local_lib_root
	export PERL_LOCAL_LIB_ROOT=""
fi
if [[ -n $PERL_MB_OPT ]]
then
	echo $PERL_MB_OPT > $CONDA_PREFIX/.perl_mb_opt
	export PERL_MB_OPT='--install_base ""'
fi
if [ -n $PERL_MM_OPT ]
then
	echo $PERL_MM_OPT > $CONDA_PREFIX/.perl_mm_opt
	export PERL_MM_OPT=INSTALL_BASE=
fi

echo $PATH > $CONDA_PREFIX/.path
