#!/bin/bash

echo -e "\n### Begin installing mintia...\n"

# Add the launcher and the scripts to the PATH by .bashrc
HERE=$(dirname $(realpath $0))

# Export the conda functions in the subshell
source activate >& /dev/null
source $CONDA_PREFIX/etc/profile.d/conda.sh

# Conda environment installing
echo -e "### Conda environment installing... \n"
if [ $(conda env list | grep 'mintia' | wc -l) = 0 ]
then
    conda env create --file $HERE/environment.yaml
    
    if [ $? = 0 ] 
    then 
        echo -e "### Conda environment installed. \n"
    
        # Conda environment refinement
        echo -e "### Conda environment configuration...\n"
    
        conda activate mintia
    else
        code=$?
        echo "FAIL : Conda environment installation FAILED."
        exit $code
    fi
    
else
    echo -e "Conda environment already installed. \n"
    echo -e "### Conda environment configuration...\n"
    
    conda activate mintia
fi

if [ $? = 0 ]
then
    ln -s $CONDA_PREFIX/lib/libcrypto.so.1.1 $CONDA_PREFIX/lib/libcrypto.so.1.0.0
    
    mkdir -p $CONDA_PREFIX/etc/conda/activate.d/
    cp set_env_config.sh $CONDA_PREFIX/etc/conda/activate.d/
    echo -e "export PATH=$HERE:$HERE/scripts/:\$PATH" >> $CONDA_PREFIX/etc/conda/activate.d/set_env_config.sh
    
    mkdir -p $CONDA_PREFIX/etc/conda/deactivate.d/
    cp unset_env_config.sh $CONDA_PREFIX/etc/conda/deactivate.d/
    
    conda deactivate
else
    code=$?
    echo "FAIL : Conda environment activation FAILED."
    exit $code
fi

if [ $? = 0 ]
then
    echo -e "### Conda environment configured.\n"
    
    conda activate mintia
else
    code=$?
    echo "FAIL : Conda environment configuration FAILED."
    exit $code
fi

echo -e "### Installation done.\n"
