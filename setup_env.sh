if mamba info --envs | grep -q coo; then
	mamba env update -f environment.yml
else 
	mamba env create -f environment.yml
fi
mamba activate coo
Rscript post_installation.R
