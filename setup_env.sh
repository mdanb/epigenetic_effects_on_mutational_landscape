if conda info --envs | grep -q coo; then
	conda env update -f environment_coo.yml
else 
	conda env create -f environment_coo.yml
fi
