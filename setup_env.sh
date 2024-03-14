if conda info --envs | grep -q coo; then
	conda env update -f environment.yml
else 
	conda env create -f environment.yml
fi
