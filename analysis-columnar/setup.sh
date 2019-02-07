conda create --copy --name coffea-env python=3.6
conda install --name coffea-env --file coffea-env.conda

source activate coffea-env
pip install fnal-column-analysis-tools
pip install jupyter