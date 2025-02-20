# Install mantid and pychop

Install mamba first per, make a "neutron" environment: https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html
Then install mantid in neutron environment per: https://www.mantidproject.org/installation/index.html

This will have to be done everytime. Install `PythonCall` (not PyCall). Open julia and set the following environment variables

```
ENV["JULIA_CONDAPKG_BACKEND"] = "Null"
ENV["JULIA_PYTHONCALL_EXE"] = "/[local_path]/miniforge3/envs/neutron/bin/python"
```

# Shiver

In your chosen environment, call

```
mamba install oauthlib
mamba install -c oncat pyoncat
mamba install -c neutrons shiver
```