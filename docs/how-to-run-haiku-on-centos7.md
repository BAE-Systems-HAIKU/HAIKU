# HOW TO RUN ON CENTOS 7

## OVERVIEW

CentOS 7 offers fewer required dependencies in its upstream repositories than Ubuntu,
so we need to build some things from scratch. The following steps were followed and
verified on a barebones CentOS 7 system.

## STEPS

### Installing YUM Dependencies
```
yum install -y epel-release
yum update -y
yum groupinstall -y "Development Tools"
yum install -y gcc openssl-devel bzip2-devel libffi-devel zlib-devel netcdf netcdf-devel
```

### Installing Python 3.8.10
Note: If you already have Python 3.8.10 installed (whether it's through conda or other means), then skip this step.
This step is simply for getting Python running on your system.

```
cd /opt
wget https://www.python.org/ftp/python/3.8.10/Python-3.8.10.tgz
tar xzf Python-3.8.10.tgz
cd Python-3.8.10
./configure --with-optimizations
make altinstall
ln -s /usr/local/bin/python3.8 /usr/local/bin/python3
```

You should now confirm your `python` version with the command: `python3 --version`

### Installing CDO from Source
The CDO binary is required for Haiku to run. On CentOS 7, we need to build this library from source.

We are using CDO version 1.9.9

```
cd /opt
wget https://code.mpimet.mpg.de/attachments/download/23323/cdo-1.9.9.tar.gz
tar xzf cdo-1.9.9.tar.gz
cd cdo-1.9.9
./configure --with-netcdf
make
make install
```

You should now confirm your `cdo` version with the command: `cdo --version`

### Installing Python Dependencies

For now, we will be using `pip` to install `python` dependencies (Poetry support can be added later)

If you are not using `conda`, then create a `python` virtual environment and activate it.
If you are using `conda`, you can most likely skip this step (although I have not tested with `conda`)

```
python3 -m venv ./venv
source ./venv/bin/activate
```

Next, install the `pip` dependencies using the `requirements.txt` file in this repo's base directory:

```
pip install -r requirements.txt
```

### Configuring PYTHONPATH

You need to make sure the project's repository base directory is located on the `PYTHONPATH` environment variable.

For example, if I checkout this repository to `/home/test/haiku`, then I can accomplish this with the following:

```
export PYTHONPATH=$PYTHONPATH:/home/test
```

I recommend pasting this line at the end of the `~/.bashrc` file, so it is run automatically with each new terminal window.
Otherwise, you'll need to run it yourself every time you open a new terminal window.

### Running the Code

With the above completed, you should now be able to run the code.

First, copy `configs/example_config.yml` and update the fields accordingly.

Then, run training with `python scripts/train.py configs/your_config.yml`

Execution output can be found in the log file specified in the configuration file.
