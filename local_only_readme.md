# Initial setup
## Git initialization

### Preferred initial setup with ssh public key in bitbucket account:

`git clone ssh://git@git.dmz.devlnk.net/darpaactm/core.git HAIKU`

### Alternative approach using https:

`git clone https://git.dmz.devlnk.net/scm/darpaactm/core.git HAIKU`

*note, if this fails, you may need to temporarily disable the git
 sslVerification as below:*

`git config --global http.sslVerify "false"`

*and also ignore a corporate proxy when connecting to the sdmz's bitbucket
 repo:*

`export no_proxy=git.dmz.devlnk.net`

*(the above line can likely go in your ~/.bashrc file so you don't need to run it
 every time.)*

``$cd HAIKU``

## Poetry (python env) setup

Let's use poetry for our python package manager.  Since this machine is only for
HAIKU, we'd probably be okay to add packages, etc to the main python3
environment, but it's generally much safer if we each have our own user
environment and can pass around configuration files to make sure any setup is
identical once we've determined which packages to use.



1. verify if you need access to poetry locally:

``poetry --version``

2. if poetry is uninstalled, issue the following command:

`curl -sSL https://install.python-poetry.org | python3 -`

3. use the current poetry environment from your git working branch:

`poetry init`

4. activate the environment so you can use the correct python version and any
installed packages

`poetry shell`

5. install the required packages listed in the pyproject.toml:

`poetry install`

6. add any new packages to your local poetry environment that isn't already
present in the system

`poetry add new_python_package`

# Continuing a new session with already set up environment

## Loading poetry environment that's already set up

1. from the project directory run:

`poetry shell`

2. you should be good to go after that, but in case pyproject.toml has changed
since you last ran python install (for instance you pulled down a new commit
from git) you will want to update your poetry build:

`poetry install`

# Poetry Alternative

As an alternative to Poetry, we can also use python virtual environments with the supplied `requirements.txt` file.

```
python3 -m venv ./venv
source ./venv/bin/activate
pip install -r requirements.txt
```

This virtual environment can now be activated at any time, and the installed libraries are kept separate
from the global libraries.

# Other Dependencies

`cdo` is required for running this software. On Ubuntu, the binaries can be easily installed with: `apt install cdo`

Specific releases can be downloaded from the website: https://code.mpimet.mpg.de/projects/cdo/files

Building CDO on CentOS 7:

```
yum install epel-release
yum install netcdf netcdf-devel
wget https://code.mpimet.mpg.de/attachments/download/23323/cdo-1.9.9.tar.gz
tar xzf cdo-1.9.9.tar.gz
cd cdo-1.9.9
./configure --with-netcdf
make
make install
```

# HAIKU software structure

 - **data**: contains code for accessing and pre-processing climate data
   for the HAIKU system

 - **climate**: contains code to run or interface directly with climate
   models (CESM1 components)

 - **koopman**: contains code to generate and process Koopman operator based models

 - **analytics**: contains code to process models or models outputs to
   generate various analyses.

 - **TODO**: do we want a separate directory for the graphical model structure?
   Or should that live in the analytics directory?

# Documentation for HAIKU

The full documentation is contained on the
'documentation' branch on our internal bitbucket repository on VC3.
That is generate using mkdocs to translate the markdown into HTML.

If you want to modify the documentation, please do so on the
documentation branch for now (that won't contain any working software
until we are ready to release the software publicly).

## Public version of documentation

### You can checkout the latest documentation this way:

`git fetch origin`

`git checkout -b documentation origin/documentation`

*If you already have a local branch for documentation:*

`git pull origin documentation:documentation`

### Viewing documentation locally:

Once on the documentation branch, you can serve up the documentation
directly with mkdocs:

`mkdocs serve`

You can then view the full site rendered as html on a web browser at
http://127.0.0.1:8000/

## Private documentation

This README.md file is only to be shared with HAIKU team members and
focusses on what is specific to our version of the system. Any generic
documentation should be contained in the publice documentation. We'd
prefer to have all results/background on the public only.  But initial
versions of code documentation can live here and will be migrated to
the public documentation when the associated code is released. And of
course, documentation specific to our VC3 operation or internal notes
will remain on this private readme only.

### Modifying VC3/internal only documentation

This README.md file is for VC3 internal use only. Feel free to modify
it directly with any updates or additional notes of use.

### Viewing documentation without committing

If updating the readme, `grip -b README.md` will serve the markdown
document and render it as it will appear once pushed to bitbucket.

# License

For now, we're including the MIT license along with this package, but
I think we need to figure out how to fold in AIMdyn's license and
discuss with legal what license they prefer.