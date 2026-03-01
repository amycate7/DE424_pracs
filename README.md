# EMDW description

A probabilistic graphical model (PGM) toolbox with the following
characteristics:

* Flexible in terms of PGM type:

	+ Currently primarily focusses on cluster graphs

	+ We might consider region graphs sometime in the future

* Flexible in terms of factor type. The existing ones are:

	+ Arbitrary sparse discrete factor types e.g. int, unsigned, string etc.
	- Also busy with implementing fast beamed versions

	+ Gaussians including:

	- Linear Gaussians

	- Sigma Point non-linear transformations

	+ On the wish list:

	- Polya

	- Conditional Gaussian

	- Conditional Linear Gaussian

	- Gaussian Wishart

	- GMMs :)

* Flexible in terms of inference algorithm

	+ loopy Shafer-Shenoy

	+ loopy Lauritzen-Spiegelhalter

	+ On the wish list: Flexible interaction between parametric and
	  sampled representations

* Automatic graph configuration

* Efficient C++ implementation

	+ Can handle large networks

	- large numbers of random variables

	- grouped in a large number of factors

	+ On the wish list: parallelization

* On the wish list: Automatic structure learning with particular focus on the order of the model.

# Installing EMDW

This guide assumes you are running an up-to-date Linux
distribution. However, if you manage to install EMDW on another
operating system, please document your process for us to update this
guide. Since EMDW depends on `prlite` you need to install it
first. Follow the instructions in its `README.md` file.

Note we use the following convention: our development environment is located at
`~/devel/` and EMDW gets installed under `~/devel/emdw/`:
```
‌‌
	mkdir ~/devel  # Create devel directory
```

Feel free to adapt our instructions to your need.


## Upgrade CMake > 3.13

Your CMake version must be 3.13 or higher, which might not yet be available on apt-get.
To get CMake > 3.13 up and running, we recommend:

1. make sure apt-get hasn't caught up yet,
2. otherwise download a `cmake-<version>-Linux-x86_64.tar.gz` binary distribution from [cmake.org/download](https://cmake.org/download/),
3. extract the `cmake-<version>-Linux-x86_64` directory within the `.tar.gz` file to `~/devel/cmake-<version>-Linux-x86_64`,
4. add it to your system `PATH` by running:
```
‌‌
	echo "export PATH=$(readlink -f $HOME/devel/cmake-*-Linux-x86_64/bin):\$PATH" >> ~/.bashrc
	source ~/.bashrc
```
5. and then test if it is working with:
```
‌‌
	cmake --version
```
## Install EMDW dependancies

EMDW makes use of the following dependancies:
```
‌‌
	# Install dependancies
	sudo apt-get install build-essential liblapack-dev libblas-dev libboost-all-dev libssl-dev git
```

## Clone EMDW from bitbucket

We recommend registering your system on Bitbucket before you continue, for an easier git workflow:
```
‌‌
	ssh-keygen             # Generate ssh-pair; empty passphrase are allowed for convenience
	cat ~/.ssh/id_rsa.pub  # Show your public key
						   # Now register your key at https://bitbucket.org/account/user/<your-username>/ssh-keys/
```

You can now clone EMDW into `~/devel/emdw`:
```
‌‌
	cd devel
	git clone git@bitbucket.org:jadupreez/cutting_emdw.git emdw  # Assuming the cutting_emdw repo
```

## Building the library

Make sure you you have a build directory:
```
‌‌
	mkdir ~/devel/emdw/build     # Create build directory
```

Then to build only the library and tstemdw, we suggest:
```
‌‌
	cd ~/devel/emdw/build        # Go to build directory
	cmake ../; make tstemdw -j4  # Build the library and tstemdw
```

Or to build everything:
```
‌‌
	cd ~/devel/emdw/build        # Go to build directory
	cmake ../; make -j4          # Build all executables
```

You can test if EMDW is functional by running tstemdw:
```
‌‌
	~/devel/emdw/build/src/bin/tstemdw
```

Note that if your system complains about any uninstalled packages, install then and try again.

## Register EMDW system wide

If you want the EMDW **library** accessible system wide, add the `.so` location to your `PATH`:
```
‌‌
	echo "export PATH=$HOME/devel/emdw/build/src:\$PATH" >> ~/.bashrc
	source ~/.bashrc
```

If you want the EMDW **binaries** (like `tstemdw`) accessible system wide, add the `bin` location to your `PATH`:
```
‌‌
	echo "export PATH=$HOME/devel/emdw/build/src/bin:\$PATH" >> ~/.bashrc
	source ~/.bashrc
```


# Developing EMDW code

In the emdw root directory you will find a `userguide.tex` Latex
(install it if you don't already have it) file. Compile it with
```
‌‌
	pdflatex userguide
```

It gives background information on the details of the whole emdw
system -- read it first.

## Expanding EMDW itself

The main branch is called `master`. Create your own branch from this
for your edits/extensions. In the `emdw/src` directory you will
find directories for all the EMDW components, as well as `bin` for EMDW
executables and `test` for tests.

Please follow the project's layout and logic when you expand or add any
emdw features or applications. To do this properly, you will have to
familiarise yourself with the `emdw/src/CMakeLists.txt` file
and the subdirectory `CMakeLists.txt` files. If you want to use
external libraries in your project, please take a look at the examples
for Google Test, xtensor, and more in `emdw/CMakeLists.txt`.

For a quick and dirty start, we suggest you make your own
work-in-progress subdirectory `emdw/src/wip-your-features-name`
and start from there.

When extending EMDW with new objects (such as factors and their operators
etc), it is a good idea to duplicate and rename existing components.
Use appropriate subdirectories and add the names of the new files
to the appropriate `CMakeLists.txt` file.

### Testing and contributing

In the `emdw/src/test` directory are all the unit tests. It provides
typical examples of how to use the testing suite and how to write tests.
Similarly, it is good practice to write unit tests for your own
code. Once your new code passed its appropriate unit tests, make a
request via the bitbucket `Issues` facility for it to be code reviewed
and then merged back to `master`.

### Editor settings

Please set your editor so that your code will:

* contain no trailing whitespace,
* uses spaces instead of tabs,
* uses unix lf and not windows crlf.

These whitespace characters only complicate file comparisons.

## External applications linking to EMDW

If your new application just uses EMDW for some purpose without
extending its functionality, you should consider to instead create a
separate/external project that just links to the EMDW library. The
`my_emdw_app.tgz` tarball in the `emdw/extras` directory provides a handy
template for doing this.

To create your project:
```
‌‌
	cd ~/devel                                     # Go to development directory
	tar -xvzf ~/devel/emdw/extras/my_emdw_app.tgz  # Extract my_emdw_app
	mv my_emdw_app <your_project_name>             # Rename my_emdw_app
```

You can now edit your `CMakeLists.txt` files according to your needs and replace
my_emdw_app" with your project's name.

You can compile and run your project as follows:
```
‌‌
	cd ~/devel/<your_project_name>           # Go to app directory
	cmake ../; make -j4 <your_project_name>  # Compile your app
	src/<your_app_name>                      # Run your app
```

### If your system for some reason still can not find the EMDW library

Do the following ONLY if at this stage your system complains that it still can't find `libemdw.so` (horrible hack)
```
‌‌
	cd /usr/local/lib
	sudo ln -sf $HOME/devel/emdw/build/src/libemdw.so
	sudo ln -sf $HOME/devel/gLinear/build/Linux-Release-x86_64/libgLinear.so
	sudo ln -sf $HOME/devel/prlite/build/src/libprlite.so
	sudo ldconfig -v
```

## More git info

Also, you need to know how to work with git. If you don't, read
`http://stackoverflow.com/questions/315911/git-for-beginners-the-definitive-practical-guide`

## Devguide

### Compiling the Devguide

The devguide LaTex document uses a combination of [GraphViz](https://graphviz.org/) and 
[PStricks](https://tug.org/PSTricks/main.cgi/) to create clean looking graphs with LaTex notation.
In order to compile this document install the following packages:

```
	sudo apt-get install texlive-base texlive-latex-extra texlive-science texlive-pstricks
```
Then compile it with:
```
	cd docs;
	pdflatex -shell-escape devguide;
	pdflatex -shell-escape devguide;
```
The `-shell-escape` argument is necessary to compile the graphs.
### Editing the Devguide with GraphViz and PStricks

If you plan on adding a figure to the devguide you can install GraphViz as well.

```
	sudo apt-get install graphviz
```

You can then create a `figure.dot` file (refer to the existing `.dot` files in `docs` for examples) and compile it inside of the `docs` directory to a `.eps` file with:
```
	dot figure.dot -Teps -o figure.eps
```

If you want to add LaTex notation create a `figure.tex` file in `docs` and use `\psfrag`. 
Refer to one of the other `.tex` files in `docs` for examples.

To use this in a separate LaTex document include the following in your main LaTex file:
```
	\usepackage{psfrag,pstool}
```