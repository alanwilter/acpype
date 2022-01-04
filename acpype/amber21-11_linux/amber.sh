# Source this script to add the variables necessary to use Amber to your shell.
# This script must be located in the Amber root folder!

# Amber was configured on 2021-10-25 at 21:04:22

# determine file path of this script (credit http://unix.stackexchange.com/questions/96203/find-location-of-sourced-shell-script)
if [ -n "$BASH_SOURCE" ]; then
    this_script="$BASH_SOURCE"
elif [ -n "$DASH_SOURCE" ]; then
    this_script="$DASH_SOURCE"
elif [ -n "$ZSH_VERSION" ]; then
    setopt function_argzero
    this_script="$0"
elif eval '[[ -n ${.sh.file} ]]' 2>/dev/null; then
    eval 'this_script=${.sh.file}'
else
    echo 1>&2 "Unsupported shell. Please use bash, dash, ksh93 or zsh."
    exit 2
fi

export AMBERHOME=$(cd "$(dirname "$this_script")"; pwd)
export PATH="$AMBERHOME/bin:$PATH"

# Add Amber lib folder to LD_LIBRARY_PATH (if your platform supports it)
# Note that LD_LIBRARY_PATH is only necessary to help Amber's Python programs find their dynamic libraries,
# unless Amber has been moved from where it was installed.
if [ 1 = 1 ]; then
	if [ -z "$LD_LIBRARY_PATH" ]; then
		export LD_LIBRARY_PATH="$AMBERHOME/lib"
	else
		export LD_LIBRARY_PATH="$AMBERHOME/lib:$LD_LIBRARY_PATH"
	fi
fi

# Add location of Amber Perl modules to default Perl search path (if your platform supports it)
if [ 1 = 1 ]; then
	if [ -z "$PERL5LIB" ]; then
		export PERL5LIB="$AMBERHOME/lib/perl"
	else
		export PERL5LIB="$AMBERHOME/lib/perl:$PERL5LIB"
	fi
fi

# Add location of Amber Python modules to default Python search path (if your platform supports it)
if [ 1 = 1 ]; then
	if [ -z "$PYTHONPATH" ]; then
		export PYTHONPATH="$AMBERHOME/lib/python3.7/site-packages"
	else
		export PYTHONPATH="$AMBERHOME/lib/python3.7/site-packages:$PYTHONPATH"
	fi
fi

# Tell QUICK where to find its data files
if [ 0 = 1 ]; then
	export QUICK_BASIS="$AMBERHOME/AmberTools/src/quick/basis"
fi
