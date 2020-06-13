 #!/bin/sh
export CPATH=
$PYTHON setup.py build
$PYTHON -m pip install . -vv
$PYTHON setup.py run_make --target 'clean'