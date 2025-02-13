#!/bin/bash

set -x
if  [[ "$1" =~ ^(1|main)$ ]]; then
  app=main
elif   [[ "$1" =~ ^(2|coder)$ ]]; then
  app=coder
else
  echo "wrong argument"
  exit
fi

[[ -f $app ]] && rm -f $app
gcc src/gflib/gflib.c src/gflib/polinom.c src/gflib/gen_polinom.c src/$app.c -o $app 
>>>>>>> test
