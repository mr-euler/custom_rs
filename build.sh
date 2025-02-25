#!/bin/bash

# set -x
if  [[ "$1" =~ ^(1|main)$ ]]; then
  app=main
elif   [[ "$1" =~ ^(2|coder)$ ]]; then
  app=coder
elif   [[ "$1" =~ ^(3|tests)$ ]]; then
  app=tests
  src="gflib"
elif   [[ "$1" =~ ^(4|fsearch)$ ]]; then
  app=forming_search
elif   [[ "$1" =~ ^(5|decoder)$ ]]; then
  app=decoder
else
  echo "wrong argument"
  exit
fi

gcc src/gflib/gflib.c src/gflib/polinom.c src/gflib/gen_polinom.c src/$src/$app.c -o $app
