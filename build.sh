#!/bin/bash

if  [[ "$1" =~ ^(1|main)$ ]]; then
  app=main
elif   [[ "$1" =~ ^(2|coder)$ ]]; then
  app=coder
else
  echo "wrong argument"
  exit
fi

[[ -f coder ]] && rm -f coder
[[ -f main ]] && rm -f main

gcc src/gflib/gflib.c src/gflib/polinom.c src/gflib/gen_polinom.c src/main.c -o $app 
