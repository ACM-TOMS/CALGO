#!/bin/bash
git clone https://github.com/avventi/Slycot.git
cat Slycot/slycot/src/{sb03md,mb01rd,sb03mx,sb03my,sb03mv,sb03mw,sb04px,select}.f > sb03md.f
rm -rf Slycot
