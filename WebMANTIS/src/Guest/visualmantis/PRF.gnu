# Generate the point response function image using myimage_***.dat file
#
# Replace 'filename' with the output myimage_.. name.
#
#	@file    PRG.gnu
#       @author  Diksha Sharma (Diksha.Sharma@fda.hhs.gov)
#       @date    Apr 13, 2012
#

reset
clear

set pm3d map
set log cb
set log z
set si sq
set cbtics format "1e%T"
unset key
set palette rgbformulae 33,13,10

splo 'myimage1.dat'

set terminal png enh font "Arial" 36 size 1028,1028 crop
set out 'PRF.png'
rep

