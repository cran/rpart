#
# Do the full shooting match of tests
#
 cat setup.s > temp
 cat stagec.s >> temp
 cat state.s >> temp
 cat gini.s >> temp
 cat cars.s >> temp
 cat treble.s treble2.s treble3.s treble4.s >> temp
 echo 'q()' >> temp
 R BATCH temp testall.out
 rm temp