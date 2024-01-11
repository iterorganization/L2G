##########################
Frequently Asked Questions
##########################

What is the difference between cutoff connection length and maximum connection length?
======================================================================================

There are two parameters called cutoff_conlen and maximum_connection_length
that a user can set for a FLT case.

The maximum connection lengths, means to what length are the fieldlines
followed before stopping. With this you reduce computation time as some
fieldlines can be kilometers long in certain cases. This also causes confusion,
if we folow a fieldline and it spans to a certain maximum reachable length, set
by us, then it cannot really be treated as a connection fieldline as it is not
really connecting to any physical geometry.

The cutoff connection length is used to determine, conservatively, the wetted
areas on your analyzed geometry. In order for us to determine whether an area
is wetted, we take a fieldine from that area, follow it and see if the length
is larger than this cutoff connection length value. If it is, the area is
considered wetted.

I changed my target meshes in YAML case file but the results do not change
==========================================================================

The program checks if a resulting file exists. It does not check the contents
entirely as that is too complex to perform. Therefore when in doubt that the
resulting file does not match your expectation, delete the resulting MED file
so that it is created from scratch based on the YAML case file.
