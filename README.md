# d2-hex

Description of the files:

k.py is the file where the pressure drops around the whole thermosyphon are written out.

new-helix.py is the file containing the code where the helical correlations were used. 

helical-hex.py is the file containing the code where the straight groove correlations were used as an approximation for the helical groove 

backwards-hex-turbulent-tube.py is the file Jeff started to follow the example from https://uspas.fnal.gov/materials/10MIT/Lecture_3.2.pdf of an approximation of the helical tube as a straight tube. I added the optimized diameter for dp=10 Pa

backwards-hex-turbulent.py is the file Jeff made to add turbulent and laminar design tools which assume finned annular geometry. 

backwards-hex-laminar.py is the file with the laminar version (so the Reynolds number is in the laminar range) of the above turbulent files.

backwards-forwards.py is the file where Kiera tried to work the backwards-hex-turbulent-tube.py code backwards (so starting with optimized diameter get dp,etc) to make sure I was maybe doing it right.

d2-hex.py is the file where the groove is just a rectangle. 

d2-hex-drawing.py is the d2-hex.py file but with the drawing of the cross section of the HEX and expanded upon. 

d2-hex-drawing-laminar.py is the d2-hex-drawing.py file but with the dimensions used to have the Reynolds number in the laminar range.
