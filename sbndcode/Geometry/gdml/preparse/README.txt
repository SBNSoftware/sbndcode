% GDML preparser for SBND
% Gustavo Valdiviesso (gustavo.valdiviesso@unifal-mg.edu.br)
% June 16, 2017


This document is written in markdown format. To render it into HTML,
use `pandoc`:
    
    pandoc -o README.html README.txt
    
To render it as PDF or another format, replace "html" with the proper suffix.



Compilation
============

The preparser is compiled and installed in `sbndcode` using the usual means.

To compile the preparser outside of a LArSoft development environment:
    
    g++ -o preparse preparse.cpp `root-config --cflags --glibs`
    


Usage
======


To see if it works, run
    
    ./preparse -help
    
and a short help screen will appear.

Some basic examples:
    
    ./preparse sbnd_v00_08_base.gdml
    
makes `sbnd_v00_08.gdml` with all the wires and default setup. precision
might be a little lacky at first, so I always specify the amount of decimal
places in scientific notation, like this:
    
    ./preparse sbnd_v00_08_base.gdml -sci 10
    
If I also want a file with no wires (which is always the case), then I use:
    
    ./preparse sbnd_v00_08_base.gdml -nowires
    
and it will make `sbnd_v00_08.gdml` and `sbnd_v00_08_nowires.gdml`.

Naming convention follows exactly what `sbndcode` asks for, as long as the
base file as a _base in the name.

If I want a different setup, then:
    
    ./preparse sbnd_v00_08_base.gdml -setup DetectorOnly
    
Of course, the setup must exist or you will get an error message.

Some more advanced options: you may not want to abide by the rules and name
your files whatever you want. In that cause, you must be more specific:
    
    ./preparse myBaseGeometry.gdml -o WithWires.gdml -nowires NoWires.gdml
    
All the options can be used in any combination. As a bonus, I am also
sending my root macro for visualization. You can used it for adding color
and transparency making for a good geometry presentation. Just edit the
root file.


Contact the author
===================

Any problems, questions or suggestions, please let me know.

Best
Gustavo.

