#!/usr/bin/perl
# Heavily revised Jan-2010 <seligman@nevis.columbia.edu>

#
# Build a description of the detectors from gdml fragments
#

# Use the command line arguments to give us the list of file fragments
# that we'll zip together.  The man page for 'GetOptions' can be read
# with: perldoc Getopt::Long

use Getopt::Long;
GetOptions( "input|i:s" => \$input,
	    "help|h" => \$help,
	    "output|o:s" => \$output);

if ( defined $help )
{
    # If the user requested help, print the usage notes and exit.
    usage();
    exit;
}

if ( ! defined $input )
{
    # If the user has not used "-i", then just take the first
    # non-option argument.
    if ( $#ARGV > -1 )
    {
	$input = $ARGV[0];
	if ( ! -r $input )
	{
	    print "Input file $input not found\n";
	    usage();
	    exit;
	}
    }
    else
    {
	# Read from STDIN
	$input = "-"; 
    }
}
else
{
    # The "-i" option was provided; check that the input file exists.
    if ( ! -r $input )
    {
	print "Input file $input not found or cannot be read\n";
	usage();
	exit;
    }
}

if ( ! defined $output )
{
    $output = "-"; # write to STDOUT 
}

# Read the list of file fragments from an XML-formatted file. 
# (Type "perldoc XML::LibXML" for how to use this package.)

use XML::LibXML;
# Create an XML parser.
$parser = new XML::LibXML;

# Read the XML file. The entire contents are slurped into a DOM
# structure.
$dom = $parser->parse_file($input);

# Note that @defFiles and @gdmlFiles are both arrays of "nodes".
@defFiles = $dom->findnodes('/config/constantfiles/filename');
@gdmlFiles = $dom->findnodes('/config/gdmlfiles/filename');

# Keep track of how many constants we read in.
$numberConstants = -1;

# Developer note: If I were slicker, I'd also use XML::LibXML to read
# in these GDML fragments. However, the files would have to be edited
# to be strict XML, and that's more effort than it's worth. The script
# generate_gdml.pl writes strict XML documents, but there's no
# guarantee that this program will be used just with those
# files.

#
# Build the table of variable replacements
#

# For each node that represents a file of constants... 
foreach $filename (@defFiles)
{
    $CONSTANTS = $filename->to_literal;
    open(CONSTANTS) or die("Could not open file $CONSTANTS");

    # For each line in the file...
    foreach $line (<CONSTANTS>) 
    {
	# Find the name and value as defined in a <constant> block; e.g.:
	#   <constant name="kInch"	value="2.54" />
	# will be parsed as $name="kInch" and value = "2.54"

	# If the line begins with "<constant"... 
	if ( $line =~ /^\s*\<constant / )
	{
	    # (See "perldoc perlfaq6" for the definition of a non-greedy search.)
	    # Do a "non-greedy" search to get the text assigned to 'name'
	    $line =~ /name="(.*?)"/;
	    $name = $1;
	    # Do a "non-greedy" search to get the text assigned to 'value'
	    $line =~ /value="(.*?)"/;
	    $value = $1;

	    # Append to the list of constants.  Put parenthesis around
	    # the value, to avoid potential problems with arithmetic.
	    ++$numberConstants;
	    $name[$numberConstants] = $name;
	    $value[$numberConstants] = "($value)";
	    # print "$numberConstants: $name = $value[$numberConstants]\n"; # debug
	}
    }
    close(CONSTANTS);
}

# The main GDML keywords, used in tags. The order here is important:
# No matter what order these blocks are in the sub-files, this is the
# order in which they must be written in the final GDML output.

# If I were being super-slick, I'd use the features of LibXML to read
# each of these GDML sub-files, then fiddle with the elements using
# XML manipulation routines. However, it's much simpler to treat them
# as big globs of text.

@keywords = qw( define materials solids structure );

# For each node that's a file of GDML statements...
foreach $filename ( @gdmlFiles )
{

    # Open the file.
    $FILE = $filename->to_literal;
    open(FILE) or die("Could not open file $FILE for read.");

    # Read the entire file into an array.
    @file = <FILE>;
    close(FILE);

    # Slurp the array into a single variable.
    $file = join("",@file);

    foreach $keyword ( @keywords )
    {
	# Search the file for a keyword block. For example, if
	# $keyword were "blorg", the following would search for
	# <blorg>some text</blorg> and snip it out of $file.
	while ( $file =~ s/(.*)\<$keyword\>(.*?)\<\/$keyword\>(.*)/\1\3/s )
	{
	    # Following the above example, this would append "some
	    # text" to $keyhash{blorg}

	    $keyhash{$keyword} = $keyhash{$keyword} . $2;
	}
    }
}

# Write the final GDML file.

$OUTPUT = ">" . $output;
open(OUTPUT) or die("Could not open $output for writing");

# The preliminary material for the GDML file. This defines the GDML
# schema and namespaces.

print OUTPUT <<EOF;
<?xml version="1.0" encoding="UTF-8" ?>
<gdml xmlns:gdml="http://cern.ch/2001/Schemas/GDML"
      xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
      xsi:noNamespaceSchemaLocation="GDMLSchema/gdml.xsd">
EOF

# Print to OUTPUT all the GDML sections, with each set of tags
# unified, in the appropriate order.

foreach $keyword ( @keywords )
{
    print OUTPUT "<$keyword>";

    # Substitute any constants in the GDML with their numeric
    # equivalents. 

    # In strict GDML, this should not be necessary; according to the
    # specification, you should be able to define a constant and use
    # it in subsequent GDML statements. If you use the standard GDML
    # parser, this not a problem. ROOT, however, uses its own parser,
    # and (as of ROOT 5.16) it could not handle variables in GDML
    # expressions. So here we side-step the problem by performing our
    # own substitutions.

    # Step through this list of constants in reverse order, to resolve
    # any dependencies.

    for ( $i = $numberConstants; $i >= 0; --$i ) 
    {
	$keyhash{$keyword} =~ s/$name[$i]/$value[$i]/g;
    }

    # Print the edited section. 
    print OUTPUT $keyhash{$keyword};

    print OUTPUT "</$keyword>\n";
}

# The final section of a GDML file is the <setup></setup> section. At
# present, we don't define alternate setups in the same GDML file, so
# these lines are constant.

print OUTPUT <<EOF;

<setup name="Default" version="1.0">
  <world ref="volWorld" />
</setup>

</gdml>
EOF

close(OUTPUT);


sub usage()
{
    print "Usage: $0 [-h|--help] [-i|--input <xml-fragments-file>] [-o|--output <output-file>]\n";
    print "       -i/--input can be omitted; if no input file, defaults to STDIN\n";
    print "       if -o is omitted, output goes to STDOUT\n";
    print "       -h prints this message, then quits\n";
}
