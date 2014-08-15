#!/usr/bin/perl

# This program creates GDML sub-files, with values supplied by user
# parameters.  Geometry/gdml/make_gdml.pl "zips" together those
# sub-files to make a single detector description.

# Packages
use Math::Trig;
use Math::BigFloat;
Math::BigFloat->precision(-10);
use XML::LibXML;
use Getopt::Long;

# Get the input parameters from an XML file. Optionally append a
# suffix to the GDML sub-files we create.

GetOptions( "input|i:s" => \$input,
	    "help|h" => \$help,
	    "suffix|s:s" => \$suffix,
	    "output|o:s" => \$output,
	    "wires|w:s" => \$wires,
	    "cryostat|c:s" => \$cryostat);

if ( defined $help )
{
    # If the user requested help, print the usage notes and exit.
    usage();
    exit;
}

if ( ! defined $suffix )
{
    # The user didn't supply a suffix, so append nothing to the file
    # names.
    $suffix = "";
}
else
{
    # Otherwise, stick a "-" before the suffix, so that a suffix of
    # "test" applied to filename.gdml becomes "filename-test.gdml".
    $suffix = "-" . $suffix;
}


# Create an XML parser.
$parser = new XML::LibXML;

# Read in the parameters from an XML file. The following command
# slurps the entire file into a DOM data structure.
$xmldata = $parser->parse_file($input);

# Go through each parameter in the DOM data structure:
foreach $parameter ( $xmldata->findnodes('/parameters/geometry/parameter') )
{
    # Get the name and value attributes for that parameter:
    $name = $parameter->getAttribute("name");
    $value = $parameter->getAttribute("value");

    # Here's the clever part: The following eval creates a variable
    # with the same name as $name. For example, if $name eq "TPCWidth",
    # then the following statement assigns the value to $TPCWidth. The
    # value is in quotes, because some of the parameters have text
    # strings in them (like "$inch").

    eval "\$$name = '$value'";
}

# Our calculations and constants depend on the geometry of the wires.
my $SinUVAngle = sin( deg2rad($UVAngle) );
my $CosUVAngle = cos( deg2rad($UVAngle) );
my $TanUVAngle = tan( deg2rad($UVAngle) );

my $inch=2.54;
my $wires_on=1; 			# turn wires on=1 or off=0

if ( defined $wires )
{
#The user supplied the wires on parameter, so using that. Otherwise the default=1 is used.
$wires_on = $wires;
}

if (defined $cryostat)
{
#The user supplied the cryo on parameter; default is set to 1 (cryostat included)
}


my $WireInterval=10;
my $NumberOfTPCPlanes=3;
my $pmt_switch="on";		#turn on or off depending on pmts wanted
my $NumberOfTestBoxes=30;
my $enclosureExtras="on";       #turn on or off depending on whether you'd like to generate the external things around the cryostat (ie. insulation, platform, stands, etc.) in the gdml file
my $vetoWall_switch="off";  #turn on or off a proposed scintillator wall infront of the cryostat


# The routines that create the GDML sub-files. Most of the explanatory
# comments are in gen_defs().
gen_defs();
gen_rotations();
gen_materials();

gen_LAr1plane();
gen_LAr1vertplane();

gen_cathode();		# physical volumes defined in gen_tpc()
gen_fieldcage();	# physical volumes defined in gen_tpc()
gen_tpc();

if ( $pmt_switch eq "on" ) {  gen_pmt();	}	# physical volumes defined in gen_cryostat()
#gen_testbox();
if ( $enclosureExtras eq "on" ) {  gen_enclosureExtras(); } #generation of insulation, etc. will happen if specified
gen_cryostat();
if ( $vetoWall_switch eq "on" ) {  gen_vetoWall();  } # physical volumes defined in gen_vetoWall()

gen_enclosure();
gen_world();
write_fragments();

exit;



sub usage()
{
    print "Usage: $0 [-h|--help] -i|--input <parameters-file> [-o|--output <fragments-file>] [-s|--suffix <string>]\n";
    print "	  [-w wires] [-c cryostat]\n";
    print "       -i/--input can be omitted; <parameters-file> contains geometry and material parameters\n";
    print "       if -o is omitted, output goes to STDOUT; <fragments-file> is input to make_gdml.pl\n";
    print "       -s <string> appends the string to the file names; useful for multiple detector versions\n";
    print " 	  -w wires excludes wires from gdml when wires=0, default to wires=1; -c cryostat excludes\n";
    print "	  cryostat, default to cryostat=1\n"; 
    print "       -h prints this message, then quits\n";
}



# Create the detector constant file. This file is actually temporary,
# since the make_gdml.pl program will interpret its contents rather
# than include it in the final GDML file.
# The detector information below is for LAr1-ND


sub gen_defs()
{
    $pi   = pi;
    $inch = 2.54;
    
    #TPCWirePlaneLength is the size in the z direction    
 	#TPCWirePlaneWidth is the size in the y direction
    $TPCWirePlaneLengthZ	=	365;
    $TPCWirePlaneWidthY	=	480;

    $WorldLength	=	100.0*$DetEnclosureLength;
    $WorldWidth		=	100.0*$DetEnclosureWidth;
    $WorldHeight	=	100.0*$DetEnclosureHeight;

    $CathodePlateDepth	=	0.1;
    $CathodeLengthZ	=	365;
    $CathodeWidthX	=	5;
    $CathodeHeightY	=	400;

	$AnodeLengthZ  	=   365;
    $AnodeWidthX 	=   20;
    $AnodeHeightY 	=   400;	

	#Below are numbers for uboone geometry
    $tpc_neg_length	=	1000;
    $tpc_neg_width	=	200;
    $tpc_neg_height	=	200;
    $wire_frame_width	=	9.5;
    $wire_plane_height	=	240;
    $wire_plane_length	=	1042;
    $wires_plength		=	($wire_plane_length - 2*$wire_frame_width) ;
    $wires_pwidth		=	($wire_plane_height - 2*$wire_frame_width) ;

    $field_cage_width			=	256;
    $field_cage_height			=	207.6;  #increased size to avoid overlap of TubeBottom and Top with TPCActive
    $field_cage_cross_length	=	sqrt(($field_cage_width)**2+($field_cage_height-50)**2);
    $field_cage_length			=	1011.6; #increased size to avoid overlap of TubeFront and Back with TPCActive
    $field_cage_loop_interval	=	1; 	# =1 is normal, =4 skips 3/4

    $FieldCageTPCClearance	=	5*$inch;
    $FieldCageTubeRadius	=	0.5*$inch;
    $FieldCageTubeThickness	=	0.25*$inch;
    $FieldCageBeamDepth		=	12.5;
    $FieldCageBeamWidth		=	2.5;
    $FieldCageCrossDepth	=	2.5;
    $FieldCageCrossWidth	=	4;
    $FieldCageCrossLength	=	$field_cage_cross_length;

    $TPCTotalLength	=	$field_cage_length;
    $TPCTotalWidth	=	$field_cage_width;
    $TPCTotalHeight	=	$field_cage_height;

    $FieldCageLoopLength	=	$TPCTotalLength+2*($FieldCageTPCClearance+2*$FieldCageTubeRadius);
    $FieldCageLoopWidth		=	$field_cage_width;
    $FieldCageLoopHeight	=	$TPCTotalHeight+2*($FieldCageTPCClearance+2*$FieldCageTubeRadius);

    $FieldCageCornerRadius	=	0.5*$FieldCageTPCClearance;
    $FieldCageCornerY		=	(0.5*$FieldCageLoopHeight)-$FieldCageCornerRadius-$FieldCageTubeRadius;
    $FieldCageCornerZ		=	(0.5*$FieldCageLoopLength)-$FieldCageCornerRadius-$FieldCageTubeRadius;

    $FieldCageHeight	=	$FieldCageLoopHeight+2*($FieldCageBeamDepth+$FieldCageCrossDepth);
    $FieldCageLength	=	$FieldCageLoopLength+2*($FieldCageBeamDepth+$FieldCageCrossDepth);
    $FieldCageWidth	=	$field_cage_width;

}


sub gen_rotations()
{
    my $WirePlusRotation = $UVAngle + 90;
    my $WireMinusRotation = $UVAngle - 90;

    $ROTATIONS = "LAr1-rotations" . $suffix . ".gdml";
    push (@gdmlFiles, $ROTATIONS); # Add file to list of GDML fragments
    $ROTATIONS = ">" . $ROTATIONS;
    open(ROTATIONS) or die("Could not open file $ROTATIONS for writing");

    print ROTATIONS <<EOF;
<?xml version='1.0'?>
<define>
   <rotation name="rPlus30AboutX"  unit="deg" x="30"  y="0"   z="0"/>
   <rotation name="rPlus60AboutX"  unit="deg" x="60"  y="0"   z="0"/>
   <rotation name="rPlus90AboutX"  unit="deg" x="90"  y="0"   z="0"/>
   <rotation name="rMinus90AboutX"  unit="deg" x="-90"  y="0"   z="0"/>
   <rotation name="rPlusUVAngleAboutX"  unit="deg" x="90+$UVAngle" y="0"   z="0"/>
   <rotation name="rPlus150AboutX"      unit="deg" x="150" y="0"   z="0"/>
   <rotation name="rPlus180AboutX"      unit="deg" x="180" y="0"   z="0"/>
   <rotation name="rMinusUVAngleAboutX" unit="deg" x="-30" y="0"   z="0"/>
   <rotation name="rPlus30AboutY"  unit="deg" x="0"   y="30"  z="0"/>
   <rotation name="rPlus60AboutY"  unit="deg" x="0"   y="60"  z="0"/>
   <rotation name="rPlus90AboutY"  unit="deg" x="0"   y="90"  z="0"/>
   <rotation name="rPlus180AboutY" unit="deg" x="0"   y="180" z="0"/>
   <rotation name="rMinus90AboutY" unit="deg" x="0"   y="-90" z="0"/>
   <rotation name="rPlus90AboutZ"  unit="deg" x="0"   y="0"   z="90"/>
   <rotation name="rMinus90AboutZ"  unit="deg" x="0"   y="0"   z="-90"/>
   <rotation name="rPlus180AboutZ"      unit="deg" x="0"   y="0"   z="180"/>
   <rotation name="rMinus180AboutZ"     unit="deg" x="0"   y="0"   z="-180"/>
   <rotation name="rMinus90AboutYPlus180AboutZ" unit="deg" x="0" y="-90" z="180"/>
   <rotation name="rMinus90AboutYMinus90AboutZ" unit="deg" x="0" y="-90" z="-90"/>
   <rotation name="rPlus90AboutYPlus180AboutZ" unit="deg" x="0" y="90" z="180"/>
   <rotation name="rMinus90AboutYPlus90AboutZ" unit="deg" x="0" y="-90" z="90"/>
   <rotation name="rPlus90AboutYMinus90AboutZ" unit="deg" x="0" y="90" z="-90"/>
   <rotation name="rPlus90AboutXPlus90AboutZ"  unit="deg" x="90" y="0"   z="90"/>
   <rotation name="rPlus90AboutXPlus180AboutZ" unit="deg" x="90" y="0"   z="180"/>
   <rotation name="rPlus90AboutXMinus90AboutY" unit="deg" x="90" y="-90" z="0"/>
   <rotation name="rPlus90AboutXMinus90AboutZ" unit="deg" x="90" y="0"   z="-90"/>
   <rotation name="rPlus90AboutXPlus90AboutY"  unit="deg"  x="90" y="90" z="0"/>
   <rotation name="rPMTRotation1"  unit="deg" x="90"  y="270"   z="0"/>
   <position name="posCenter" unit="mm" x="0" y="0" z="0"/>

</define>
EOF
    close (ROTATIONS);
}


sub gen_materials()
{
    # Create the materials file name and open it.
    $MATERIALS = "materials" . $suffix . ".gdml";
    push (@gdmlFiles, $MATERIALS); # Add file to list of GDML fragments
    $MATERIALS = ">" . $MATERIALS;
    open(MATERIALS) or die("Could not open file $MATERIALS for writing");

    # Write the standard XML prefix.
    print MATERIALS <<EOF;
<?xml version='1.0'?>
EOF

    # Go back the DOM structure read in near the beginning of the
    # program. For each <materials /> element (and there'll probably
    # be only one):
    foreach $materials ( $xmldata->findnodes('/parameters/materials') )
    {
	# Convert that element back to text, and write it out.
	print MATERIALS $materials->toString;
    }

    close (MATERIALS);
}


# This is a re-write of Brian Rebel's gen_microvertplane.C into
# Perl. It contructs the TPC wire plane for the Y view.

sub gen_LAr1vertplane()
{
  my $NumberWires = 0;
  
  if ( $wires_on == 1 ) 
    { $NumberWires = int( $TPCWirePlaneLengthZ / $TPCWirePitch ); }
  
    

    $GDML = "LAr1-vertplane" . $suffix . ".gdml";
    push (@gdmlFiles, $GDML); # Add file to list of GDML fragments
    $GDML = ">" . $GDML;
    open(GDML) or die("Could not open file $GDML for writing");

    # Define the solids and structures: the wires and the TPC wire plane.
    print GDML <<EOF;
<?xml version='1.0'?>
<gdml>
<solids>
<tube name="TPCWireVert"
  rmax="0.5*$TPCWireThickness"
  z="$TPCWirePlaneWidthY"
  deltaphi="360"
  aunit="deg"
  lunit="cm"/>
<box name="TPCPlaneVert"
  x="$TPCWirePlaneThickness" 
  y="$TPCWirePlaneWidthY" 
  z="$TPCWirePlaneLengthZ"
  lunit="cm"/>
</solids>
<structure>
  <volume name="volTPCWireVert">
    <materialref ref="Titanium"/>
    <solidref ref="TPCWireVert"/>
  </volume>
  <volume name="volTPCPlaneVert">
    <materialref ref="LAr"/>       
    <solidref ref="TPCPlaneVert"/>
EOF

    # the wires 
    for ( $i = 0; $i < $NumberWires; ++$i)
    {
	print GDML <<EOF;
 	    <physvol>
	     <volumeref ref="volTPCWireVert"/>
	     <position name="posTPCWireVert$i" unit="cm" z="-0.5*$TPCWirePlaneLengthZ+$TPCWirePitch*($i+1)" x="$TPCWidth/2" y="0"/>
	     <rotationref ref="rPlus90AboutX"/>
	    </physvol>
EOF
    }

    print GDML <<EOF;
  </volume>
</structure>
</gdml>
EOF

    close(GDML);
}


# This is a re-write of Brian Rebel's gen_microplane.C into Perl. It
# constructs the TPC wire plane for the U or V view.

sub gen_LAr1plane()
{

    $TPCWireXPitch=$TPCWirePitch/$CosUVAngle;

    $GDML = "LAr1-plane" . $suffix . ".gdml";
    push (@gdmlFiles, $GDML); # Add file to list of GDML fragments
    $GDML = ">" . $GDML;
    open(GDML) or die("Could not open file $GDML for writing");

    # Calculate the number of wire ends on a given y-edge of the plane.
    my $TPCYWirePitch = $TPCWirePitch / $CosUVAngle;

   my $NumberWiresPerEdge = 0;
   if ( $wires_on == 1 ) 
    {  $NumberWiresPerEdge = int( $TPCWirePlaneLengthZ / $TPCYWirePitch );}
  
    # How many side wires will be "cut off" by the lower or higher
    # z-edge?
   my $NumberSideWires = 0;
    if ( $wires_on == 1 ) 
    { $NumberSideWires = int( $TanUVAngle * $TPCWirePlaneWidthY / $TPCYWirePitch ); }

    # The number of full-length "center" wires.
  
   my $NumberCenterWires = 0.;
   if ( $wires_on == 1 ) 
    { $NumberCenterWires = $NumberWiresPerEdge - $NumberSideWires; }

    # define the solids
    print GDML <<EOF;
<?xml version='1.0'?>
<gdml>
<solids>
EOF

    # wires on either end of the tpc #subtraction added to avoid overlap for all wires
    for($i = 0; $i < $NumberSideWires; ++$i)
    {
	print GDML <<EOF;
<tube name="TPCWire$i"
  rmax="0.5*$TPCWireThickness"
  z="$TPCWireXPitch*($i+1)/$SinUVAngle-0.01498648"    
  deltaphi="360"
  aunit="deg"
  lunit="cm"/> 
EOF
    }

    # The solids for the middle wire and the TPC wire plane, and start off the structures.
    #subtraction added to avoid overlap for all wires
    print GDML <<EOF;
<tube name="TPCWireCommon"
  rmax="0.5*$TPCWireThickness"
  z="$TPCWirePlaneWidthY/$CosUVAngle-0.0259573526"
  deltaphi="360"
  aunit="deg"
  lunit="cm"/>
<box name="TPCPlane"
  x="$TPCWirePlaneThickness"
  y="$TPCWirePlaneWidthY"
  z="$TPCWirePlaneLengthZ"
  lunit="cm"/>
</solids>
<structure>
EOF
 
    # the wires at either end of the plane
    for ($i = 0; $i < $NumberSideWires; ++$i)
    {
	print GDML <<EOF;
    <volume name="volTPCWire$i">
	<materialref ref="Titanium"/>
	<solidref ref="TPCWire$i"/>
    </volume>
EOF
    }

  
    # The wires in the middle of the plane, and the plane itself.
    print GDML <<EOF;
    <volume name="volTPCWireCommon">
      <materialref ref="Titanium"/>
      <solidref ref="TPCWireCommon"/>
    </volume>
    <volume name="volTPCPlane">
      <materialref ref="LAr"/>
      <solidref ref="TPCPlane"/>
EOF

  # the wires at the -z end
  for ($i = 0; $i < $NumberSideWires; ++$i)
  {
    $j=$NumberSideWires-$i-1;

    print GDML <<EOF;
    <physvol>
     <volumeref ref="volTPCWire$j"/> 
     <position name="posTPCWireF$j" unit="cm" y="-0.5*($i)*$TPCWireXPitch/$TanUVAngle" z="-0.5*$TPCWirePlaneLengthZ+0.5*$TPCWireXPitch*($j+1)" x="0"/>
     <rotationref ref="rPlusUVAngleAboutX"/>
    </physvol>
EOF
  $ypos=-0.5*($i)*$TPCWireXPitch/$TanUVAngle;
  $zpos=-0.5*$TPCWirePlaneLengthZ+0.5*$TPCWireXPitch*($j);
  open (MYFILE, '>>data.txt');
  print MYFILE "TPCWire$j y=$ypos z=$zpos\n";
  }

  # the wires at the +z end
  for ($i = 0; $i < $NumberSideWires; ++$i)
  {
      my $j = $NumberSideWires-$i-1;
      my $k = $NumberCenterWires+$NumberSideWires+$i;
      $zposlast=-0.5*$TPCWirePlaneLengthZ+$TPCWireXPitch*(0.5*$NumberSideWires+$NumberCenterWires);
      $zpos=$zposlast+0.5*$TPCWireXPitch*($i);

    print GDML <<EOF;
    <physvol>
     <volumeref ref="volTPCWire$j"/> 
     <position name="posTPCWireB$j" unit="cm" y="0.5*($i)*$TPCWireXPitch/$TanUVAngle" z="$zpos" x="0"/>
     <rotationref ref="rPlusUVAngleAboutX"/>
    </physvol>
EOF


  }

  # The wires in the middle.
  for ($i = 0; $i < $NumberCenterWires - 1; ++$i)
  {
      my $j = $NumberSideWires+$i;
      $zpos=-0.5*$TPCWirePlaneLengthZ+$TPCWireXPitch*(0.5*$NumberSideWires+$i+1);

      print GDML <<EOF;
    <physvol>
     <volumeref ref="volTPCWireCommon"/>
     <position name="posTPCWire$j" unit="cm" y="0" z="$zpos" x="0"/>
     <rotationref ref="rPlusUVAngleAboutX"/>
    </physvol>
EOF
  }


      print GDML <<EOF;
  </volume>
</structure>
</gdml>
EOF

  close(GDML);
}


# subdirectory to write field cage
sub gen_fieldcage() {

    # Set up the output file.
    $FIELDCAGE = "LAr1-fieldcage.gdml";
    push (@gdmlFiles, $FIELDCAGE); # Add file to list of constant files
    $FIELDCAGE = ">" . $FIELDCAGE;
    open(FIELDCAGE) or die("Could not open file $FIELDCAGE for writing");

    # set the field cage constants

    # Prints Field Cage solids
    print FIELDCAGE <<EOF;
<solids>
 <tube name="FieldCageTubeZ"
  rmin="$FieldCageTubeRadius-$FieldCageTubeThickness"
  rmax="$FieldCageTubeRadius"  
  z="$FieldCageLoopLength-2*($FieldCageCornerRadius+$FieldCageTubeRadius)" 
  deltaphi="360" 
  aunit="deg" 
  lunit="cm"/> 
 <tube name="FieldCageTubeY" 
  rmin="$FieldCageTubeRadius-$FieldCageTubeThickness"
  rmax="$FieldCageTubeRadius"  
  z="$FieldCageLoopHeight-2*($FieldCageCornerRadius+$FieldCageTubeRadius)"  
  deltaphi="360" 
  aunit="deg" 
  lunit="cm"/> 
</solids> 
EOF

    # Prints Field Cage tube loop sub-structure
    print FIELDCAGE <<EOF;
<structure> 
 <volume name="volFieldCageTubeTop"> 
  <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/> 
  <solidref ref="FieldCageTubeZ"/> 
 </volume> 
 <volume name="volFieldCageTubeBot"> 
  <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/> 
  <solidref ref="FieldCageTubeZ"/> 
 </volume> 
 <volume name="volFieldCageTubeFront"> 
  <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/> 
  <solidref ref="FieldCageTubeY"/> 
 </volume> 
 <volume name="volFieldCageTubeBack"> 
  <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/> 
  <solidref ref="FieldCageTubeY"/> 
 </volume> 
</structure>
EOF

close(FIELDCAGE); 

}

# Cathode solids and volumes
sub gen_cathode() {

    $CATHODE = "LAr1-cathode" . $suffix . ".gdml";
    push (@gdmlFiles, $CATHODE); # Add file to list of GDML fragments
    $CATHODE = ">" . $CATHODE;
    open(CATHODE) or die("Could not open file $CATHODE for writing");

    print CATHODE <<EOF;

<?xml version='1.0'?>
<gdml>
<solids>
 <box name="CathodePlate" lunit="cm" x="$CathodeWidthX" y="$CathodeHeightY" z="$CathodeLengthZ"/>
</solids>
<structure>
 <volume name="volCathodePlate">
  <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
  <solidref ref="CathodePlate"/>
 </volume>
</structure>
</gdml>
EOF

}


# Parameterize the TPC and the planes within it.
sub gen_tpc()
{
    # Set up the output file.
    $GDML = "LAr1-tpc" . $suffix . ".gdml";
    push (@gdmlFiles, $GDML); # Add file to list of GDML fragments
    $GDML = ">" . $GDML;
    open(GDML) or die("Could not open file $GDML for writing");


    # Size info for active TPC volume (LAr inside)
    $aTPC_xos_cathode = 0.5*$CathodeWidthX + 0.5*$CathodePlateDepth;
    $aTPC_xos_wires = 3*$TPCWirePlaneSpacing + $TPCWirePlaneThickness;
    $aTPC_xoffset = $aTPC_xos_cathode - $aTPC_xos_wires ; 
    $TPCActiveDepth  = $TPCWidth - $aTPC_xos_cathode - $aTPC_xos_wires ; 
    $TPCActiveHeight = $TPCWirePlaneWidthY;   #  
    $TPCActiveLength = $TPCWirePlaneLengthZ-0.2;   # extra subtraction to arrive at TPCActive values in the TDR

    print GDML <<EOF;
<?xml version='1.0'?>
<gdml>
<solids>
 <box name="TPC" lunit="cm" x="$TPCWidth" y="$TPCHeight" z="$TPCLength"/>
 <box name="TPCActive" lunit="cm" x="$TPCActiveDepth" y="$TPCActiveHeight" z="$TPCActiveLength"/>
 <box name="TPCFrameA" lunit="cm" x="$AnodeWidthX" y="$AnodeHeightY" z="$AnodeLengthZ"/>
 <box name="TPCFrameB" lunit="cm" x="$AnodeWidthX+ 0.1" y="$AnodeHeightY-20" z="$AnodeLengthZ -20"/>
 <box name="TPCSideBeam" lunit="cm" x="$TPCWidth+$AnodeWidthX" y="5" z="10"/>
 
 <subtraction name="TPCFrame">
    <first ref="TPCFrameA"/> <second ref="TPCFrameB"/>
    <position name="posTPCSubtraction" x="0" y="0" z="0"/>
  </subtraction> 

</solids>
<structure>
 <volume name="volTPCActive">
  <materialref ref="LAr"/>
  <solidref ref="TPCActive"/>
 </volume>
 <volume name="volTPCFrame">
  <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
  <solidref ref="TPCFrame"/>
 </volume>
 <volume name="volTPCSideBeam">
  <materialref ref="G10"/>
  <solidref ref="TPCSideBeam"/>
 </volume>
 <volume name="volTPC">
  <materialref ref="LAr"/>
  <solidref ref="TPC"/>
EOF


    # Cathode Plane physical volumes
    # Center = (0.5*($TPCWidth-$CathodeWidthX),0,0)
    $Cathode_X=0.5*($TPCWidth-$CathodeWidthX);
    $Cathode_plate_offset=0.5*$CathodePlateDepth;
   # $Cathode_frame_position=$Cathode_X+$Cathode_plate_offset;

    print GDML <<EOF;
  <physvol>
   <volumeref ref="volCathodePlate"/>
   <position name="posCathodePlate" unit="cm" x="0" y="0" z="0" />
  </physvol>
  <physvol>
   <volumeref ref="volTPCFrame"/>
   <position name="posTPCFrame" unit="cm" x="$TPCWidth/2 +10" y="0" z="0"/>
  </physvol>
  <physvol>
	 <volumeref ref="volTPCFrame"/>
     <position name="posTPCFrame1" unit="cm" x="-$TPCWidth/2 - 10" y="0" z="0"/>
  </physvol>
  <physvol>
	 <volumeref ref="volTPCSideBeam"/>
     <position name="posTPCSideBeam0" unit="cm" x="0" y="80" z="$TPCLength/2"/>
  </physvol>
  <physvol>
	 <volumeref ref="volTPCSideBeam"/>
     <position name="posTPCSideBeam1" unit="cm" x="0" y="80*2" z="$TPCLength/2"/>
  </physvol>
  <physvol>
	 <volumeref ref="volTPCSideBeam"/>
     <position name="posTPCSideBeam2" unit="cm" x="0" y="-80" z="$TPCLength/2"/>
  </physvol>
  <physvol>
	 <volumeref ref="volTPCSideBeam"/>
     <position name="posTPCSideBeam3" unit="cm" x="0" y="-80*2" z="$TPCLength/2"/>
  </physvol>
  <physvol>
	 <volumeref ref="volTPCSideBeam"/>
     <position name="posTPCSideBeam4" unit="cm" x="0" y="0" z="$TPCLength/2"/>
  </physvol>
  <physvol>
	 <volumeref ref="volTPCSideBeam"/>
     <position name="posTPCSideBeam0" unit="cm" x="0" y="80" z="-$TPCLength/2"/>
  </physvol>
  <physvol>
	 <volumeref ref="volTPCSideBeam"/>
     <position name="posTPCSideBeam1" unit="cm" x="0" y="80*2" z="-$TPCLength/2"/>
  </physvol>
  <physvol>
	 <volumeref ref="volTPCSideBeam"/>
     <position name="posTPCSideBeam2" unit="cm" x="0" y="-80" z="-$TPCLength/2"/>
  </physvol>
  <physvol>
	 <volumeref ref="volTPCSideBeam"/>
     <position name="posTPCSideBeam3" unit="cm" x="0" y="-80*2" z="-$TPCLength/2"/>
  </physvol>
  <physvol>
	 <volumeref ref="volTPCSideBeam"/>
     <position name="posTPCSideBeam4" unit="cm" x="0" y="0" z="-$TPCLength/2"/>
  </physvol>
  <physvol>
	 <volumeref ref="volTPCSideBeam"/>
     <position name="posTPCSideBeam0" unit="cm" x="0" y="$TPCHeight/2" z="0"/>
  </physvol>
  <physvol>
	 <volumeref ref="volTPCSideBeam"/>
     <position name="posTPCSideBeam1" unit="cm" x="0" y="$TPCHeight/2" z="80"/>
  </physvol>
  <physvol>
	 <volumeref ref="volTPCSideBeam"/>
     <position name="posTPCSideBeam2" unit="cm" x="0" y="$TPCHeight/2" z="80*2"/>
  </physvol>
  <physvol>
	 <volumeref ref="volTPCSideBeam"/>
     <position name="posTPCSideBeam3" unit="cm" x="0" y="$TPCHeight/2" z="-80"/>
  </physvol>
  <physvol>
	 <volumeref ref="volTPCSideBeam"/>
     <position name="posTPCSideBeam4" unit="cm" x="0" y="$TPCHeight/2" z="-80*2"/>
  </physvol>
  <physvol>
	 <volumeref ref="volTPCSideBeam"/>
     <position name="posTPCSideBeam0" unit="cm" x="0" y="-$TPCHeight/2" z="0"/>
  </physvol>
  <physvol>
	 <volumeref ref="volTPCSideBeam"/>
     <position name="posTPCSideBeam1" unit="cm" x="0" y="-$TPCHeight/2" z="80"/>
  </physvol>
  <physvol>
	 <volumeref ref="volTPCSideBeam"/>
     <position name="posTPCSideBeam2" unit="cm" x="0" y="-$TPCHeight/2" z="80*2"/>
  </physvol>
  <physvol>
	 <volumeref ref="volTPCSideBeam"/>
     <position name="posTPCSideBeam3" unit="cm" x="0" y="-$TPCHeight/2" z="-80"/>
  </physvol>
  <physvol>
	 <volumeref ref="volTPCSideBeam"/>
     <position name="posTPCSideBeam4" unit="cm" x="0" y="-$TPCHeight/2" z="-80*2"/>
  </physvol>



EOF


    # Wire Plane physical volumes 
    # Center = ( -0.5*($TPCHeight-kWireFrameWidth) , 0 , 0 )
    print GDML <<EOF;
<!--  <physvol>
   <volumeref ref="volTPCPlaneVert"/>
   <position name="posTPCPlaneVert" unit="cm" x="-0.5*$TPCActiveDepth-2*$TPCWirePlaneSpacing" y="0" z="0"/>
  </physvol>
  <physvol>
   <volumeref ref="volTPCPlane"/>
   <position name="posTPCPlane" unit="cm" x="-0.5*$TPCActiveDepth-$TPCWirePlaneSpacing" y="0" z="0"/>
  </physvol>
  <physvol>
   <volumeref ref="volTPCPlane"/>
   <position name="posTPCPlane2" unit="cm" x="-0.5*$TPCActiveDepth" y="0" z="0"/>
   <rotationref ref="rPlus180AboutY"/>
  </physvol> -->
EOF


$space=0;
$i=1;
while ( $space < ( $field_cage_width / 2 ) ) {
	$xPos=$space+2;
	print GDML <<EOF;
<!--  <physvol>
   <volumeref ref="volFieldCageTubeTop"/>
   <position name="posFieldCageTubeTopA$i" unit="cm" x="$xPos" y="0.5*($FieldCageLoopHeight-2*$FieldCageTubeRadius)" z="0"/>
  </physvol>
  <physvol>
   <volumeref ref="volFieldCageTubeBot"/>
   <position name="posFieldCageTubeBotA$i" unit="cm" x="$xPos" y="-0.5*($FieldCageLoopHeight-2*$FieldCageTubeRadius)" z="0"/>
  </physvol>
  <physvol>
   <volumeref ref="volFieldCageTubeFront"/>
   <position name="posFieldCageTubeFrontA$i" unit="cm" x="$xPos" y="0" z="0.5*($FieldCageLoopLength-2*$FieldCageTubeRadius)"/>
   <rotation name="rFieldCageVertPlusA$i" unit="deg" x="90" y="0" z="0"/>
  </physvol>
  <physvol>
   <volumeref ref="volFieldCageTubeBack"/>
   <position name="posFieldCageTubeBackA$i" unit="cm" x="$xPos" y="0" z="-0.5*($FieldCageLoopLength-2*$FieldCageTubeRadius)"/>
   <rotation name="rFieldCageVertMinusA$i" unit="deg" x="-90" y="0" z="0"/>
  </physvol>  -->
EOF
	print GDML <<EOF;
<!--  <physvol>
   <volumeref ref="volFieldCageTubeTop"/>
   <position name="posFieldCageTubeTopB$i" unit="cm" x="-$xPos" y="0.5*($FieldCageLoopHeight-2*$FieldCageTubeRadius)" z="0"/>
  </physvol>
  <physvol>
   <volumeref ref="volFieldCageTubeBot"/>
   <position name="posFieldCageTubeBotB$i" unit="cm" x="-$xPos" y="-0.5*($FieldCageLoopHeight-2*$FieldCageTubeRadius)" z="0"/>
  </physvol>
  <physvol>
   <volumeref ref="volFieldCageTubeFront"/>
   <position name="posFieldCageTubeFrontB$i" unit="cm" x="-$xPos" y="0" z="0.5*($FieldCageLoopLength-2*$FieldCageTubeRadius)"/>
   <rotation name="rFieldCageVertFrontB$i" unit="deg" x="90" y="0" z="0"/>
  </physvol>
  <physvol>
   <volumeref ref="volFieldCageTubeBack"/>
   <position name="posFieldCageTubeBackB$i" unit="cm" x="-$xPos" y="0" z="-0.5*($FieldCageLoopLength-2*$FieldCageTubeRadius)"/>
   <rotation name="rFieldCageVertBackB$i" unit="deg" x="-90" y="0" z="0"/>
  </physvol>  -->
EOF
	$space+=4*$field_cage_loop_interval;
	$i++;
}


    print GDML <<EOF;
<!--   <physvol>
    <volumeref ref="volTPCActive"/>
    <position name="posTPCActive" unit="cm" x="-$aTPC_xoffset" y="0" z="0"/>
   </physvol> -->
EOF

    # Closes TPC volume definition space
    print GDML <<EOF;
 </volume> 
</structure>
</gdml>
EOF

   close(GDML);
}


# Generates Ben Jones's PMT micro-pmtdef (with temporary edit to ellipsoid shapes
sub gen_pmt {

    $PMT = "LAr1-pmtdef" . $suffix . ".gdml";
    push (@gdmlFiles, $PMT); # Add file to list of GDML fragments
    $PMT = ">" . $PMT;
    open(PMT) or die("Could not open file $PMT for writing");

	print PMT <<EOF;
<solids>
 <tube name="PMTVolume"
  rmax="(6.1*2.54)"
  z="(11.1*2.54)+0.005"
  deltaphi="360"
  aunit="deg"
  lunit="cm"/>

 <tube name="PMT_AcrylicPlate"
  rmax="(6.0*2.54)"
  z="(0.2)"
  deltaphi="360"
  aunit="deg"
  lunit="cm"/>
 <tube name="PMT_Stalk"
  rmax="(1.25*2.54)"
  z="(3.0*2.54)"
  deltaphi="360"
  aunit="deg"
  lunit="cm"/>
 <tube name="PMT_SteelBase"
  rmax="(6.0*2.54)"
  z="(1.5*2.54)"
  deltaphi="360"
  aunit="deg"
  lunit="cm"/>
 <tube name="PMT_Underside"
  rmax="2.54*4.0"
  z="2.54*2.5"
  deltaphi="360"
  aunit="deg"
  lunit="cm"/>
EOF
	print PMT <<EOF;
 <tube name="PMT_Lens"
  rmax="2.54*4.0"
  z="2.54*2.5"
  deltaphi="360"
  aunit="deg"
  lunit="cm"/>
EOF

	print PMT <<EOF;
</solids>
<structure>
 <volume name="volOpDetSensitive">
  <materialref ref="LAr"/>
  <solidref ref="PMT_AcrylicPlate"/>
 </volume>
 <volume name="vol_PMT_AcrylicPlate">
  <materialref ref="Acrylic"/>
  <solidref ref="PMT_AcrylicPlate"/>
 </volume>
 <volume name="vol_PMT_Stalk">
  <materialref ref="Glass"/>
  <solidref ref="PMT_Stalk"/>
 </volume>
 <volume name="vol_PMT_SteelBase">
  <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
  <solidref ref="PMT_SteelBase"/>
 </volume>
 <volume name="vol_PMT_Underside">
  <materialref ref="Glass"/>
  <solidref ref="PMT_Underside"/>
 </volume>
EOF
	print PMT <<EOF;
 <volume name="vol_PMT_Lens">
  <materialref ref="LAr"/>
  <solidref ref="PMT_Lens"/>
 </volume>
EOF
	print PMT <<EOF;
 <volume name="volPMT">
  <materialref ref="LAr"/>
  <solidref ref="PMTVolume"/>
  <physvol>
   <volumeref ref="volOpDetSensitive"/>
   <position name="posOpDetSensitive" unit="cm" x="0" y="0" z="(5.5 * 2.54) - 0.1"/>
  </physvol>
  <physvol>
   <volumeref ref="vol_PMT_AcrylicPlate"/>
   <position name="pos_PMT_AcrylicPlate" unit="cm" x="0" y="0" z="(5.5 * 2.54) - 0.3"/>
  </physvol>
  <physvol>
   <volumeref ref="vol_PMT_Stalk"/>
   <position name="pos_PMT_Stalk" unit="cm" x="0" y="0" z="(3.0 * 2.54)-(5.5 * 2.54)"/>
  </physvol>
  <physvol>
   <volumeref ref="vol_PMT_SteelBase"/>
   <position name="pos_PMT_SteelBase" unit="cm" x="0" y="0" z="(0.75 * 2.54)-(5.5 * 2.54)"/>
  </physvol>
  <physvol>
   <volumeref ref="vol_PMT_Lens"/>
   <position name="pos_PMT_Lens" unit="cm" x="0" y="0" z="(7.0 * 2.54)-(5.5 * 2.54)"/>
  </physvol>
  <physvol>
   <volumeref ref="vol_PMT_Underside"/>
   <position name="pos_PMT_Underside" unit="cm" x="0" y="0" z="(7.0 * 2.54)-(5.5 * 2.54)"/>
  </physvol>
EOF

	print PMT <<EOF;
 </volume>
</structure>
EOF
}


#Parameterize the steel cryostat that encloses the TPC.
sub gen_cryostat()
{

    # Set up the output file.
    $CRYOSTAT = "LAr1-cryostat" . $suffix . ".gdml";
    push (@gdmlFiles, $CRYOSTAT); # Add file to list of GDML fragments
    $CRYOSTAT = ">" . $CRYOSTAT;
    open(CRYOSTAT) or die("Could not open file $CRYOSTAT for writing");

    print CRYOSTAT <<EOF;
<?xml version='1.0'?>
<gdml>
<solids>
 <tube name="Cryostat" rmax="93" z="2000" deltaphi="360" aunit="deg" lunit="cm"/>
<sphere name="EndCap" rmin="144*2.54" rmax="144.5*2.54" deltaphi="360" deltatheta="31.3822" aunit="deg" lunit="cm"/>
EOF
	
  print CRYOSTAT <<EOF;
</solids>

<structure>
<volume name="volEndCap">
   <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
   <solidref ref="EndCap"/>
  </volume>
 <volume name="volSteelTube">
  <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
  <solidref ref="SteelTube"/>
 </volume>
 <volume name="volCryostat">
  <materialref ref="LAr"/>
  <solidref ref="Cryostat"/>
<physvol>
   <volumeref ref="volEndCap"/>
   <position name="posEndCap1" unit="cm" x="0" y="0" z="427.75*2.54/2- 2.54*sqrt(144.5^2-75.5^2)"/>
   </physvol> 
  <physvol>
   <volumeref ref="volTPC"/>
   <position name="posTPC" unit="cm" x="0.0" y="0" z="0"/>
  </physvol> 
EOF


  @pmt_pos = ( ' x="-147.8"  y="3.21654"  z="-472"',
               ' x="-147.8"  y="3.21654"  z="472"' );

  if ( $pmt_switch eq "on" ) {
    for ( $i=0; $i<2; ++$i ){
      print CRYOSTAT <<EOF;
<!--  <physvol>
   <volumeref ref="volPMT"/>
   <position name="posPMT$i" unit="cm" @pmt_pos[$i]/>
   <rotationref ref="rPMTRotation1"/>
  </physvol> -->
EOF
    }
  }

	print CRYOSTAT <<EOF;
 </volume>
</structure>
</gdml>
EOF

   close(CRYOSTAT);
}

#Generates veto wall
sub gen_vetoWall()
{
  # Set up the output file.
  $VW = "LAr1-vetoWall" . $suffix . ".gdml";
  push (@gdmlFiles, $VW); # Add file to list of GDML fragments
  $VW = ">" . $VW;
  open(VW) or die("Could not open file $VW for writing");
  
  print VW <<EOF;
  <?xml version='1.0'?>
  <gdml>
    <solids>
      <box lunit="cm" name="AuxDet0"    x="10"    y="20"    z="2" />
      <box lunit="cm" name="AuxDet1"    x="10"    y="20"    z="1" />
  
  </solids>
  <structure>
    <volume name="volAuxDet0">
      <materialref ref="Polystyrene"/>
      <solidref ref="AuxDet0"/>
    </volume>
    <volume name="volAuxDet1">
      <materialref ref="Polystyrene"/>
      <solidref ref="AuxDet1"/>
    </volume>
  </structure>
EOF
}

# Parameterize the cryostat's surroundings.
sub gen_enclosure()
{
    # Set up the output file.
    $GDML = "LAr1-enclosure" . $suffix . ".gdml";
    push (@gdmlFiles, $GDML); # Add file to list of GDML fragments
    $GDML = ">" . $GDML;
    open(GDML) or die("Could not open file $GDML for writing");

    print GDML <<EOF;
<?xml version='1.0'?>
<gdml>
<solids>
 <box name="DetEnclosure" lunit="cm" x="$DetEnclosureWidth" y="$DetEnclosureHeight" z="$DetEnclosureLength" />
</solids>

<structure>
 <volume name="volDetEnclosure">
  <materialref ref="Air"/>
  <solidref ref="DetEnclosure"/>
  <physvol>
   <volumeref ref="volCryostat"/>
   <position name="posCryostat" unit="cm" x="0" y="0" z="0"/>
  </physvol>
EOF
#if the extra pieces in the enclosure are desired, put them in the correct places within volEnclosure
  if ( $enclosureExtras eq "on" ) {
    print GDML <<EOF;
EOF
}

  
  if ( $vetoWall_switch eq "on" ) {
    print GDML <<EOF;
    
    <physvol>
    <volumeref ref="volAuxDet0"/>
    <position name="posAuxDet0" unit="cm" x="0" y="0" z="-800"/>
    <rotation name="rAuxDet0" unit="deg" x="0" y="0" z="0"/>
    </physvol>
    <physvol>
    <volumeref ref="volAuxDet1"/>
    <position name="posAuxDet1" unit="cm" x="0" y="0" z="-805"/>
    <rotation name="rAuxDet1" unit="deg" x="0" y="0" z="0"/>
    </physvol>
    
EOF
  }

    print GDML <<EOF;
 </volume>
</structure>
</gdml>
EOF

   close(GDML);
}


# Parameterize the dirt mound that surrounds the enclosure.
sub gen_world()
{
    # Set up the output file.
    $GDML = "LAr1-world" . $suffix . ".gdml";
    push (@gdmlFiles, $GDML); # Add file to list of GDML fragments
    $GDML = ">" . $GDML;
    open(GDML) or die("Could not open file $GDML for writing");

    print GDML <<EOF;
<?xml version='1.0'?>
<gdml>
<solids>
  <box name="World" lunit="cm" x="$WorldWidth" y="$WorldHeight" z="$WorldLength"/>
  <tube name="Ground" rmin="312*2.54" rmax="((50*12)+310)*2.54" z="41*12*2.54" deltaphi="360" lunit="cm" aunit="deg"/>
</solids>

<structure>
  <volume name="volGround" >
    <materialref ref="Dirt" />
    <solidref ref="Ground" />
  </volume>
  <volume name="volWorld" >
    <materialref ref="Air"/> 
    <solidref ref="World"/>
    <physvol>
      <volumeref ref="volDetEnclosure"/>
      <position name="posDetEnclosure" unit="cm" x="0" y="0" z="0"/>
    </physvol>
      <!--<position name="posDetEnclosure" unit="cm" x="0.5*$TPCActiveDepth" y="0" z="0.5*$TPCWirePlaneLengthZ"/>-->
  </volume> 
</structure>
</gdml>
EOF

   close(GDML);
}

####################################
sub gen_enclosureExtras()
{
    $GDML = "LAr1-enclosureExtras" . $suffix . ".gdml";
    push (@gdmlFiles, $GDML); # Add file to list of GDML fragments
    $GDML = ">" . $GDML;
    open(GDML) or die("Could not open file $GDML for writing");

    # Define the solids and structures associated with the "enclosureExtras" (platform, two walkways, insulation, cement stands and foam saddles, tanks on the floor  (7 - 5 big, 2, small) and platform (2), computer racks on top of the platform (10)). Names of particular volumes are explanatory
    print GDML <<EOF;
<?xml version='1.0'?>
<gdml>

</gdml>
EOF

    close(GDML);
}
#####################################3


sub write_fragments()
{
    # The output file is a list of the GDML sub-files created by this
    # script.

    if ( ! defined $output )
    {
	$output = "-"; # write to STDOUT 
    }

    # Set up the output file.
    $OUTPUT = ">" . $output;
    open(OUTPUT) or die("Could not open file $OUTPUT");

    print OUTPUT <<EOF;
<?xml version='1.0'?>

<!-- Input to Geometry/gdml/make_gdml.pl; define the GDML fragments
     that will be zipped together to create a detector description. 
     -->

<config>

   <constantfiles>

      <!-- These files contain GDML <constant></constant>
           blocks. They are read in separately, so they can be
           interpreted into the remaining GDML. See make_gdml.pl for
           more information. 
	   -->
	   
EOF

    foreach $filename (@defFiles)
    {
	print OUTPUT <<EOF;
      <filename> $filename </filename>
EOF
    }

    print OUTPUT <<EOF;

   </constantfiles>

   <gdmlfiles>

      <!-- The GDML file fragments to be zipped together. -->

EOF

    foreach $filename (@gdmlFiles)
    {
	print OUTPUT <<EOF;
      <filename> $filename </filename>
EOF
    }

    print OUTPUT <<EOF;

   </gdmlfiles>

</config>
EOF

    close(OUTPUT);
}
