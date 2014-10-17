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
		"tpb|t:s" => \$tpb) ;

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
    # with the same name as $name. For example, if $name eq "TPCDepth",
    # then the following statement assigns the value to $TPCDepth. The
    # value is in quotes, because some of the parameters have text
    # strings in them (like "$inch").

    eval "\$$name = '$value'";
}

# Our calculations and constants depend on the geometry of the wires.
my $SinUVAngle = sin( deg2rad($UVAngle) );
my $CosUVAngle = cos( deg2rad($UVAngle) );
my $TanUVAngle = tan( deg2rad($UVAngle) );

my $inch=2.54;
my $wires_on=0; 			# turn wires on=1 or off=0
my $tpb_coverage = 0 ;      # multiplier to determine amount fo tpb coverage--default is 0

if ( defined $wires )
{
#The user supplied the wires on parameter, so using that. Otherwise the default=1 is used.
$wires_on = $wires;
}
if ( defined $tpb )
{
$tpb_coverage = $tpb ;
}


my $NumberOfTPCPlanes=3;
my $pmt_switch="on";		#turn on or off depending on pmts wanted
my $enclosureExtras="on"; 

# The routines that create the GDML sub-files. Most of the explanatory
# comments are in gen_defs().
gen_defs();
gen_rotations();
gen_materials();


#if ( $wires eq "on"){
gen_wirevertplane();
gen_wireplane(); 
#}

gen_cathode();		# physical volumes defined in gen_tpc()
gen_tpc();

if ( $pmt_switch eq "on" ) {  gen_pmt();	}	# physical volumes defined in gen_cryostat()
if ( $enclosureExtras eq "on" ) {  gen_enclosureExtras(); } #generation of insulation, etc. will happen if specified
gen_cryostat();

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
    print " 	  -w wires excludes wires from gdml when wires=0, default to wires=1\n";
    print "       -h prints this message, then quits\n";
}


# Create the detector constant file. This file is actually temporary,
# since the make_gdml.pl program will interpret its contents rather
# than include it in the final GDML file.

sub gen_defs()
{
    #TPCWirePlaneLength is the size in the z direction
    #TPCWirePlaneWidth is the size in the y direction
    $TPCWirePlaneLengthZ	=	365  ;
    $TPCWirePlaneWidthY		=	400  ;

    $pi   = pi;
    $inch = 2.54;

    $WorldWidth		=	100.0*$DetEnclosureWidth;
    $WorldHeight	=	100.0*$DetEnclosureHeight;
    $WorldLength	=	100.0*$DetEnclosureLength;

    $CathodePlateDepth	=	0.1;
    $CathodeLengthZ	=	365;
    $CathodeWidthX	=	5;
    $CathodeHeightY	=	400;


    $AnodeLengthZ   =   365;
    $AnodeWidthX    =   20;
    $AnodeHeightY   =   400;

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
   <rotation name="rPlus45AboutX" unit="deg" x="45" y="0" z="0"/>
   <rotation name="rPlus90AboutX" unit="deg" x="90" y="0" z="0"/>
   <rotation name="rPlusUVAngleAboutX" unit="deg" x="30" y="0" z="0"/>
   <rotation name="rPlus180AboutY" unit="deg" x="0" y="180" z="0"/>
   <rotation name="rPMTRotation"  unit="deg" x="90"  y="270"   z="0"/>

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

#=begincomment
sub gen_wireplane()
{
    $GDML = "LAr1-wireplane" . $suffix . ".gdml";
    push (@gdmlFiles, $GDML); # Add file to list of GDML fragments
    $GDML = ">" . $GDML;
    open(GDML) or die("Could not open file $GDML for writing");

    my $TPCYWirePitch = $TPCWirePitch / $SinUVAngle;
    my $TPCZWirePitch = $TPCWirePitch / $CosUVAngle;

  # Calculate the number of wire ends on a given z-edge of the plane.
    my $NumberWiresPerEdge = 0;
    if ( $wires_on == 1 )
      {  $NumberWiresPerEdge = int( $TPCWirePlaneLengthZ / $TPCZWirePitch );}

  # The number of full-length "center" wires.
   my $NumberCenterWires = 0.;
   if ( $wires_on == 1 )
	  { $NumberCenterWires = $UVWireCount - 2*$NumberWiresPerEdge ; }

 #####################
 # Define the solids
  print GDML <<EOF;
<?xml version='1.0'?>
<gdml>
<solids>
EOF


  # End wires 
  for($i = 0; $i < $NumberWiresPerEdge; $i++)
{
  print GDML <<EOF;
<tube name="TPCWire$i"
  rmax="0.5*$TPCWireThickness"
  z="$TPCYWirePitch*($i+1) * 2  - 0.075" 
  deltaphi="360"
  aunit="deg"
  lunit="cm"/> 
EOF
}

   
 print GDML <<EOF;
<tube name="TPCWireCommon"
  rmax="0.5*$TPCWireThickness"
  z="$TPCWirePlaneLengthZ/$SinUVAngle - 0.075 "
  deltaphi="360"
  aunit="deg"
  lunit="cm"/>
<box name="TPCPlane"
  x="$TPCWirePlaneThickness"
  y="$TPCWirePlaneWidthY"
  z="$TPCWirePlaneLengthZ + 0.35"
  lunit="cm"/>
</solids>
<structure>
EOF

    # Wires of varying lengths 
  for ($i = 0; $i < $NumberWiresPerEdge; $i++)
    {
    print GDML <<EOF;
    <volume name="volTPCWire$i">
    <materialref ref="Titanium"/>
    <solidref ref="TPCWire$i"/>
    </volume>
EOF
    }

    # Center wires and plane
   print GDML <<EOF;
    <volume name="volTPCWireCommon">
      <materialref ref="Titanium"/>
      <solidref ref="TPCWireCommon"/>
    </volume>
    <volume name="volTPCPlane">
      <materialref ref="LAr"/>
      <solidref ref="TPCPlane"/>
EOF

  # The wires at the -z, +y end (For +60 deg-- can rotate by 180 later for -60)
  for ($i = 0; $i < $NumberWiresPerEdge ; $i++)
  {

    print GDML <<EOF;
  <physvol>
     <volumeref ref="volTPCWire$i"/> 
     <position name="posTPCWire$i" unit="cm" y="0.5*$TPCWirePlaneWidthY - 0.5*$TPCYWirePitch*($i+1)" z="-0.5*$TPCWirePlaneLengthZ+0.5*$TPCZWirePitch*($i+1)" x="0"/>
     <rotationref ref="rPlusUVAngleAboutX"/> 
    </physvol>  
EOF
  $ypos=0.5*$TPCWirePlaneWidthY - 0.5*$TPCYWirePitch*($i+1);
  $zpos=-0.5*$TPCWirePlaneLengthZ+0.5*$TPCZWirePitch*($i+1);
  open (MYFILE, '>>data.txt');
  print MYFILE "TPCWire$i y=$ypos z=$zpos\n";
  }

  # The wires in the middle.
for ($i = 0; $i < $NumberCenterWires ; $i++)
  {
      my $j = $NumberWiresPerEdge  +$i;
      $ypos=0.5*$TPCWirePlaneWidthY - $TPCYWirePitch*(0.5*$NumberWiresPerEdge + $i+1) ; 

      print GDML <<EOF;
   <physvol>
     <volumeref ref="volTPCWireCommon"/>
     <position name="posTPCWire$j" unit="cm" y="$ypos" z="0" x="0"/>
     <rotationref ref="rPlusUVAngleAboutX"/>
    </physvol> 
EOF
  }

  # The wires at the +z end
  for ($i = 0; $i < $NumberWiresPerEdge; $i++)
  {

	  my $j = $NumberWiresPerEdge + $NumberCenterWires + $i ;
	  my $k = $NumberWiresPerEdge - $i - 1 ;
      $ypos =0.5*$TPCWirePlaneWidthY - 0.5*$TPCYWirePitch*($NumberWiresPerEdge + 2*$NumberCenterWires + $i +1 ) ; 
	
    print GDML <<EOF;
   <physvol>
     <volumeref ref="volTPCWire$k"/> 
     <position name="posTPCWireB$j" unit="cm" y="$ypos" z="0.5*$TPCZWirePitch*($i+1)" x="0"/>
     <rotationref ref="rPlusUVAngleAboutX"/>
    </physvol> 
EOF
  }

      print GDML <<EOF;
  </volume>

    <volume name="volTPCPlane2">
      <materialref ref="LAr"/>
      <solidref ref="TPCPlane"/>
EOF

  # The wires at the -z, +y end (For +60 deg-- can rotate by 180 later for -60)
  for ($i = 0; $i < $NumberWiresPerEdge ; $i++)
  {

    print GDML <<EOF;
  <physvol>
     <volumeref ref="volTPCWire$i"/> 
     <position name="posTPCWireVol2$i" unit="cm" y="0.5*$TPCWirePlaneWidthY - 0.5*$TPCYWirePitch*($i+1)" z="-0.5*$TPCWirePlaneLengthZ+0.5*$TPCZWirePitch*($i+1)" x="0"/>
     <rotationref ref="rPlusUVAngleAboutX"/> 
    </physvol>  
EOF
  $ypos=0.5*$TPCWirePlaneWidthY - 0.5*$TPCYWirePitch*($i+1);
  $zpos=-0.5*$TPCWirePlaneLengthZ+0.5*$TPCZWirePitch*($i+1);
  open (MYFILE, '>>data.txt');
  print MYFILE "TPCWire$i y=$ypos z=$zpos\n";
  }

  # The wires in the middle.
for ($i = 0; $i < $NumberCenterWires ; $i++)
  {
      my $j = $NumberWiresPerEdge  +$i;
      $ypos=0.5*$TPCWirePlaneWidthY - $TPCYWirePitch*(0.5*$NumberWiresPerEdge + $i+1) ; 

      print GDML <<EOF;
   <physvol>
     <volumeref ref="volTPCWireCommon"/>
     <position name="posTPCWireVol2$j" unit="cm" y="$ypos" z="0" x="0"/>
     <rotationref ref="rPlusUVAngleAboutX"/>
    </physvol> 
EOF
  }

  # The wires at the +z end
  for ($i = 0; $i < $NumberWiresPerEdge; $i++)
  {

	  my $j = $NumberWiresPerEdge + $NumberCenterWires + $i ;
	  my $k = $NumberWiresPerEdge - $i - 1 ;
      $ypos =0.5*$TPCWirePlaneWidthY - 0.5*$TPCYWirePitch*($NumberWiresPerEdge + 2*$NumberCenterWires + $i +1 ) ; 
	
    print GDML <<EOF;
   <physvol>
     <volumeref ref="volTPCWire$k"/> 
     <position name="posTPCWireBVol2$j" unit="cm" y="$ypos" z="0.5*$TPCZWirePitch*($i+1)" x="0"/>
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
#=endcomment
#=cut


#=begin comment
sub gen_wirevertplane()
{
  my $NumberWires = 0; 

 if ( $wires_on == 1 )  
    { $NumberWires = int( $TPCWirePlaneLengthZ / $TPCWirePitch );} 


    $GDML = "LAr1-wirevertplane" . $suffix . ".gdml";
    push (@gdmlFiles, $GDML); # Add file to list of GDML fragments
    $GDML = ">" . $GDML;
    open(GDML) or die("Could not open file $GDML for writing");

    #Define the solids and structures: the wires and the TPC wire plane.
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
    for ( $i = 0; $i < $NumberWires; ++$i) {    
    print GDML <<EOF;
      <physvol>
        <volumeref ref="volTPCWireVert"/>
        <position name="posTPCWireVert$i" unit="cm" z="-0.5*$TPCWirePlaneLengthZ+$TPCWirePitch*($i+1)" x="0" y="0"/>
        <rotationref ref="rPlus90AboutX"/>
      </physvol>  
EOF
	}

	print GDML <<EOF ;
  </volume>
  <volume name="volTPCPlaneVert2">
    <materialref ref="LAr"/>       
    <solidref ref="TPCPlaneVert"/>
EOF
    for ( $i = 0; $i < $NumberWires; ++$i){    
    print GDML <<EOF;
      <physvol>
        <volumeref ref="volTPCWireVert"/>
        <position name="posTPCWireVertVol2$i" unit="cm" z="-0.5*$TPCWirePlaneLengthZ+$TPCWirePitch*($i+1)" x="0" y="0"/>
        <rotationref ref="rPlus90AboutX"/>
      </physvol>  
EOF
	}

	print GDML <<EOF ;
  </volume>

</structure>
</gdml>
EOF
	close(GDML);
}

#=end comment
#=cut

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

close(GDML);

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
    $TPCActiveHeight = $TPCWirePlaneWidthY-0.2;   #  
    $TPCActiveLength = $TPCWirePlaneLengthZ-0.2;   # extra subtraction to arrive at TPCActive values in the TDR

#The addition of 0.35 below is to avoid overlap extrusions
    print GDML <<EOF;
<?xml version='1.0'?>
<gdml>
<solids>
 <box name="TPC" lunit="cm" x="$TPCWidth" y="$TPCHeight" z="$TPCLength+0.35"/>
 <box name="TPCActive" lunit="cm" x="$TPCActiveDepth" y="$TPCActiveHeight" z="$TPCActiveLength"/>

 <box name="TPBLayerXY" lunit="cm" x="$tpb_coverage*$TPCActiveDepth" y="$tpb_coverage*$TPCActiveHeight" z="0.01"/>
 <box name="TPBLayerXZ" lunit="cm" x="$tpb_coverage*$TPCActiveDepth" y="0.01" z="$tpb_coverage*$TPCActiveLength"/>
 <box name="TPBLayerCathode" lunit="cm" x="0.01" y="$tpb_coverage*$TPCActiveHeight" z="$tpb_coverage*$TPCActiveLength"/>
 
</solids>

<structure>
 <volume name="volTPCActive">
   <materialref ref="LAr"/>
   <solidref ref="TPCActive"/>
 </volume>
 <volume name="volTPCActive2">
   <materialref ref="LAr"/>
   <solidref ref="TPCActive"/>
 </volume>
 <volume name="volTPBLayerXY">
    <materialref ref="LAr"/>
    <solidref ref="TPBLayerXY"/>
 </volume>
 <volume name="volTPBLayerXZ">
    <materialref ref="LAr"/>
    <solidref ref="TPBLayerXZ"/>
 </volume>
 <volume name="volTPBLayerCathode">
    <materialref ref="LAr"/>
    <solidref ref="TPBLayerCathode"/>
 </volume>
 <volume name="volTPC">
   <materialref ref="LAr"/>
   <solidref ref="TPC"/>
EOF

     print GDML <<EOF;
	<physvol>
		 <volumeref ref="volTPCPlaneVert"/>
		 <position name="posTPCPlaneVert2" unit="cm" x="$TPCWidth/2-0.1" y="0" z="0" />
	 </physvol>
	<physvol>
		<volumeref ref="volTPCPlane"/>
		<position name="posTPCPlane2" unit="cm" x="$TPCWidth/2 - 0.6 -0.1" y="0" z="0" />
     	<rotationref ref="rPlus180AboutY"/> 
	</physvol>
	<physvol>
		<volumeref ref="volTPCPlane"/>
		<position name="posTPCPlane3" unit="cm" x="$TPCWidth/2 - 0.3 -0.1" y="0" z="0" />
	</physvol>

EOF

if( $tpb_coverage != 0 ){ 
    print GDML <<EOF ;
    <physvol>
        <volumeref ref="volTPBLayerXY"/>
        <position name="posTPBLayerXY0" unit="cm" x="0" y="0" z="$TPCActiveLength/2-0.005"/>
    </physvol>
    <physvol>
        <volumeref ref="volTPBLayerXY"/>
        <position name="posTPBLayerXY1" unit="cm" x="0" y="0" z="-$TPCActiveLength/2 + 0.005"/>
    </physvol>
    <physvol>
        <volumeref ref="volTPBLayerXZ"/>
        <position name="posTPBLayerXZ0" unit="cm" x="0" y="$TPCActiveHeight/2 - 0.005" z="0"/>
    </physvol>
    <physvol>
        <volumeref ref="volTPBLayerXZ"/>
        <position name="posTPBLayerXZ1" unit="cm" x="0" y="-$TPCActiveHeight/2+0.005" z="0"/>
    </physvol>

	<physvol>
		<volumeref ref="volTPBLayerCathode"/>
		<position name="posTPBLayerCathode" unit="cm" x="-$TPCActiveDepth/2+0.005" y="0" z="0"/>
	</physvol>

EOF
}

	print GDML <<EOF;
	<physvol>
	 	 <volumeref ref="volTPCActive"/>
	     <position name="posTPCActive1" unit="cm" x="0" y="0" z="0"/>
	</physvol>
EOF

    # Closes TPC volume definition space
    print GDML <<EOF;
 </volume> 

<!-- <volume name="volTPC2">
   <materialref ref="LAr"/>
   <solidref ref="TPC"/>
EOF

     print GDML <<EOF;
	<physvol>
		 <volumeref ref="volTPCPlaneVert2"/>
		 <position name="posTPCPlaneVert3" unit="cm" x="-$TPCWidth/2 +0.1" y="0" z="0" />
	 </physvol>
	<physvol>
		<volumeref ref="volTPCPlane2"/>
		<position name="posTPCPlane4" unit="cm" x="-$TPCWidth/2 + 0.3 +0.1" y="0" z="0" />
     	<rotationref ref="rPlus180AboutY"/> 
	</physvol>
	<physvol>
		<volumeref ref="volTPCPlane2"/>
		<position name="posTPCPlane5" unit="cm" x="-$TPCWidth/2 + 0.6 +0.1" y="0" z="0" />
	</physvol>

EOF
	print GDML <<EOF;
	<physvol>
	 	 <volumeref ref="volTPCActive2"/>
	     <position name="posTPCActive2" unit="cm" x="0" y="0" z="0"/>
	</physvol>
</volume>-->

</structure>
</gdml>
EOF

   close(GDML);
}

sub make_APA()
{

#Place only physical volumes here
#This goes inside volCryostat
	
	print CRYOSTAT <<EOF;
 	 <physvol>
   		 <volumeref ref="volAPAFrame"/>
    	 <position name="posAPAFrame" unit="cm" x="$TPCWidth +10" y="0" z="0"/>
   	 </physvol>
   	 <physvol>
    	 <volumeref ref="volAPAFrame"/>
         <position name="posAPAFrame1" unit="cm" x="-$TPCWidth - 10" y="0" z="0"/>
   	 </physvol>
   	 <physvol>
    	<volumeref ref="volHorizontalBeam"/>
        <position name="posHorizontalBeam0" unit="cm" x="0" y="80" z="$TPCLength/2+5"/>
     </physvol>
     <physvol>
     	 <volumeref ref="volHorizontalBeam"/>
         <position name="posHorizontalBeam1" unit="cm" x="0" y="80*2" z="$TPCLength/2+5"/>
     </physvol>
     <physvol>
    	 <volumeref ref="volHorizontalBeam"/>
         <position name="posHorizontalBeam2" unit="cm" x="0" y="-80" z="$TPCLength/2+ 5"/>
     </physvol>
     <physvol>
     	<volumeref ref="volHorizontalBeam"/>
        <position name="posHorizontalBeam3" unit="cm" x="0" y="-80*2" z="$TPCLength/2+5"/>
     </physvol>
     <physvol>
     	<volumeref ref="volHorizontalBeam"/>
        <position name="posHorizontalBeam4" unit="cm" x="0" y="0" z="$TPCLength/2 + 5"/>
     </physvol>
     <physvol>
    	 <volumeref ref="volHorizontalBeam"/>
         <position name="posHorizontalBeam5" unit="cm" x="0" y="80" z="-$TPCLength/2 - 5"/>
     </physvol>
     <physvol>
         <volumeref ref="volHorizontalBeam"/>
         <position name="posHorizontalBeam6" unit="cm" x="0" y="80*2" z="-$TPCLength/2 - 5"/>
     </physvol>
     <physvol>
     	 <volumeref ref="volHorizontalBeam"/>
         <position name="posHorizontalBeam7" unit="cm" x="0" y="-80" z="-$TPCLength/2 - 5"/>
   	 </physvol>
   	 <physvol>
      	 <volumeref ref="volHorizontalBeam"/>
         <position name="posHorizontalBeam8" unit="cm" x="0" y="-80*2" z="-$TPCLength/2 - 5"/>
     </physvol>
     <physvol>
     	<volumeref ref="volHorizontalBeam"/>
        <position name="posHorizontalBeam9" unit="cm" x="0" y="0" z="-$TPCLength/2 - 5 "/>
     </physvol>
     <physvol>
     	<volumeref ref="volHorizontalBeam"/>
        <position name="posHorizontalBeam10" unit="cm" x="0" y="$TPCHeight/2 + 5" z="0"/>
     </physvol>
     <physvol>
     	<volumeref ref="volHorizontalBeam"/>
        <position name="posHorizontalBeam11" unit="cm" x="0" y="$TPCHeight/2 +5" z="80"/>
     </physvol> 
     <physvol>
     	<volumeref ref="volHorizontalBeam"/>
        <position name="posHorizontalBeam12" unit="cm" x="0" y="$TPCHeight/2+5" z="80*2"/>
     </physvol>
    <physvol>
        <volumeref ref="volHorizontalBeam"/>
        <position name="posHorizontalBeam13" unit="cm" x="0" y="$TPCHeight/2+5" z="-80"/>
     </physvol> 
     <physvol>
    	 <volumeref ref="volHorizontalBeam"/>
      	 <position name="posHorizontalBeam14" unit="cm" x="0" y="$TPCHeight/2+5" z="-80*2"/>
     </physvol>
     <physvol>
    	 <volumeref ref="volHorizontalBeam"/>
      	 <position name="posHorizontalBeam15" unit="cm" x="0" y="-$TPCHeight/2-5" z="0"/>
     </physvol>
     <physvol>
      	 <volumeref ref="volHorizontalBeam"/>
         <position name="posHorizontalBeam16" unit="cm" x="0" y="-$TPCHeight/2-5" z="80"/>
     </physvol>
     <physvol>
     	 <volumeref ref="volHorizontalBeam"/>
    	 <position name="posHorizontalBeam17" unit="cm" x="0" y="-$TPCHeight/2-5" z="80*2"/>
  	 </physvol>
     <physvol>
     	 <volumeref ref="volHorizontalBeam"/>
         <position name="posHorizontalBeam18" unit="cm" x="0" y="-$TPCHeight/2-5" z="-80"/>
   	 </physvol>
     <physvol>
     	 <volumeref ref="volHorizontalBeam"/>
         <position name="posHorizontalBeam19" unit="cm" x="0" y="-$TPCHeight/2-5" z="-80*2"/>
     </physvol>   
	<physvol>
		<volumeref ref="volTPC"/>
		<position name="posTPC1" unit="cm" x="$TPCWidth/2+2.5" y="0" z="0" />
	</physvol>
	<physvol>
		<volumeref ref="volTPC"/>
		<position name="posTPC2" unit="cm" x="-$TPCWidth/2-2.5" y="0" z="0" />
		<rotationref ref="rPlus180AboutY"/> 
	</physvol>  

EOF

}

# Generates Ben Jones's PMT micro-pmtdef (with temporary edit to ellipsoid shapes
sub gen_pmt {

    $PMT = "lar1nd-lightguidedef" . $suffix . ".gdml";
    push (@gdmlFiles, $PMT); # Add file to list of GDML fragments
    $PMT = ">" . $PMT;
    open(PMT) or die("Could not open file $PMT for writing");

  #The below pmt geo is a placeholder--"PMTVolume" contains pmt info from microboone--ahack
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


	print PMT <<EOF;
<solids>
 <box name="lightguidebar"
  lunit="cm"
  x="0.25*2.54"
  y="100"
  z="2.54"/>
 <box name="lightguidebarcoating"
  lunit="cm"
  x=".05"
  y="100"
  z="2.54"/>
 <box name="lightguidedetector"
  lunit="cm"
  x="0.25*2.54+0.1"
  y="100"
  z="2.54"/>
</solids>
<structure>
 <volume name="volOpDetSensitive">
  <materialref ref="LAr"/>
  <solidref ref="lightguidebarcoating"/>
 </volume>
 <volume name="vollightguidebar">
  <materialref ref="Acrylic"/>
  <solidref ref="lightguidebar"/>
 </volume>
 <volume name="vollightguidedetector">
  <materialref ref="LAr"/>
  <solidref ref="lightguidedetector"/>
  <physvol>
   <volumeref ref="volOpDetSensitive"/>
   <position name="posOpDetSensitive" unit="cm" x="(0.125*2.54)+0.025" y="0" z="0"/>
  </physvol>
  <physvol>
   <volumeref ref="vollightguidebar"/>
   <position name="poslightguidebar" unit="cm" x="0" y="0" z="0"/>
  </physvol>
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
 <box name="Cryostat" lunit="cm" x="$CryostatWidth+0.1" y="$CryostatHeight+0.1" z="$CryostatLength+0.1" /> 
 <box name="SteelBoxA" lunit="cm" x="$CryostatWidth" y="$CryostatHeight" z="$CryostatLength"/>
 <box name="SteelBoxB" lunit="cm" x="$CryostatWidth-5" y="$CryostatHeight-5" z="$CryostatLength-5"/>


 <box name="APAFrameA" lunit="cm" x="$AnodeWidthX" y="$AnodeHeightY" z="$AnodeLengthZ"/>
 <box name="APAFrameB" lunit="cm" x="$AnodeWidthX+ 0.1" y="$AnodeHeightY-20" z="$AnodeLengthZ -20"/>
 <box name="APAFrameC" lunit="cm" x="$AnodeWidthX/2" y="$AnodeHeightY-20" z="$AnodeWidthX/2"/>
 <box name="APAFrameD" lunit="cm" x="$AnodeWidthX/2" y="$AnodeHeightY-35" z="$AnodeWidthX/2"/>

 <box name="HorizontalBeam" lunit="cm" x="$TPCWidth*2" y="5" z="5"/>

 <box name="APASideCrossA" lunit="cm" x="$AnodeWidthX/2" y="127.28" z="127.28"/> 
 <box name="APASideCrossB" lunit="cm" x="$AnodeWidthX/2+ 0.1" y="127.28-$AnodeWidthX/2" z="127.28-$AnodeWidthX/2"/> 


  <subtraction name="SteelBox">
	<first ref="SteelBoxA"/> <second ref="SteelBoxB"/>
	<position name="posSteelBoxSubtraction" x="0" y="0" z="0"/>
  </subtraction>

   <subtraction name="APAFrame0">
     <first ref="APAFrameA"/> <second ref="APAFrameB"/>
     <position name="posTPCSubtraction" x="0" y="0" z="0"/>
   </subtraction>

	<union name="APAFrame1">
	  <first ref="APAFrame0"/> <second ref="APAFrameC"/>
	  <position name="posTPCUnion0" x="0" y="0" z="0"/>
	</union>

	<union name="APAFrame2">
	   <first ref="APAFrame1"/> <second ref="APAFrameD"/>
	   <position name="posTPCUnion1" x="0" y="0" z="0"/>
	   <rotationref ref="rPlus90AboutX"/>
	</union>
 
   <subtraction name="APASideCross">
	 <first ref="APASideCrossA"/> <second ref="APASideCrossB"/>
	 <position name="posAPASideCross" x="0" y="0" z="0"/>
   </subtraction>

	<union name="APAFrame3">
		<first ref="APAFrame2"/> <second ref="APASideCross"/>
		 <position name="posAPASideCross0" unit="cm" x="0" y="($AnodeHeightY-20)/4" z="($AnodeLengthZ-20)/4"/>
 		 <rotationref ref="rPlus45AboutX"/>
      </union>
	<union name="APAFrame4">
		<first ref="APAFrame3"/> <second ref="APASideCross"/>
		 <position name="posAPASideCross1" unit="cm" x="0" y="($AnodeHeightY-20)/4" z="-($AnodeLengthZ-20)/4"/>
 		 <rotationref ref="rPlus45AboutX"/>
      </union>
	<union name="APAFrame5">
		<first ref="APAFrame4"/> <second ref="APASideCross"/>
		 <position name="posAPASideCross2" unit="cm" x="0" y="-($AnodeHeightY-20)/4" z="-($AnodeLengthZ-20)/4"/>
 		 <rotationref ref="rPlus45AboutX"/>
      </union>
	<union name="APAFrame">
		<first ref="APAFrame5"/> <second ref="APASideCross"/>
		 <position name="posAPASideCross3" unit="cm" x="0" y="-($AnodeHeightY-20)/4" z="($AnodeLengthZ-20)/4"/>
 		 <rotationref ref="rPlus45AboutX"/>
      </union>

</solids>

<structure>
 <volume name="volAPAFrame">
   <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
   <solidref ref="APAFrame"/>
 </volume>
 <volume name="volHorizontalBeam">
   <materialref ref="G10"/>
   <solidref ref="HorizontalBeam"/>
 </volume>
 <volume name="volSteelBox">
   <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
   <solidref ref="SteelBox"/>
 </volume>

 <volume name="volCryostat">
   <materialref ref="LAr"/>
   <solidref ref="Cryostat"/>
	<physvol>
	  <volumeref ref="volSteelBox"/>
	  <position name="posSteelBox" unit="cm" x="0" y="0" z="0"/>
	</physvol>  
    <physvol> 
       <volumeref ref="volCathodePlate"/>
 	   <position name="posCathodePlate" unit="cm" x="0" y="0" z="0"/>
	 </physvol> 
EOF
	make_APA();


	@pmt_pos = ( ' x="-220" y="0 " z="0" ',
				 ' x="-220" y="-(400-20)/4 " z="(365 - 20)/4" ',
				 ' x="-220" y="-(400-20)/4 " z="(365 - 20)/4" ',
				 ' x="-220" y="(400-20)/4 " z="-(365 - 20)/4" ',
				 ' x="-220" y="-(400-20)/4 " z="-(365 - 20)/4" ' );

#		 <position name="posAPASideCross0" unit="cm" x="0" y="($AnodeHeightY-20)/4" z="($AnodeLengthZ-20)/4"/>

    if ( $pmt_switch eq "on" ) {
	  for( $ii=0; $ii < 3; ++$ii){
 	 print CRYOSTAT <<EOF ;
  <!--	  <physvol>
	  <volumeref ref="volPMT"/>
	  <position name="posPMT$ii" unit="cm" @pmt_pos[$ii]/> 
	  <rotationref ref="rPMTRotation"/>
	  </physvol> -->
EOF
	}

}


=begin
#Commented out temporarily to do some pmt positioning

	for ( $i=0; $i<125; ++$i ){
	    for ($j=0; $j<4; ++$j ){
		for ($k=0; $k<2; ++$k){  
		    print CRYOSTAT <<EOF;
		    <physvol>
			<volumeref ref="vollightguidedetector"/>
			<position name="lightguidedetector-$i-$j-$k" unit="cm" x="(-1)**$k*($TPCWidth+20+(0.25*2.54+0.1)/2)" y="(2*$j+1)*$TPCHeight/8-$TPCHeight/2" z="-($TPCLength/2-20)+2.6*($i+1)" />
EOF
			if ( $k==0 ) {
                            print CRYOSTAT <<EOF;
                            <rotationref ref="rPlus180AboutY"/>
EOF
                        }  
                   print CRYOSTAT<<EOF;
                   </physvol>
EOF
        }
      }
    }
=end
=cut
 
# }

	print CRYOSTAT <<EOF;
 </volume>
</structure>
</gdml>
EOF

   close(CRYOSTAT);
}


# Parameterize the cryostat's surroundings.
sub gen_enclosure()
{
    # Set up the output file.
    $GDML = "LAr1-enclosure" . $suffix . ".gdml";
    push (@gdmlFiles, $GDML); # Add file to list of GDML fragments
    $GDML = ">" . $GDML;
    open(GDML) or die("Could not open file $GDML for writing");


#Note for below--To understand coordinate placements for Cryostat etc, 
#see page 20&21 of Dec 31, 2013 LAr1ND proposal
    print GDML <<EOF;
<?xml version='1.0'?>
<gdml>
<solids>
 <box name="DetEnclosure" lunit="cm" x="$DetEnclosureWidth+0.1" y="$DetEnclosureHeight+0.1" z="$DetEnclosureLength+0.1" />
</solids>

<structure>
 <volume name="volDetEnclosure">
  <materialref ref="Air"/>
  <solidref ref="DetEnclosure"/>
  <physvol>
   <volumeref ref="volCryostat"/>
   <position name="posCryostat" unit="cm" x="40" y="0" z="0"/>
  </physvol>
EOF
  if ( $enclosureExtras eq "on" ) {
    print GDML <<EOF;
     <physvol>
        <volumeref ref="volInsulation"/>
        <position name="posInsulation" unit="cm" x="40" y="2" z="0"/>
      </physvol> 
EOF
	}

#if the extra pieces in the enclosure are desired, place within volEnclosure(insulation, racks, etc)

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
	#Ground info from Fig4.1(Right) from LAr1ND Dec 2013 proposal(for now)
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
  <box name="ConcreteEnclosureOuter" lunit="cm" x="$DetEnclosureWidth+8" y="$DetEnclosureHeight+8" z="$DetEnclosureLength+8"/>
  <box name="ConcreteEnclosureInner" lunit="cm" x="$DetEnclosureWidth+0.1" y="$DetEnclosureHeight+0.1" z="$DetEnclosureLength+0.1"/>
  <box name="GroundOuterTemp" lunit="cm" x="$DetEnclosureWidth + 8 + 100" y="(28*12)*2.54" z="$DetEnclosureLength+8+100"/>
  <box name="GroundOuter" lunit="cm" x="$DetEnclosureWidth + 8 + 1000" y="(28*12)*2.54+1000" z="$DetEnclosureLength+ 8 + 2500" />
  <box name="GroundInner" lunit="cm" x="$DetEnclosureWidth + 8.1" y="(28*12)*2.54 +5" z="$DetEnclosureLength+8.1"/>


  <subtraction name="ConcreteEnclosure">
	<first ref="ConcreteEnclosureOuter"/> <second ref="ConcreteEnclosureInner"/>
	<position name="posConcreteEnclosureSubtraction" unit="cm" x="0" y="0" z="0"/>
  </subtraction>

  <subtraction name="GroundTemp">
	<first ref="GroundOuterTemp"/> <second ref="GroundInner"/>
	<position name="posGroundTempSubtraction" unit="cm" x="0" y="0" z="0"/>
  </subtraction>

  <subtraction name="Ground">
	<first ref="GroundOuter"/> <second ref="GroundInner"/>
	<position name="posGroundSubtraction" unit="cm" x="0" y="0" z="859.05"/>
  </subtraction>

</solids>

<structure>
  <volume name="volConcreteEnclosure">
	<materialref ref="Concrete"/>
	<solidref ref="ConcreteEnclosure"/>
  </volume>
  <volume name="volGround">
	<materialref ref="Dirt"/>
	<solidref ref="Ground"/>
  </volume>
  <volume name="volTempGround">
	<materialref ref="Dirt"/>
	<solidref ref="GroundTemp"/>
  </volume>
  <volume name="volWorld" >
    <materialref ref="Air"/> 
    <solidref ref="World"/>
    <physvol>
      <volumeref ref="volDetEnclosure"/>
      <position name="posDetEnclosure" unit="cm" x="-40" y="0" z="($TPCLength)/2"/>
    </physvol> 
	<physvol>
	  <volumeref ref="volConcreteEnclosure"/>
	  <position name="posConcreteEnclosure" unit="cm" x="-40" y="0" z="($TPCLength)/2"/>
	 </physvol>
	<physvol>
	  <volumeref ref="volGround"/>
	  <position name="posGround" unit="cm" x="-40" y="0" z="($TPCLength)/2-859.05"/>
	</physvol>

<!--	<physvol>
	  <volumeref ref="volTempGround"/>
	  <position name="posGroundTemp" unit="cm" x="-40" y="0" z="($TPCLength)/2"/>
	</physvol> -->
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

# Define the solids and structures associated with the "enclosureExtras" 
 print GDML <<EOF;
<?xml version='1.0'?>
<gdml>
<solids>
  <box name="InsulationOuter" lunit="cm" x="600" y="550" z="490"/>
  <box name="InsulationInner" lunit="cm" x="$CryostatWidth+0.1" y="$CryostatHeight+0.1" z="$CryostatLength+0.1"/>

   <subtraction name="Insulation">
     <first ref="InsulationOuter"/> <second ref="InsulationInner"/>
     <position name="posInsulationSubtraction" x="0" y="-20" z="0"/>
   </subtraction>

</solids>

<structure>
  <volume name="volInsulation">
	<materialref ref="PU_foam_light"/>
	<solidref ref="Insulation"/>
  </volume>
</structure>
</gdml>
EOF

	close(GDML);
}




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
