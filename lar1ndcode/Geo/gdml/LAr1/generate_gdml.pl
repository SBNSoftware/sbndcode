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
	    "wires|w:s" => \$wires);

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

if ( defined $wires )
{
#The user supplied the wires on parameter, so using that. Otherwise the default=1 is used.
$wires_on = $wires;
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
    for ( $i = 0; $i < $NumberWires; ++$i){    
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
</solids>

<structure>
 <volume name="volTPCActive">
   <materialref ref="LAr"/>
   <solidref ref="TPCActive"/>
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
		<position name="posTPCPlane2" unit="cm" x="$TPCWidth/2 - 0.3 -0.1" y="0" z="0" />
     	<rotationref ref="rPlus180AboutY"/> 
	</physvol>
	<physvol>
		<volumeref ref="volTPCPlane"/>
		<position name="posTPCPlane3" unit="cm" x="$TPCWidth/2 - 0.6 -0.1" y="0" z="0" />
	</physvol>

EOF
	print GDML <<EOF;
	<physvol>
	 	 <volumeref ref="volTPCActive"/>
	     <position name="posTPCActive" unit="cm" x="0" y="0" z="0"/>
	</physvol>
EOF

    # Closes TPC volume definition space
    print GDML <<EOF;
 </volume> 
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
   		 <volumeref ref="volTPCFrame"/>
    	 <position name="posTPCFrame" unit="cm" x="$TPCWidth +10" y="0" z="0"/>
   	 </physvol>
   	 <physvol>
    	 <volumeref ref="volTPCFrame"/>
         <position name="posTPCFrame1" unit="cm" x="-$TPCWidth - 10" y="0" z="0"/>
   	 </physvol>
   	 <physvol>
    	<volumeref ref="volTPCHorizontalBeam"/>
        <position name="posTPCHorizontalBeam0" unit="cm" x="0" y="80" z="$TPCLength/2"/>
     </physvol>
     <physvol>
     	 <volumeref ref="volTPCHorizontalBeam"/>
         <position name="posTPCHorizontalBeam1" unit="cm" x="0" y="80*2" z="$TPCLength/2"/>
     </physvol>
     <physvol>
    	 <volumeref ref="volTPCHorizontalBeam"/>
         <position name="posTPCHorizontalBeam2" unit="cm" x="0" y="-80" z="$TPCLength/2"/>
     </physvol>
     <physvol>
     	<volumeref ref="volTPCHorizontalBeam"/>
        <position name="posTPCHorizontalBeam3" unit="cm" x="0" y="-80*2" z="$TPCLength/2"/>
     </physvol>
     <physvol>
     	<volumeref ref="volTPCHorizontalBeam"/>
        <position name="posTPCHorizontalBeam4" unit="cm" x="0" y="0" z="$TPCLength/2"/>
     </physvol>
     <physvol>
    	 <volumeref ref="volTPCHorizontalBeam"/>
         <position name="posTPCHorizontalBeam5" unit="cm" x="0" y="80" z="-$TPCLength/2"/>
     </physvol>
     <physvol>
         <volumeref ref="volTPCHorizontalBeam"/>
         <position name="posTPCHorizontalBeam6" unit="cm" x="0" y="80*2" z="-$TPCLength/2"/>
     </physvol>
     <physvol>
     	 <volumeref ref="volTPCHorizontalBeam"/>
         <position name="posTPCHorizontalBeam7" unit="cm" x="0" y="-80" z="-$TPCLength/2"/>
   	 </physvol>
   	 <physvol>
      	 <volumeref ref="volTPCHorizontalBeam"/>
         <position name="posTPCHorizontalBeam8" unit="cm" x="0" y="-80*2" z="-$TPCLength/2"/>
     </physvol>
     <physvol>
     	<volumeref ref="volTPCHorizontalBeam"/>
        <position name="posTPCHorizontalBeam9" unit="cm" x="0" y="0" z="-$TPCLength/2"/>
     </physvol>
     <physvol>
     	<volumeref ref="volTPCHorizontalBeam"/>
        <position name="posTPCHorizontalBeam10" unit="cm" x="0" y="$TPCHeight/2" z="0"/>
     </physvol>
     <physvol>
     	<volumeref ref="volTPCHorizontalBeam"/>
        <position name="posTPCHorizontalBeam11" unit="cm" x="0" y="$TPCHeight/2" z="80"/>
     </physvol>
     <physvol>
     	<volumeref ref="volTPCHorizontalBeam"/>
        <position name="posTPCHorizontalBeam12" unit="cm" x="0" y="$TPCHeight/2" z="80*2"/>
     </physvol>
     <physvol>
        <volumeref ref="volTPCHorizontalBeam"/>
        <position name="posTPCHorizontalBeam13" unit="cm" x="0" y="$TPCHeight/2" z="-80"/>
     </physvol>
     <physvol>
    	 <volumeref ref="volTPCHorizontalBeam"/>
      	 <position name="posTPCHorizontalBeam14" unit="cm" x="0" y="$TPCHeight/2" z="-80*2"/>
     </physvol>
     <physvol>
    	 <volumeref ref="volTPCHorizontalBeam"/>
      	 <position name="posTPCHorizontalBeam15" unit="cm" x="0" y="-$TPCHeight/2" z="0"/>
     </physvol>
     <physvol>
      	 <volumeref ref="volTPCHorizontalBeam"/>
         <position name="posTPCHorizontalBeam16" unit="cm" x="0" y="-$TPCHeight/2" z="80"/>
     </physvol>
     <physvol>
     	 <volumeref ref="volTPCHorizontalBeam"/>
    	 <position name="posTPCHorizontalBeam17" unit="cm" x="0" y="-$TPCHeight/2" z="80*2"/>
  	 </physvol>
     <physvol>
     	 <volumeref ref="volTPCHorizontalBeam"/>
         <position name="posTPCHorizontalBeam18" unit="cm" x="0" y="-$TPCHeight/2" z="-80"/>
   	 </physvol>
     <physvol>
     	 <volumeref ref="volTPCHorizontalBeam"/>
         <position name="posTPCHorizontalBeam19" unit="cm" x="0" y="-$TPCHeight/2" z="-80*2"/>
     </physvol> 
	<physvol>
		<volumeref ref="volTPC"/>
		<position name="posTPC1" unit="cm" x="$TPCWidth/2" y="0" z="0" />
	</physvol>
	<physvol>
		<volumeref ref="volTPC"/>
		<position name="posTPC2" unit="cm" x="-$TPCWidth/2" y="0" z="0" />
		<rotationref ref="rPlus180AboutY"/>
	</physvol>  

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
 
 <box name="TPCFrameA" lunit="cm" x="$AnodeWidthX" y="$AnodeHeightY" z="$AnodeLengthZ"/>
 <box name="TPCFrameB" lunit="cm" x="$AnodeWidthX+ 0.1" y="$AnodeHeightY-20" z="$AnodeLengthZ -20"/>
 <box name="TPCFrameC" lunit="cm" x="$AnodeWidthX/2" y="$AnodeHeightY-20" z="$AnodeWidthX/2"/>
 <box name="TPCFrameD" lunit="cm" x="$AnodeWidthX/2" y="$AnodeHeightY-35" z="$AnodeWidthX/2"/>

 <box name="TPCHorizontalBeam" lunit="cm" x="$TPCWidth*2" y="5" z="10"/>

 <box name="TPCSideCrossA" lunit="cm" x="$AnodeWidthX/2" y="127.28" z="127.28"/> 
 <box name="TPCSideCrossB" lunit="cm" x="$AnodeWidthX/2+ 0.1" y="127.28-$AnodeWidthX/2" z="127.28-$AnodeWidthX/2"/> 


  <subtraction name="SteelBox">
	<first ref="SteelBoxA"/> <second ref="SteelBoxB"/>
	<position name="posSteelBoxSubtraction" x="0" y="0" z="0"/>
  </subtraction>

   <subtraction name="TPCFrame0">
     <first ref="TPCFrameA"/> <second ref="TPCFrameB"/>
     <position name="posTPCSubtraction" x="0" y="0" z="0"/>
   </subtraction>

	<union name="TPCFrame1">
	  <first ref="TPCFrame0"/> <second ref="TPCFrameC"/>
	  <position name="posTPCUnion0" x="0" y="0" z="0"/>
	</union>

	<union name="TPCFrame2">
	   <first ref="TPCFrame1"/> <second ref="TPCFrameD"/>
	   <position name="posTPCUnion1" x="0" y="0" z="0"/>
	   <rotationref ref="rPlus90AboutX"/>
	</union>
 
   <subtraction name="TPCSideCross">
	 <first ref="TPCSideCrossA"/> <second ref="TPCSideCrossB"/>
	 <position name="posTPCSideCross" x="0" y="0" z="0"/>
   </subtraction>

	<union name="TPCFrame3">
		<first ref="TPCFrame2"/> <second ref="TPCSideCross"/>
		 <position name="posTPCSideCross0" unit="cm" x="0" y="($AnodeHeightY-20)/4" z="($AnodeLengthZ-20)/4"/>
 		 <rotationref ref="rPlus45AboutX"/>
      </union>
	<union name="TPCFrame4">
		<first ref="TPCFrame3"/> <second ref="TPCSideCross"/>
		 <position name="posTPCSideCross1" unit="cm" x="0" y="($AnodeHeightY-20)/4" z="-($AnodeLengthZ-20)/4"/>
 		 <rotationref ref="rPlus45AboutX"/>
      </union>
	<union name="TPCFrame5">
		<first ref="TPCFrame4"/> <second ref="TPCSideCross"/>
		 <position name="posTPCSideCross2" unit="cm" x="0" y="-($AnodeHeightY-20)/4" z="-($AnodeLengthZ-20)/4"/>
 		 <rotationref ref="rPlus45AboutX"/>
      </union>
	<union name="TPCFrame">
		<first ref="TPCFrame5"/> <second ref="TPCSideCross"/>
		 <position name="posTPCSideCross3" unit="cm" x="0" y="-($AnodeHeightY-20)/4" z="($AnodeLengthZ-20)/4"/>
 		 <rotationref ref="rPlus45AboutX"/>
      </union>

</solids>

<structure>
 <volume name="volTPCFrame">
   <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
   <solidref ref="TPCFrame"/>
 </volume>
 <volume name="volTPCHorizontalBeam">
   <materialref ref="G10"/>
   <solidref ref="TPCHorizontalBeam"/>
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
   <position name="posCryostat" unit="cm" x="0" y="0" z="0"/>
  </physvol>
EOF
  if ( $enclosureExtras eq "on" ) {
    print GDML <<EOF;
     <physvol>
        <volumeref ref="volInsulation"/>
        <position name="posInsulation" unit="cm" x="0" y="2" z="0"/>
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

  <subtraction name="ConcreteEnclosure">
	<first ref="ConcreteEnclosureOuter"/> <second ref="ConcreteEnclosureInner"/>
	<position name="posConcreteEnclosureSubtraction" x="0" y="0" z="0"/>
  </subtraction>
</solids>

<structure>
  <volume name="volConcreteEnclosure">
	<materialref ref="Concrete"/>
	<solidref ref="ConcreteEnclosure"/>
  </volume>
  <volume name="volWorld" >
    <materialref ref="Air"/> 
    <solidref ref="World"/>
    <physvol>
      <volumeref ref="volDetEnclosure"/>
      <position name="posDetEnclosure" unit="cm" x="0" y="0" z="0"/>
    </physvol> 
	<physvol>
	  <volumeref ref="volConcreteEnclosure"/>
	  <position name="posConcreteEnclosure" unit="cm" x="-40" y="0" z="0"/>
	 </physvol>


	
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
