#include <iostream>
#include "TList.h"
#include "TString.h"
#include "TFormula.h"
#include <vector>
#include <map>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <limits> // std::numeric_limits<>
#include <getopt.h> // getopt_long(), option
#include <stdexcept> // std::runtime_error

bool sci = false;
bool setupChoice = false;
int debug = 0;
TString mySetup;
int prec = std::numeric_limits<double>::max_digits10;

using namespace std;

struct aLoop{
	std::string	varName;
	double varTo;
	double varStep;
};	

struct aSetup{
	std::string	stpName;
	std::string	stpWorldVolume;
	std::string	stpVersion;
};	

typedef std::map<std::string, double> Variables_t;

bool isComment(TString key) {

	int pos1 = key.Index("<!--");
	int pos2 = key.Index("-->");
	bool k=0;
	TString s1="X", s2="X";
	if(pos1>=0 && pos2>pos1){
		s1=key;
		s1.Remove(pos1);
		s2=key;
		s2.Replace(0,pos2+3,"");
	}
	for(int i=0; i<s1.Length(); i++) if(s1[i]==char(9)) {s1.Replace(i,1," ");}
	for(int i=0; i<s2.Length(); i++) if(s2[i]==char(9)) {s2.Replace(i,1," ");}
	k = s1.IsWhitespace() && s2.IsWhitespace();
	return k;
}

void replaceVariable(TString &key, Variables_t const& variables ){

	Variables_t::const_iterator it;
	std::stringstream stream;
	TString varValue;
	int i;
	TString temp;
	stream.precision(prec);
	for (it=variables.begin(); it!=variables.end(); ++it) { 
		TString const& varName = it->first;
		double value = it->second;
		if (!key.Contains("=") || !key.Contains(varName) ) continue;
		if (debug >= 3) {
			std::cout.precision(prec);
			std::cout << "REPL '" << varName << "' => " << prec << std::endl;
		}
		stream.str("");
		stream.clear();	
		if(0) { stream << scientific << value;} else {stream << value;}
		stream >> varValue;
		do {
			i = key.Index(varName);
			if (i>0) {
				temp = key;
				temp.Remove(i);
				temp += varValue;
				key.Replace(0,i+varName.Length(),"");
				temp += key;
				key = temp;
			}
		} while (i>0);
	}
}

void replaceKeyword(TString &key, TString word, TString keyword){

	int pos1 = key.Index(word);
	if(pos1>=1){
		int pos2 = word.Length();
		TString delimiter1 = key[pos1-1];
		TString delimiter2 = key[pos1+pos2];
		if(delimiter1=="\"" && delimiter2=="\"") key.Replace(pos1,pos2,keyword);
	}

}

void getVariable(TString key, Variables_t& variables){

	TString name;

	replaceVariable(key,variables);

	int pos = key.First("\"")+1;
	key.Replace(0,pos,"");
	pos = key.First("\"");
	TString varName = key;
	varName.Remove(pos);
	pos = key.First("\"")+1;
	key.Replace(0,pos,"");
	pos = key.First("\"")+1;
	key.Replace(0,pos,"");
	pos = key.First("\"");
	TString varValue = key;
	varValue.Remove(pos);

	TFormula formula("formula",varValue);
	name = varName;
	double value = formula.Eval(0);

	if (debug >= 2) {
		std::cout.precision(prec);
		std::cout << "VAR " << name << "=\"" << varValue << "\" => " << value << std::endl;
	}
	variables[name.Data()]=value;

//	if(varName.Length()==1) cout << varName << "\t"<< value<< endl;
}

void evalFormulas(TString &key){

	TString storeKey;
	int pos1,pos2;
	const int n=66;
	TString lookFor[n]={"x=\"","y=\"","z=\"","r=\"","rmax=\"","rmin=\"","rmax1=\"","rmin1=\"","rmax2=\"","rmin2=\"",
	"startphi=\"","deltaphi=\"","starttheta=\"","deltatheta=\"","ax=\"","by=\"","cz=\"","dx=\"","dy=\"","dz=\"",
	"zcut1=\"","zcut2=\"","rlo=\"","rhi=\"","alpha=\"","theta=\"","phi=\"","numsides=\"","x1=\"","x2=\"","y1=\"","y2=\"","x3=\"","x4=\"",
	"alpha1=\"","alpha2=\"","inst=\"","outst=\"","lowX=\"","lowY=\"","lowZ=\"","highX=\"","highY=\"","highZ=\"","zOrder=\"","zPosition=\"",
	"xOffset=\"","yOffset=\"","scalingFactor=\"","v1x=\"","v1y=\"","v2x=\"","v2y=\"","v3x=\"","v3y=\"","v4x=\"","v4y=\"","v5x=\"","v5y=\"",
	"v6x=\"","v6y=\"","v7x=\"","v7y=\"","v8x=\"","v8y=\"","dz=\""};



	for(int i=0; i<n; i++)
	if(key.Contains(lookFor[i])){
		storeKey=key;
		pos1 = key.Index(lookFor[i])+lookFor[i].Length();
		key.Replace(0,pos1,"");
		pos2 = key.First("\"");
		TString varFormula = key;
		varFormula.Remove(pos2);
		TFormula formula("",varFormula);
		double value = formula.Eval(0);
		if (debug >= 2) {
			std::cout.precision(prec);
			std::cout << "FML " << lookFor[i] << "='" << varFormula << "' => " << value << std::endl;
		}
		TString varValue;
		std::stringstream stream;
		stream.precision(prec);
		if(sci) { stream << scientific << value;} else {stream << value;}
		stream << "\"";
		stream >> varValue;
		key = storeKey;
		key.Replace(pos1,varFormula.Length()+1,varValue);
	}


}

void makeLoop(TString &key, aLoop &theLoop){
	int pos = key.Index("for")+3;
	key.Replace(0,pos,"");
	pos = key.First("\"")+1;
	key.Replace(0,pos,"");
	pos = key.First("\"");
	TString name = key;
	name.Remove(pos);

	pos = key.Index("to")+2;
	key.Replace(0,pos,"");
	pos = key.First("\"")+1;
	key.Replace(0,pos,"");
	pos = key.First("\"");
	TString value = key;
	value.Remove(pos);

	pos = key.Index("step")+1;
	key.Replace(0,pos,"");
	pos = key.First("\"")+1;
	key.Replace(0,pos,"");
	pos = key.First("\"");
	TString step = key;
	step.Remove(pos);

	theLoop.varName = name;
	theLoop.varTo = value.Atof();
	theLoop.varStep = step.Atof();
}

void getFileNames(TString inName, TString &outName1, TString &outName2){
	int pos = inName.Index("_base");
	TString temp = inName;
	temp.Remove(pos);
	
	if (outName1.IsNull()) outName1 = temp+".gdml";
	if (outName2.IsNull()) outName2 = temp+"_nowires.gdml";
}

void SetDefaultVariables(Variables_t& variables) {
	
	variables["degree"] = M_PI / 180.0;
	
} // SetDefaultVariables()


void getSetup(TString key, aSetup &theSetup){
	int pos;
	TString s;
	pos = key.Index("name=");
	if(pos>=0){
		pos+=6;
		s = key;
		s.Replace(0,pos,"");
		pos = s.Index("\"");
		s.Replace(pos,s.Length()-pos,"");
		theSetup.stpName = s;
	}
	
	pos = key.Index("version=\"");
	if(pos>=0){
		pos+=9;
		s = key;
		s.Replace(0,pos,"");
		pos = s.Index("\"");
		s.Replace(pos,s.Length()-pos,"");
		theSetup.stpVersion = s;
	}

	pos = key.Index("ref=\"");
	if(pos>=0 && key.Contains("world")){
		pos+=5;
		s = key;
		s.Replace(0,pos,"");
		pos = s.Index("\"");
		s.Replace(pos,s.Length()-pos,"");
		theSetup.stpWorldVolume = s;
	}

}

int preparse(TString inName="sbnd_base.gdml", TString outName1="sbnd.gdml", TString outName2=""){
	
	std::cout << "Preparsing '" << inName << "'"
		"\n => '" << outName1 << "' (complete description)";
  if (!outName2.IsNull())
		std::cout << "\n => '" << outName2 << "' (without TPC wires)";
	std::cout << "\n" << std::endl;
  
//	cout << sci << "\t" << prec << endl;

	bool noWiresFiles = 1;
	if(outName2!="") noWiresFiles = 0;

	ofstream output1, output2;
	output1.open(outName1);
	if(!noWiresFiles) output2.open(outName2);
	

	ifstream input(inName);
	TString key;
	Variables_t variables;
	SetDefaultVariables(variables);
	
	std::vector<std::string> loopLine;
	aLoop theLoop;

	std::vector<aSetup> stpList;
	aSetup theSetup;	

	TString name, loopkey;
	bool inLoop = 0;
	bool inSetup = 0;
	bool isWire = 0;
	bool keepLoop = 1;
	do{
		key.ReadLine(input,0);
//		cout << i << endl; i++;
		if ((key.Contains("<variable") && !key.Contains("<!--")) || (key.Contains("<constant") && !key.Contains("<!--"))) {
			getVariable(key,variables);
		} else if (key.Contains("<loop")) {
			inLoop = 1;
			if(key.Contains("#wire")) isWire=1;
			makeLoop(key,theLoop);
		} else if (key.Contains("</loop")) {
			do{
				for (std::vector<std::string>::iterator it = loopLine.begin() ; it != loopLine.end(); ++it){
					loopkey = *it;
					replaceVariable(loopkey,variables);
					if(loopkey.Contains("--")) {
					  if (debug) {
					    std::cout << "checking double minus sign " << std::endl;
					    std::cout << loopkey << std::endl;
					  }
					  loopkey.ReplaceAll("--","+");
					  if (debug) std::cout << loopkey << std::endl; 
					}
					evalFormulas(loopkey);
					output1 << loopkey << endl;
					if(!isWire && !noWiresFiles) output2 << loopkey << endl;
				}			
				variables[theLoop.varName] += theLoop.varStep;
				if(theLoop.varStep>0) {keepLoop = (variables[theLoop.varName] <= theLoop.varTo);}
				else {keepLoop = (variables[theLoop.varName] >= theLoop.varTo);}
			} while (keepLoop);
			inLoop = 0;
			isWire = 0;
			loopLine.clear();
		} else if (inLoop) {
			loopLine.push_back(key.Data());
		} else {
			if( !( isComment(key) || key.IsWhitespace() ) ){
				if(key.Contains("<setup")){
					inSetup = 1;
					getSetup(key,theSetup);
				}
				if(inSetup) getSetup(key,theSetup);
				if(key.Contains("</setup")) { 
					stpList.push_back(theSetup);
					theSetup.stpName="";
					theSetup.stpVersion="";
					theSetup.stpWorldVolume="";
					inSetup = 0;	
				}
				replaceVariable(key,variables);
				evalFormulas(key);
				output1 << key << endl;
				if(!noWiresFiles) output2 << key << endl;
			}
		}

	}while(!input.eof());

	
	output1.close();
	if(!noWiresFiles) output2.close();
	input.close();

	bool goAhead = setupChoice;
	if(setupChoice)	{

		TString ver="1.0";
		int pos=mySetup.Index(":");
		if(pos>=0){
			ver=mySetup;
			ver.Replace(0,pos+1,"");
			mySetup.Replace(pos,mySetup.Length()-pos,"");
		}

		// is the indicated setup among the ones available?
		for (std::vector<aSetup>::iterator itx = stpList.begin() ; itx != stpList.end(); ++itx){

			theSetup=(*itx);
			if( theSetup.stpName==mySetup && theSetup.stpVersion==ver ){ goAhead=1; break;} 
			else {goAhead=0;}
		}
		
		if(goAhead){

			TString inName1, inName2, comm;
			inName1 = outName1+"~";
			inName2 = outName2+"~";
			TString inNameVec[]={inName1,inName2};
			TString outNameVec[]={outName1,outName2};

			const int kmax = 1+(!noWiresFiles);
	
			for(int k=0; k<kmax; k++){
				errno = 0;
				rename(outNameVec[k], inNameVec[k]);
				if (errno != 0) {
				  std::cerr << "Error renaming '" << outNameVec[k] << "' into '" << inNameVec[k] << "': "
				    << strerror(errno) << std::endl;
          return errno;
				}
				output1.open(outNameVec[k]);
				ifstream input1(inNameVec[k]);

				inSetup=0;
				do{
					key.ReadLine(input1,0);
					if((!key.Contains("<setup") && !inSetup) ){
						replaceKeyword(key,"volWorld","volIgnoredOnThisSetup");
						if(!key.Contains("volumeref"))replaceKeyword(key,theSetup.stpWorldVolume,"volWorld");
						output1 << key << endl;
					}else {inSetup=1;}
				}while(!input1.eof());
	
				output1 << "<setup name=\"Default\" version=\"1.0\">" << endl;
				output1 << "\t<world ref=\"volWorld\" />" << endl;
				output1 << "</setup>" << endl;
				output1 << "</gdml>" << endl;
				output1.close();
				input1.close();
				remove(inNameVec[k]);
			}
		}else { cout << "Setup or its version not found. " << endl;}
	}
	return 0;
} // preparse()


//---- Main program ------------------------------------------------------------
# ifndef __CINT__

bool genNoWires = false, printHelp = false;
TString inFile="sbnd_base.gdml", outFile, outFile2, noWires, withWires;

void parseArguments(unsigned int argc, char** argv) {

  static const option longopts[] = {
    { "output",    required_argument, NULL, 'o' },
    { "nowires",   optional_argument, NULL, 'w' },
    { "sci",       optional_argument, NULL, 's' },
    { "prec",      required_argument, NULL, 'p' },
    { "setup",     required_argument, NULL, 'S' },
    { "debug",     optional_argument, NULL, 'd' },
    { "help",      no_argument,       NULL, 'h' },
    { NULL,        0,                 NULL,  0  }
  }; // longopts
  
  // automatically build the short option string from the long one
  std::string shortopts = ":"; // this means no error printout
  for (auto const& longopt: longopts) {
    if (!longopt.name) break;
    if (longopt.val == 0) continue;
    shortopts += (char) longopt.val;
    if (longopt.has_arg != no_argument) shortopts += ':';
    if (longopt.has_arg == optional_argument) shortopts += ':';
  } // for
  
  // ----------------------------------------------------------------------
  // options
  char ch;
  optind = 1;
  while
    ((ch = getopt_long(argc, argv, shortopts.c_str(), longopts, NULL)) != -1)
  {
    switch (ch) {
      case 'o': // -o, --output
        outFile = optarg;
        continue;
      case 'w': // -w, --nowires
        genNoWires = true;
        if (optarg) outFile2 = optarg;
        continue;
      case 'S': // -S, --setup
        setupChoice = true;
        mySetup = optarg;
        continue;
      case 's': // -s, --sci
        sci = true;
        if (optarg) {
          TString par = optarg;
          if (!par.IsDigit())
            throw std::runtime_error("Invalid precision in -s option.");
          prec = par.Atoi();
        }
        continue;
      case 'd': // -d, --debug
        if (optarg) {
          TString par = optarg;
          if (!par.IsDigit())
            throw std::runtime_error("Invalid debug level in -d option.");
          debug = par.Atoi();
        }
        else debug = 1;
        continue;
      case 'p': // -p, --prec
        {
          TString par = optarg;
          if (!par.IsDigit())
            throw std::runtime_error("Invalid precision in -s option.");
          prec = par.Atoi();
        }
        continue;
      case 'h': // -h, --help
        printHelp = true;
        continue;
     } // switch
  } // while
  
  // ----------------------------------------------------------------------
  // arguments
  if (optind < (int) argc) inFile = argv[optind++];
  
  if (optind < (int) argc) {
    std::cerr << "Spurious arguments: '" << argv[optind] << "'";
    if (optind + 1 < (int) argc)
      std::cerr << " and " << (argc - optind - 1) << " more";
    std::cerr << "." << std::endl;
    throw std::runtime_error("Too many arguments on the command line!");
  }
  
} // parseArguments()



int main(int argc, char **argv)
{

	parseArguments(argc, argv);
	
	if (printHelp) {
		std::string programName = argv[0];
		if (programName.rfind('/') != std::string::npos) programName.erase(programName.rfind('/'));
		
		cout <<
		  "\nGDML Preparser v 1.0 (gustavo.valdiviesso@unifal-mg.edu.br)"
		  "\nBasic usage:"
		  "\n  " << programName << " [file_base.gdml] [flags]"
		  "\n"
		  "\nGenerate file.gdml with the following configuration flags:"
		  "\n-w , --nowires[=nowirename]"
		  "\n    Generate file_nowires.gdml, ignoring all <loop> with <!--wire--> comment."
		  "\n-p , --prec=precision"
		  "\n    Set numerical output precision (current precision: " << prec << " digits)"
		  "\n-s , --sci[=precision]"
		  "\n    Uses scientific notation for numbers (optionally alter numerical output precision)"
		  "\n-o , --output=outputname"
		  "\n    Uses this file name as main output name"
		  "\n-S , --setup=setupname"
		  "\n    Uses the specified setup name instead of default one"
		  "\n-d , --debug"
		  "\n    Enable debugging messages"
		  "\n"
		  "\nExamples:"
		  "\n"
		  "\n  " << programName << " geometry_base.gdml --nowires"
		  "\n    Outputs geometry.gdml and geometry_nowires.gdml with standard numeric notation."
		  "\n"
		  "\n  " << programName << " geometry_base.gdml -setup Cryostat:1.0"
		  "\n    Outputs geometry.gdml making setup Cryostat (version 1.0) the volWorld required by LArSoft."
		  "\n"
		  << endl;
		return 0;
	}

	ifstream input(inFile);
	if (!input.is_open()) { 
		std::cerr << "Error opening '" << inFile << "': " << strerror(errno) << std::endl;
		return errno;
	}

	if(inFile.Contains("_base.gdml")) {getFileNames(inFile,withWires,noWires);}
	else {
		withWires = inFile+"_preparsed.gdml";
		noWires = inFile+"_nowires.gdml";
	}
	if(outFile!="") {
		withWires = outFile;
		noWires = outFile;
		int pos = outFile.Index(".gdml");
		noWires.Replace(pos,5,"_nowires.gdml");
	}
	if(outFile2!="") { noWires = outFile2; }
	if(!genNoWires) noWires="";
	
	return preparse(inFile,withWires,noWires);
	
} // main()
#endif // __CINT__
