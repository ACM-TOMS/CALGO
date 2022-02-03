/*=============================================================================
author        :Walter Schreppers
filename      :ParseCmdLine.h
created       :/
modified      :25/02/2002
version       :4
copyright     :Walter Schreppers
bugreport(log):/
=============================================================================*/


#ifndef IO_MODE
#define IO_MODE
/*=============================================================================
  enumeration to set the IO Mode for Mpieee and Rational
=============================================================================*/
enum ioMode{PARAM,DECIMAL,BINARY,BINREP,HEXREP,FLAGS,RATIONAL};
#endif


#ifndef PARSECMDLINE_H
#define PARSECMDLINE_H

#include "Parser.h" //For getting the enumeration roundType.
                    //Redefining the enumeration would be errorprone (due to redundancy)
class ParseCmdLine{
  public:

    //constructor/destructor
    //======================
    ParseCmdLine();                 //constructor
    ~ParseCmdLine(){};              //destructor


    //members
    //=======
		bool parse(int argsize,char* argtable[]);
		
    bool getNoFor() const;
    bool getConfigFile() const;
    bool getConvertFile() const;
    bool getForceChanges() const;
    bool getDefault() const;
    bool getInit() const;
    bool getNoPrintf() const;
    
    bool getPreParse() const;
    bool getConstants() const;
    
    int getRadix() const;
    int getPrecision() const;
    long getExprangeLow() const;
    long getExprangeUp() const;
    roundType getRounding() const;
    
    char* getInFileName() const;
    char* getOutFileName() const;
    char* getConfFileName() const;
    char* getConvertFileName() const;
    char* getInitFileName() const;
    
    string getUsage() const;
    string getIoModeStr() const;

  private:
  
		//private locals
		//==============
    bool bNoFor, bConf, bInit, bForceChanges, bVerbose, bDefault, bConvert;
    bool bPreParse, bNoPrintf, bConstants;
    char *infile, *outfile, *conffile, *progname, *convertfile, *initfile;
    int fRadix, fPrecision;
    long fExprangeLow, fExprangeUp;
		roundType fRounding;
    vector<ioMode> fMpieee, fRational;

		//private members
		//===============
	  void init();
    bool setRoundType(char);
    bool addIOMpieee(char);
		bool addIORational(char);
    bool setIOMpieee(string);
    bool setIORational(string);

    bool setExponentBits(int);
    long lpow(long,long);
};

#endif

